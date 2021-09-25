 function [time_vec,activity_data,spike_data,other_data]...
      = simGksSFNeuronalNet(gks,neuron_graph,varargin)
% SIMGKSSFNEURONALNET Run a full simluation of a network of Hodgkin-Huxley 
% like neuron conductance based model with an M-type potassium current.
% Current model allows for hub removal (top node degree neurons) and spike
% timing dependent plasticity (STDP).
%
% INPUTS:
%   gks : [0:.05:1.5] (mS/cm2) Slow potassium conductance where 0
%   corresponds to High ACh and 1.5 corresponds to Low ACh.
%
%   neruon_graph : [n x n] directed matrix of connections between neurons
%
%   ei_labels : [1 x n] categorical vectors denoting whether a node is an
%   excitatory neuron ("e") or an inhibitory neuron ("i")
%
% VARARGIN:
%   'run_time' : [num] (ms) Total runtime of the simulation.
%
%   'gsyn' : [1 x 4] (mS/cm2) Synaptic conductance between types of neurons
%   [E-E,E-I,I-E,I-I].
%
%   'noise_amp' : [num] (microA) Amplitude of a random current pulse.
%
%   'noise_duration' : [num] (ms) Duration of a random current pulse.
%
%   'noise_probability' : [0:1] Probability of initiating a random current
%   pulse at any given time point.
%
%   'hub_del_ct' : [0:n] Number of top node degree neurons to ablate
%
%   'hub_del_start_time' : [num] (ms) Time at which the number of specified
%   hub neurons will be removed.
%
%   'STDP_start_time' : [num] (ms) Time at which spike timing dependent
%   plasticity will begin. By default, STDP is turned off (-1).
%
% OUTPUTS:
%   time_vec : [1 x m] (ms) Vector of time steps of the simulation
%
%   activity_data : [n x m] (mV) neuron x time point voltage acitvity
%
%   spike_data : [n x m] neuron x time point of spike detections (voltage
%   above a specified threshold)
%
%   other_data : [struct]
%       dt : time step
%       neuron_graph0 : neuron_graph before simulation
%       neuron_graph  : neuron_graph after simulation
%       gsyn_matrix0  : synaptic conducatance graph before simulation
%       gysn_matrix   : synaptic conductance graph after simulation

%% AUTHOR       : Jack Lin
%% VERSION      : 1.0
%% TESTED       : (R2019a)

%%
neuron_ct = size(neuron_graph,1);
               
defaultEILabels = repmat(categorical("e"),1,neuron_ct);
defaultRuntime = 2000;  %ms
defaultGsyn = [.04,.04,.04,.01];
defaultNoiseAmp = .7;   %(mA/s2)
defaultNoiseDur = 2;    %ms
defaultNoiseProb = .02;
defaultHubDel = 0;    
defaultHubDelStartTime = defaultRuntime+1;  %ms
defaultSTDPStartTime = -1;  %ms

validEILabels = @(x) isvector(x) && length(x) == size(neuron_graph,1) &&...
    length(x) == sum(ismember(x,{'e','i'}));

p = inputParser; 
p.addRequired('gks',@(x) ismember(x,[0:.05:1.5]));
p.addRequired('neuron_graph',@(x) (isnumeric(x) || islogical(x)) && size(x,1) == size(x,2));
p.addOptional('ei_labels',defaultEILabels,validEILabels);
p.addParameter('run_time',defaultRuntime,@(x) x>0);
p.addParameter('gsyn',defaultGsyn,@(x) isvector(x) && isnumeric(x));
p.addParameter('noise_amp',defaultNoiseAmp,@(x) isscalar(x) && isnumeric(x));
p.addParameter('noise_duration',defaultNoiseDur,@(x) isscalar(x) && x>=0);
p.addParameter('noise_probability',defaultNoiseProb,@(x) isscalar(x) && x>=0);
p.addParameter('hub_del_ct',defaultHubDel,@(x) isscalar(x) && x>=0 && x<=size(neuron_graph,1));
p.addParameter('hub_del_start_time',defaultHubDelStartTime,@(x) isscalar(x) && x>=0);
p.addParameter('STDP_start_time',defaultSTDPStartTime,@(x) isscalar(x) && x>=0);
p.parse(gks,neuron_graph,varargin{:});

ei_labels = p.Results.ei_labels;
run_time = p.Results.run_time;
gsyn = p.Results.gsyn;
noise = p.Results.noise_amp;
noise_duration = p.Results.noise_duration;
noise_probability = p.Results.noise_probability;
hub_del_ct = p.Results.hub_del_ct;
hub_del_start_time = p.Results.hub_del_start_time;
STDP_start_time = p.Results.STDP_start_time;
STDP_on = false;
if(STDP_start_time >= 0 && STDP_start_time <= run_time); STDP_on = true; end

%%%OTHER DEFAULT PARAMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = 0;  %firing threshold to qualify spiking
dt = .1;        %step time (ms)
%%%SDTP PARMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STDP_rate = 20 * 10^-4;     %learning rate (Al) (mS/cm^2)
tau_STDP = 10;              %learning delay constant (ms)
Lmaxwindow = 40;            %learning delay cutoff window (ms)
Lminwindow = 0;

%%%GRAPH PARAMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_idxs = find(ei_labels == 'e');
i_idxs = find(ei_labels == 'i');
neuron_graph(neuron_graph > 0) = 1;     %strip off multiple similar connections

%%%GRAPH SORT
graph_sort_indeg = false;   %determine if you want to sort graph by # of in degrees
if(graph_sort_indeg)
    [~,graph_mapping] = sort(sum(neuron_graph(:,e_idxs),1),'descend');
else
    [~,graph_mapping] = sort(sum(neuron_graph(:,e_idxs),1)+ sum(neuron_graph(e_idxs,:),2)','descend');
end
graph_mapping = [graph_mapping,i_idxs];
neuron_graph = neuron_graph(graph_mapping,graph_mapping);

[~,totdeg_mapping] = sort(sum(neuron_graph(:,e_idxs),1)+ sum(neuron_graph(e_idxs,:),2)','descend');
hub_neuron_idxs = totdeg_mapping(1:hub_del_ct);

%%%Gsyn synaptic strength mapping
gsyn_matrix = ones(neuron_ct);
gsyn_matrix(e_idxs,e_idxs) = gsyn_matrix(e_idxs,e_idxs) .* gsyn(1,1);
gsyn_matrix(e_idxs,i_idxs) = gsyn_matrix(e_idxs,i_idxs) .* gsyn(1,2);
gsyn_matrix(i_idxs,e_idxs) = gsyn_matrix(i_idxs,e_idxs) .* gsyn(1,3);
gsyn_matrix(i_idxs,i_idxs) = gsyn_matrix(i_idxs,i_idxs) .* gsyn(1,4);

gsyn_matrix0 = gsyn_matrix;
neuron_graph0 = neuron_graph;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_vec = [0:dt:run_time];
activity_data = zeros(neuron_ct,length(time_vec));
spike_data = activity_data;     %Binary form of activity data only focusing on spiking

%%%NEURON PROPERTIES & STATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*Matrix is faster than table
neuron_properties = [];
%starting voltage parameters
Vo = [-60,-50];              %[min_vo,max_vo] Starting mem. potential for all neurons (mV)
min_Vo = Vo(1);
max_Vo = Vo(2);
if(min_Vo>max_Vo); warning('min_Vo cannot be greater than max_Vo'); end
activity_data(:,1) = (rand(neuron_ct,1)*(max_Vo-min_Vo))+min_Vo;

%input current parameters (1)
Ic_idx = 1;
load('IF_responses_05res.mat','gks_IF_responses');
Driving_freq = [0 0];       %[min_ic,max_ic] Base DC to all neurons (mA/s2)
min_Ic = Iapp_by_freq(gks_IF_responses,Driving_freq(1),gks);
max_Ic = Iapp_by_freq(gks_IF_responses,Driving_freq(2),gks);
if(min_Ic>max_Ic); warning('min_ic cannot be greater than max_ic'); end
neuron_properties(:,Ic_idx) = (rand(neuron_ct,1)*(max_Ic-min_Ic))+min_Ic;

%starting gating variables (h=2,n=3,z=4)
h_idx = 2;
n_idx = 3;
z_idx = 4;
vo_vec = activity_data(:,1);
neuron_properties(:,h_idx) = 1./(1+exp((vo_vec+53)./7));        %h steady state
neuron_properties(:,n_idx) = 1./(1+exp((-vo_vec-30)./10));      %n steady state
neuron_properties(:,z_idx) = 1./(1+exp((-vo_vec-39)./5));       %z steady state

%gKs value (5)
gks_idx = 5;
neuron_properties(:,gks_idx) = gks * ones(neuron_ct,1);

%Esyn synaptic reversal potential (6);
esyn_idx = 6;
e_reversal = 0;         %reversal potential for excitatory neuron (mv)
i_reversal = -75;       %reversal potential for inhibitory neuron (mv)
neuron_properties(e_idxs,esyn_idx) = e_reversal * ones(length(e_idxs),1);
neuron_properties(i_idxs,esyn_idx) = i_reversal * ones(length(i_idxs),1);

%neuron type (7)
is_excitatory_idx = 7;
neuron_properties(e_idxs,is_excitatory_idx) = ones(length(e_idxs),1);
neuron_properties(i_idxs,is_excitatory_idx) = zeros(length(i_idxs),1);

%last spike time paramters (last_spike=8,spike_flag=9)
last_spike_idx = 8;
spike_flag_idx = 9;
neuron_properties(:,last_spike_idx) = -1000*ones(neuron_ct,1);
neuron_properties(:,spike_flag_idx) = zeros(neuron_ct,1);

%noise current parameters
last_noise_idx = 10;
neuron_properties(:,last_noise_idx) = -1*ones(neuron_ct,1);

%%%MODEL PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ionic channel parameters
gNa = 24;
gKd = 3;
gL = .02;
ENa = 55;
EK = -90;
EL = -60;
tauZ = 75;

tauD = .5; %synaptic decay constant (ms)
tauR = .2; %synaptic rise constant (ms)

%%%ODE for conductance based neuron and gating function
dX_vectorized = @(v,h,n,z,gks,cur_input,Isyn, Inoise) [...
     ((-gNa.*((1./(1+exp((-v-30)./9.5))).^3).*h.*(v-ENa))- (gKd.*(n.^4).*(v-EK))-...   %dV
     (gks.*z.*(v-EK))- (gL.*(v-EL)) + cur_input - Isyn + Inoise),...
     (((1./(1+exp((v+53)./7)))-h)./(.37+(2.78./(1+exp((v+40.5)./6))))),...             %dXh
     (((1./(1+exp((-v-30)./10)))-n)./(.37+(1.85./(1+exp((v+27)./15))))),...            %dXn
     (((1./(1+exp((-v-39)./5)))-z)./tauZ)];                                            %dXz
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%RUN SIMULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rep_esyn = repmat(neuron_properties(:,esyn_idx),[1,neuron_ct]);     %used for vectorized Isyn calc
one_vec = ones(1,neuron_ct);        %used for vectorized Isyn calc
one_mat = ones(neuron_ct);

hubdel_flag = 0;
if(hub_del_ct > 0); hubdel_flag = 1; end

noise_activity = activity_data;
for idx = 1:length(time_vec)-1
    cur_t = dt*idx;
    v_vec = activity_data(:,idx);
    updated_spike_time = neuron_properties(:,last_spike_idx);
    updated_spike_flag = neuron_properties(:,spike_flag_idx);
    
    %%%handle hub deletion
    if(hubdel_flag && cur_t > hub_del_start_time)
        hubdel_flag = 0;
        neuron_graph(hub_neuron_idxs,:) = 0;
        neuron_graph(:,hub_neuron_idxs) = 0;
        gsyn_matrix(hub_neuron_idxs,:) = 0;
        gsyn_matrix(:,hub_neuron_idxs) = 0;
    end
        
    %%%handle noise
    cur_noise_vec = zeros(neuron_ct,1);
    random_noise_vec = rand(neuron_ct,1);
    last_spike_vec = cur_t - neuron_properties(:,last_noise_idx);
    
    cur_noise_vec(last_spike_vec < noise_duration,1) = noise;
    
    newly_init_noise_vec = last_spike_vec >= noise_duration &...
        random_noise_vec <= noise_probability;
    cur_noise_vec(newly_init_noise_vec,1) = noise;
    neuron_properties(newly_init_noise_vec,last_noise_idx) = cur_t;  

    noise_activity(:,idx) = cur_noise_vec;
    
    %%%handle Isyn calculations
    Isyn_vec = zeros(neuron_ct,1);
    cur_t_vec = one_vec .* cur_t;
    driving_forces = v_vec' - rep_esyn;         %column is driving_forces to a given neuron
    
    adjust_delays = exp(-(cur_t_vec-neuron_properties(:,last_spike_idx))./tauD)-...
        exp(-(cur_t_vec-neuron_properties(:,last_spike_idx))./tauR);
    
    Isyn_vec = sum(driving_forces .* gsyn_matrix .* adjust_delays .* neuron_graph,1);
    Isyn_vec = Isyn_vec';
    
    %%%Runge Kutta 4th order for one step (vectorized across neurons)
    gKs_vec = neuron_properties(:,gks_idx);
    h_vec = neuron_properties(:,h_idx);
    n_vec = neuron_properties(:,n_idx);
    z_vec = neuron_properties(:,z_idx);
    cur_input_vec = neuron_properties(:,Ic_idx);
    
    dx1 = dX_vectorized(v_vec,h_vec,n_vec,z_vec,gKs_vec,cur_input_vec,Isyn_vec,cur_noise_vec);
    K1 =dx1(:,1);
    hk1=dx1(:,2);
    nk1=dx1(:,3);
    zk1=dx1(:,4);

    dx2 = dX_vectorized(v_vec+((dt/2).*K1), h_vec+((dt/2).*hk1), n_vec+((dt/2).*nk1), z_vec+((dt/2).*zk1)...
        ,gKs_vec,cur_input_vec,Isyn_vec,cur_noise_vec);
    K2 =dx2(:,1);
    hk2=dx2(:,2);
    nk2=dx2(:,3);
    zk2=dx2(:,4);

    dx3 = dX_vectorized(v_vec+((dt/2).*K2), h_vec+((dt/2).*hk2), n_vec+((dt/2).*nk2), z_vec+((dt/2).*zk2)...
        ,gKs_vec,cur_input_vec,Isyn_vec,cur_noise_vec);
    K3 =dx3(:,1);
    hk3=dx3(:,2);
    nk3=dx3(:,3);
    zk3=dx3(:,4);

    dx4 = dX_vectorized(v_vec+dt.*K3, h_vec+dt.*hk3, n_vec+dt.*nk3, z_vec+dt.*zk3...
        ,gKs_vec,cur_input_vec,Isyn_vec,cur_noise_vec);
    K4 =dx4(:,1);
    hk4=dx4(:,2);
    nk4=dx4(:,3);
    zk4=dx4(:,4);

    newV_vec = v_vec + (1/6)*dt.*(K1+(2.*K2)+(2.*K3)+K4);
    neuron_properties(:,h_idx) = h_vec + (1/6)*dt.*(hk1+(2.*hk2)+(2.*hk3)+hk4);
    neuron_properties(:,n_idx) = n_vec + (1/6)*dt.*(nk1+(2.*nk2)+(2.*nk3)+nk4);
    neuron_properties(:,z_idx) = z_vec + (1/6)*dt.*(zk1+(2.*zk2)+(2.*zk3)+zk4);
    activity_data(:,idx+1) = newV_vec;
    
    %update spiking time
    neuron_properties(newV_vec < threshold & updated_spike_flag == 2,spike_flag_idx) = 0;
    neuron_properties(newV_vec >= threshold & updated_spike_flag == 1,spike_flag_idx) = 2;
    newly_spike_vec = newV_vec >= threshold & updated_spike_flag == 0;
    neuron_properties(newly_spike_vec,spike_flag_idx) = 1;
    neuron_properties(newly_spike_vec,last_spike_idx) = cur_t;
    spike_data(newly_spike_vec,idx+1) = 1;    
    
    
    if(STDP_on && cur_t >= STDP_start_time && sum(newly_spike_vec)>0)
        delta_spike_times = neuron_properties(:,last_spike_idx) - neuron_properties(:,last_spike_idx)';
        sign_mat = one_mat;
        sign_mat(delta_spike_times<0) = -1;
        
        sign_mat(~newly_spike_vec,~newly_spike_vec) = 0;        %ignore all neurons either not spiking nor connected with spiking neurons
        sign_mat(abs(delta_spike_times) > Lmaxwindow) = 0;
        sign_mat(abs(delta_spike_times) < Lminwindow) = 0;
        
        delta_gsyn = neuron_graph .* sign_mat .* (STDP_rate.*exp(-abs(delta_spike_times)./tau_STDP));
        gsyn_matrix(e_idxs,e_idxs) = gsyn_matrix(e_idxs,e_idxs) + delta_gsyn(e_idxs,e_idxs);
        gsyn_matrix(gsyn_matrix < 0) = 0;
    end
end

%Other return data
other_data = struct();
other_data.dt = dt;
other_data.e_idxs = e_idxs;
other_data.i_idxs = i_idxs;
other_data.neuron_graph0 = neuron_graph0;
other_data.neuron_graph = neuron_graph;
other_data.gsyn_matrix0 = gsyn_matrix0;
other_data.gsyn_matrix = gsyn_matrix;
