%% Example Simulation

neuron_ct = 250;        %# of total neurons
in_percentage = .5;     %In degree perentage (0 to 1)
inhibitory_percentage = 0;  %Set % neurons to inhibitory (0 to 1)

%Generate Scale Free Network
[ei_graph,ei_labels] = genEIScaleFreeGraph(neuron_ct,...
    'inhibitory_per',inhibitory_percentage,...
    'in_deg_per',in_percentage);

%Plotting Node Degree Information
% figure();
% subplot(1,2,1);
% plot(sort(sum(ei_graph,2)' + sum(ei_graph,1),'descend'));
% xlim([0, neuron_ct]);
% title('Node Degree');
% ylabel('Degree','FontSize',12);
% xlabel('Neuron ID','FontSize',12);
% subplot(1,2,2);
% plot((sum(ei_graph,1)./(sum(ei_graph,2)' + sum(ei_graph,1))));
% xlim([0, neuron_ct]);
% title('In Degree Percentage');
% ylabel('In Degree Percentage','FontSize',12);
% xlabel('Neuron ID','FontSize',12);

%Run Simulation 
gks = 0;                  %mS/cm^2 Limited to values [0:.05:1.5]
[time_vec,activity_data,spike_data,other_data] = simGksSFNeuronalNet(gks,ei_graph,ei_labels);

%% Example Simulation with Hub ablation
neuron_ct = 250;        %# of total neurons
in_percentage = .5;     %In degree perentage (0 to 1)
inhibitory_percentage = 0;  %Set % neurons to inhibitory (0 to 1)
hub_del_ct = 25;
hub_del_start_time = 1000;
run_time = 2000;

%Generate Scale Free Network
[ei_graph,ei_labels] = genEIScaleFreeGraph(neuron_ct,...
    'inhibitory_per',inhibitory_percentage,...
    'in_deg_per',in_percentage);

%Run Simulation 
gks = 0;                  %mS/cm^2 Limited to values [0:.05:1.5]
[time_vec,activity_data,spike_data,other_data] = simGksSFNeuronalNet(gks,ei_graph,...
    ei_labels,'hub_del_ct',hub_del_ct,'run_time',run_time,...
    'hub_del_start_time',hub_del_start_time);
