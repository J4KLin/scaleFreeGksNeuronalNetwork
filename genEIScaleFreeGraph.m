function [eiMat,eiLabels]= genEIScaleFreeGraph(node_ct,varargin)
% GENEISCALEFREEGRAPH create a n x n directed weighted adjacency matrix 
% where 'to be excitatory' nodes branch out following a scale free rule of 
% connections while 'to be inhibitory' nodes branch out randomly.  The scale 
% free rule is an implementation of the Barabasi-Albert Linearized Chord 
% Diagram. *Note that this function only generates a network graph so the
% nodes do not inherently possess the ability to be 'excitatory' or
% 'inhibitory.'
%
% INPUTS:
%   node_ct: [int] number of nodes
%
% VARARGIN:
%   'inhibitory_per' : [0:1] Percentage of the nodes to be assigned as
%   'inhibitory' to follow a random branching pattern.  (Default) is 0.
%
%   'in_deg_per': [0:1] Setting the "in_degree" percentage.  The effects
%   will be most prominent in high degree nodes ('hubs').  (Default) is 0.5.
%
%   'iterations' : [int] Number of passthroughs of the Linearlized Chord Diagram
%   algorithm with each passthrough generating an extra branch for every
%   node in the existing network. *Must be greater than 1. (Default) is 15.
%
%   'self_connect' : [boolean] Allow for self connecting branches.
%   (Default) is false.
%
% OUTPUTS:
%   'eiMat' : n x n directed weighted adjacency matrix of connections.
%   Weighted since multi-connections are allowed.  y-axis represents from
%   node i and x-axis represents to node j.
%
%   'eiLables' : n x 1 label for excitatory "e" and inhibitory "i" nodes

%% AUTHOR       : Jack Lin
%% VERSION      : 1.0
%% TESTED       : (R2019a)

%%
p = inputParser; 
p.addRequired('node_ct',@(x) floor(x) == x);
p.addParameter('inhibitory_per',0,@(x)(x>= 0)&&(x<=1));
p.addParameter('in_deg_per',.5,@(x)(x>=0)&&(x<=1));
p.addParameter('iterations',15,@(x) isinteger(x) && x > 1);
p.addParameter('self_connect',false,@isboolean);
p.parse(node_ct,varargin{:});

node_ct = p.Results.node_ct;
inhibitory_per = p.Results.inhibitory_per;
in_deg_per = p.Results.in_deg_per;
iterations = p.Results.iterations;
self_connect = p.Results.self_connect;

eiMat = zeros(node_ct);
eiMat(1,1) = 1;
excitatory_nct = node_ct - floor(inhibitory_per*node_ct);   %# of excitatory neurons

eiLabels = repmat(categorical("e"),1,node_ct);
eiLabels(excitatory_nct+1:end) = categorical("i");

node_grab_bag = [1,1];      %choose which node to connect to 

%First passthrough of LCD with each newly generated node getting a new
%outward connection
for new_node = 2:node_ct
    temp_grab_bag = [node_grab_bag, new_node];
    
    if(new_node > excitatory_nct)
        connecting_node = randperm(new_node,1);           %inhibitory cells branches out randomly
    else
        random_choice_idx = randperm(length(temp_grab_bag),1);  %excitatory cells branches out following sf law
        connecting_node = temp_grab_bag(1,random_choice_idx);   %node to which the new node will connect to
        eiMat(connecting_node, new_node) = eiMat(connecting_node,new_node) + 1; %e-e is matrix for percentage adjustment
    end
    
    eiMat(new_node,connecting_node) = eiMat(new_node,connecting_node) + 1;
    node_grab_bag = [node_grab_bag, connecting_node, new_node];
end

%Rerun LCD through each of the existing nodes for 'iterations' number of
%passthroughs
%IMPORTANT: Node A is allowed to connect to Node B multiple times
for it = 2:iterations
    for node = 1:node_ct
        if(node > excitatory_nct)
            connecting_node = randperm(node_ct,1);
        else
            randperms = randperm(length(node_grab_bag));
            random_choice_idx = randperms(1);
            connecting_node = node_grab_bag(1,random_choice_idx);
            if(connecting_node <= excitatory_nct)
                eiMat(connecting_node,node) = eiMat(connecting_node,...
                node) + 1;
            end
        end
        eiMat(node,connecting_node) = eiMat(node,...
            connecting_node) + 1;
        node_grab_bag = [node_grab_bag,connecting_node,node];
    end
end

%remove connections to self
if(~self_connect)
    eiMat = eiMat- diag(diag(eiMat));
end

%Scale Free In Degree Adjustment
e_idxs = find(eiLabels == "e");
eMat = eiMat(e_idxs,e_idxs);
for row = 1:size(eMat,1)    %working with upper diagonal
    for col = row:size(eMat,2)
        link_ct = eMat(row,col);
        link_prob = rand(1,link_ct);
        eMat(row,col) = sum(link_prob >= in_deg_per);
        eMat(col,row) = sum(link_prob < in_deg_per);
    end
end
eiMat(e_idxs,e_idxs) = eMat;
