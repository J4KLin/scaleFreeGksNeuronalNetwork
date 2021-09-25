%Given a set of IF curves of varying gKs parameters, estimate the
%Input currents for each IF curve respectively necessary to 
%simulate a defined firing rate
% INPUTS:
%   gks_IF_responses: [m x 3] Columns respectively gKs, input current and
%   spiking frequency.
%
%   desired_freq: [num] Specified neuron spiking frequency (Hz)
%
%   desired_gks: [num] Specified gKs
%
% OUTPUTS:
%   Iapp_results : Resulted input current

%% AUTHOR       : Jack Lin
%% VERSION      : 1.0
%% TESTED       : (R2019a)

function [Iapp_results] = Iapp_by_freq(gks_IF_responses,desired_freq,desired_gks)

gks_input_vals = round(unique(gks_IF_responses(:,1)),4);
dc_input_vals = unique(gks_IF_responses(:,2));
input_dt = dc_input_vals(2) - dc_input_vals(1);

IF_responses = reshape(gks_IF_responses(:,3),[length(dc_input_vals),length(gks_input_vals)]);
    
if(exist('desired_gks'))
    IF_responses = IF_responses(:,find(gks_input_vals == desired_gks));
end

Iapp_results = zeros(size(IF_responses,2),1);

for IF = 1:length(Iapp_results)
    cur_if = IF_responses(:,IF);
    lower_bound_idx = 1;
    for i = 1:length(cur_if)
        if(cur_if(i) > desired_freq)
            break;
        else
            lower_bound_idx = i;
        end
    end
    lower_bound = cur_if(lower_bound_idx);
    upper_bound = cur_if(find(cur_if >= desired_freq,1));
    %pause;
    if(lower_bound == upper_bound)
        Iapp_results(IF,1) = dc_input_vals(lower_bound_idx);
    else
        slope = (upper_bound-lower_bound)/input_dt;
        intercept = lower_bound-(slope*dc_input_vals(lower_bound_idx));
        Iapp_results(IF,1) = (desired_freq-intercept)/slope;
    end
end
