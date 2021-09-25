Code for publication (undergoing review): "Dynamical Mechanism Underlying Scale-Free Network Reorganization in Low Acetylcholine States Corresponding with Slow Wave Sleep".

*Matlab codes tested in version R2019a

CODE(s):
	testSimulation : Example script for running a full single simulation of a network of G_ks modulated HH-like neurons.
	
	genEIScaleFreeGraph : Generate a connectivity/adjacency matrix of scale free excitatory connections and random inhibitory connections (if specified).  This function also allows for the adjustment of the network's 'in-degree percentage.'
	
	simGksSFNeuronalNet : Run a full single simulation of G_ks modulated HH-like neurons given a defined G_ks value and a connectivity graph.  This function also allows for hub neuron (highest node degree neuron) removal and activaition of spike timing dependent plasticity (STDP).
	
	Iapp_by_freq : Utility/helper function for finding input current associated with a specified G_ks value and a single neuron spiking freuqency given all the sets of input-response (IF) curves.
	

DATAFILE(s):
	gks_IF_responses.mat :  Current-frequency response curve data of single G_ks modulated neurons for varying values of G_ks in the range of [0:.05:1.5].  This is used by the 'simGksSFNeuronalNet' function in conjunction with the 'Iapp_by_freq' function to find the specific input current associated with a desired firing frequency.
		'dc_input_vals': All the tested current inputs
		'gks_input_vals': All the tested G_ks values
		'gks_IF_responses': [m x 3] Each column respectively representing the G_ks, current input and spiking frequency.