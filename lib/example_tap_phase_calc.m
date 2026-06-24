%An example of running the tap changer and phase shifter analysis  on the 
%IEEE 300 bus system with two tap changers and two phase shifters


%% Define case data
file = 'case300';

%% 1. Define the transformers
%tap_changers_data = [line of insertion, control node, control voltage]
tap_changers_data = [
    1,9001,.95
    2,9005, 1
    ];

%phase_shifters_data = [line of insertion, control power]
phase_shifters_data = [
    100 , 0
    110 , 1
    ];

%% Run the analysis
[results, stevilo_iteracij, success] = taps_and_phases_analysis(file, tap_changers_data, phase_shifters_data);

