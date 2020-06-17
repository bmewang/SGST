function y = sigmaHOProjected(Parameters, Rot, experiment, stretch )
% Read out function for "sigmaHO" corresponding to a specific Experiment 
% "experiment". For definitions of the different options look at function "F".
% Allows for fibre families to vary from the principle
% direction of the experimental setup via rotation matrix "Rot". 

    [FF, Readout] = F(stretch,Rot, experiment);
    SigmaPre = sigmaHO(Parameters,FF);

y = Readout(1,:)*SigmaPre*Readout(2,:)'-Readout(3,:)*SigmaPre*Readout(3,:)';

