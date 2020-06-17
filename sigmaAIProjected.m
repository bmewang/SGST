function y=sigmaAIProjected(N, rho, PrincipleDirection, a, b, Rot, Experiment,Stretch)
% Read out function for "sigmaAI" corresponding to a specific experiment 
% "Experiment". For definitions of the different options look at function "F".
% Allowing for fibre families to vary from the principle
% direction of the experimental setup via rotation matrix "Rot".
%   "N"         = Discretization for AI in both dimensions
%   "rho"       = von Mises dispersion coefficient
%   "Stretch"   = shear/strain depending on experiment
%   "a,b"       = Holzapfel-Ogden parameters for fiber family

    [FF, Readout] = F(Stretch, Rot, Experiment);
    sigma = sigmaAI(N, rho, PrincipleDirection, a,b, FF);
    
y = Readout(1,:)*sigma*Readout(2,:)' - Readout(3,:)*sigma*Readout(3,:)';
