function y=sigmaNGSTProjected(N, StructureTensorList, PrincipleDirection, a, b, Rot, Experiment,Stretch)
% Takes function "sigmaNGST" and projects it corresponding to a given
% experiment "Experiment", allowing for fibre families to vary from the principle
% direction "PrincipleDirection" of the experimental setup via rotation matrix "Rot".  
% "N" denotes to which order the nGST-model is evaluated, demanding an equal amount 
% of structure tensors "StructureTensorList" which have to be calculated beforehand 
% using function "MN". Material parameters "a" and "b" resemble fibre 
% stiffness as within the HO-model. "Stretch" refers to biaxial stretch or shear
% rate depending on the experiment.

    [FF, Readout] = F(Stretch, Rot, Experiment);
    sigma = sigmaNGST(N, StructureTensorList, PrincipleDirection, a, b, FF);
    
y = Readout(1,:)*sigma*Readout(2,:)' - Readout(3,:)*sigma*Readout(3,:)';
