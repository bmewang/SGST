function [FF, Readout] = F(shear, Rot, experiment)
% F returns the deformation gradient belonging to a given experimental
% protocol plus a corresponding readout matrix for later projection. 
% Possible experimental protocol "Exp" are as follows:
%       Simple Shear
%           1:   (fs)       2:  (sn)        3:  (fn)
%           4:   (sf)       5:  (ns)        6:  (nf)
%       Biaxial Stretch(direction of readout is capitalized)
%           7:  (Fn)        8:  (fN)        9:  (Fs)
%           10: (fS)        11: (Sn)        12: (sN)
%   Biaxial Tests have to be supplemented with a stretch ratio in the form 
%   "experiment=[8 0.5]" which,  multiplied with the actual stretch, yields
%   stretch for the second index.

FF = eye(3); Readout = zeros(3);

switch experiment(1)
    case 1
        FF(2,1) = shear; Readout(1,1) = 1; Readout(2,2)=1;
    case 2
        FF(3,2) = shear; Readout(1,2) = 1; Readout(2,3)=1;
    case 3
        FF(3,1) = shear; Readout(1,1) = 1; Readout(2,3)=1;
    case 4
        FF(1,2) = shear; Readout(1,2) = 1; Readout(2,1)=1;
    case 5
        FF(2,3) = shear; Readout(1,3) = 1; Readout(2,2)=1;
    case 6
        FF(1,3) = shear; Readout(1,3) = 1; Readout(2,1)=1;
    case 7
        FF(1,1) = 1+shear;        
        FF(3,3) = 1+shear*experiment(2);
        FF(2,2) = 1/(FF(1,1)*FF(3,3));
        Readout(1,1)= 1; Readout(2,1)=1; Readout(3,2)=1;
    case 8
        FF(1,1) = 1+shear;        
        FF(3,3) = 1+shear*experiment(2);
        FF(2,2) = 1/(FF(1,1)*FF(3,3));
        Readout(1,3)= 1; Readout(2,3)=1; Readout(3,2)=1;
    case 9
        FF(1,1) = 1+shear;        
        FF(2,2) = 1+shear*experiment(2);
        FF(3,3) = 1/(FF(1,1)*FF(2,2));
        Readout(1,1)= 1; Readout(2,1)=1; Readout(3,3)=1;
    case 10
        FF(1,1) = 1+shear;        
        FF(2,2) = 1+shear*experiment(2);
        FF(3,3) = 1/(FF(1,1)*FF(2,2));
        Readout(1,2)= 1; Readout(2,2)=1; Readout(3,3)=1;
    case 11
        FF(2,2) = 1+shear;        
        FF(3,3) = 1+shear*experiment(2);
        FF(1,1) = 1/(FF(2,2)*FF(3,3));
        Readout(1,2)= 1; Readout(2,2)=1; Readout(3,1)=1;
    case 12
        FF(2,2) = 1+shear;        
        FF(3,3) = 1+shear*experiment(2);
        FF(1,1) = 1/(FF(2,2)*FF(3,3));
        Readout(1,3)= 1; Readout(2,3)=1; Readout(3,1)=1;
end

FF  = Rot'*FF*Rot;
Readout = Readout*Rot;
end
