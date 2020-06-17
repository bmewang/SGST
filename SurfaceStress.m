function y = SurfaceStress( Model, Parameters, Exp, Stretch, Rotation)
% SurfaceStress( Model, Parameters, Exp, Stretch, Rotation) returns real 
% Cauchy Stress for a given experimental protocol "Exp" projected onto the 
% respective readout-surface. The protocol may be altered by a rotation
% "Rotation" of the prepared tissue.
%   Implemented Models "Model" are
%       1. Holzapfel-Ogden Model ("HO") –– 8 Constitutive parameters
%       2. Generalized Structure Tensor Model of different orders
%       ("0GST"/"1GST"/"2GST"/"3GST"/"4GST")
%       3. Squared Generalized Structure Tensor Model of different
%       orders ("0SGST"/"1SGST"/"2SGST"/"3SGST")
%       4. Angular Integration Model ("AI"): This scheme is a numerical approximation
%       using primitive trapezoidal rule and obligates you to support the 
%       amount of integration steps N.
%   Please note that all dispersion based models expect 8 constitutive
%   parameters + fibre and sheet (optional) distribution parameters = 9/10 (or max 11 for AI).
%   
%   Possible eperimental protocol "Exp" are as follows:
%       simple shear
%           1:   (fs)       2:  (sn)        3:  (fn)
%           4:   (sf)       5:  (ns)        6:  (nf)
%       Biaxial (direction of readout is capitalized)
%           7:  (Fn)        8:  (fN)        9:  (Fs)
%           10: (fS)        11: (Sn)        12: (sN)
%   Biaxial Tests have to be supplemented with a stretch ratio in the form 
%   "Exp=[8 0.5]" which,  multiplied with the actual stretch, yields
%   stretch for the second index.
%   
% SurfaceStress expects "Stretch" to be a (column-) vector valued entity
% which is matched to the "Exp" input. If they are of same length, they
% are evaluated in pairs and if "Exp" has a length of 1, it is assumed to
% account for all choices of "Stretch".
% "Rotation" represents an optional rotation matrix which should contain the
% pictures of base vectors in order f-s-n.
%
%   Example: 
%       SurfaceStress("HO", [1 7 1 10 0 1 0 1], [7 0.5], ...
%       [0:0.01:0.1.]', Rod(pi/2,pi/2,3*pi/180))


if nargin <5 || isempty(Rotation)
    Rotation = eye(3);
end
sizeE=size(Exp);
sizeS=size(Stretch);
Experiment = ones(sizeS(1),2);


if sizeE(2) == 2
    Experiment(1,2)= Exp(1,2);
end

if sizeE(1) == 1
    Experiment(:,1) = Exp(1,1)*ones(sizeS(1),1);
    Experiment(:,2) = Experiment(1,2)*ones(sizeS(1),1);
else
    Experiment(:,1) = Exp(:,1);
    if sizeE(2) == 1
        Experiment(:,2) = Experiment(1,2)*ones(sizeS(1),1);
    else
        Experiment(:,2) = Exp(:,2);
    end
end

if sizeE(2) == 1
    Experiment(:,2)=ones(sizeS(1),1);
end


switch Model
    case "HO"
        if length(Parameters) ~= 8
            error("The Holzapfel-Ogden model requires at 8 parameters. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        y = arrayfun(@(gamma, protocol, ratio) sigmaHOProjected(Parameters, Rotation, [protocol ratio], gamma ), Stretch, Experiment(:,1), Experiment(:,2));
        
    case "AI"
        if length(Parameters) == 10
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma); 
            Sigmas = @(gamma,protocol,ratio) 0;
        elseif length(Parameters) ==11
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1,Parameters(7:8)], Rotation, [protocol ratio], gamma); 
            Sigmas = @(gamma,protocol,ratio) sigmaAIProjected(Parameters(end), Parameters(10), 2, Parameters(5),Parameters(6), Rotation, [protocol ratio],gamma)  ;
        else
            error("The Angular Integration model requires 10 or 11 Parameters. The first 8 are constitutive Parameters, followed by fibre and sheet(optional) dispersion coefficient and the amount of Integration steps N");
        end
        Sigmaf = @(gamma,protocol,ratio) sigmaAIProjected(Parameters(end), Parameters(9), 1, Parameters(3),Parameters(4), Rotation, [protocol ratio],gamma)  ;
        y = arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));
        
    case "0GST"
        if length(Parameters) == 9
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) 0;
        elseif length(Parameters) == 10
            rhos = Parameters(10);
            Structs{1} = MN(1,rhos);
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1, Parameters(7:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) sigmaNGSTProjected(0, Structs, 2, Parameters(5), Parameters(6), Rotation, [protocol ratio], gamma);           
        else
            error("The Generalized structure tensor models implemented need 9 or 10 parameters in total of which the last two parameters describe dispersion of fibres and sheets if needed. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        rhof = Parameters(9);
        Structf{1} = MN(1,rhof);
        Sigmaf= @(gamma, protocol, ratio) sigmaNGSTProjected(0, Structf, 1, Parameters(3), Parameters(4), Rotation, [protocol ratio], gamma);
        
        y= arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));
        
    case "1GST"
        disp("1GST is equivalent to 0GST.")
        y= SurfaceStress( "0GST", Parameters, Exp, Stretch, Rotation);
        
    case "2GST"
        if length(Parameters) == 9
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) 0;
        elseif length(Parameters) == 10
            rhos = Parameters(10);
            Structs{1} = MN(1,rhos);
            Structs{2} = MN(2,rhos);
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1, Parameters(7:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) sigmaNGSTProjected(2, Structs, 2, Parameters(5), Parameters(6), Rotation, [protocol ratio], gamma);           
        else
            error("The Generalized structure tensor models implemented need 9 or 10 parameters in total of which the last two parameters describe dispersion of fibres and sheets if needed. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        rhof = Parameters(9);
        Structf{1} = MN(1,rhof);
        Structf{2} = MN(2,rhof);
        Sigmaf= @(gamma, protocol, ratio) sigmaNGSTProjected(2, Structf, 1, Parameters(3), Parameters(4), Rotation, [protocol ratio], gamma);
        
        y= arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));
        
    case "3GST"
        if length(Parameters) == 9
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) 0;
        elseif length(Parameters) == 10
            rhos = Parameters(10);
            Structs{1} = MN(1,rhos);
            Structs{2} = MN(2,rhos);
            Structs{3} = MN(3,rhos);
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1, Parameters(7:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) sigmaNGSTProjected(3, Structs, 2, Parameters(5), Parameters(6), Rotation, [protocol ratio], gamma);           
        else
            error("The Generalized structure tensor models implemented need 9 or 10 parameters in total of which the last two parameters describe dispersion of fibres and sheets if needed. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        rhof = Parameters(9); 
        Structf{1} = MN(1,rhof);
        Structf{2} = MN(2,rhof);
        Structf{3} = MN(3,rhof);
        

        Sigmaf= @(gamma, protocol, ratio) sigmaNGSTProjected(3, Structf, 1, Parameters(3), Parameters(4), Rotation, [protocol ratio], gamma);
        y= arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));

    case "4GST"
        if length(Parameters) == 9
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) 0;
        elseif length(Parameters) == 10
            rhos = Parameters(10);
            Structs{1} = MN(1,rhos);
            Structs{2} = MN(2,rhos);
            Structs{3} = MN(3,rhos);
            Structs{4} = MN(4,rhos);
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1, Parameters(7:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) sigmaNGSTProjected(4, Structs, 2, Parameters(5), Parameters(6), Rotation, [protocol ratio], gamma);           
        else
            error("The Generalized structure tensor models implemented need 9 or 10 parameters in total of which the last two parameters describe dispersion of fibres and sheets if needed. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        rhof = Parameters(9);
        Structf{1} = MN(1,rhof);
        Structf{2} = MN(2,rhof);
        Structf{3} = MN(3,rhof);
        Structf{4} = MN(4,rhof);
        
        
        Sigmaf= @(gamma, protocol, ratio) sigmaNGSTProjected(4, Structf, 1, Parameters(3), Parameters(4), Rotation, [protocol ratio], gamma);
        y= arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));
        
    case "0SGST"
        if length(Parameters) == 9
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) 0;
        elseif length(Parameters) == 10
            rhos = Parameters(10);
            Structs{1} = MN(2,rhos);
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1, Parameters(7:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) sigmaNSGSTProjected(0, Structs, 2, Parameters(5), Parameters(6), Rotation, [protocol ratio], gamma);           
        else
            error("The Generalized structure tensor models implemented need 9 or 10 parameters in total of which the last two parameters describe dispersion of fibres and sheets if needed. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        rhof = Parameters(9);
        Structf{1} = MN(2,rhof);
 
        Sigmaf= @(gamma, protocol, ratio) sigmaNSGSTProjected(0, Structf, 1, Parameters(3), Parameters(4), Rotation, [protocol ratio], gamma);
        
        y= arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));
        
    case "1SGST"
        disp("1SGST is equivalent to 0SGST.")
        y=SurfaceStress( "0SGST", Parameters, Exp, Stretch, Rotation);
    case "2SGST"
        if length(Parameters) == 9
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) 0;
        elseif length(Parameters) == 10
            rhos = Parameters(10);
            Structs{1} = MN(2,rhos);
            Structs{2} = MN(4,rhos);
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1, Parameters(7:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) sigmaNSGSTProjected(2, Structs, 2, Parameters(5), Parameters(6), Rotation, [protocol ratio], gamma);           
        else
            error("The Generalized structure tensor models implemented need 9 or 10 parameters in total of which the last two parameters describe dispersion of fibres and sheets if needed. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        rhof = Parameters(9);
        Structf{1} = MN(2,rhof);
        Structf{2} = MN(4,rhof);
        Sigmaf= @(gamma, protocol, ratio) sigmaNSGSTProjected(2, Structf, 1, Parameters(3), Parameters(4), Rotation, [protocol ratio], gamma);
        
        y= arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));
        
    case "3SGST"
        if length(Parameters) == 9
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, Parameters(5:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) 0;
        elseif length(Parameters) == 10
            rhos = Parameters(10);
            Structs{1} = MN(2,rhos);
            Structs{2} = MN(4,rhos);
            Structs{3} = MN(6,rhos);
            SigmaIsofs = @(gamma, protocol, ratio) sigmaHOProjected([Parameters(1:2), 0, 1, 0, 1, Parameters(7:8)], Rotation, [protocol ratio], gamma);
            Sigmas= @(gamma, protocol, ratio) sigmaNSGSTProjected(3, Structs, 2, Parameters(5), Parameters(6), Rotation, [protocol ratio], gamma);           
        else
            error("The Generalized structure tensor models implemented need 9 or 10 parameters in total of which the last two parameters describe dispersion of fibres and sheets if needed. If you do not intend to make use of all terms, please set the prefactor to zero.");
        end
        rhof = Parameters(9); 
        Structf{1} = MN(2,rhof);
        Structf{2} = MN(4,rhof);
        Structf{3} = MN(6,rhof);
        

        Sigmaf= @(gamma, protocol, ratio) sigmaNSGSTProjected(3, Structf, 1, Parameters(3), Parameters(4), Rotation, [protocol ratio], gamma);
        y= arrayfun(@(gamma, protocol, ratio) SigmaIsofs(gamma, protocol, ratio)+Sigmaf(gamma, protocol, ratio)+Sigmas(gamma, protocol, ratio)   ,Stretch, Experiment(:,1), Experiment(:,2));

end

end
