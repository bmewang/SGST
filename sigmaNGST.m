function StressTensor= sigmaNGST(N, StructureTensorList, PrincipleDirection, a, b, FF)
% For a given deformation gradient "FF" sigmaNGST returns the true Cauchy 
% stress tensor where "PrincipleDirection" signifies which index corresponds to 
% the mean fibre direction. Thus, you have to make sure that "FF" is defined
% in the same base system. "N" denotes to which order the nGST-model is
% evaluated, demanding an equal amount of structure tensors which have to be
% calculated beforehand using function "MN". Models available range from 
% zeroth to fourth order. Material parameters "a" and "b" resemble fibre 
% stiffness as within the HO-model.

if N>4
    error("sigmaNGST only supports expansion up to order 4")
elseif N> length(StructureTensorList)
    error("The order of the Taylor-Series must not exceed the amount of Structure Tensors handed to sigmaNGST")
end

if floor(PrincipleDirection)~=PrincipleDirection || PrincipleDirection>3
    error("PrincipleDirection mus take values 1,2, or 3 to determine which is the mean fibre direction")
end

switch PrincipleDirection
    case 1
        main = 1; side1 = 2; side2 = 3;
    case 2
        main = 2; side1 = 1; side2 = 3;
    case 3
        main = 3; side1 = 2; side2 = 1;
end
C = FF'*FF;


kappa = StructureTensorList{1};

M = kappa*eye(3);
M(main, main) = 1-2*kappa;

B = FF*M*FF';

I4 = trace(B);
xi = I4 -1;
Prefactor = 2*a*exp(b*xi^2);

Accumulant = xi*M;

if N>=2
    
    M2 = StructureTensorList{2};
    C2MM = C(main,main)^2*M2(1)+ M2(3)*( C(side1,side1)^2+C(side2,side2)^2)+...
        2*( M2(2)*(C(main,main)*( C(side1,side1)+C(side2,side2) )+2*C(main,side2)^2+2*C(main,side1)^2 )+ M2(4)*(C(side1,side1)*C(side2,side2)+2*C(side1,side2)^2));    

    MM2 = C2MM - I4^2;

    MMD=zeros(3,3);
    MMD(main,side1) = 2*M2(2)*C(main,side1);
    MMD(main,side2) = 2*M2(2)*C(main,side2);
    MMD(side1,side2) = 2*M2(4)*C(side1,side2);
    MMD = MMD+MMD';
    MMD(main,main) = M2(1)*C(main,main) + M2(2)*(C(side1,side1) + C(side2,side2));
    MMD(side1,side1) = M2(2)*C(main,main) + M2(3)*C(side1,side1) + M2(4)*C(side2,side2);
    MMD(side2,side2) = M2(2)*C(main,main) + M2(4)*C(side1,side1) + M2(3)*C(side2,side2);

    MMD_next = 2*(MMD-M*I4);

    Accumulant = Accumulant + 1/2* (M*MM2* (4*b^2*xi^3+ 6*b*xi) +MMD_next*(2*b*xi^2+1));
    
    if N>=3
            M3 = StructureTensorList{3};

            C3M3= C(main,main).^3.*M3(1) ...
                +(12.*C(main,main).*C(main,side1).^2+12.*C(main,main).*C(main,side2).^2+3.*C(main,main).^2.*C(side1,side1)+3.*C(main,main).^2.*C(side2,side2)).*M3(2) ...
                +(12.*C(main,side1).^2.*C(side1,side1)+3.*C(main,main).*C(side1,side1).^2+12.*C(main,side2).^2.*C(side2,side2)+3.*C(main,main).*C(side2,side2).^2).*M3(3) ...
                +(12.*C(main,side2).^2.*C(side1,side1)+48.*C(main,side1).*C(main,side2).*C(side1,side2)+12.*C(main,main).*C(side1,side2).^2+12.*C(main, side1).^2.*C(side2,side2)+6.*C(main,main).*C(side1,side1).*C(side2,side2)).*M3(4) ...
                +(C(side1,side1).^3+C(side2,side2).^3).*M3(5) ...
                +(12.*C(side1,side1).*C(side1,side2).^2+3.*C(side1,side1).^2.*C(side2,side2)+12.*C(side1,side2).^2.*C(side2,side2)+3.*C(side1,side1).*C(side2,side2).^2).*M3(6);

            M3D = zeros(3);
            M3D(main,side1) = 12.*C(main,main).*C(main,side1).*M3(2) ...
            +12.*C(main,side1).*C(side1,side1).*M3(3) ...
            +(1/2).*(48.*C(main,side2).*C(side1,side2)+24.*C(main,side1).*C(side2,side2)).*M3(4);
            M3D(main,side2) = 12.*C(main,main).*C(main,side2).*M3(2) ...
            +12.*C(main,side2).*C(side2,side2).*M3(3) ...
            +(1/2).*(24.*C(main,side2).*C(side1,side1)+48.*C(main,side1).*C(side1,side2)).*M3(4);
            M3D(side1,side2) =  (1/2).*(48.*C(main,side1).*C(main,side2)+24.*C(main,main).*C(side1,side2)).*M3(4)...
            +(1/2).*(24.*C(side1,side1).*C(side1,side2)+24.*C(side1,side2).*C(side2,side2)).*M3(6);
            M3D = M3D+M3D';
            M3D(main,main) = 3.*C(main,main).^2.*M3(1) ...
            +(12.*C(main,side1).^2+12.*C(main,side2).^2+6.*C(main,main).*C(side1,side1)+6.*C(main,main).*C(side2,side2)).*M3(2) ...
            +(3.*C(side1,side1).^2+3.*C(side2,side2).^2).*M3(3) ...
            +(12.*C(side1,side2).^2+6.*C(side1,side1).*C(side2,side2)).*M3(4);
            M3D(side1,side1) = 3.*C(main,main).^2.*M3(2) ...
            +(12.*C(main,side1).^2+6.*C(main,main).*C(side1,side1)).*M3(3) ...
            +(12.*C(main,side2).^2+6.*C(main,main).*C(side2,side2)).*M3(4) ...
            +3.*C(side1,side1).^2.*M3(5) ...
            +(12.*C(side1,side2).^2+6.*C(side1,side1).*C(side2,side2)+3.*C(side2,side2).^2).*M3(6);
            M3D(side2,side2) = 3.*C(main,main).^2.*M3(2) ...
            +(12.*C(main,side2).^2+6.*C(main,main).*C(side2,side2)).*M3(3) ...
            +(12.*C(main,side1).^2+6.*C(main,main).*C(side1,side1)).*M3(4) ...
            +3.*C(side2,side2).^2.*M3(5) ...
            +(3.*C(side1,side1).^2+12.*C(side1,side2).^2+6.*C(side1,side1).*C(side2,side2)).*M3(6);
            M3D = M3D/3;


            M3D_next = 3*(M3D + 2* (M*I4^2 - MMD*I4) -M*C2MM); 

            Accumulant = Accumulant +b* ( (1+4*b*xi^2+4/3*b^2*xi^4)* (C3M3-3*C2MM*I4 +2*I4^3) *M +xi*(1+2/3*b*xi^2)*M3D_next);
        if N>=4
            M4 = StructureTensorList{4};

            C4M4 = C(main,main).^4.*M4(1)+(24.*C(main,main).^2.*C(main,side1).^2+24.*C(main,main).^2.*C(main,side2).^2+4.*C(main,main).^3.*C(side1,side1)+4.*C(main,main).^3.*C(side2,side2)).*M4(2)+(16.*C(main,side1).^4+16.*C(main,side2).^4+48.*C(main,main).*C(main,side1).^2.*C(side1,side1)+6.*C(main,main).^2.*C(side1,side1).^2+48.*C(main,main).*C(main,side2).^2.*C(side2,side2)+6.*C(main,main).^2.*C(side2,side2).^2).*M4(3)+(96.*C(main,side1).^2.*C(main,side2).^2+48.*C(main,main).*C(main,side2).^2.*C(side1,side1)+192.*C(main,main).*C(main,side1).*C(main,side2).*C(side1,side2)+24.*C(main,main).^2.*C(side1,side2).^2+48.*C(main,main).*C(main,side1).^2.*C(side2,side2)+12.*C(main,main).^2.*C(side1,side1).*C(side2,side2)).*M4(4)+(24.*C(main,side1).^2.*C(side1,side1).^2+4.*C(main,main).*C(side1,side1).^3+24.*C(main,side2).^2.*C(side2,side2).^2+4.*C(main,main).*C(side2,side2).^3).*M4(5)+(24.*C(main,side2).^2.*C(side1,side1).^2+192.*C(main,side1).*C(main,side2).*C(side1,side1).*C(side1,side2)+96.*C(main,side1).^2.*C(side1,side2).^2+96.*C(main,side2).^2.*C(side1,side2).^2+48.*C(main,main).*C(side1,side1).*C(side1,side2).^2+48.*C(main,side1).^2.*C(side1,side1).*C(side2,side2)+48.*C(main,side2).^2.*C(side1,side1).*C(side2,side2)+12.*C(main,main).*C(side1,side1).^2.*C(side2,side2)+192.*C(main,side1).*C(main,side2).*C(side1,side2).*C(side2,side2)+48.*C(main,main).*C(side1,side2).^2.*C(side2,side2)+24.*C(main,side1).^2.*C(side2,side2).^2+12.*C(main,main).*C(side1,side1).*C(side2,side2).^2).*M4(6)+(C(side1,side1).^4+C(side2,side2).^4).*M4(7)+(24.*C(side1,side1).^2.*C(side1,side2).^2+4.*C(side1,side1).^3.*C(side2,side2)+24.*C(side1,side2).^2.*C(side2,side2).^2+4.*C(side1,side1).*C(side2,side2).^3).*M4(8)+(16.*C(side1,side2).^4+48.*C(side1,side1).*C(side1,side2).^2.*C(side2,side2)+6.*C(side1,side1).^2.*C(side2,side2).^2).*M4(9);

            
            M4D = zeros(3,3);
            M4D(main,side1)=24.*C(main,main).^2.*C(main,side1).*M4(2)+(1/2).*(64.*C(main,side1).^3+96.*C(main,main).*C(main,side1).*C(side1,side1)).*M4(3)+(1/2).*(192.*C(main,side1).*C(main,side2).^2+192.*C(main,main).*C(main,side2).*C(side1,side2)+96.*C(main,main).*C(main,side1).*C(side2,side2)).*M4(4)+24.*C(main,side1).*C(side1,side1).^2.*M4(5)+(1/2).*(192.*C(main,side2).*C(side1,side1).*C(side1,side2)+192.*C(main,side1).*C(side1,side2).^2+96.*C(main,side1).*C(side1,side1).*C(side2,side2)+192.*C(main,side2).*C(side1,side2).*C(side2,side2)+48.*C(main,side1).*C(side2,side2).^2).*M4(6);
            M4D(main,side2)=24.*C(main,main).^2.*C(main,side2).*M4(2)+(1/2).*(64.*C(main,side2).^3+96.*C(main,main).*C(main,side2).*C(side2,side2)).*M4(3)+(1/2).*(192.*C(main,side1).^2.*C(main,side2)+96.*C(main,main).*C(main,side2).*C(side1,side1)+192.*C(main,main).*C(main,side1).*C(side1,side2)).*M4(4)+24.*C(main,side2).*C(side2,side2).^2.*M4(5)+(1/2).*(48.*C(main,side2).*C(side1,side1).^2+192.*C(main,side1).*C(side1,side1).*C(side1,side2)+192.*C(main,side2).*C(side1,side2).^2+96.*C(main,side2).*C(side1,side1).*C(side2,side2)+192.*C(main,side1).*C(side1,side2).*C(side2,side2)).*M4(6);
            M4D(side1,side2)=(1/2).*(192.*C(main,main).*C(main,side1).*C(main,side2)+48.*C(main,main).^2.*C(side1,side2)).*M4(4)+(1/2).*(192.*C(main,side1).*C(main,side2).*C(side1,side1)+192.*C(main,side1).^2.*C(side1,side2)+192.*C(main,side2).^2.*C(side1,side2)+96.*C(main,main).*C(side1,side1).*C(side1,side2)+192.*C(main,side1).*C(main,side2).*C(side2,side2)+96.*C(main,main).*C(side1,side2).*C(side2,side2)).*M4(6)+(1/2).*(48.*C(side1,side1).^2.*C(side1,side2)+48.*C(side1,side2).*C(side2,side2).^2).*M4(8)+(1/2).*(64.*C(side1,side2).^3+96.*C(side1,side1).*C(side1,side2).*C(side2,side2)).*M4(9);
            M4D = M4D+M4D';
            M4D(main,main)=4.*C(main,main).^3.*M4(1)+(48.*C(main,main).*C(main,side1).^2+48.*C(main,main).*C(main,side2).^2+12.*C(main,main).^2.*C(side1,side1)+12.*C(main,main).^2.*C(side2,side2)).*M4(2)+(48.*C(main,side1).^2.*C(side1,side1)+12.*C(main,main).*C(side1,side1).^2+48.*C(main,side2).^2.*C(side2,side2)+12.*C(main,main).*C(side2,side2).^2).*M4(3)+(48.*C(main,side2).^2.*C(side1,side1)+192.*C(main,side1).*C(main,side2).*C(side1,side2)+48.*C(main,main).*C(side1,side2).^2+48.*C(main,side1).^2.*C(side2,side2)+24.*C(main,main).*C(side1,side1).*C(side2,side2)).*M4(4)+(4.*C(side1,side1).^3+4.*C(side2,side2).^3).*M4(5)+(48.*C(side1,side1).*C(side1,side2).^2+12.*C(side1,side1).^2.*C(side2,side2)+48.*C(side1,side2).^2.*C(side2,side2)+12.*C(side1,side1).*C(side2,side2).^2).*M4(6);
            M4D(side1,side1)=4.*C(main,main).^3.*M4(2)+(48.*C(main,main).*C(main,side1).^2+12.*C(main,main).^2.*C(side1,side1)).*M4(3)+(48.*C(main,main).*C(main,side2).^2+12.*C(main,main).^2.*C(side2,side2)).*M4(4)+(48.*C(main,side1).^2.*C(side1,side1)+12.*C(main,main).*C(side1,side1).^2).*M4(5)+(48.*C(main,side2).^2.*C(side1,side1)+192.*C(main,side1).*C(main,side2).*C(side1,side2)+48.*C(main,main).*C(side1,side2).^2+48.*C(main,side1).^2.*C(side2,side2)+48.*C(main,side2).^2.*C(side2,side2)+24.*C(main,main).*C(side1,side1).*C(side2,side2)+12.*C(main,main).*C(side2,side2).^2).*M4(6)+4.*C(side1,side1).^3.*M4(7)+(48.*C(side1,side1).*C(side1,side2).^2+12.*C(side1,side1).^2.*C(side2,side2)+4.*C(side2,side2).^3).*M4(8)+(48.*C(side1,side2).^2.*C(side2,side2)+12.*C(side1,side1).*C(side2,side2).^2).*M4(9);
            M4D(side2,side2)=4.*C(main,main).^3.*M4(2)+(48.*C(main,main).*C(main,side2).^2+12.*C(main,main).^2.*C(side2,side2)).*M4(3)+(48.*C(main,main).*C(main,side1).^2+12.*C(main,main).^2.*C(side1,side1)).*M4(4)+(48.*C(main,side2).^2.*C(side2,side2)+12.*C(main,main).*C(side2,side2).^2).*M4(5)+(48.*C(main,side1).^2.*C(side1,side1)+48.*C(main,side2).^2.*C(side1,side1)+12.*C(main,main).*C(side1,side1).^2+192.*C(main,side1).*C(main,side2).*C(side1,side2)+48.*C(main,main).*C(side1,side2).^2+48.*C(main,side1).^2.*C(side2,side2)+24.*C(main,main).*C(side1,side1).*C(side2,side2)).*M4(6)+4.*C(side2,side2).^3.*M4(7)+(4.*C(side1,side1).^3+48.*C(side1,side2).^2.*C(side2,side2)+12.*C(side1,side1).*C(side2,side2).^2).*M4(8)+(48.*C(side1,side1).*C(side1,side2).^2+12.*C(side1,side1).^2.*C(side2,side2)).*M4(9);
            M4D=M4D/4;
            

            M4D_next = 4*(M4D - 3*I4*M3D +3*I4^2*MMD) +4*(3*I4*C2MM - C3M3 -3*I4^3)*M;
            
            Accumulant = Accumulant +b^2*xi/2*(5+20/3*b*xi^2+4/3*b^2*xi^4)*(C4M4-4*I4*C3M3+6*I4^2*C2MM-3*I4^4)  *M + b/2*(0.5+2*b*xi^2+2/3*b^2*xi^4)*M4D_next;
        end
    end
end
StressTensor = Prefactor * FF*Accumulant*FF';
