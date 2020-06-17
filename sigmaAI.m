function y= sigmaAI(N, rho, PrincipleDirection, a,b, FF)
% For a given deformation gradient "FF" sigmaAI returns the true Cauchy 
% stress where "PrincipleDirection" signifies which index corresponds to 
% the mean fibre direction. For fiber arrangements that do not comply with
% the principle directions, apply corresponding base transformation to "FF".
% Material parameters "a" and "b" resemble fibre stiffness as within the HO-model.
% As MATLAB's inbuilt integration schemes struggle with vector valued output
% for double integrations, a simplistic trapezoidal rule is employed,
% where paramter "N" defines the amount of trapezoids per dimension (2). 
% This scheme converges slowly and is to be used with care.
% Currently, exclusion of compressed fibres is not accounted for but can be
% achieved easily via inclusion of logic expressions in function
% "sigmaAI_raw".


dtheta = pi /N;
theta = 0:dtheta:pi;
dphi = dtheta;
phi = 0:dphi:(2*pi);
dtheta = ones(N+1,1) * dtheta; dtheta(1) = dtheta(1)/2; dtheta(end) = dtheta(end)/2;
dphi = ones(2*N+1,1) * dphi; dphi(1) = dphi(1)/2; dphi(end) = dphi(end)/2;

[s11, s22, s33, s12, s13, s23]=sigmaAI_raw(theta ,phi, rho, PrincipleDirection, a, b,  FF);

Psi_avg11 = dphi'* s11 *dtheta ;
Psi_avg22 = dphi'* s22 *dtheta ;
Psi_avg33 = dphi'* s33 *dtheta ;
Psi_avg12 = dphi'* s12 *dtheta ;
Psi_avg13 = dphi'* s13 *dtheta ;
Psi_avg23 = dphi'* s23 *dtheta ;


y = [Psi_avg11, Psi_avg12, Psi_avg13;...
    Psi_avg12, Psi_avg22, Psi_avg23;...
    Psi_avg13, Psi_avg23, Psi_avg33]*sqrt(rho/2/pi)/erfi(sqrt(2*rho))/pi;
