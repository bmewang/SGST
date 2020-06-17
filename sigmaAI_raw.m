function [s11, s22, s33, s12, s13, s23]=sigmaAI_raw(theta ,phi, rho, PrincipleDirection, a, b,  FF)
% Psi4fs_raw is the first iterate to calculate the AI (sigmaAI) fiber stress (PrincipleDirection=1)
% or sheet stress (PrincipleDirection=2). This 
% function is not meant to be used on its own. Mean fiber
% direction is aligned with first (x) axis while sheet direction points in
% second (y) axis direction. If experimental protocols involve rotated
% fiber and sheet direction it is captured via a rotated deformation
% protocol FF.
%   "rho"       = von Mises dispersion coefficient
%   "theta/phi" = classical spherical coordinates.
%   "a,b"       = Holzapfel-Ogden parameters for fiber family

    Psi4 = @(f) a*(norm(f)^2-1)*exp(b*(norm(f)^2-1)^2);


    switch PrincipleDirection % Switches between Fibre, Sheet & Normal direction
        case 1 
            v0 = @(theta,phi) [cos(theta),sin(theta).*cos(phi),sin(theta).*sin(phi)]';
        case 2
            v0 = @(theta,phi) [cos(phi).*sin(theta),cos(theta),sin(theta).*sin(phi)]';
        case 3
            v0 = @(theta,phi) [cos(phi).*sin(theta),sin(theta).*sin(phi),cos(theta)]';
    end
    v1 = @(theta, phi) FF(1,:)*v0(theta,phi);
    v2 = @(theta, phi) FF(2,:)*v0(theta,phi);
    v3 = @(theta, phi) FF(3,:)*v0(theta,phi);
    
    % Include (norm(v(theta,phi))>1) to discard Myofibers under Compression
    Psi4_raw11 = @(theta,phi) arrayfun(@(theta,phi) 2.*Psi4([v1(theta,phi); v2(theta,phi); v3(theta,phi)]) .*v1(theta,phi)*v1(theta,phi)'.*sin(theta)*exp(rho*(cos(2*theta)+1)) , theta,phi);
    Psi4_raw22 = @(theta,phi) arrayfun(@(theta,phi) 2.*Psi4([v1(theta,phi); v2(theta,phi); v3(theta,phi)]) .*v2(theta,phi)*v2(theta,phi)'.*sin(theta)*exp(rho*(cos(2*theta)+1)) , theta,phi);
    Psi4_raw33 = @(theta,phi) arrayfun(@(theta,phi) 2.*Psi4([v1(theta,phi); v2(theta,phi); v3(theta,phi)]) .*v3(theta,phi)*v3(theta,phi)'.*sin(theta)*exp(rho*(cos(2*theta)+1)) , theta,phi);
    Psi4_raw12 = @(theta,phi) arrayfun(@(theta,phi) 2.*Psi4([v1(theta,phi); v2(theta,phi); v3(theta,phi)]) .*v1(theta,phi)*v2(theta,phi)'.*sin(theta)*exp(rho*(cos(2*theta)+1)) , theta,phi);
    Psi4_raw13 = @(theta,phi) arrayfun(@(theta,phi) 2.*Psi4([v1(theta,phi); v2(theta,phi); v3(theta,phi)]) .*v1(theta,phi)*v3(theta,phi)'.*sin(theta)*exp(rho*(cos(2*theta)+1)) , theta,phi);
    Psi4_raw23 = @(theta,phi) arrayfun(@(theta,phi) 2.*Psi4([v1(theta,phi); v2(theta,phi); v3(theta,phi)]) .*v2(theta,phi)*v3(theta,phi)'.*sin(theta)*exp(rho*(cos(2*theta)+1)) , theta,phi);

    [X,Y]=meshgrid(theta,phi);
    s11 = Psi4_raw11(X,Y);
    s22 = Psi4_raw22(X,Y);
    s33 = Psi4_raw33(X,Y);
    s12 = Psi4_raw12(X,Y);
    s13 = Psi4_raw13(X,Y);
    s23 = Psi4_raw23(X,Y);
end
