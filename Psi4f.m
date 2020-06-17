function y = Psi4f(f,Theta)
% Psi4f returns the fiber stress of a deformed fiber "f" given the model
% parameters "Theta(3)" and "Theta(4)" corresponding to the exponential 
% Holzapfel-Ogden model
    I4 = norm(f).^2;
y= Theta(3)*(I4-1).*exp(Theta(4)*(I4-1)^2);
