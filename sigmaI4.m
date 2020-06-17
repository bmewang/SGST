function y = sigmaI4(f,Theta)
% Second and third term of the HO-model corresponding to pseudo-invariant
% I4. Expects model Parameters "Theta(1)" and "Theta(2)" ~ [af,bf] or [as,bs]
    I4 = norm(f).^2;
y= Theta(1)*(I4-1).*exp(Theta(2)*(I4-1)^2);
