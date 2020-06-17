function y= sigmaI1(C, Theta)
% First term of HO-model corresponding to invariant I1.
% Expects model parameters "Theta(1)" and "Theta(2)" ~ [a,b] from
% Holzapfel-Ogden model
y = Theta(1) *exp(Theta(2)*(trace(C)-3))/2;
