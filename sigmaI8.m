function y=sigmaI8(f,s, Theta)
% Fourth term of the HO-model corresponding to pseudo-invariant I8.    
% Expects model Parameters "Theta(1)" and "Theta(2)" ~ [afs,bfs]
    I8 = s'*f;
y = Theta(1)*I8*exp(Theta(2)*I8^2);
