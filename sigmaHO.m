function y = sigmaHO(Parameters, FF)
% For a given deformation gradient "FF" sigmaHO returns the true Cauchy stress tensor.
% Principle direction are such that index 1 corresponds to fibres, 2 to
% sheets and ultimately 3 to the normal direction. For alternative use,
% apply corresponding base transformation to "FF". sigmaHO expects a list
% of 8 "Parameters" in order [a,b,af,bf,as,bs,afs,bfs].


    B=FF*FF';
    f=FF(:,1);
    s=FF(:,2);

y = 2*sigmaI1(FF'*FF, Parameters)*B+ 2*sigmaI4(f, Parameters(3:4))*f*f' + 2*sigmaI4(s, Parameters(5:6))*s*s'+sigmaI8(f,s, Parameters(7:8))*(f*s'+s*f');

