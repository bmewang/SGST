function y =  Rod(theta, psi, omega)
% Rod employs the Rodriguez formula to rotate about vector "v" by an angle
% "omega". v on the other hand is defined in polar coordinates "theta" and
% "psi". "theta" is measured with respect to the first axis and "psi" with
% repsect to the second.

v=[ cos(theta), sin(theta)*cos(psi), sin(theta)*sin(psi) ];


y=[[cos(omega)+v(1)^2*(1-cos(omega)) ,          v(1)*v(2)*(1-cos(omega))-v(3)*sin(omega) ,  v(2)*sin(omega)+v(1)*v(3)*(1-cos(omega))];
   [v(3)*sin(omega)+v(1)*v(2)*(1-cos(omega)) ,  cos(omega)+v(2)^2*(1-cos(omega)) ,          -v(1)*sin(omega)+v(2)*v(3)*(1-cos(omega))];
   [-v(2)*sin(omega)+v(1)*v(3)*(1-cos(omega)) , v(1)*sin(omega)+v(2)*v(3)*(1-cos(omega)) ,  cos(omega)+v(3)^2*(1-cos(omega))]];