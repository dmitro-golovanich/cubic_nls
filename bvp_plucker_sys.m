%% Dynamical system formulation of Plucker problem (without singular term)
% for the bvp solver bvp4c
function dydx = bvp_plucker_sys(c,ep,lambda,r,u,x,y)
dydx = zeros(5,1);
dydx(1) = y(2) - y(3); 
dydx(2) = Wm(r,u,x)*y(1) + Vep(c,ep,x)*y(1) + y(5);
dydx(3) = -Wp(r,u,x)*y(1) - Vep(c,ep,x)*y(1) - y(5); 
dydx(4) = lambda*y(1);
dydx(5) = Wp(r,u,x)*y(2) + Vep(c,ep,x)*y(2) - Wm(r,u,x)*y(3) - Vep(c,ep,x)*y(3) - 2*lambda*y(4);
end