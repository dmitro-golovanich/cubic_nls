%% Analytic Jacobian of the Dynamical system bvp_plucker_sys.m
%
function J = bvp_plucker_jac(c,ep,lambda,r,u,x,y)
J = zeros(5,5);
J(1,2) = 1;
J(1,3) = -1;
J(2,1) = Wm(r,u,x) + Vep(c,ep,x);
J(2,5) = 1;
J(3,1) =  -Wp(r,u,x)- Vep(c,ep,x);
J(3,5) = -1;
J(4,1) = lambda;
J(5,2) = Wp(r,u,x)+ Vep(c,ep,x);
J(5,3) = -Wm(r,u,x) - Vep(c,ep,x);
J(5,4) = -2*lambda;
end