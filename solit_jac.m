%% Analytical Jacobian of the soliton problem
%
function J = solit_jac(c,ep,x,y)
J = zeros(2,2);
J(1,2) = 1;
J(2,1) = (1-ep^2*Vep(c,ep,0))-3*(y(1)^2) + Vep(c,ep,x);
end