%% Analytical Jacobian of boundary conditions for Plucker problem 
%
function [dBCdya, dBCdyb] = plucker_bc_jac(ya,yb)
    dBCdya = eye(5);
    dBCdyb = zeros(5,5);
end