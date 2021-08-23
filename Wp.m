%% Create Wp
%
function f = Wp(r,u,z)
f = 1-3*soliton(r,u,z)*soliton(r,u,z);
end