%% Clear previous values and clear screen
%
clear;
clc;

%% Set the Potential magnitude "c" and epsilon scaling "ep" for potential Vep.m
%
c = -10;
ep = 0.1;

%% Set coupling constant lambda for Plucker system
%
lambda = 1.5;

%% BEGIN COMPUTE SOLITON
%
%
%
%% Set soliton amplitude guess "amp0", interval [x0,xmax0], and initial # of points in mesh "npts0"
%
amp0 = 4.3373876;

x0 = 0;
xmax0 = 15;
npts0 = 20;

%% Set further parameters needed for extending soliton to [xmax0,xmax1]
%
deltaxmax = 5;
xmax1  = 850;
nmax = floor((xmax1-xmax0)/deltaxmax);
xmax = xmax0;

%% Create initial guess
%
guessinit = @(x)solit_guess(amp0,x);
xmesh = linspace(x0,xmax0,npts0);
sol = bvpinit(xmesh, guessinit);

%% Set singular term "S" and call bvp4c to calculate soliton on [x0,xmax0]
% 
S = [0 0; 0 -2];

solit_sys_handle = @(x,y)solit_sys(c,ep,x,y);
solit_jac_handle = @(x,y)solit_jac(c,ep,x,y);

solit_bc_handle = @(ya,yb)solit_bc(xmax0,ya,yb);
solit_bc_jac_handle = @(ya,yb)solit_bc_jac(xmax0,ya,yb);

options = bvpset('SingularTerm',S,'FJacobian',solit_jac_handle,'BCJacobian',solit_bc_jac_handle);    
sol = bvp4c(solit_sys_handle, solit_bc_handle, sol, options);

%% Iterate & extend soliton to [x0,xmax1]
%
for i=1:nmax
    xmax = xmax + deltaxmax;
    sol = bvpinit(sol,[x0,xmax]);
    
    solit_sys_handle = @(x,y)solit_sys(c,ep,x,y);
    solit_jac_handle = @(x,y)solit_jac(c,ep,x,y);

    solit_bc_handle = @(ya,yb)solit_bc(xmax,ya,yb);
    solit_bc_jac_handle = @(ya,yb)solit_bc_jac(xmax,ya,yb);

    options = bvpset('SingularTerm',S,'FJacobian',solit_jac_handle,'BCJacobian',solit_bc_jac_handle);    
    sol = bvp4c(solit_sys_handle, solit_bc_handle, sol, options);
end

%% Plot soliton on [0,10]
%
plot(sol.x, sol.y(1,:));
axis([0 10 0 inf]);

%% Vectorize the soliton so that it may be made into a function.
%
r = sol.x';
u = sol.y(1,:)';

%% Clear all un-necessary values
%
clear amp0 solit_bc_jac_handle deltaxmax guessinit i;
clear solit_jac_handle nmax npts0 options S solit_bc_handle;
clear solit_sys_handle x0 xmax0 xmax1 xmax xmesh sol;

%% END COMPUTE SOLITON
%
%
%
%% BEGIN COMPUTE PLUCKER SYSTEM
%
%
%
%% Set initial interval of solution [x0,xmax0], set initial number of mesh points npts0
%
x0=0;
xmax0 = 15;
npts0 = 20;

%% Create initial guess
%
guessinit = @(x)plucker_guess(x);
xmesh = linspace(x0,xmax0,npts0);
sol = bvpinit(xmesh, guessinit);

%% Set singular term S
%
S = zeros(5,5);
S(2,2) = -2;
S(3,3) = -2;
S(4,4) = -2;
S(5,5) = -4;

%% Call bvp4c to calculate plucker system on [x0,xmax0]
%
bvp_plucker_sys_handle = @(x,y)bvp_plucker_sys(c,ep,lambda,r,u,x,y);
bvp_plucker_jac_handle = @(x,y)bvp_plucker_jac(c,ep,lambda,r,u,x,y);

plucker_bc_handle = @(ya,yb)plucker_bc(ya,yb);
plucker_bc_jac_handle = @(ya,yb)plucker_bc_jac(ya,yb);

options = bvpset('SingularTerm',S,'FJacobian',bvp_plucker_jac_handle,'BCJacobian',plucker_bc_jac_handle);    
sol = bvp4c(bvp_plucker_sys_handle, plucker_bc_handle, sol, options);

%% Get the components of the 5D Plucker system on [x0, xmax0]
%
x = sol.x;
s1 = sol.y(1,:);
s2 = sol.y(2,:);
s3 = sol.y(3,:);
s4 = sol.y(4,:);
s5 = sol.y(5,:);

%% Calculate the angle tht of det^2 on [x0, xmax0]
%
re = s1-s5;
im = s2-s3;
top = complex(re,im);
bottom = complex(re,-im);
det2 = top./bottom;
tht = angle(det2); 

%% Plot det2 angle on [x0, xmax0]
%
subplot(1,1,1);
plot(x, tht, x, x*0, 'LineWidth',2);
xlim([0,15]);
xlabel('x');
legend('Angle');

%% Clear all un-necessary variables
%
clear bottom det2 guessinit IM im1 npts0 options top;
clear plucker_bc_handle plucker_bc_jac_handle bvp_plucker_jac_handle;
clear bvp_plucker_sys_handle re RE S sol x xmesh tht;

%% Set initial condition Y0 at xmax0
%
Y0 = zeros(5,1);
Y0(1,1) = s1(end);
Y0(2,1) = s2(end);
Y0(3,1) = s3(end);
Y0(4,1) = s4(end);
Y0(5,1) = s5(end);

%% Set right endpoint xmax1, and set xspan = [xmax0, xmax1]
%
xmax1 = 500;
xspan = [xmax0 xmax1];

%% Clear s1-s5
%
clear s1 s2 s3 s4 s5;

%% Solve the system on [xmax0, xmax1] using ode45
%
ode_plucker_sys_handle = @(x,y)ode_plucker_sys(c,ep,lambda,r,u,x,y);
ode_plucker_jac_handle = @(x,y)ode_plucker_jac(c,ep,lambda,r,u,x,y);

options = odeset('Jacobian',ode_plucker_jac_handle);

[x,y] = ode45(ode_plucker_sys_handle,xspan,Y0);

%% Extract Plucker system on [xmax0, xmax1]
%
y = y';
s1 = y(1,:);
s2 = y(2,:);
s3 = y(3,:);
s4 = y(4,:);
s5 = y(5,:);

%% Calculate the angle tht of det^2 on [xmax0, xmax1]
%
re = s1-s5;
im = s2-s3;
top = complex(re,im);
bottom = complex(re,-im);
det2 = top./bottom;
tht = angle(det2); 

%% Plot det2 angle on [xmax0, xmax1]
%
subplot(1,1,1);
plot(x, tht, x, x*0, 'LineWidth',2);
xlim([xmax0,xmax1]);
xlabel('x');
legend('Angle');

%% END COMPUTE PLUCKER SYSTEM
%
%
%
%% BEGIN COMPUTE NUMBER OF L+/- EIGENSTATES
%












