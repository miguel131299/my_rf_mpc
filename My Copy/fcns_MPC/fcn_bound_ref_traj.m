function [p,Xt,Ut] = fcn_bound_ref_traj(p)
% This function finds the initial condition for periodic bounding
% The calculation is based on the paper (citation):
% Park, Hae-Won, Patrick M. Wensing, and Sangbae Kim. 
% "High-speed bounding with the MIT Cheetah 2: Control design and experiments."
% The International Journal of Robotics Research 36, no. 2 (2017): 167-192.

%% Variable Initialization
[mass,J,g,Tst,Tsw] = deal(p.mass,p.J,p.g,p.Tst,p.Tsw);
T = Tst + Tsw;
Tair = 1/2 * (Tsw - Tst);

b_co = [0 0.8 1 1 0.8 0];
b_ = mean(b_co);

%% Fz - vertical ground reaction force
% 2 * alpha * b_ * Tst = mass * g * T
alpha_z = (mass * g * T) / (2 * b_ * Tst);

Fz_co = alpha_z * b_co; % coefficients of Fz


dz_co = bz_int(Fz_co/mass-g,0,Tst);     % vertial velocity
z_co = bz_int(dz_co,0,Tst);             % vertical position

% first principle: integration
dz0 = -1/(Tst+Tair)*(z_co(end) + Tair*(dz_co(end)+g*Tst)-1/2*g*((Tst+Tair)^2-Tst^2));

dz_co = bz_int(Fz_co/mass-g,dz0,Tst);
z_co = bz_int(dz_co,p.z0,Tst);

%% theta - pitch angle during bounding motion
alpha_th = 140 * J(2,2);
tau_co = -alpha_th * b_co;              % coefficients of torque

dth_co = bz_int(tau_co/J(2,2),0,Tst);
% by symmetry
dth0 = -1/2 * dth_co(end);

th0 = dth0*Tair/2;
dth_co = bz_int(tau_co/J(2,2),dth0,Tst);    % angular velocity
th_co = bz_int(dth_co,th0,Tst);             % pitch angle

%% output B-spline coefficient
p.Fz_co = Fz_co;            % vertical ground reaction force
p.dz_co = dz_co;            % vertical velocity
p.z_co = z_co;              % vertical position

p.tau_co = tau_co;          % coefficients of torque
p.dth_co = dth_co;          % angular velocity
p.th_co = th_co;            % pitch angle theta

%% intial condition
% expm - matrix exponential (using Taylor Expansion)
R0 = expm(hatMap([0;th0;0]));               % Rotation Matrix
Xt = [0;0;p.z0;0;0;dz0;R0(:);0;dth0;0];     % inital state vector
Xt(19:30) = p.pf34(:);
Ut = repmat([0;0;1/4*p.mass*p.g],[4,1]);    % control input matrix

end