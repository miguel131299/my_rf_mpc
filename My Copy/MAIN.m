%% initialization
clear all; close all; clc
addpath fcns fcns_MPC fcns_animation\


%% --- parameters ---
% ---- gait ----
% 0-trot; 1-bound; 2-pacing; 3-gallop; 4-trot run; 5-crawl

gait = 2;
p = get_params(gait);
p.playSpeed = 1;
p.flag_movie = 1;       % 1 - make movie
use_qpSWIFT = 0;        % 0 - quadprog, 1 - qpSWIFT (external)

dt_sim = p.simTimeStep;
SimTimeDuration = 1.5; % [sec]
MAX_ITER = floor(SimTimeDuration/p.simTimeStep);

% desired trajectory
p.acc_d = 1;
p.vel_d = [0.5;0];
p.yaw_d = 0.5;

%% Model Predictive Control
% --- initial condition ---
% Xt = [pc dpc vR wb pf]': [30,1]

if gait == 1
    [p, Xt, Ut] = fcn_bound_ref_traj(p);
else
    [Xt, Ut] = fcn_gen_XdUd(0, [], [1;1;1;1], p);
end

% --- logging ---
tstart = 0;
tend = dt_sim;

[tout, Xout, Uout, Xdout, Udout, Uext, FSMout] = deal([]);

% --- simulation ---
h_waitbar = waitbar(0, 'Calculating...');

tic % start timer

for ii = 1:MAX_ITER

    % --- time vector ---
    % calculated based on current iteration and prediction horizon
    t_ = dt_sim * (ii - 1) + p.Tmpc * (0:p.predHorizon - 1);

    % --- FSM ---
    % determine current state and desired state and control input
    if gait == 1
        [FSM, Xd, Ud, Xt] = fcn_FSM_bound(t_, Xt, p);
    else 
        [FSM, Xd, Ud, Xt] = fcn_FSM(t_, Xt, p);
    end

    % --- MPC ---
    % formulate as Quadratic Problem (QP)
    [H, g, Aineq, bineq, Aeq, beq] = fcn_get_QP_form_eta(Xt, Ut, Xd, Ud, p);

    if ~use_qpSWIFT
    % solve QP using quadprog
    [zval] = quadprog(H,g,Aineq,bineq,Aeq,beq,[],[]);
    else
        % interface with the QP solver qpSWIFT
    end
    
    % Update control input based on QP solution (only first step)
    Ut = Ut + zval(1:12);

    % --- external disturbance ---
    [u_ext, p_ext] = fcn_get_disturbance(tstart, p);
    p.p_ext = p_ext;        % position of external force
    u_ext = 0 * u_ext;      % why 0?

    % --- simulate ---
    % system dynamics solved using the ode45 solver and the updated control
    % input and disturbances

    % ode45(system dynamics, time interval, initial condition)

    % @(t, X) creates a function handle that represents the dynamics_SRB
    % function with fixed inouts, while allowing t and X to vary during the
    % integration process
    [t, X] = ode45(@(t, X)dynamics_SRB(t, X, Ut, Xd, 0*u_ext, p), [tstart, tend], Xt);

    % --- update ---
    % set current state to last state of the simulation
    Xt = X(end, :)';
    tstart = tend;
    tend = tstart + dt_sim;

    % --- log ---
    lent = length(t(2:end));
    tout = [tout;t(2:end)];
    Xout = [Xout;X(2:end,:)];
    Uout = [Uout;repmat(Ut',[lent,1])];
    Xdout = [Xdout;repmat(Xd(:,1)',[lent,1])];
    Udout = [Udout;repmat(Ud(:,1)',[lent,1])];
    Uext = [Uext;repmat(u_ext',[lent,1])];
    FSMout = [FSMout;repmat(FSM',[lent,1])];

    waitbar(ii/MAX_ITER, h_waitbar, "Calculating...");
end

close(h_waitbar);
fprintf('Calculation Complete!\n')
toc

%% Animation
[t,EA,EAd] = fig_animate(tout,Xout,Uout,Xdout,Udout,Uext,p);

