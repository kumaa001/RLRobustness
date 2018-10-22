%% DoylesExampleRL
% This file studies the example from  Doyle ('78 TAC) using the
% reinforcement learning framework.  The original example is in
% continuous-time but this file uses a discretized version of the example.
% An abstraction of random search is used to search for an dynamic
% output feedback controller.
%
% Ref: J.C. Doyle, Guaranteed Margins for LQG Regulators, IEEE TAC,
% vol. 23, no. 4, p.756-757, 1978.

%% Set random seed generator
% This is used to obtain repeatable results. Comment out the next
% line to generate (new) random results each time the file is run.
rng(0)

%% (Discrete-Time) Plant and Cost Data
% Plant is in the form:
%    x(k+1) = A x(k) + Bu u(k) + Bw w(k)
%      y(k) = C x(k) + Du u(k) + v(k)
Ts = 0.1;
[A,Bu,Bw,C,Du,Q,R,W,V] = DataForDoylesExample(Ts);

%% Optimal LQG Controller
% Use lqr/kalman to solve for the state feedback and estimator gains

% LQR Gain: u = -Klqr*x
[Klqr,Pc] = dlqr(A,Bu,Q,R);

% Kalman Filter Gain
[Lt,Pe] = dlqr(A',C',Bw*W*Bw',V);
Lkf = Lt';

% Form LQG Controller (Input is y and output is u)
Alqg = A-Bu*Klqr-Lkf*C;
Blqg = Lkf;
Clqg = -Klqr;
Klqg = ss(Alqg,Blqg,Clqg,0,Ts);

% Compute Infinite-Horizon Cost
Jlqg = computeIHCost(A,Bu,Bw,C,Du,Alqg,Blqg,Clqg,Q,R,W,V);

% Gain/Phase/Disk Margins
[GMlqg,PMlqg,DMlqg] = computeMargins(A,Bu,C,Du,Alqg,Blqg,Clqg,Ts);

% LQG Controller in Companion Form
% plqg is the 4-by-1 parameter vector corresponding to the LQG controller
Kcan = canon( Klqg, 'companion');
plqg = [Kcan.A(1,2); Kcan.A(2,2); Kcan.C'];

%% Random Search
% Nrollout = # of rollouts
% Npert = # of perturbed plants per roll-out
% PertLevel = Relative perturbation mag for each entry of plant matrices
% Hypercube = lower and upper bounds on controller parameters

% Perform Random Search
Opt.Nrollout = 1e5; 
Opt.Npert = 1;
Opt.PertLevel = 0;
Opt.Hypercube = [-0.2 0; -0.2 0; -40 0; 0 40];
Kfh = @Kcompanion;
[Jrs,prs] = randomSearch(Kfh,A,Bu,Bw,C,Du,Q,R,W,V,Ts,Opt);

% Gain/Phase/Disk Margins
[Ars,Brs,Crs] = Kfh(prs);
[GMrs,PMrs,DMrs] = computeMargins(A,Bu,C,Du,Ars,Brs,Crs,Ts);

% Compare LQG/RS Costs (= -Reward )
fprintf('\n Cost (=-Reward): LQG = %4.4e \t RS = %4.4e \t Ratio = %4.4f',...
    Jlqg,Jrs,Jrs/Jlqg);

% Compare LQG/RS Margins
fprintf('\n Lower Gain Margins: LQG = %4.4f \t RS = %4.4f',GMlqg(1),GMrs(1));
fprintf('\n Upper Gain Margins: LQG = %4.4f \t RS = %4.4f',GMlqg(2),GMrs(2));
fprintf('\n Phase Margins: LQG = %4.4f \t RS = %4.4f',PMlqg,PMrs);
fprintf('\n Disk Margins: LQG = %4.4f \t RS = %4.4f\n',DMlqg,DMrs);
fprintf('\n\n');

% Compare LQG / RS controller companion form parameters
[plqg prs]

% Results for RNG(0)
% Cost (-Reward): LQG = 1.3728e+05 	 RS = 1.9197e+05 	 Ratio = 1.3983
% Lower Gain Margins: LQG = 0.9802 	 RS = 0.8929
% Upper Gain Margins: LQG = 1.0007 	 RS = 1.0433
% Phase Margins: LQG = 0.0702 	 RS = 2.7427
% Disk Margins: LQG = 1.0007 	 RS = 1.0433    
    