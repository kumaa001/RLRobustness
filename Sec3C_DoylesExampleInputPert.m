%% DoylesExampleInputPert
% This file studies the example from  Doyle ('78 TAC) using the
% reinforcement learning framework.  The original example is in
% continuous-time but this file uses a discretized version of the example.
% An abstraction of random search is used to search for an dynamic
% output feedback controller. The random search evaluates the cost by
% averaging the performance over systems with random input gain
% perturbations.  These random perturbations sacrifice some performance
% (reward) but have the benefit of increasing the robustness the
% robustness margins at the input.
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
Opt.Npert = 10; 
Opt.PertLevel = 0;
Opt.Hypercube = [-0.2 0; -0.2 0; -40 0; 0 40];
Kfh = @Kcompanion;

% Repeat search for many trials and perturbation levels
Ntrial = 100;
PertPercent = [0 10 20 30 40 50];
Nlev = numel(PertPercent);

Jrs = zeros(Ntrial,Nlev);
DMrs = zeros(Ntrial,Nlev);
prs = zeros(4,Ntrial,Nlev);
rng(0)
for i1=1:Nlev
    fprintf('\n Level = %d of %d',i1,Nlev)
    
    % Repeat many trials at each level
    for i2=1:Ntrial
        % Random Search
        Opt.PertLevel = PertPercent(i1)/100;
        [J,p] = randomSearch(Kfh,A,Bu,Bw,C,Du,Q,R,W,V,Ts,Opt);
        Jrs(i2,i1) = J;
        if ~isempty(p)
            prs(:,i2,i1) = p;
            
            
            % Gain/Phase/Disk Margins
            [Ars,Brs,Crs] = Kfh(p);
            [GM,PM,DM] = computeMargins(A,Bu,C,Du,Ars,Brs,Crs,Ts);
        else
            DM=inf;
        end
        DMrs(i2,i1) = DM;
    end        
end
fprintf('\n');

%% Plot LQG Results and average results
% Uncomment next line to load the results used to generate plot in paper:
%  load Sec3C_DoylesExampleInputPert_Results

Javg = zeros(1,Nlev);
Jstd = zeros(1,Nlev);
DMavg = zeros(1,Nlev);
DMstd = zeros(1,Nlev);
for i1=1:Nlev
    % Compute average and standard deviations for cost and DM     
    % (After removing infinite entries)
    Jtmp = Jrs(:,i1);
    Jtmp( isinf(Jtmp) ) = [];
    Javg(i1) = mean(Jtmp);
    Jstd(i1) = std(Jtmp);
    
    DMtmp = DMrs(:,i1);
    DMtmp( isinf(DMtmp) ) = [];
    DMavg(i1) = mean(DMtmp);
    DMstd(i1) = std(DMtmp);
    
    % Plot Results
    figure(3)
    ph1 = plot(PertPercent(i1),Jrs(:,i1),'bx');
    xlabel('Pert Percent')
    ylabel('LQG Cost')
    hold on;
    
    figure(4)
    ph2 = plot(PertPercent(i1),DMrs(:,i1),'bx');
    xlabel('Pert Percent')
    ylabel('Disk Margin')
    hold on;
    drawnow
end

% Plot LQG results and mean/standard deviation of RL results
lw = 2;
figure(3)
ph3 = plot(PertPercent,Jlqg*ones(1,Nlev),'r--','LineWidth',lw);
ph4 = plot(PertPercent,Javg,'c--',PertPercent,Javg+Jstd,'c--',...
    PertPercent,Javg-Jstd,'c--','LineWidth',lw);
xlabel('Pert Percent')
ylabel('LQG Cost')
grid on
legend([ph3; ph1(1); ph4(1)],'Optimal LQG','Rand. Search',...
    'Rand. Search Mean+/-Std','Location','Northwest')
hold off;

figure(4)
ph5 = plot(PertPercent,DMlqg*ones(1,Nlev),'r--','LineWidth',lw);
ph6 = plot(PertPercent,DMavg,'c--',PertPercent,DMavg+DMstd,'c--',...
    PertPercent,DMavg-DMstd,'c--','LineWidth',lw);
xlabel('Pert Percent')
ylabel('Disk Margin')
grid on
legend([ph5; ph2(1); ph6(1)],'Optimal LQG','Rand. Search',...
    'Rand. Search Mean+/-Std','Location','Northwest')
hold off;
ax = axis; ax(3) = 1; axis(ax);

