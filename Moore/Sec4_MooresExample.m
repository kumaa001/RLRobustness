%% MooresExampleInputPert
% 
% Ref 1: J.B. Moore, D. Gangsaas, and J. Blight, Performance and 
% Robustness Trades in LQG Regulator Design, IEEE CDC, 1981.
%
% Ref 2: B.D.O. Anderson and J.B. Moore, Optimal Control: Linear Quadratic
% Methods, Prentice Hall, 1989.

%% Set random seed generator
% This is used to obtain repeatable results. Comment out the next
% line to generate (new) random results each time the file is run.
rng(0)

%% (Discrete-Time) Plant and Cost Data
% Plant is in the form:
%    x(k+1) = A x(k) + Bu u(k) + Bw w(k)
%      y(k) = C x(k) + Du u(k) + v(k)

% Continuous-time data
% A = [-1 0 0 1;0 0 1 0;0 -100 -0.2 0;0 0 0 -1];
wn = 10;
zeta = 1e-2;  
A = [-1 0 0 1;0 0 1 0;0 -wn^2 -2*zeta*wn 0;0 0 0 -1];
Bu = [1; 0; 100; 0];
Bw = [0; 0; 0; 0.45];
C = [1 10 0 1];
Du  = 0;

% Variance and cost data
W = 1;
V = 0.01;
Q = [4 0 0 0;zeros(3,4)];
R = 1;

% Discretize dynamics
Ts = 0.09;
Gc = ss(A,[Bu Bw],C,0);
Gd = c2d(Gc,Ts);
[A,B] = ssdata(Gd);
Bu = B(:,1);
Bw = B(:,2);

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

%% Random Search
% Nrollout = # of rollouts
% Npert = # of perturbed plants per roll-out
% PertLevel = Relative perturbation mag for each entry of plant matrices
% Hypercube = lower and upper bounds on controller parameters

% Perform Random Search agumaneted with Gradient Descent
Opt.Nrollout = 3e4;
Opt.Npert = 10; 
Opt.PertLevel = 0;
Opt.Hypercube = [0 1; -2 0; 0 2; -0.1 0;0 0.3; -0.3 0];
Kfh = @Kcontrollable3;

% Hyper params for Gradient descent
% The Step size of each descent and number of descents are stored bellow 
Opt.StepSize = 1;
Opt.DescentCount = 100;

% Repeat search for many trials and perturbation levels
Ntrial = 40;

%PertPercent = [0];

PertPercent = [0 10 20 30 40 50];
Nlev = numel(PertPercent);

Jrs = zeros(Ntrial,Nlev);
DMrs = zeros(Ntrial,Nlev);
prs = zeros(6,Ntrial,Nlev);
rng(0)
for i1=1:Nlev
    fprintf('\n Level = %d of %d',i1,Nlev)
    
    % Repeat many trials at each level
    for i2=1:Ntrial
        % Random Search
        Opt.PertLevel = PertPercent(i1)/100;
        [J,p] = gdRandomSearch(Kfh,A,Bu,Bw,C,Du,Q,R,W,V,Ts,Opt);
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
%  load Sec4_MooresExample_Results
Javg = zeros(1,Nlev);
Jstd = zeros(1,Nlev);
DMavg = zeros(1,Nlev);
DMstd = zeros(1,Nlev);
for i1=1:5
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
ph3 = plot(PertPercent(1:5),Jlqg*ones(1,5),'r--','LineWidth',lw);
ph4 = plot(PertPercent(1:5),Javg(1,1:5),'c--',PertPercent(1:5),Javg(1,1:5)+Jstd(1,1:5),'c--',...
    PertPercent(1:5),Javg(1,1:5)-Jstd(1,1:5),'c--','LineWidth',lw);
xlabel('Pert Percent')
ylabel('LQG Cost')
grid on
legend([ph3; ph1(1); ph4(1)],'Optimal LQG','GD with RS',...
    'GD with RS Mean+/-Std','Location','Northwest')
hold off;

figure(4)
ph5 = plot(PertPercent(1:5),DMlqg*ones(1,5),'r--','LineWidth',lw);
ph6 = plot(PertPercent(1:5),DMavg(1,1:5),'c--',PertPercent(1:5),DMavg(1,1:5)+DMstd(1,1:5),'c--',...
    PertPercent(1:5),DMavg(1,1:5)-DMstd(1,1:5),'c--','LineWidth',lw);
xlabel('Pert Percent')
ylabel('Disk Margin')
grid on
legend([ph5; ph2(1); ph6(1)],'Optimal LQG','GD with RS',...
    'GD with RS Mean+/-Std','Location','Northwest')
hold off;
ax = axis; ax(3) = 1; axis(ax);

