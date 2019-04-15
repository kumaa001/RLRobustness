%% DoylesExampleLQG
% This file studies the example from  Doyle ('78 TAC) showing that LQG
% regulators can have poor margins. The original example is in
% continuous-time. This file performs an analysis on a discretized
% version of the example.
%
% Ref: J.C. Doyle, Guaranteed Margins for LQG Regulators, IEEE TAC,
% vol. 23, no. 4, p.756-757, 1978.

%% (Discrete-Time) Plant and Cost Data
% Plant is in the form:
%    x(k+1) = A x(k) + Bu u(k) + Bw w(k)
%      y(k) = C x(k) + Du u(k) + v(k)
Ts = 0.1;
[A,Bu,Bw,C,Du,Q,R,W,V] = DataForDoylesExample(Ts);

%% Solve Discrete-Time LQG Problem
% Use lqr/kalman to solve for the state feedback and estimator gains

% LQR Gain: u = -Klqr*x
[Klqr,Pc] = dlqr(A,Bu,Q,R);
%[Klqr2,Pc2] = lqr( ss(A,Bu,eye(2),0,1), Q,R);

% Compare with analytical (continuous-time) solution Ka from Doyle
% Note: The CT and DT state-feedback gains should be similar if the
% sample time is "small", e.g. if Ts=1e-3.
q = Q(1,1);
f = 2+sqrt(4+q);
Ka = [1 1]*f;

fprintf('\n State Feedback Gains:')
fprintf('\n Disc-Time: K = [%4.4f, %4.4f]',Klqr(1),Klqr(2))
fprintf('\n Cont-Time: K = [%4.4f, %4.4f]\n',Ka(1),Ka(2))

% Kalman Filter Gain
[Lt,Pe] = dlqr(A',C',Bw*W*Bw',V);
Lkf = Lt';
%[~,Lkf2,Pe2] = kalman( ss(A,[Bu Bw],C,[Du Dw],1),W,V);

% Check DARE solutions and LQR/Kalman gains.
% The DARE results should be = 0.  The LQG/Kalman gains should match
% the results given by the expressions involving the DARE solutions.
if false
    A'*Pc*A-A'*Pc*Bu*inv(R+Bu'*Pc*Bu)*Bu'*Pc*A+Q - Pc
    [Klqr; inv(Bu'*Pc*Bu+R)*Bu'*Pc*A]
    
    A*Pe*A'-A*Pe*C'*inv(V+C*Pe*C')*C*Pe*A'+Bw*W*Bw' - Pe
    [Lkf, A*Pe*C'*inv(C*Pe*C'+V)]
end

% Compare with analytical (continuous-time) solution La from Doyle
% Note: The CT and DT Kalman filter gains should be similar if the
% sample time is "small", e.g. if Ts=1e-3, and the DT gain is
% appropriately scaled by Ts.
sig = W;
d = 2+sqrt(4+sig);
La = [1; 1]*d;

fprintf('\n Kalman Filter Gains:')
fprintf('\n Disc-Time: L = [%4.4f, %4.4f]',Lkf(1),Lkf(2))
fprintf('\n Cont-Time: L*Ts = [%4.4f, %4.4f]\n\n',La(1)*Ts, La(2)*Ts)

%% Compute Margins for DT Example

% LQG Controller
Alqg = A-Bu*Klqr-Lkf*C;
Blqg = Lkf;
Clqg = -Klqr;
Klqg= ss(Alqg,Blqg,Clqg,0,Ts);

% Gain/Phase/Disk Margins
[GM,PM,DM] = computeMargins(A,Bu,C,Du,Alqg,Blqg,Clqg,Ts)

%% Infinite-Horizon Cost

% Compute steady-state solution using Riccati solutions
% (This holds form only holds for the optimal LQG controller)
Jss = trace( Pc*(Bw*W*Bw') ) + trace( (Q+A'*Pc*A-Pc)*Pe)

% Compute steady-state cost using Lyapunov solution
% The solution of a Lyapunov equation yields the steady-state variance of
% the closed-loop state. This is used to compute the steady-state cost.  
Jss2 = computeIHCost(A,Bu,Bw,C,Du,Alqg,Blqg,Clqg,Q,R,W,V)

%% Additional Check: Solve LQG Problem Using LQG Command
% Note: In R2018b the LQG command has options to construct the optimal
% controller in either "current" or "delayed" form. The default is
% "delayed" and this uses x[n|n-1] as the state estimate.  This introduces
% a one step delay from measurement y[n] to control input u[n].  The
% alternative is "current" which uses x[n|n] as the state estimate. The
% code below uses the default "delayed" form which is consistent with the
% form used above.  These options don't exist in R2016b. Moreover, the
% code in R2016b does not return the second output (with gains, etc) and
% hence the code below will error out.  I'm not sure if the updated LQG
% code exists in R2017a-R2018a.
if false
    G = ss(A,Bu,C,Du,Ts);
    Qxu = blkdiag(Q,R);
    Qwv = blkdiag(Bw*W*Bw',V);
    [Klqg2,lqginfo] = lqg( G, Qxu, Qwv );
    
    % Compare lqg controllers (Norm should be close to zero)
    norm(Klqg-Klqg2,inf)
    
    % Compare Kalman gains (Norm should be close to zero)
    norm(Lkf-lqginfo.L)
    
    % Compare LQR gains (Norm should be close to zero)
    norm(Klqr-lqginfo.Kx)
end

%% Additional Check: Solve LQG Problem Using H2SYN
% Use H2SYN to compute the optimal controller. This is simply
% for comparison. The construction of the "error" assumes that w
% and v are scalar signals.
%
% Note: The optimal H2 controller corresponds to the "current" form
% returned by the LQG command. In Matlab R2016b, H2SYN actually returned
% the controller in "delayed" form. Thus the results below (e.g. cost Jss3
% and controller Kh2) will agree with all other results in the file if
% run in R2016b.  However, H2SYN was updated by R2018b to correctly return
% the optimal controller in "current" form.  This yields a slightly lower
% cost for Jss3 and slightly different controller as compared to the
% results given above.  I'm not sure which form (current or delayed) is
% returned by H2SYN from R2017a - R2018a.
if false
    Bh2 = [Bw*sqrt(W) zeros(2,1) Bu];
    Ch2 = [sqrt(q)*[1 1]; 0 0; C];
    D = [0 0 0; 0 0 sqrt(R); 0 sqrt(V) Du];
    P = ss(A,Bh2,Ch2,D,Ts);
    [Kh2,CLh2,GAMh2,INFOh2] = h2syn(P,1,1);
    
    % The optimal H2 cost squared should equal the steady-state LQG cost.
    Jss3 = GAMh2^2
    
    % Compare optimal LQG and H2SYN controllers (They should match)
    figure(1);
    bode(Klqg,'b',Kh2,'r--')
end
