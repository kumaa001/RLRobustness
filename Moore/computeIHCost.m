function [J,Xcl,Acl,Mcl,Bcl,Wcl] = ...
    computeIHCost(A,Bu,Bw,C,Du,Ak,Bk,Ck,Q,R,W,V)
%% computeIHCost
% This function computes the infinite-horizon LQG quadratic cost J for a 
% plant and controller with the following assumptions:
%
% Plant Dynamics:
%    x(k+1) = A x(k) + Bu u(k) + Bw w(k)
%      y(k) = C x(k) + Du u(k) + v(k)
% Process & sensor noise are white noise, zero mean, Gaussian with 
% variances W = E[w^2] and V = E[v^2].
%
% Controller Dynamics:
%    z(k+1) = Ak z(k) + Bk y(k)
%      u(k) = Ck z(k)
% Note that it is assumed the controller has no feedthrough (Dk=0).
% The infinite horizon cost is:
%     E[ sum_{k=0}^infty x(k)' Q x(k) + u(k)' R u(k) ]

%% Signal Dimensions
Nx = size(A,1);
Nv = size(C,1);
Nw = size(Bw,2);
Nk = size(Ak,1);

%% Form Closed Loop from [w;v] to xcl = [x; z] 
%  xcl(k+1) = Acl xcl(k) + Bcl [w(k); v(k)]
Acl = [A Bu*Ck; Bk*C (Ak+Bk*Du*Ck)];
Bcl = [Bw zeros(Nx,Nv); zeros(Nk,Nw) Bk];

% Check Closed-Loop Stability
ecl = eig(Acl);
if max(abs(ecl))>=1
    % Closed-loop is unstable
    J = inf;
    Xcl = inf(7);
    Mcl = blkdiag(Q, Ck'*R*Ck);
    Wcl = blkdiag(W,V);
    return
end



%% Variance and Cost Matrices for Closed-Loop 
% Wcl is the variance of [w(k); v(k)]
% Cost per iteration is 
%   E[ x(k)' Q x(k) + u(k)' R u(k)] = E[xcl(k)' M xcl(k)]
Wcl = blkdiag(W,V);
Mcl = blkdiag(Q, Ck'*R*Ck);

%% Solve Lyapunov Equation for Steady-State Variance of xcl
Xcl = dlyap(Acl,Bcl*Wcl*Bcl');

%% Steady-State Quadratic Cost
J = trace(Mcl*Xcl);

% function J = computeIHCost(A,B,C,D,Ak,Bk,Ck,Q,R,Qw,Qv)
% %%
% % This function computes the infinite-horizon quadratic cost J for a plant 
% % (A,B,C,D) and controller (Ak,Bk,Ck,Dk) given LQG iteration cost 
% % matrices (Q,R) and process/sensor noise variances (Qw,Qv).  All 
% % dimensions are assumed to be consistent with Doyle's example.
% %
% % The plant inputs are assumed to be ordered as [w; v; u] where w is the
% % process noise, v is the sensor noise, and u is the control input. It
% % is also assumed that Dk=0 (no controller feedthrough) and Du = 0 (no 
% % plant feedthrough from u to y).
% 
% %% Form Closed Loop from [w;n] to xcl = [xG; xC] 
% %  xcl(k+1) = Acl xcl(k) + Bcl [w(k); v(k)]
% Bu = B(:,3);
% Bwv = B(:,1:2);
% Dwv = D(:,1:2);
% 
% Acl = [A Bu*Ck; Bk*C Ak];
% Bcl = [Bwv; Bk*Dwv];
% 
% %% Check Stability
% ecl = eig(Acl);
% if max(abs(ecl))>=1
%     % Closed-loop is unstable
%     J = inf;
%     return
% end
% 
% %% Cost Matrices for Closed-Loop 
% % Cost per iteration is 
% % E[ x(k)' Q x(k) + u(k)' R u(k)] = E[xcl(k)' M xcl(k)] = 
% Qcl = blkdiag(Qw,Qv);
% Mcl = blkdiag(Q, Ck'*R*Ck);
% 
% %% Solve Lyapunov Equation for Steady-State Variance of xcl
% S = dlyap(Acl,Bcl*Qcl*Bcl');
% 
% %% Steady-State Quadratic Cost
% J = trace(Mcl*S);
% 
