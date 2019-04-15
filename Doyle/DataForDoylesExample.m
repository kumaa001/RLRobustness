function [A,Bu,Bw,C,Du,Q,R,W,V] = DataForDoylesExample(Ts)
%% DataForDoylesExample
% This file generates plant/cost data for the example from  Doyle 
% ('78 TAC). The original example is in continuous-time and shows that  
% LQG regulators can have poor robustness margins. This file constructs
% a discrete-time version of the example.
%
% The generated data is for the plant:
%  x(k+1) = A x(k) + Bu u(k) + Bw w(k) 
%    y(k) = C x(k) + Du u(k) + v(k)
% These dynamics are obtained by discretizing the continuous-time data 
% from Doyle's example with sample time Ts. The process & sensor noise 
% are white noise, zero mean, Gaussian with variances W = E[w^2] and 
% V = E[v^2].
% The infinite horizon cost is:
%     E[ sum_{k=0}^infty x(k)' Q x(k) + u(k)' R u(k) ]
%
% Ref: J.C. Doyle, Guaranteed Margins for LQG Regulators, IEEE TAC,
% vol. 23, no. 4, p.756-757, 1978.

%% Plant Data (Continuous-Time)
%  dx/dt = A x + Bu u + Bw w + Bv v
%      y = C x + Du u + Dw w + Dv v;
A = [1 1; 0 1];

Bu = [0; 1];
Bw = [1; 1];
% Bv = [0; 0];  % Not used given the form of the output DT system

C = [1 0];

Du = 0;
% Dw = 0;   % Not used given the form of the output DT system
% Dv = 1;   % Not used given the form of the output DT system

%% LQG Cost Data
% Continuous Time cost is E[ int_0^inf x'Qx + u'Ru dt ]
% Process and sensor noise are white noise, zero mean, Gaussian
% with variances W = E[w^2] and V = E[v^2].
q = 1000;
Q = q*[1;1]*[1 1];
R = 1;

sig = q;
W = sig;
V = 1;

%% Discrete-Time Plant Dynamics

% Discretize dynamics
Gc = ss(A,[Bu Bw],eye(2),0);
Gd = c2d(Gc,Ts);

% Overwrite state matrices (A,B) with discrete-time values
[A,B] = ssdata(Gd);
Bu = B(:,1);
Bw = B(:,2);
