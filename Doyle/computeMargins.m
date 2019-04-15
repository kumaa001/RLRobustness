function [GM,PM,DM] = computeMargins(A,Bu,C,Du,Ak,Bk,Ck,Ts)
%% computeMargins
% This function computes the classical gain/phase margins and the
% symmetric disk gain margins for a plant (A,Bu,C,Du) in *positive 
% feedback* with a controller (Ak,Bk,Ck,0). 
%
% The plant dynamics are expressed without the process/sensor noise:
%    x(k+1) = A x(k) + Bu u(k) 
%      y(k) = C x(k) + Du u(k) 
%
% The controller dynamics assume no feedthrough (Dk=0):
%    z(k+1) = Ak z(k) + Bk y(k)
%      u(k) = Ck z(k)


% Form plant (without noises, i.e. input u and output y)
G = ss(A,Bu,C,Du,Ts);

% Form Controller (input y and output u)
K = ss(Ak,Bk,Ck,0,Ts);

% Check stability
% The syntax below assumes K in positive feedback
CLstable = isstable( feedback(G,K,+1) );
if ~CLstable
    % No margin for nominally unstable feedback systems
    GM = [1 1]; PM = 0; DM = 1;
    return
end

% Classical Gain, Phase, and Symmetric Disk Margins
% The loopmargin syntax below assumes K in positive feedback
[cm,dm] = loopmargin(G,-K);

% Extract most restrictive gain margins
% (loop margin returns gain margins associated with all phase crossover
%  frequencies. The code below extracts the most restrictive margins)
AllGM = cm.GainMargin;
GMlb = max( AllGM( AllGM<1 ) );
if isempty(GMlb)
    GMlb = 0;
end
GMub = min( AllGM( AllGM>1 ) );
if isempty(GMub)
    GMub = inf;
end
GM =  [GMlb GMub];

% Extract most restrictive phase margin
AllPM = cm.PhaseMargin;
PM = min( abs(AllPM) );
if isempty(PM)
    PM = inf;
end

% Extract (upper) disk margin
DM = dm.GainMargin(2);


