function [Jrs,prs] = gdRandomSearch(K,A,Bu,Bw,C,Du,Q,R,W,V,Ts,Opt)
%% randomSearch
% Add documentation

%% Options for random search
Nro = Opt.Nrollout;
Npert = Opt.Npert;
PertLevel = Opt.PertLevel;
H = Opt.Hypercube;

%% Create Perturbed Systems
% The input to the plant from the controller is perturbed by a random,
% multiplicative factor of the control input.
if PertLevel==0
    % No perturbations
    Bp = Bu; 
    Dp = Du;
    Npert = 1;
else
    % Randomly Generated Perturbations
    % XXX - Always include Nominal as one case (?)
    % Perts to the inputs up = u (1+/-delta); 
    % XXX - Assumes single input....
    Pert = 1+PertLevel*(2*rand(1,Npert)-1);
    Pert(1) = 1;
    Bp = Bu*Pert;
    Dp = Du*Pert;
end

%% Random Search Loop

% Pre-compute random samples of controller parameters
Hlb = H(:,1);
Hdiff = H(:,2) - H(:,1);
Np = size(H,1);
p = Hlb + Hdiff.*rand(Np,Nro);

% Perform Gradient descent with random search
Jrs = inf;
Jpert = zeros(Npert,1);
prs = [];
for i1 = 1:Nro
    % Form controller from i^th parameter sample
    pi = p(:,i1);    
    [Ak,Bk,Ck] = K( pi );
    
    % GD being applied on the nominal plant to reduce cost and get the best possible
    % controller cost
    [Jn,pi] = gdSteps(A,Bu,Bw,C,Du,Ak,Bk,Ck,Q,R,W,V,Opt);
    [Ak,Bk,Ck] = K( pi );    
                   
    % Cost for each perturbed plant
    for i2=1:Npert        
        Bui = Bp(:,i2);
        Dui = Dp(:,i2);
        Jpert(i2) = computeIHCost(A,Bui,Bw,C,Dui,Ak,Bk,Ck,Q,R,W,V);
        %[Jpert(i2),pi] = gdSteps(A,Bui,Bw,C,Dui,Ak,Bk,Ck,Q,R,W,V,Opt);
    end
    Ji = sum(Jpert)/Npert;
    
    % Compare current sample against best sample thus far
    if Ji < Jrs
        Jrs = Ji;
        prs = pi;
    end
end

% if nargout==3
%     % Construct Optimal Controller
%     [Ak,Bk,Ck] = K(prs);
%     Krs = ss(Ak,Bk,Ck,0,Ts);
% end

