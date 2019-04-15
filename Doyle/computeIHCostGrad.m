function [JG] = computeIHCostGrad(A,Bui,Bw,C,Dui,Q,R,W,V,p,pc,S,Mcl,Acl)
% Specific to Doyle's Example
% This function computes the Gradient of the cost function of the closed
% loop system at the point on parameter space p along the direction defined 
% by the parameter index pc  

% Del_J = Trace(Mcl*Del_S + S*Del_Mcl) where Del is the Gradient operator 

%% Compute the gradient of the Mcl 
if(pc<=2)
    Del_Mcl = zeros(4);
elseif(pc == 3)
    Del_Mcl = [zeros(2,4);zeros(2) R.*[2*p(3) p(4);p(4) 0]];
else
    Del_Mcl = [zeros(2,4);zeros(2) R.*[0 p(3);p(3) 2*p(4)]];    
end



%% Compute the Gradient of S
% Del_S = Acl*Del_S*Acl' + Del_Acl*S*Acl' + Acl*S*Del_Acl'. The above Del_S can be computed
% using iteration as the solution is convergent or can be solved as a
% Discrete lyapanov equation just as done in the compute cost function

% Del_Acl is assigned in one of the following cases
if(pc == 1)
    Del_Acl = [zeros(2,4);zeros(2) [0 1;0 0]];
elseif(pc == 2)
    Del_Acl = [zeros(2,4);zeros(2) [0 0;0 1]];
elseif(pc == 3)
    Del_Acl = [zeros(2) Bui zeros(2,1);zeros(2,4)];
else
    Del_Acl = [zeros(2,3) Bui;zeros(2,4)];
end

Del_S = dlyap(Acl,(Del_Acl*S*Acl'+Acl*S*Del_Acl'));
    
JG = trace(Mcl*Del_S + Del_Mcl*S);

end

