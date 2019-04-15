function [JG] = computeIHCostGrad(A,Bui,Bw,C,Dui,Q,R,W,V,p,pc,S,Mcl,Acl)

% This function computes the Gradient of the cost function of the closed
% loop system at the point on parameter space p along the direction defined 
% by the parameter index pc  

% Del_J = Trace(Mcl*Del_S + S*Del_Mcl) where Del is the Gradient operator 

%% Compute the gradient of the Mcl 
if(pc<=3)
    Del_Mcl = zeros(7);
elseif(pc == 4)
    Del_Mcl = [zeros(4,7);zeros(3,4) R.*[2*p(4) p(5) p(6);p(5) 0 0;p(6) 0 0]];
elseif(pc == 5)
    Del_Mcl = [zeros(4,7);zeros(3,4) R.*[0 p(4) 0;p(4) 2*p(5) p(6);0 p(6) 0]];
else
    Del_Mcl = [zeros(4,7);zeros(3,4) R.*[0 0 p(4);0  0 p(5);p(4) p(5) 2*p(6)]];
end


%% Compute the Gradient of S
% Del_S = Acl*Del_S*Acl' + Del_Acl*S*Acl' + Acl*S*Del_Acl'. The above Del_S can be computed
% using iteration as the solution is convergent or can be solved as a
% Discrete lyapanov equation just as done in the compute cost function

% Del_Acl is assigned in one of the following cases
if(pc == 1)
    Del_Acl = [zeros(4,7);zeros(3,4) [0 0 0;0 0 0;1 0 0]];
elseif(pc == 2)
    Del_Acl = [zeros(4,7);zeros(3,4) [0 0 0;0 0 0;0 1 0]];
elseif(pc == 3)
    Del_Acl = [zeros(4,7);zeros(3,4) [0 0 0;0 0 0;0 0 1]];
elseif(pc == 4)
    Del_Acl = [zeros(4,4) Bui zeros(4,2);zeros(3,7)];
elseif(pc == 5)
    Del_Acl = [zeros(4,5) Bui zeros(4,1);zeros(3,7)];
else
    Del_Acl = [zeros(4,6) Bui;zeros(3,7)];
end

Del_S = dlyap(Acl,(Del_Acl*S*Acl'+Acl*S*Del_Acl'));
    
JG = trace(Mcl*Del_S + Del_Mcl*S);

end

