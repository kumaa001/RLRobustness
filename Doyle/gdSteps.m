function [Jgd,Kp] = gdSteps(A,Bui,Bw,C,Dui,Aki,Bki,Cki,Q,R,W,V,Opt)

% Specific for Doyle's Example
% This function takes in the intial controller checks for stability and then 
% uses gradient descent to reduce the cost for a fixed nuumber of
% descent iterations

[Jk,S,Acl,Mcl]  = computeIHCost(A,Bui,Bw,C,Dui,Aki,Bki,Cki,Q,R,W,V);

if (Jk == inf)
    %disp('unstable control')
    Jgd = inf;
    Kp = [Aki(:,2)' Cki]';
    return;
end


%% Gradient descent loop
stepCount = Opt.DescentCount;
eta = Opt.StepSize;
Jgd = Jk;
Ak = Aki; Bk = Bki; Ck = Cki;
p = [Ak(:,2)' Ck]';
count = 0;


while(count <= stepCount)
   % Compute the cost function gradient at the current controller
   JG1 = computeIHCostGrad(A,Bui,Bw,C,Dui,Q,R,W,V,p,1,S,Mcl,Acl); 
   JG2 = computeIHCostGrad(A,Bui,Bw,C,Dui,Q,R,W,V,p,2,S,Mcl,Acl); 
   JG3 = computeIHCostGrad(A,Bui,Bw,C,Dui,Q,R,W,V,p,3,S,Mcl,Acl); 
   JG4 = computeIHCostGrad(A,Bui,Bw,C,Dui,Q,R,W,V,p,4,S,Mcl,Acl);
   
   count = count + 1;
   
   % Storing previous controller parameters
   pPrv = p;
   MclPrv = Mcl;
   AclPrv = Acl;
   SPrv = S;
   JkPrv = Jk;

   % Gradient vector
   JG = [JG1;JG2;JG3;JG4];
   p = p - eta.*JG;
   
   % Controller with the new given parameter
   [Ak,Bk,Ck] = Kcompanion(p);
   
   % Compute the new controller cost
   
   [Jk,S,Acl,Mcl]  = computeIHCost(A,Bui,Bw,C,Dui,Ak,Bk,Ck,Q,R,W,V);
   
   % Comparing the new cost Jk with JkPrv to decide if the step size has to
   % be reduced
   if(Jk > JkPrv && eta >= 1e-12)
       % disp(['trying smaller step size of :',num2str(eta/10)]);
       eta = eta/10;
       p = pPrv;
       Jk = JkPrv;
       Mcl = MclPrv;
       Acl = AclPrv;
       S = SPrv;
   else
       if(Jk > JkPrv && eta < 1e-12)
          Jk = JkPrv;
          [Ak,Bk,Ck] = Kcompanion(pPrv);
          break;
       elseif(eta < 1e-10)
          %disp(['reverting to original stepsize :',num2str(1e-11)]);
          eta = Opt.StepSize;
       end
   end
   
   
end

Jgd = Jk;
Kp = [Aki(:,2)' Cki]';

end


