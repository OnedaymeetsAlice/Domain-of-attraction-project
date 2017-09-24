% Define variables
x = sdpvar(2,1);
% Define constraints and objective
Constraints = [sum(x) <= 1, x(1)==0, x(2) >= 0.5];
Objective = x'*x+norm(x,1);
% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
% Solve the problem
sol = optimize(Constraints,Objective,options);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution = value(x)
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end