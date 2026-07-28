function data = client_ob_acc(n, k, Sigma, V, A, C, options, setting)
%% Manifold factory.
data = NaN(3, 5);
problem.M = obliquefactory(n,k);

%% Set-up objective function.
VVt=V*V';
problem.cost = @(X) trace(X'*C*X)+Sigma*(norm(X*V,'fro')^2-1);
problem.egrad = @(X) 2*C*X+X*(VVt);
problem.ehess = @(X,dX) 2*C*dX+dX*(VVt);


%% Set-up constraints.
ineqnum = n*k;
constraints_cost = cell(ineqnum, 1);
constraints_grad = cell(ineqnum, 1);
constraints_hess = cell(ineqnum, 1);

nnconst_idx = 1;
for row = 1: n
    for col = 1: k
        % cost
        constraints_cost{nnconst_idx} = @(Y) -Y(row, col);
        % egrad
        constraintgrad = zeros(n, k);
        constraintgrad(row, col) = -1;
        constraints_grad{nnconst_idx} = @(U) constraintgrad;
        % ehess
        constrainthess = zeros(n, k);
        constraints_hess{nnconst_idx} = @(X, U) constrainthess;
        nnconst_idx = nnconst_idx + 1;
    end
end

problem.ineq_constraint_cost = constraints_cost;
problem.ineq_constraint_grad = constraints_grad;
problem.ineq_constraint_hess = constraints_hess;

%%%% setting.initialpoint by LQH
if strcmp(setting.initialpoint, "feasible_region")
    fprintf('Starting LQH to calculate a strictly feasible point\n');
    feasible_problem.M = problem.M;
    feasible_problem.cost = problem.cost;
    feasible_problem.egrad = problem.egrad;
    feasible_problem.ehess = problem.ehess;
    feasible_problem.ineq_constraint_cost = constraints_cost;
    feasible_problem.ineq_constraint_grad = constraints_grad;
    feasible_problem.ineq_constraint_hess = constraints_hess;

    feasible_options.tolKKTres = 10^(-1); 
    feasible_options.maxOuterIter = 30; 
    feasible_options.maxiter = 30;  
    feasible_options.maxtime = 30;  
    feasible_options.verbosity = 1;  %            
    feasible_options.KrylovIterMethod = 0;       
    feasible_options.startingtolgradnorm = max(1e-3, 10^(-1 + 3)); 
    feasible_options.endingtolgradnorm = 10^(-1);
    feasible_options.outerverbosity = options.verbosity;
    feasible_x0 = problem.M.rand();
    [x0, ~, ~] = exactpenaltyViaSmoothinglqh(feasible_problem, feasible_x0, feasible_options);
    x0=abs(x0);
    
else
    x0=problem.M.rand();
    x0=abs(x0);
end


%% Calculating by solvers

fprintf('Starting QP_free \n');
[xfinal, costfinal, residual, stopreason] = QP_free_Oblique(n, k, A, C, x0, options);

% fprintf('Starting QP_free_mex \n');
% [xfinal, costfinal, residual, stopreason] = QP_free_Oblique_mex(n, k, A, C, x0, options);

fprintf('Starting Riemannian SQP \n');
[xfinal, costfinal, residual, info,~, stopreason] = SQP(problem, x0, options);

fprintf('Starting ALM \n');
[xfinal, info, residual, stopreason] = almbddmultiplier(problem, x0, options);

fprintf('Starting LQH \n');
[xfinal, info, residual, stopreason] = exactpenaltyViaSmoothinglqh(problem, x0, options);

fprintf('Starting LSE \n');
[xfinal, info, residual, stopreason] = exactpenaltyViaSmoothinglse(problem, x0, options);

end
