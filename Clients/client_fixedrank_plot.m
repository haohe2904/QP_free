function data = client_fixedrank_plot(m, n, k, P, A, eqindices, ineqindices, methodoptions, setting)
%% Manifold factory

M = fixedrankembeddedfactory(m, n, k); % m: # of rows, n: # of cols, k: rank
problem.M = M;

%% Set-up Constraints

%%% nonnegativity (ineq) and equality (eq)
[eqnum,~] = size(eqindices);
[ineqnum,~] = size(ineqindices);

eq_constraints_cost = cell(eqnum,1);
eq_constraints_grad = cell(eqnum,1);
eq_constraints_hess = cell(eqnum,1);
nn_constraints_cost = cell(ineqnum, 1);
nn_constraints_grad = cell(ineqnum, 1);
nn_constraints_hess = cell(ineqnum, 1);

for i = 1 : eqnum
    [row,col] = ind2sub([m,n],eqindices(i));
    eq_constraints_cost{i} = @(Y) eqcostfun(Y, row, col);

    constraintgrad = zeros(m, n);
    constraintgrad(row, col) = 1;
    eq_constraints_grad{i} = @(U) constraintgrad;
            
    constrainthess = zeros(m, n);
    eq_constraints_hess{i} = @(X, U) constrainthess;
end

for j = 1 : ineqnum
    [row,col] = ind2sub([m,n],ineqindices(j));
    nn_constraints_cost{j} = @(Y) nncostfun(Y, row, col);

    constraintgrad = zeros(m, n);
    constraintgrad(row, col) = -1;
    nn_constraints_grad{j} = @(U) constraintgrad;
            
    constrainthess = zeros(m, n);
    nn_constraints_hess{j} = @(X, U) constrainthess;
end


function val = nncostfun(Y, row, col)
    Vt = Y.V.';
    val = - Y.U(row,:) * Y.S * Vt(:,col);
end

function val = eqcostfun(Y, row, col)
    Vt = Y.V.';
    val = Y.U(row,:) * Y.S * Vt(:,col) - A(row,col);
end

ineq_constraints_cost = nn_constraints_cost;
ineq_constraints_grad = nn_constraints_grad;
ineq_constraints_hess = nn_constraints_hess;

% constraints setting
problem.ineq_constraint_cost = ineq_constraints_cost;
problem.ineq_constraint_grad = ineq_constraints_grad;
problem.ineq_constraint_hess = ineq_constraints_hess;

problem.eq_constraint_cost = eq_constraints_cost;
problem.eq_constraint_grad = eq_constraints_grad;
problem.eq_constraint_hess = eq_constraints_hess;

%     Debug Only
%     checkconstraints_upto2ndorder(problem) 
     
condet = constraintsdetail(problem);

%% Setting a objective function
% Note that the observed elements (used as equality constraints) are
% eliminated when constructing the objective function.

% Define the problem cost function. The input X is a structure with
% fields U, S, V representing a rank k matrix as U*S*V'.
% f(X) = 1/2 * || P.*(X-A) ||^2
PA = P .* A;

problem.cost = @objcost;
function f = objcost(X)
    % Note that it is very much inefficient to explicitly construct the
    % matrix X in this way. Seen as we only need to know the entries
    % of Xmat corresponding to the mask P, it would be far more
    % efficient to compute those only.
    Xmat = X.U*X.S*X.V';
    f = .5*norm( P.*Xmat - PA , 'fro')^2;
end

% Define the Euclidean gradient of the cost function, that is, the
% gradient of f(X) seen as a standard function of X.
% nabla f(X) = P.*(X-A)
problem.egrad = @eobjgrad;
function G = eobjgrad(X)
    % Same comment here about Xmat.
    Xmat = X.U*X.S*X.V';
    G = P.*Xmat - PA;
end

% Define the Euclidean Hessian of the cost at X, along H, where H is
% represented as a tangent vector: a structure with fields Up, Vp, M.
% This is the directional derivative of nabla f(X) at X along Xdot:
% nabla^2 f(X)[Xdot] = P.*Xdot
problem.ehess = @euclidean_objhessian;
function ehess = euclidean_objhessian(X, H)
    % The function tangent2ambient transforms H (a tangent vector) into
    % its equivalent ambient vector representation. The output is a
    % structure with fields U, S, V such that U*S*V' is an mxn matrix
    % corresponding to the tangent vector H. Note that there are no
    % additional guarantees about U, S and V. In particular, U and V
    % are not orthonormal.
    ambient_H = problem.M.tangent2ambient(X, H);
    Xdot = ambient_H.U*ambient_H.S*ambient_H.V';
    % Same comment here about explicitly constructing the ambient
    % vector as an mxn matrix Xdot: we only need its entries
    % corresponding to the mask P, and this could be computed
    % efficiently.
    ehess = P.*Xdot;
end

%     DEBUG only
%      figure;
%      checkgradient(problem);
%      figure;
%      checkhessian(problem); 


%% Generating x0 such that it is strictly feasible.
if strcmp(setting.initialpoint, "eye")
    x0 = struct();
    x0.U = [eye(k);zeros(m-k,k)];
    x0.S = eye(k);
    x0.V = [eye(k);zeros(n-k,k)];
    
elseif strcmp(setting.initialpoint, "feasible_region")
    feasible_problem.M = M;
    feasible_problem.cost = @zerofun;
    feasible_problem.egrad = @egradzerofun;
    feasible_problem.ehess = @ehesszerofun;
    feasible_problem.ineq_constraint_cost = ineq_constraints_cost;
    feasible_problem.ineq_constraint_grad = ineq_constraints_grad;
    feasible_problem.ineq_constraint_hess = ineq_constraints_hess;
    feasible_problem.eq_constraint_cost = eq_constraints_cost;
    feasible_problem.eq_constraint_grad = eq_constraints_grad;
    feasible_problem.eq_constraint_hess = eq_constraints_hess; 
  %  feasible_x0 = M.rand();
    feasible_options.maxOuterIter = 200;
    feasible_options.maxtime = 40;  % 40
    feasible_options.outerverbosity = 1;  % 1
    feasible_options.tolKKTres = 10^(-2);  
    feasible_options.startingtolgradnorm = 1;
    feasible_options.endingtolgradnorm = feasible_options.tolKKTres;
    
    %     Debug Only
    %     checkconstraints_upto2ndorder(problem) 
    
    fprintf('Starting LQH to calculate a strictly feasible point\n');
    B=A;
    A=0.5*A;
    l = 0;
    while l < 10
    feasible_x0 = M.rand();
    [x0, ~, ~] = exactpenaltyViaSmoothinglqh(feasible_problem, feasible_x0, feasible_options);
    A=B;
    tag=0;
    if condet.has_ineq_cost
        for i=1:condet.n_ineq_constraint_cost
            ineqcosthandle = problem.ineq_constraint_cost{i};
            ineqcost=ineqcosthandle(x0);
            if ineqcost>=0
                tag = tag + 1; l = l + 1;
                continue;
            end
        end
    end
    if condet.has_eq_cost
        for i=1:condet.n_eq_constraint_cost
            eqcosthandle = problem.eq_constraint_cost{i};
            eqcost=eqcosthandle(x0);
            if eqcost>=0
                tag = tag + 1; l = l + 1;
                continue;
            end
        end
    end
    if tag==0
        break;
    end
    end 
    
  %  filename = sprintf('RC_nnlc_LQH_feasible_initial_point_%s.csv',filepath);
   % struct2csv(info, filename);      
else
    x0 = M.rand();
end
setting.x0 = x0.U * x0.S * x0.V';
%A=B;
if l >= 10
    [x0.U, x0.S, x0.V] = svds(0.1*A, k);
end

%% Calculating by solvers

options = methodoptions;

fprintf('Starting QP_free \n');
 [xfinal, costfinal, residual, KKTres, Time, X, stopreason] = QP_free_rank2_plot(m, n, k, P, A, eqindices, ineqindices, x0, options);
 data.KQPfree = KKTres; data.TQPfree = Time; data.XQPfree = X;

fprintf('Starting ALM \n');
[xfinal, info, residual, stopreason] = almbddmultiplier(problem, x0, options);
data.KALM = [info.KKT_residual]; data.TALM = [info.time];

fprintf('Starting Riemannian SQP \n');
[xfinal, costfinal, residual, info,~, ~] = SQP(problem, x0, options);
data.KSQP = [info.KKT_residual]; data.TSQP = [info.time];

fprintf('Starting LQH \n');
[xfinal, info, residual, stopreason] = exactpenaltyViaSmoothinglqh(problem, x0, options);
data.KLQH = [info.KKT_residual]; data.TLQH = [info.time];

fprintf('Starting LSE \n');
[xfinal, info, residual, stopreason] = exactpenaltyViaSmoothinglse(problem, x0, options);
data.KLSE = [info.KKT_residual]; data.TLSE = [info.time];

%     
%     filename = sprintf('RC_nnlc_Info_%s.csv',filepath);
%     struct2csv(setting, filename);
    
%% Sub functions

    function val = zerofun(X)
        val = 0;
    end

    function val = egradzerofun(X)
        val = zeros(m,n);
    end

    function val = ehesszerofun(X,H)
        val = zeros(m,n);
    end
end