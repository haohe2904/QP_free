%-------------------------- Solve this problem on Onlique Manifold---------
%-------------------------- Measurement of Accuracy -----------------------
%close all; clear all; clc;
%%

% randstate = 11100+1000444433;%
% randn('state',double(randstate));
% rand('state',double(randstate));
randstate = 111+1000444433;%
randn('state',double(randstate));
rand('state',double(randstate));

% Set problem size
d = 50;            % [20 30 40 50 60]
s = 15;
Upsilon = 0.5;     % penalty parameter.

% Set Object matrix
A = rand(d,s); C = -A*A';
V = ones(s,1); V = V/norm(V,"fro");

% KKT residual
tolKKT = 16; %%%%%%%%%%%%%% log10 scale, i.e., 1e-* tolerance

%________Experiment_____
% common options
options.tolKKTres = 10^(-tolKKT); % tolKKTres = 16
options.maxOuterIter = 5000; % for RALM, REPMs.
options.maxiter = 50;  % for RSQP, QP-free.
options.maxtime = 300;  %%%%%% 60
options.verbosity = 3;  %          
% for RIPM 
options.KrylovIterMethod = 0;   
% for RALM, REPMs.
options.startingtolgradnorm = max(1e-3, 10^(-tolKKT + 3)); %%%%%% 1e-2
options.endingtolgradnorm = 10^(-16);
options.outerverbosity = options.verbosity;
% for RSQP, QP-free
options.mineigval_correction = 1e-5;  % 1e-5

%________for initial point_____
setting.initialpoint="feasible_region";

%________Setting________
setting.repeat = 1;
setting.row_dim = d;
setting.col_dim = s;

setting.tolKKTres =  options.tolKKTres;
setting.maxOuterIter = options.maxOuterIter;
setting.maxiter = options.maxiter;
setting.maxtime = options.maxtime;
setting.verbosity = options.verbosity;

setting.C = C;

setting.filepath = sprintf('nrep%d_Row%d_Col%d_KKTtol%.1e',...
    setting.repeat,setting.row_dim,setting.col_dim,setting.tolKKTres);

result = client_ob_acc(d, s, Upsilon, V, A, C, options, setting);


