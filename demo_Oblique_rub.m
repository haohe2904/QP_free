%-------------------------- Solve this problem on Onlique Manifold---------
%-------------------------- Measurement of rubost -----------------------
close all; clear all; clc;
%%

restoredefaultpath;
addpath(genpath('clients'));
addpath(genpath('connectors'));
addpath(genpath('fix_rank_functions'));
addpath(genpath('manopt'));
addpath(genpath('solvers'));
addpath(genpath('supportingfunctions'));

% Set problem size
d = 70;            % [20 30 40 50 60]
s = 10;
Upsilon = 0.5;     % penalty parameter.

repeat = 20;
QPfree = zeros(1,repeat); QPfree_time = zeros(1,repeat);
SQP = zeros(1,repeat); SQP_time = zeros(1,repeat);
ALM = zeros(1,repeat); ALM_time = zeros(1,repeat);
LQH = zeros(1,repeat); LQH_time = zeros(1,repeat);
LSE = zeros(1,repeat); LSE_time = zeros(1,repeat);

for i = 1:repeat
    
    randstate = 111*i+1000;
    randn('state',double(randstate));
    rand('state',double(randstate));
    
    % Set Object matrix
    A = rand(d,s); C = -A*A';
    V = ones(s,1); V = V/norm(V,"fro");
    
    % KKT residual
    tolKKT = 9; %%%%%%%%%%%%%% log10 scale, i.e., 1e-* tolerance
    
    %________Experiment_____
    % common options
    options.tolKKTres = 10^(-tolKKT); % tolKKTres = 16
    options.maxOuterIter = 5000; % for RALM, REPMs.
    options.maxiter = 1000;  % for RSQP, QP-free.
    options.maxtime = 300;  %%%%%% 60
    options.verbosity = 2;  %          
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
    setting.repeat = repeat;
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
    
    result = client_ob_rub(d, s, Upsilon, V, A, C, options, setting);

    %% experimental results
    QPfree(i) = result.QPfree; QPfree_time(i) = result.QPfree_time;
    SQP(i) = result.SQP; SQP_time(i) = result.SQP_time;
    ALM(i) = result.ALM; ALM_time(i) = result.ALM_time;
    LQH(i) = result.LQH; LQH_time(i) = result.LQH_time;
    LSE(i) = result.LSE; LSE_time(i) = result.LSE_time;
end

%% Print results
fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'QP_free', sum(QPfree), sum(QPfree)/repeat, sum(QPfree_time)/sum(QPfree));
fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'SQP', sum(SQP), sum(SQP)/repeat, sum(SQP_time)/sum(SQP));
fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'ALM', sum(ALM), sum(ALM)/repeat, sum(ALM_time)/sum(ALM));
fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'LQH', sum(LQH), sum(LQH)/repeat, sum(LQH_time)/sum(LQH));
fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'LSE', sum(LSE), sum(LSE)/repeat, sum(LSE_time)/sum(LSE));