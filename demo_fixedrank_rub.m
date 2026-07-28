%%
%-------------------------- Fixed rank with nonnegativity and equality ----
%-------------------------- Measurement of Robustness -----------------------
% This code comes from Obara et al. https://github.com/shirokumakur0/Sequential-quadratic-programming-on-manifold.
% We made some minor changes.
close all; clear all; clc;

restoredefaultpath;
addpath(genpath('clients'));
addpath(genpath('connectors'));
addpath(genpath('fix_rank_functions'));
addpath(genpath('manopt'));
addpath(genpath('solvers'));
addpath(genpath('supportingfunctions'));

tolKKTres = 7;

for rdim = [10, 20]
    for cdim = [rdim 2*rdim]
        for rankval = [0.3*rdim, 0.4*rdim]
            eqconstratio = 0.25;  % [0.5]
            maskratio = 0.25;  % [0.5]
            %_______Set up data______
            repeat = 20;
            QPfree = zeros(1,repeat); QPfree_time = zeros(1,repeat);
            SQP = zeros(1,repeat); SQP_time = zeros(1,repeat);
            ALM = zeros(1,repeat); ALM_time = zeros(1,repeat);
            LQH = zeros(1,repeat); LQH_time = zeros(1,repeat);
            LSE = zeros(1,repeat); LSE_time = zeros(1,repeat);
            
            for i = 1:repeat
                
                randstate = 11111*i+10000;%
                randn('state',double(randstate));
                rand('state',double(randstate));

%                  randstate = 1111111*i+10000;%尝试跑
%                  randn('state',double(randstate));
%                  rand('state',double(randstate));
                %______Set Object matrix_____
                % Generate a random mxn matrix A of rank k
                while true
                    L = rand(rdim, rankval);
                    R = rand(cdim, rankval);
                    A =  L*R';
                    rankA = rank(A);
                    if rankA == rankval
                        break
                    end
                end
                % Generate a random mask for observed entries: P(i, j) = 1 if the entry
                % (i, j) of A is observed, and 0 otherwise.
                initP = zeros(rdim, cdim);
                num_N = ceil(0.8*rdim*cdim); num_J = ceil(0.25*num_N);
                index_N = randperm(rdim*cdim, num_N);
                index_J = index_N(1:num_J);
                ineqindices = index_N(num_J+1:end)';
                for obs = index_J
                    initP(obs) = 1;
                end
                
                %%% For nonnegativity inequality and equality
                nonzero_num = nnz(initP);
                eqnum = ceil(nonzero_num * eqconstratio);
                nonzeroidcs = find(initP);
                s = RandStream('mlfg6331_64');
                eqindices = sort(randsample(s,nonzeroidcs,eqnum));
                for eqindex = eqindices'
                    initP(eqindex) = 0;
                end
                P = initP;
                
                %________Experiment_____
                options.maxOuterIter = 5000;  % for Riemannian methods
                options.maxiter = 2000;  % for RSQO (RSQP), QP-free
                options.maxtime = 600;  % 600
                options.verbosity = 0;  % 1
                options.tolKKTres = 10^(-tolKKTres);%5*
                options.startingtolgradnorm = max(1e-3, 10^(-tolKKTres + 3));
                options.endingtolgradnorm = 10^(-16); %for RALM, REPMs.
                options.outerverbosity = options.verbosity;
                options.mineigval_correction = 1e-5;  % 1e-5
                
                %________for initial point_____
                setting.initialpoint =  'feasible_region';  % other cndidates: 'eye', 'random'
                
                %________Setting________
                setting.repeat = repeat;
                setting.row_dim = rdim;
                setting.col_dim = cdim;
                setting.rank = rankval;
                setting.mineigval_correction = options.mineigval_correction;
                setting.tolKKTres =  tolKKTres;
                setting.maxOuterIter = options.maxOuterIter;
                setting.maxtime = options.maxtime;
                setting.verbosity = options.verbosity;
                setting.eqconstratio = eqconstratio;
                setting.maskratio = maskratio;
                setting.filepath = sprintf('nrep%dRowdim%dColdim%dRank%dTol%dEqratio%0.1eMaskratio%0.1e',...
                    setting.repeat, setting.row_dim, setting.col_dim, setting.rank, setting.tolKKTres,setting.eqconstratio, setting.maskratio);
                setting.P = P;
                setting.A = A;
                
                
                result = client_fixedrank_rub(rdim, cdim, rankval, P, A, eqindices, ineqindices, options, setting);
                
                %% experimental results
                QPfree(i) = result.QPfree; QPfree_time(i) = result.QPfree_time;
                SQP(i) = result.SQP; SQP_time(i) = result.SQP_time;
                ALM(i) = result.ALM; ALM_time(i) = result.ALM_time;
                LQH(i) = result.LQH; LQH_time(i) = result.LQH_time;
                LSE(i) = result.LSE; LSE_time(i) = result.LSE_time;
            end
            
            %% Print results
            fprintf('rdim: %d, cdim: %d, rank: %d \n', rdim, cdim, rankval);
            fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'QP_free', sum(QPfree), sum(QPfree)/repeat, sum(QPfree_time)/sum(QPfree));
            fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'SQP', sum(SQP), sum(SQP)/repeat, sum(SQP_time)/sum(SQP));
            fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'ALM', sum(ALM), sum(ALM)/repeat, sum(ALM_time)/sum(ALM));
            fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'LQH', sum(LQH), sum(LQH)/repeat, sum(LQH_time)/sum(LQH));
            fprintf('Method: %s, Num_success: %d, Rate_success: %f, Ave_time: %f \n', 'LSE', sum(LSE), sum(LSE)/repeat, sum(LSE_time)/sum(LSE));
            fprintf('***************************************************************************\n');
            fprintf('***************************************************************************\n');
        end
    end
end 

