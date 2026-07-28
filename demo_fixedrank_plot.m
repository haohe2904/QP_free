%%
%-------------------------- Fixed rank with nonnegativity and equality ----
%-------------------------- Plot figure -----------------------
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

randstate = 11111;
randn('state',double(randstate));
rand('state',double(randstate));

specifier.matlabversion = 0; % 0 if older than 2015 1 otherwise

repeat = 1;
rdim = 10;  % [6,8,10]
cdim = 15;
tolKKTres = 16; % log10 scale, i.e., 1e-* tolerance
rankval =0.4*rdim;
eqconstratio = 0.25;  % [0.5]
maskratio = 0.25;  % [0.5]
%_______Set up data______

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
options.maxOuterIter = 1000;  % for Riemannian methods
options.maxiter = 400;  % for RSQO (RSQP), QP-free
options.maxtime = 600;  % 600
options.verbosity = 2;  % 1
options.tolKKTres = 10^(-tolKKTres);
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


result = client_fixedrank_plot(rdim, cdim, rankval, P, A, eqindices, ineqindices, options,setting);

Xfinal = result.XQPfree{end}; Dist = zeros(length(result.XQPfree),1);
for ii = 1:length(result.XQPfree)
    Dist(ii) = norm(result.XQPfree{ii}-Xfinal,'fro');
end

%% Plot the figures
figure(1)
set(gcf,'position',[50 50 1000 300])
subplot(1,2,1);     % subplot(x,y,n)
plot(0:length(result.KQPfree)-1,result.KQPfree, 'r-','linewidth',2);
xlabel('Iter');
ylabel('KKT residual');
axis([0 16  1e-12 100]);
set(gca, 'YScale', 'log');
grid on   
grid minor

subplot(1,2,2);     % subplot(x,y,n)
plot(0:length(result.KQPfree)-1,Dist, 'r-','linewidth',2);
xlabel('Iter');
ylabel('$\|X^k-X^*\|_F$','Interpreter', 'latex');
axis([0 16  1e-10 100]);
set(gca, 'YScale', 'log');
grid on   
grid minor

figure(2)
set(gcf,'position',[50 50 1000 300])
subplot(1,2,1);     % subplot(x,y,n)
plot(0:length(result.KQPfree)-1,result.KQPfree, 'r-','linewidth',2);
hold on
plot(0:length(result.KALM)-1,result.KALM, 'b-','linewidth',2);
hold on
plot(0:length(result.KSQP)-1,result.KSQP, 'g-','linewidth',2);
hold on
plot(0:length(result.KLQH)-1,result.KLQH, 'k-','linewidth',2);
hold on
plot(0:length(result.KLSE)-1,result.KLSE, 'm-','linewidth',2);
xlabel('Iter');
ylabel('KKT residual');
axis([0 100  1e-15 100]);
set(gca, 'YScale', 'log');
grid on   
grid minor
legend('RQO-free','RALM','RSQO','REPM-LQH', 'REPM-LSE');

subplot(1,2,2);     % subplot(x,y,n)
plot(result.TQPfree,result.KQPfree, 'r-','linewidth',2);
hold on
plot(result.TALM,result.KALM, 'b-','linewidth',2);
hold on
plot(result.TSQP,result.KSQP, 'g-','linewidth',2);
hold on
plot(result.TLQH,result.KLQH, 'k-','linewidth',2);
hold on
plot(result.TLSE,result.KLSE, 'm-','linewidth',2);
axis([0 35  1e-15 100]);
xlabel('Time');
ylabel('KKT residual');
set(gca, 'YScale', 'log');
grid on   
grid minor