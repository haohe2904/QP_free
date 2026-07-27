
addpath(genpath('connectors'));


mex dgeqrf_mex.c -lmwlapack;
mex apply_q_mex.c -lmwlapack;