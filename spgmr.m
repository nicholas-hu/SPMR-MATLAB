% SPGMR

clear variables; format long;

load('test_matrices/navier_stokes_N16.mat');

% Create blocks

n = size_n;
m = size_m;

% (A is already in .mat file)
G1 = B;

tic;

% Create saddle-point matrix and RHS vector

K = [A G1'; G1 sparse(m, m)];

% Create preconditioner

setup.type = 'nofill';
[L1, U1] = ilu(A, setup);
[L2, U2] = ilu(A', setup);
A_func = transfunc_wrapper(@(x) U1 \ (L1 \ x), @(x) U2 \ (L2 \ x));

P = spmr_sc_matrix(A_func, G1, G1, n, m);
M_func = @(x) spgmr_inner(P, x, n+m);

% Run iteration

[x, flag, relres, iter, resvec] = gmres(K, [f; g], ...
                                        50, 1e-10, n+m, M_func);

toc;