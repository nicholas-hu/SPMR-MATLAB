% SPGMR

clear variables; format long;

% Create blocks

n = 2000;
m = 1000;

A = sparse(gallery('grcar', n));
F = 100 * speye(m, m);
G1 = [F F];

tic;

% Create saddle-point matrix and RHS vector

K = [A G1'; G1 sparse(m, m)];
g = ones(m, 1);

% Create preconditioner

setup.type = 'nofill';
[L1, U1] = ilu(A, setup);
[L2, U2] = ilu(A', setup);
A_func = transfunc_wrapper(@(x) U1 \ (L1 \ x), @(x) U2 \ (L2 \ x));

P = spmr_sc_matrix(A_func, G1, G1, n, m);
M_func = @(x) spgmr_inner(P, x, 10);

% Run iteration

[x, flag, relres, iter, resvec] = gmres(K, [zeros(n, 1); g], ...
                                        1, 1e-10, n, M_func);

toc;