% Grcar example using -SC family methods

n = 1000;
m = 500;

A = gallery('grcar', n);
F = 100 * speye(m, m);
G1 = [F F];

[L1, U1, P1] = lu(A);
[L2, U2, P2] = lu(A');

f = transfunc_wrapper(@(x) U1 \ (L1 \ (P1' * x)), ...
                      @(x) U2 \ (L2 \ (P2' * x)));

K = spmr_sc_matrix(f, G1, G1, n, m);

g = ones(m, 1);

tic;
result = spqmr_sc(K, g, 'tol', 1e-10, 'maxit', 2*m);
toc;