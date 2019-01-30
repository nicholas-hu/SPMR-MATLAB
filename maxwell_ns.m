% Maxwell example using -NS family methods

load('matrices_homogeneous3.mat');

G1 = B;
G2 = B;

indexat = @(expr, index) expr(index);

[L1, U1, P1] = lu([eye(n) G1'; G1 zeros(m)]);
[L2, U2, P2] = lu([eye(n) G2'; G2 zeros(m)]);

H1_proj = @(x) indexat(U1 \ (L1 \ (P1' * [x; zeros(m, 1)])), 1:n);
H2_proj = @(x) indexat(U2 \ (L2 \ (P2' * [x; zeros(m, 1)])), 1:n);

H1 = transfunc_wrapper(H1_proj, H1_proj);
H2 = transfunc_wrapper(H2_proj, H2_proj);

K = spmr_ns_matrix(A, H1, H2, n, m, n);

f = linspace(-1, 1, n)';

tic;
result = spqmr_ns(K, f, 'tol', 1e-10, 'maxit', 2*m);
toc;