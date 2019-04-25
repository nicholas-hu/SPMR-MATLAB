% Grcar example using -NS family methods

clear variables; format long;

n = 1000;
m = 500;

A = sparse(gallery('grcar', n));
F = 100 * speye(m, m);
G1 = [F F];

indexat = @(expr, index) expr(index);

[L1, U1, P1] = lu([eye(n) G1'; G1 zeros(m)]);
H1_proj = @(x) indexat(U1 \ (L1 \ (P1 * [x; zeros(m, 1)])), 1:n);
H1 = transfunc_wrapper(H1_proj, H1_proj);

K = spmr_ns_matrix(A, H1, H1, n, m, n);

f = linspace(-1, 1, n)';

tic;
result = spmr_ns(K, f, 'tol', 1e-10, 'maxit', 2*m);
toc;