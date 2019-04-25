function sol = spgmr_inner(K, rhs, maxit)
    x_hat = ldiv(K.A, rhs(1:K.n));
    
    result = spmr_sc(K, rhs(K.n + (1:K.m)) - K.G2 * x_hat, ...
                     'tol', 1e-10, 'maxit', maxit);
    
    sol = [result.x + x_hat; result.y];
end