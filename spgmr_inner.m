function sol = spgmr_inner(K, rhs, maxit)
    x_hat   = ldiv(K.A, rhs(1:K.n));
    g       = rhs(K.n + (1:K.m)) - K.G2 * x_hat;
    
    % The tolerance must be scaled as it is *relative*
    tol     = norm(rhs) / norm(g) * 1e-10;
    
    result  = spmr_sc(K, g, 'tol', tol, 'maxit', maxit);
    
    sol     = [result.x + x_hat; result.y];
end