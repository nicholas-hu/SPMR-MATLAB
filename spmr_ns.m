function result = spmr_ns(varargin)
    % Parse arguments

    inp         = inputParser;

    inp.addRequired('K');
    inp.addRequired('f');
    inp.addParameter('tol', 1e-6);
    inp.addParameter('maxit', 10);

    inp.parse(varargin{:});
    args        = inp.Results;

    % Extract arguments

    K           = args.K;
    f           = args.f;
    tol         = args.tol;
    maxit       = args.maxit;

    % Initialize SIMBA-NS variables
    
    H1t_f       = multrans(K.H1, -f);

    beta        = norm(H1t_f);
    v           = H1t_f / beta;
    delta       = beta;
    z           = v;

    u           = mul(K.H2, v);
    w           = mul(K.H1, z);
    u_hat       = mul(K.A, u);
    w_hat       = multrans(K.A, w);

    xi          = dot(u_hat, w);
    alpha       = sqrt(abs(xi));
    gamma       = alpha;

    u           = sign(xi) / alpha * u;
    w           = sign(xi) / gamma * w;
    
    % Initialize QR variables

    rho_bar     = gamma;
    phi_bar     = delta;

    % Initialize variables for x updates

    p           = zeros(K.n, 1);
    d           = u;

    % Initialize variables for residual updates

    result.resvec   = zeros(min(K.n - K.m, maxit), 1);

    % Finish setting up
    
    if abs(xi) < eps
        result.x        = -p;
        result.flag     = SpmrFlag.OTHER;
        result.resvec   = [];
        return;
    end

    result.flag = SpmrFlag.MAXIT_EXCEEDED;
    result.iter = 0;

    % SPMR-NS iteration

    while result.iter < K.m
        if result.iter > maxit
            result.x        = -p;
            result.iter     = maxit;
            result.resvec   = result.resvec(1:result.iter);
            return;
        elseif abs(xi) < eps
            result.x        = -p;
            result.flag     = SpmrFlag.OTHER;
            result.iter     = result.iter - 1;
            result.resvec   = result.resvec(1:result.iter);
            return;
        end

        % SIMBA-NS iteration

        u_hat   = sign(xi) / alpha * u_hat;
        w_hat   = sign(xi) / gamma * w_hat;
        
        v       = multrans(K.H2, w_hat) - alpha * v;
        z       = multrans(K.H1, u_hat) - gamma * z;
        
        beta    = norm(v);
        v       = normalize(v, 'norm');
        delta   = norm(z);
        z       = normalize(z, 'norm');
        
        u       = mul(K.H2, v) - sign(xi) * beta * u;
        w       = mul(K.H1, z) - sign(xi) * delta * w;
        u_hat   = mul(K.A, u);
        w_hat   = multrans(K.A, w);

        xi      = dot(u_hat, w);
        alpha   = sqrt(abs(xi));
        gamma   = alpha;

        u       = sign(xi) / alpha * u;
        w       = sign(xi) / gamma * w;

        % Update QR decomposition

        rho     = norm([rho_bar; delta]);
        c       = rho_bar / rho;
        s       = delta / rho;

        sigma   = s * gamma;
        rho_bar = -c * gamma;
        phi     = c * phi_bar;
        phi_bar = s * phi_bar;

        % Update x

        p       = p + phi/rho * d;
        d       = u - sigma/rho * d;

        % Update residual

        result.iter     = result.iter + 1;

        if result.iter == 1
            result.resvec(result.iter)  = s;
        else
            result.resvec(result.iter)  = result.resvec(result.iter - 1) * s;
        end

        if result.resvec(result.iter) < tol
            result.x        = -p;
            result.flag     = SpmrFlag.CONVERGED;
            result.resvec   = result.resvec(1:result.iter);
            return;
        end
    end
    
    result.x    = -p;
end