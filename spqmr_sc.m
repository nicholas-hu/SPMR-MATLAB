function result = spqmr_sc(varargin)
    % Parse arguments

    inp     = inputParser;

    inp.addRequired('K');
    inp.addRequired('g');
    inp.addParameter('tol', 1e-6);
    inp.addParameter('maxit', 10);
    inp.addParameter('precond', []);

    inp.parse(varargin{:});
    args    = inp.Results;

    % Extract arguments

    K       = args.K;
    g       = args.g;
    tol     = args.tol;
    maxit   = args.maxit;
    M       = args.precond;

    % Initialize SIMBO-SC variables

    v_hat   = ldivpc(M, g);
    z_hat   = v_hat;

    chi     = dot(g, v_hat);
    delta   = sqrt(abs(chi));
    beta    = sign(chi) * delta;

    v       = g / delta;
    z       = g / beta;
    v_hat   = v_hat / delta;
    z_hat   = z_hat / beta;

    u_hat   = multrans(K.G1, v_hat);
    w_hat   = multrans(K.G2, z_hat);
    u       = ldiv(K.A, u_hat);
    w       = ldivtrans(K.A, w_hat);

    xi      = dot(u_hat, w);
    alpha   = sqrt(abs(xi));
    gamma   = alpha;

    u       = sign(xi) / alpha * u;
    w       = sign(xi) / gamma * w;

    % Initialize QR variables

    rho_bar = gamma;
    phi_bar = delta;

    % Initialize variables for x updates

    result.x    = zeros(K.n, 1);
    d           = u;

    % Initialize variables for y updates

    result.y    = zeros(K.m, 1);

    T           = zeros(K.m, 3);
    T(:, 1)     = v;
    mu          = 0;
    nu          = 0;

    alpha_prev  = alpha;
    xi_prev     = xi;
    sigma_prev  = 0;

    % Initialize variables for residual updates

    result.resvec   = zeros(min(K.m, maxit), 1);

    relres          = 1;
    norm_g          = norm(g);

    % Finish setting up

    result.iter = 0;

    if abs(xi) < eps
        result.flag     = SpmrFlag.OTHER;
        result.resvec   = [];
        return;
    end

    result.flag = SpmrFlag.MAXIT_EXCEEDED;

    % SPQMR-SC iteration

    while result.iter < K.m
        if result.iter > maxit
            result.iter     = maxit;
            break;
        elseif abs(xi) < eps
            result.flag     = SpmrFlag.OTHER;
            result.iter     = result.iter - 1;
            break;
        end

        % SIMBO-SC iteration

        v       = mul(K.G2, u) - gamma * v;
        z       = mul(K.G1, w) - alpha * z;
        v_hat   = ldivpc(M, v);
        z_hat   = ldivpc(M, z);

        chi     = dot(z, v_hat);
        delta   = sqrt(abs(chi));
        beta    = sign(chi) * delta;

        v       = v / delta;
        z       = z / beta;
        v_hat   = v_hat / delta;
        z_hat   = z_hat / beta;

        u_hat   = multrans(K.G1, v_hat);
        w_hat   = multrans(K.G2, z_hat);

        u       = ldiv(K.A, u_hat) - sign(xi) * beta * u;
        w       = ldivtrans(K.A, w_hat) - sign(xi) * delta * w;

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

        result.x    = result.x + phi/rho * d;
        d           = u - sigma/rho * d;

        % Update y

        lambda      = sign(xi_prev) * rho * alpha_prev;
        T(:, 1)     = (T(:, 1) - mu * T(:, 2) - nu * T(:, 3)) / lambda;
        result.y    = result.y - phi * T(:, 1);

        T(:, 3)     = T(:, 2);
        T(:, 2)     = T(:, 1);
        T(:, 1)     = v;
        mu          = sign(xi_prev) * rho * beta + sign(xi) * sigma * alpha;
        nu          = sign(xi_prev) * sigma_prev * beta;

        alpha_prev  = alpha;
        xi_prev     = xi;
        sigma_prev  = sigma;

        % Update residual

        result.iter = result.iter + 1;

        relres                      = relres * s;
        result.resvec(result.iter)  = sqrt(result.iter) * relres;

        if result.resvec(result.iter) < tol
            norm_r  = norm(mul(K.G2, result.x) - g);

            if norm_r < tol * norm_g
                result.flag     = SpmrFlag.CONVERGED;
                break;
            end
        end
    end

    result.y        = ldivpc(M, result.y);
    result.resvec   = result.resvec(1:result.iter);
end
