function result = spmr_ns(varargin)

% SPMR_NS   Saddle-Point Minimal Residual, Null Space-based Method.
%    result = SPMR_NS(K, f) solves the saddle-point system 
%       K * [x; y] = [f; zeros(m, 1)],
%    where K is a saddle-point matrix represented by a struct constructed 
%    using spmr_ns_matrix.
%
%    Optional parameters:
%       tol:     the relative residual tolerance (default: 1e-6),
%       maxit:   the maximum number of iterations to perform (default: 10),
%       precond: a symmetric positive-definite preconditioner, given as a
%                matrix or function handle (default: [no preconditioner]).
%
%    result is a struct with fields x, flag, iter, and resvec, where
%       flag is one of
%          CONVERGED:      the relative residual or estimate thereof fell 
%                          below the prescribed tolerance,
%          MAXIT_EXCEEDED: the maximum number of iterations was performed,
%          OTHER:          some computed quantity became too small,
%       iter is the number of iterations performed, and
%       resvec is the vector of relative residuals or estimates thereof.
%    
%    Note that y must be recovered separately!
%
%    See also SPMR_NS_MATRIX.

    % Parse arguments

    inp     = inputParser;

    inp.addRequired('K');
    inp.addRequired('f');
    inp.addParameter('tol', 1e-6);
    inp.addParameter('maxit', 10);
    inp.addParameter('precond', []);

    inp.parse(varargin{:});
    args    = inp.Results;

    % Extract arguments

    K       = args.K;
    f       = args.f;
    tol     = args.tol;
    maxit   = args.maxit;
    M       = args.precond;

    % Initialize SIMBA-NS variables

    H1t_f   = multrans(K.H1, -f);

    v_hat   = ldivpc(M, H1t_f);
    beta    = sqrt(dot(v_hat, H1t_f));
    v       = H1t_f / beta;
    v_hat   = v_hat / beta;

    delta   = beta;
    z       = v;
    z_hat   = v_hat;

    u       = mul(K.H2, v_hat);
    w       = mul(K.H1, z_hat);
    u_hat   = mul(K.A, u);
    w_hat   = multrans(K.A, w);

    xi      = dot(u_hat, w);
    alpha   = sqrt(abs(xi));
    gamma   = alpha;

    u       = sign(xi) / alpha * u;
    w       = sign(xi) / gamma * w;

    % Initialize QR variables

    rho_bar = gamma;
    phi_bar = delta;

    % Initialize variables for x updates

    p       = zeros(K.n, 1);
    d       = u;

    % Initialize variables for residual updates

    result.resvec   = zeros(min(K.n - K.m, maxit), 1);

    % Finish setting up

    result.iter = 0;

    if abs(xi) < eps
        result.x        = -p;
        result.flag     = SpmrFlag.OTHER;
        result.resvec   = [];
        return;
    end

    result.flag = SpmrFlag.MAXIT_EXCEEDED;

    % SPMR-NS iteration

    while result.iter < K.m
        if result.iter > maxit
            result.iter     = maxit;
            break;
        elseif abs(xi) < eps
            result.flag     = SpmrFlag.OTHER;
            result.iter     = result.iter - 1;
            break;
        end

        % SIMBA-NS iteration

        u_hat   = sign(xi) / alpha * u_hat;
        w_hat   = sign(xi) / gamma * w_hat;

        v       = multrans(K.H2, w_hat) - alpha * v;
        z       = multrans(K.H1, u_hat) - gamma * z;
        v_hat   = ldivpc(M, v);
        z_hat   = ldivpc(M, z);

        beta    = sqrt(dot(v_hat, v));
        delta   = sqrt(dot(z_hat, z));

        v       = v / beta;
        z       = z / delta;
        v_hat   = v_hat / beta;
        z_hat   = z_hat / delta;

        u       = mul(K.H2, v_hat) - sign(xi) * beta * u;
        w       = mul(K.H1, z_hat) - sign(xi) * delta * w;
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

        result.iter = result.iter + 1;

        if result.iter == 1
            result.resvec(result.iter)  = s;
        else
            result.resvec(result.iter)  = result.resvec(result.iter - 1) * s;
        end

        if result.resvec(result.iter) < tol
            result.flag     = SpmrFlag.CONVERGED;
            break;
        end
    end

    result.x        = -p;
    result.resvec   = result.resvec(1:result.iter);
end
