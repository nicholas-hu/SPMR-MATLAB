function K = spmr_sc_matrix(A, G1, G2, n, m)

% SPMR_SC_MATRIX   Construct a saddle-point matrix for use with Schur
%                  complement-based methods.
%    K = SPMR_SC_MATRIX(A, G1, G2, n, m) 'constructs' the saddle-point 
%    matrix
%       [A G1'; G2 zeros(m, m)],
%    where A, G1, and G2 are matrices or function handles; A is n-by-n; 
%    and G1 and G2 are m-by-n.
%   
%    When a function handle f is given in place of A, f(x, 1) 
%    should return A \ x, and f(x, 2) should return A' \ x.
%   
%    When a function handle f is given in place of G1 (resp., G2), f(x, 1) 
%    should return G1 * x (resp., G2 * x), and f(x, 2) should return 
%    G1' * x (resp., G2' * x).
%    
%    See also: TRANSFUNC_WRAPPER.

    if ~(isa(A, 'numeric') || isa(A, 'function_handle'))
        error('A must be numeric or a function handle')
    elseif ~(isa(G1, 'numeric') || isa(G1, 'function_handle'))
        error('G1 must be numeric or a function handle')
    elseif ~(isa(G2, 'numeric') || isa(G2, 'function_handle'))
        error('G2 must be numeric or a function handle')
    end
    
    if isa(A, 'numeric') && ~isequal(size(A), [n, n])
        error('A should have size (%d, %d)', n, n);
    elseif isa(A, 'function_handle')
        if ~isequal(size(A(zeros(n, 1), 1)), [n, 1])
            error('A should have size (%d, %d)', n, n);
        elseif ~isequal(size(A(zeros(n, 1), 2)), [n, 1])
            error('A'' should have size (%d, %d)', n, n);
        end
    end
    
    if isa(G1, 'numeric') && ~isequal(size(G1), [m, n])
        error('G1 should have size (%d, %d)', m, n);
    elseif isa(G1, 'function_handle')
        if ~isequal(size(G1(zeros(n, 1), 1)), [m, 1])
            error('G1 should have size (%d, %d)', m, n);
        elseif ~isequal(size(G1(zeros(m, 1), 2)), [n, 1])
            error('G1'' should have size (%d, %d)', n, m);
        end
    end
    
    if isa(G2, 'numeric') && ~isequal(size(G2), [m, n])
        error('G2 should have size (%d, %d)', m, n);
    elseif isa(G2, 'function_handle')
        if ~isequal(size(G2(zeros(n, 1), 1)), [m, 1])
            error('G2 should have size (%d, %d)', m, n);
        elseif ~isequal(size(G2(zeros(m, 1), 2)), [n, 1])
            error('G2'' should have size (%d, %d)', n, m);
        end
    end

    K.A     = A;
    K.G1    = G1;
    K.G2    = G2;
    K.n     = n;
    K.m     = m;
end