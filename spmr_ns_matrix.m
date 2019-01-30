function K = spmr_ns_matrix(A, H1, H2, n, m, l)
    if ~(isa(A, 'numeric') || isa(A, 'function_handle'))
        error('A must be numeric or a function handle')
    elseif ~(isa(H1, 'numeric') || isa(H1, 'function_handle'))
        error('H1 must be numeric or a function handle')
    elseif ~(isa(H2, 'numeric') || isa(H2, 'function_handle'))
        error('H2 must be numeric or a function handle')
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
    
    if isa(H1, 'numeric') && ~isequal(size(H1), [n, l])
        error('H1 should have size (%d, %d)', n, l);
    elseif isa(H1, 'function_handle')
        if ~isequal(size(H1(zeros(l, 1), 1)), [n, 1])
            error('H1 should have size (%d, %d)', n, l);
        elseif ~isequal(size(H1(zeros(n, 1), 2)), [l, 1])
            error('H1'' should have size (%d, %d)', l, n);
        end
    end
    
    if isa(H2, 'numeric') && ~isequal(size(H2), [n, l])
        error('H2 should have size (%d, %d)', n, l);
    elseif isa(H2, 'function_handle')
        if ~isequal(size(H2(zeros(l, 1), 1)), [n, 1])
            error('H2 should have size (%d, %d)', n, l);
        elseif ~isequal(size(H2(zeros(n, 1), 2)), [l, 1])
            error('H2'' should have size (%d, %d)', l, n);
        end
    end

    K.A     = A;
    K.H1    = H1;
    K.H2    = H2;
    K.n     = n;
    K.m     = m;
    K.l     = l;
end