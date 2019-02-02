function b = ldivpc(A, x)
    if isa(A, 'numeric')
        if isempty(A)  % Indicates no preconditioning
            b = x;
        else
            b = A \ x;
        end
    else
        b = A(x);
    end
end