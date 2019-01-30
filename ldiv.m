function b = ldiv(A, x)
    if isa(A, 'numeric')
        b = A \ x;
    else
        b = A(x, 1);
    end
end