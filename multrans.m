function y = multrans(A, x)
    if isa(A, 'numeric')
        y = A' * x;
    else
        y = A(x, 2);
    end
end