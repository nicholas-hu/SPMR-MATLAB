function y = mul(A, x)
    if isa(A, 'numeric')
        y = A * x;
    else
        y = A(x, 1);
    end
end