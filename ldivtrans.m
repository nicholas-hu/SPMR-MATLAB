function b = ldivtrans(A, x)
    if isa(A, 'numeric')
        b = A' \ x;
    else
        b = A(x, 2);
    end
end