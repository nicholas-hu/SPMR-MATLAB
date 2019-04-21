function f = transfunc_wrapper(f1, f2)

% TRANSFUNC_WRAPPER   Wrap two single-parameter functions into a single two-parameter function.
%    f = TRANSFUNC_WRAPPER(f1, f2) returns a function handle f such that
%    f(x, 1) = f1(x) and f(x, 2) = f2(x) for all x.
%
%    This is useful for providing arguments to spmr_sc_matrix or
%    spmr_ns_matrix.
%
%    See also: spmr_sc_matrix, spmr_ns_matrix.

    function y = g(x, t)
        if t == 1
            y = f1(x);
        elseif t == 2
            y = f2(x);
        end
    end

    f = @g;
end