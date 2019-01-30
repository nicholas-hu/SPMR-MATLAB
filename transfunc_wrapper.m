function f = transfunc_wrapper(f1, f2)
    function y = g(x, t)
        if t == 1
            y = f1(x);
        elseif t == 2
            y = f2(x);
        end
    end

    f = @g;
end