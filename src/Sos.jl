module Sos
export issos
using Polyopt

function issos{T<:Number}(p::Poly{T},eps::Float64)
    
    order = int(p.deg / 2); 
    prob = momentprob(order,p)
    xx, yy, objval, solsta = solve_mosek(prob);
    abs(objval) < eps
    
end


function parse_dual(xx::Array{Any})

    offset = 1
    s = Array(Any, length(pineq))
    for i = 1:length(pineq)
        p = pineq[i];
        v = basis(order,x,p)
        s[i] = v'*xx[i+offset]*v
    end

    offset = length(s)
    lambda = Array(Any, length(peq))
    for i = 1:length(peq)
        p = pineq[i];
        v = basis(order,x,p)
        lambda[i] = v'*xx[i+offset]*v
    end

end



end

