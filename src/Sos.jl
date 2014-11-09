module Sos
export issos
using Polyopt

eps = 10e-9

function issos{T<:Number}(p::Poly{T})
    issos(p,eps)
end

function issos{T<:Number}(p::Poly{T},eps::Float64)
    
    order = int(p.deg / 2); 
    prob = momentprob(order,p)
    xx, yy, objval, solsta = solve_mosek(prob);
    gramMatrix = xx[1]

    cnstTerm = evalpoly(p,vec(zeros(1,p.n)));
    lowerBndSos = cnstTerm+objval
    gramMatrix.data[1,1] = cnstTerm

    (lowerBndSos > -eps, gramMatrix, Polyopt.basis(order,variables(p.syms)))
    
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

