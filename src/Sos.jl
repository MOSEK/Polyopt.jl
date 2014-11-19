module Sos
export issos, sosprog,solve,Constraint,SosProg,parse_dual
using Polyopt: moment, indexmap, basis
using Polyopt

eps = 10e-9

immutable Constraint
    p :: Poly{Number}
    mdeg :: Int
end

immutable SosProg
    prog :: MomentProb
    basis :: Array{Poly{Number},1}
    pineq :: Array{Constraint,1}
    peq :: Array{Constraint,1}
    obj :: Poly{Number}
end

function degree(constraint::Constraint)
    constraint.p.deg + 2*constraint.mdeg
end

function computebasis(sosprog::Constraint)
    basis(constraint.mdeg,variables(constraint.p.syms))
end

function computebasis(constraint::Constraint)
    basis(constraint.mdeg,variables(constraint.p.syms))
end




function sosprog{S}(obj::Poly{S}, pineq::Array{Constraint,1}, peq::Array{Constraint,1})

    #Find the degree of the certificate
    degreeCert = obj.deg;
    for k=1:length(pineq)
       degreeCert = max(degree(pineq[k]),degreeCert)
    end
    
    for k=1:length(peq)
       degreeCert = max(degree(peq[k]),degreeCert)
    end

    degreeCert = degreeCert + mod(degreeCert,2)
   
    v = monomials(degreeCert, variables(obj.syms))
    l = length(v)
    imap = indexmap(v)

    p = vectorize(obj, l, imap)
    mom = Array(Any, length(pineq)+1)
    mom[1] = vectorize(moment(int(degreeCert/2), obj.syms), l, imap)

    for k=1:length(pineq) 
       basis = computebasis(pineq[k])
       mom[k+1] = vectorize(pineq[k].p*basis*basis', l, imap)
    end

    momeq = Array(Any, length(peq))
    for k=1:length(peq)
        basis = computebasis(peq[k])
        momeq[k] = vectorize(peq[k].p*basis*basis', l, imap)
    end

    basisCert = monomials(int(degreeCert/2), variables(obj.syms))
    SosProg(MomentProb(degreeCert, v,  p, mom, momeq),basisCert,pineq,peq,obj)

end



function issos{T<:Number}(p::Poly{T})
    issos(p,eps)
end

function issos{T<:Number}(p::Poly{T},eps::Float64)
    
    order = int(p.deg / 2); 
    prob = momentprob(order,p)
    xx, yy, objval, solsta = solve_mosek(prob);
    xx = xx[1]
    gramMatrix = xx[1]

    cnstTerm = evalpoly(p,vec(zeros(1,p.n)));
    lowerBndSos = cnstTerm+objval
    gramMatrix.data[1,1] = cnstTerm

    (lowerBndSos > -eps, gramMatrix, Polyopt.basis(order,variables(p.syms)))
    
end

function svect(x)

    n = size(x,1);
    y = Array(Any, int(n*n/2+n/2))
    cnt = 1
    for j=1:size(x,1)
        for i=j:size(x,2)
            y[cnt] = x[i,j]
            cnt = cnt + 1
        end
    end

    y
    
end


function parse_dual(x,prog::SosProg)

    xpsd = x[1]
    xfree = x[2]
    println(typeof(xpsd))
    pineq = prog.pineq
    peq = prog.peq
    
    cert = prog.basis'*xpsd[1]*prog.basis
    cert = cert[1]
    
    offset = 1
    s = Array(Any, length(pineq))
    for i = 1:length(pineq) 
        v = computebasis(pineq[i])
        temp = v'*xpsd[i+offset]*v
        s[i] = temp[1]
    end


    lambda = Array(Any, length(peq))
    for i = 1:length(peq)
        v =  basis(peq[i].mdeg,variables(peq[i].p.syms))
        v = svect(v*v')
        temp = xfree'*v
        lambda[i] = temp[1]
    end

    (cert,s,lambda)

end

function solve(sosprog::SosProg)
    x,y,optval = solve_mosek(sosprog.prog)
    cert,s,lambda = parse_dual(x,sosprog)
    (optval,cert,s,lambda)
end

end

