module Polyopt

export solve_mosek, monomials
export MomentProb, momentprob, momentprob_chordal, momentprob_chordalembedding 
export BSOSProb, bsosprob_chordal, bsosprob_chordal2
export solve_mosek2

import Base.^, Base.start, Base.next, Base.done, Base.sub, Base.length

include("polynomial.jl")
include("cliques.jl")

immutable MomentProb{T<:Number}
    order :: Int
    basis :: Vector{Poly{Int}}
    obj   :: AbstractArray{T}
    mom   :: Array{Any,1}
    eq    :: Array{Any,1}
end

immutable BSOSProb{T<:Number}
    degree :: Int
    order  :: Int
    obj    :: AbstractArray{T}
    Al     :: Array{Any,1}
    El     :: Array{Any,1}
    As     :: Array{Any,1}
    lb     :: Array{Any,1}
end

include("solver_mosek.jl")
include("latex.jl")


immutable MonomialPowers
    n :: Int
    d :: Int
end

monomialpowers(n, d) = MonomialPowers(n, d)
start(m::MonomialPowers) = zeros(Int, m.n)
done(m::MonomialPowers, state::Vector{Int}) = (state[1] == m.d+1)
length(m::MonomialPowers) = binomial(m.n + m.d, m.d)

function next(m::MonomialPowers, powers::Vector{Int})

    # find index of element to increment
    k = m.n
    while (k>1)
        if sum(view(powers,1:k)) < m.d
            break
        end
        k -= 1
    end     

    state = copy(powers)    
    state[k+1:end] = 0                        
    state[k] += 1
    
    powers, state
end

function _merge_monomials{T<:Number}(a::Array{Poly{T},1}, b::Array{Poly{T},1})

    n, m = length(a), length(b)
    y = Array(Poly{T}, n + m)
    
    i, j, k = 1, 1, 1
    
    while (i <= n) && (j <= m)        
        if a[i] <= b[j]
            y[k] = a[i]
            i += 1
        else
            y[k] = b[j]
            j += 1        
        end
        k += 1
    end
    
    if i <= n
        y[k:n+m] = a[i:n]
    else
        y[k:n+m] = b[j:m]
    end
    
    y
end

# compute all monomials of degree 'deg' or less
function monomials{T<:Number}(deg::Int, vars::Array{Poly{T},1})

    if length(vars) == 1        
        monoms = Array(Poly{T}, deg+1)
        monoms[1] = vars[1]^0
        for k=1:deg            
            monoms[k+1] = monoms[k] * vars[1] 
        end
    else
        x, tail = vars[1], vars[2:end]
        m = monomials(deg, tail)
        monoms = m[:]

        t = x
        for i=1:deg
            monoms = _merge_monomials(monoms, monomials(deg-i, tail)*t)
            t = t*x
        end
    end
    monoms
end

# function monomials{T<:Number}(deg::Int, vars::Array{Poly{T},1})
# 
#     # find index of element to increment
#     function findindex(a::Array{Int,2}, deg::Int)
#         k=length(a)    
#         while (k>1)
#             if sum(a[1:k]) < deg
#                 break
#             end
#             k -= 1
#         end     
#         k
#     end
#     
#     n = length(vars)
#     m = binomial(n+deg,deg)
#     
#     a = zeros(Int, 1, n)    
#     r = Array(Poly{T}, m)            
#     
#     for i=1:m
#         r[i] = Poly{T}( vars[1].syms, copy(a), [1] )
#     
#         k = findindex(a, deg)
#         a[k+1:end] = 0                        
#         a[k] += 1
#     end
#     r
# end
# 
# function monomials{T<:Number}(deg::Int, vars::Array{Poly{T},1})
#     [ Poly{T}(vars[1].syms, a', [1]) for a in monomialpowers(vars[1].n, deg) ]    
# end

function moment(order::Int, syms::Symbols)
    v = monomials(order, variables(syms))
    v*v'
end

function moment(order::Int, syms::Symbols, I::Array{Int})
    v = monomials(order, variables(syms)[I])
    v*v'
end

function moment{T<:Number}(order::Int, syms::Symbols, p::Poly{T})
    v = monomials(order - (p.deg+1) >> 1, variables(syms))
    p*v*v'
end

function moment{T<:Number}(order::Int, syms::Symbols, p::Poly{T}, I::Array{Int})
    v = monomials(order - (p.deg+1) >> 1, variables(syms)[I])
    p*v*v'
end

function basis_index(a::Vector{Int}, degree::Int)
    n = length(a)
    idx = 1        
    for j=1:n-1                
        if a[j] > 0 
            dj = degree - sum(a[1:j-1])            
            for i=0:a[j]-1            
                nj = n-j
                di = dj - i
                idx += binomial(nj+di, di)
            end
        end
    end
    idx + a[n]
end

basis_index(p::Poly, degree::Int) = basis_index(vec(p.alpha), degree)

linear_index(p::Poly, d::Int) = Int[ basis_index(vec(p.alpha[i,:]), d) for i=1:p.m ]

vectorize{T<:Number}(p::Polyopt.Poly{T}, d::Int) = sparsevec(linear_index(p, d), p.c, binomial(p.n+d,d))

function vectorize{T<:Number}(A::AbstractArray{Poly{T}}, d::Int)
    dims = length(size(A))
    if dims==1
        m, n = length(A), 1
    elseif dims==2
        m, n = size(A)
    else
        error("dimension must be <= 2")
    end
    subi = Int[]
    subj = Int[]
    val  = T[]
    for j=1:n
        for i=1:m
            I = linear_index(A[i,j], d)
            v = A[i,j].c
            push!(subi, (i+(j-1)*m)*ones(Int, length(v))...)
            push!(subj, I...)
            push!(val, v...)
        end
    end
    sparse(subi, subj, val, m*n, binomial(A[1].n+d,d))
end

function momentprob{S,T,U}(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}, peq::Array{Poly{U},1})

    v = monomials(2*order, variables(obj.syms))
    obj.deg <= 2*order || error("obj has degree higher than 2*order")

    p = vectorize(obj, 2*order)
    mom = Array(Any, length(pineq)+1)
    mom[1] = vectorize(moment(order, obj.syms), 2*order)
    for k=1:length(pineq)
       pineq[k].deg <= 2*order || error("pineq[$(k)] has degree higher than 2*order")
       mom[k+1] = vectorize(moment(order, pineq[k].syms, pineq[k]), 2*order)
    end

    momeq = Array(Any, length(peq))
    for k=1:length(peq)
        peq[k].deg <= 2*order || error("peq[$(k)] has degree higher than 2*order")
        momeq[k] = vectorize(moment(order, peq[k].syms, peq[k]), 2*order)
    end

    MomentProb(order, v, p, mom, momeq)
end

momentprob{S}(order::Int, obj::Poly{S}) =
    momentprob(order, obj, Poly{Int}[], Poly{Int}[])

momentprob{S,T}(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}) =
    momentprob(order, obj, pineq, Poly{Int}[])

function momentprob_chordal{S,T,U}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S},
                                   pineq::Array{Poly{T},1}, peq::Array{Poly{U},1})

    v = monomials(2*order, variables(obj.syms))
    obj.deg <= 2*order || error("obj has degree higher than 2*order")
    
    p = vectorize(obj, 2*order)
    mom = []
    for k=1:length(cliques)
        push!(mom,vectorize(moment(order, obj.syms, cliques[k]), 2*order))
    end
    
    for k=1:length(pineq)
        pineq[k].deg <= 2*order || error("pineq[$(k)] has degree higher than 2*order")
        for j=clique_index(cliques, find(sum(pineq[k].alpha,1)), obj.n)
            push!(mom,vectorize(moment(order, pineq[k].syms, pineq[k], cliques[j]), 2*order))
        end
    end

    momeq = []
    for k=1:length(peq)
        peq[k].deg <= 2*order || error("peq[$(k)] has degree higher than 2*order")
        for j=clique_index(cliques, find(sum(peq[k].alpha,1)), obj.n)
            push!(momeq, vectorize(moment(order, peq[k].syms, peq[k], cliques[j]), 2*order))
        end
    end

    MomentProb(order, v,  p, mom, momeq)
end

momentprob_chordal{S,T}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S}, pineq::Array{Poly{T},1}) = 
    momentprob_chordal(order, cliques, obj, pineq, Poly{Int}[])

momentprob_chordal{S}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S}) = 
    momentprob_chordal(order, cliques, obj, Poly{Int}[], Poly{Int}[])
                                  
function correlative_sparsity{S,T}(obj::Poly{S}, p::Array{Poly{T},1})
    A = eye(Int,obj.n,obj.n)

    for k=1:obj.m
        I = find(obj.alpha[k,:])
        A[I,I] = 1
    end

    for j=1:length(p)
        I = find(sum(p[j].alpha,1))
        A[I,I] = 1
    end
    sparse(A)
end

function momentprob_chordalembedding{S,T,U}(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}, peq::Array{Poly{U},1})
    C = correlative_sparsity(obj, [pineq; peq])
    cliques = chordal_embedding(C)
    momentprob_chordal(order, cliques, obj, pineq, peq)
end

momentprob_chordalembedding{S,T}(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}) =
    momentprob_chordalembedding(order, obj, pineq, Poly{Int}[])

momentprob_chordalembedding{S}(order::Int, obj::Poly{S}) =
    momentprob_chordalembedding(order, obj, Poly{Int}[], Poly{Int}[])

function symbol_restrict{T<:Number}(p::Poly{T}, syms::Symbols, I::Array{Int,1})

    if p.n == 0
        p
    else
        mask = ones(Int, p.n)
        mask[I] = 0    
        J = find(p.alpha * mask .== 0)
        Poly{T}(syms, p.alpha[J,I], p.c[J])
    end
end

function bsosprob_chordal{S,T,V}(degree::Int, order::Int, cliques::Array{Array{Int,1}}, 
                                    obj::Poly{S}, pineq::Array{Poly{T},1}, peq::Array{Poly{V},1})

    n = obj.n
    
    allp = [pineq; peq]
    dmax = max(obj.deg, 2*order, degree*maximum([gi.deg for gi in allp]))
    m = binomial(n+dmax,dmax)

    x = variables(obj.syms)
    f = vectorize(obj, dmax)

    Ji = Vector{Int}[]
    for j=1:length(cliques)
        push!(Ji, [])
    end

    for (i, gi) in enumerate(pineq)
        for j=clique_index(cliques, find(sum(gi.alpha,1)), n)
            append!(Ji[j], i)
        end
    end

    Je = Vector{Int}[]
    for j=1:length(cliques)
        push!(Je, [])
    end

    for (i, hi) in enumerate(peq)
        for j=clique_index(cliques, find(sum(hi.alpha,1)), n)
            append!(Je[j], i)
        end
    end
        
    As, Al, El, lb = [], [], [], []
    for (j,c) in enumerate(cliques)
        symc = Symbols(obj.syms.names[c])
        xc = Poly{Int}[symbol_restrict(x[i], symc, c) for i in c ]
        
        u = monomials(order, xc)
        M = u*u'  
        push!(As, vectorize(M, dmax))

        pj = [ symbol_restrict(pineq[i], symc, c) for i=Ji[j] ]
        
        p = vcat( pj, 
                  [ (1-pji) for pji in pj ], 
                  [ symbol_restrict(peq[i], symc, c) for i=Je[j] ])
        
        ab_d = monomialpowers(length(p), degree)
        
        lbj = zeros(length(ab_d))
        ai, aj, av = Int[], Int[], Float64[]
        for (i, ab) in enumerate(ab_d)
            h_ab = prod(p .^ ab)
            #println("h($(ab)): ", h_ab)
            ak = vectorize(h_ab, dmax)

            push!(ai, i*ones(ak.nzind)...)
            push!(aj, ak.nzind...)
            push!(av, ak.nzval...)         
            
            lbj[i] = (sum(ab[2*length(pj)+1:end]) > 0 ? -Inf : 0.0)
        end            
        push!(Al, sparse(ai, aj, av, length(ab_d), binomial(length(c)+dmax,dmax)))  
        push!(lb, lbj)
           
        v = monomials(dmax, x[c])    
        k = Int[ basis_index(vi, dmax) for vi=v ]
        push!(El, SparseMatrixCSC{Int,Int}(m,length(k),collect(1:length(k)+1),k,ones(Int,length(k))))
                       
    end
    
    BSOSProb(degree, order, f, Al, El, As, lb)
end

bsosprob_chordal{S,T}(degree::Int, order::Int, cliques::Array{Array{Int,1}}, obj::Poly{S}, pineq::Array{Poly{T},1}) =
    bsosprob_chordal(degree, order, cliques, obj, pineq, Poly{Int}[])

end
