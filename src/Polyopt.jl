module Polyopt

export MomentProb, momentprob, momentprob_chordalembedding, solve_mosek, monomials

import Base.^

include("polynomial.jl")
include("cliques.jl")

immutable MomentProb{T<:Number}
    order :: Int
    basis :: Vector{Poly{Int}}
    obj   :: AbstractArray{T}
    mom   :: Array{Any,1}
    eq    :: Array{Any,1}
end

include("solver_mosek.jl")
include("latex.jl")

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

function indexmap{T<:Number}(v::Array{Poly{T},1})
    Dict( [(v[k].alpha, k) for k=1:length(v)] )
end

linear_index(imap::Dict{Array{Int,2},Int}, p::Poly) = Int[ imap[p.alpha[i,:]] for i=1:p.m ]

function vectorize{T<:Number}(p::Polyopt.Poly{T}, imap::Dict{Array{Int,2},Int})
    sparsevec(linear_index(imap,p), p.c, length(imap))
end

function inverse_indexmap{T<:Number}(v::Array{Poly{T},1}, imap::Dict{Array{Int,2},Int})

    lu = binomial(v[end].n + round(Int,v[end].deg/2), v[end].n)
    u  = v[1:lu]
    y = Array(Array{Tuple{Int64,Int64},1},length(v))
    for j=1:length(y) y[j] = [] end

    for j=1:lu
        for i=j:lu
            push!(y[ imap[(u[i]*u[j]).alpha] ], (i,j) )
        end
    end
    y
end

function vectorize{T<:Number}(A::AbstractArray{Poly{T}}, imap::Dict{Array{Int,2},Int})
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
            I = linear_index(imap,A[i,j])
            v = A[i,j].c
            push!(subi, (i+(j-1)*m)*ones(Int, length(v))...)
            push!(subj, I...)
            push!(val, v...)
        end
    end
    sparse(subi, subj, val, m*n, length(imap))
end

function vectorize{T<:Number}(A::Symmetric{Poly{T}}, imap::Dict{Array{Int,2},Int})
    m, n = size(A)
    subi = Int[]
    subj = Int[]
    val  = T[]
    for j=1:m
        v = A[j,j].c
        push!(subi, (j+(j-1)*m)*ones(Int, length(v))...)
        push!(subj, linear_index(imap,M[j,j])...)
        push!(val, v...)
        for i=j+1:m
            I = linear_index(imap,M[i,j])
            v = A[i,j].c
            push!(subi, (i+(j-1)*m)*ones(Int, length(v))...)
            push!(subj, I...)
            push!(val, v...)
        end
    end
    sparse(subi, subj, val, m*m, length(imap))
end

function momentprob{S,T,U}(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}, peq::Array{Poly{U},1})

    v = monomials(2*order, variables(obj.syms))
    imap = indexmap(v)

    obj.deg <= 2*order || error("obj has degree higher than 2*order")

    p = vectorize(obj, imap)
    mom = Array(Any, length(pineq)+1)
    mom[1] = vectorize(moment(order, obj.syms), imap)
    for k=1:length(pineq)
       pineq[k].deg <= 2*order || error("pineq[$(k)] has degree higher than 2*order")
       mom[k+1] = vectorize(moment(order, pineq[k].syms, pineq[k]), imap)
    end

    momeq = Array(Any, length(peq))
    for k=1:length(peq)
        peq[k].deg <= 2*order || error("peq[$(k)] has degree higher than 2*order")
        momeq[k] = vectorize(moment(order, peq[k].syms, peq[k]), imap)
    end

    MomentProb(order, v, p, mom, momeq)
end

momentprob{S}(order::Int, obj::Poly{S}) =
    momentprob(order, obj, Poly{Int}[], Poly{Int}[])

momentprob{S,T}(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}) =
    momentprob(order, obj, pineq, Poly{Int}[])

function momentprob_chordal{S,T,U}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S},
                                       pineq::Array{Poly{T},1}, pineq_index::Array{Int,1},
                                       peq::Array{Poly{U},1}, peq_index::Array{Int,1})

    v = monomials(2*order, variables(obj.syms))
    imap = indexmap(v)

    obj.deg <= 2*order || error("obj has degree higher than 2*order")

    p = vectorize(obj, imap)
    mom = Array(Any, length(pineq) + length(cliques))
    for k=1:length(cliques)
        mom[k] = vectorize(moment(order, obj.syms, cliques[k]), imap)
    end

    for k=1:length(pineq)
        pineq[k].deg <= 2*order || error("pineq[$(k)] has degree higher than 2*order")
        mom[length(cliques)+k] = vectorize(moment(order, pineq[k].syms, pineq[k], cliques[pineq_index[k]]), imap)
    end

    momeq = Array(Any, length(peq))
    for k=1:length(peq)
        peq[k].deg <= 2*order || error("peq[$(k)] has degree higher than 2*order")
        momeq[k] = vectorize(moment(order, peq[k].syms, peq[k], cliques[peq_index[k]]), imap)
    end

    MomentProb(order, v,  p, mom, momeq)
end

momentprob_chordal{S,T}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S}, pineq::Array{Poly{T},1}, pineq_index::Array{Int,1}) = 
    momentprob_chordal(order, cliques, obj, pineq, pineq_index, Poly{Int}[], Int[])

momentprob_chordal{S}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S}) = 
    momentprob_chordal(order, cliques, obj, Poly{Int}[], Int[], Poly{Int}[], Int[])
                                  
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

    momentprob_chordal(order, cliques, obj,
                       pineq, [ clique_index(cliques, find(sum(p.alpha,1))) for p = pineq ],
                       peq,   [ clique_index(cliques, find(sum(p.alpha,1))) for p = peq   ])
end

momentprob_chordalembedding{S,T}(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}) =
    momentprob_chordalembedding(order, obj, pineq, Poly{Int}[])

momentprob_chordalembedding{S}(order::Int, obj::Poly{S}) =
    momentprob_chordalembedding(order, obj, Poly{Int}[], Poly{Int}[])

end

