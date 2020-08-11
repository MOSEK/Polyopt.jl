module Polyopt

export solve_mosek, monomials
export MomentProb, momentprob, momentprob_chordal, momentprob_chordalembedding 
export BSOSProb, bsosprob_chordal

import Base.^, Base.length, Base.iterate, Base.eltype

include("polynomial.jl")
include("cliques.jl")

struct MomentProb{T<:Number}
    order :: Int
    basis :: Vector{Poly{Int}}
    obj   :: AbstractArray{T}
    mom   :: Array{Any,1}
    eq    :: Array{Any,1}
end

struct BSOSProb{T<:Number}
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


struct MonomialPowers
    n :: Int
    d :: Int
end

monomialpowers(n, d) = MonomialPowers(n, d)
start(m::MonomialPowers) = zeros(Int, m.n)
done(m::MonomialPowers, state::Vector{Int}) = (state[1] == m.d+1)
length(m::MonomialPowers) = binomial(m.n + m.d, m.d)
eltype(m::MonomialPowers) = Vector{Int}

function iterate(m::MonomialPowers, powers::Vector{Int}=zeros(Int, m.n))

    if powers[1] > m.d return nothing end
    
    # find index of element to increment
    k = m.n
    while (k>1)
        if sum(view(powers,1:k)) < m.d
            break
        end
        k -= 1
    end     

    state = copy(powers)    
    state[k+1:end] .= 0                        
    state[k] += 1
    
    powers, state
end

function monomials(deg::Int, vars::Array{Poly{T},1}) where {T<:Number}
    n = length(vars)
    I = [ findfirst(x -> x != 0, vi.alpha) for vi=vars ]
    
    v = Array{Poly{Int}}(undef, binomial(n+deg, deg))
    alpha = zeros(Int, 1, vars[1].n)
    for (i,a)=enumerate(monomialpowers(n, deg))  
        alpha[I] = a
        v[i] = Poly{Int}(vars[1].syms, copy(alpha), [1])
        alpha[I] .= 0
    end
    v        
end

function moment(order::Int, syms::Symbols)
    v = monomials(order, variables(syms))
    v*v'
end

function moment(order::Int, syms::Symbols, I::Array{Int})
    v = monomials(order, variables(syms)[I])
    v*v'
end

function moment(order::Int, syms::Symbols, p::Poly{T}) where {T<:Number}
    v = monomials(order - (p.deg+1)>>1, variables(syms))
    p*v*v'
end

function moment(order::Int, syms::Symbols, p::Poly{T}, I::Array{Int}) where {T<:Number}
    v = monomials(order - (p.deg+1)>>1, variables(syms)[I])
    p*v*v'
end

function basis_index(a::Vector{Int}, degree::Int)
    n = length(a)
    idx = 1      
    
    @inbounds begin  
    for j=1:n-1                
        if a[j] > 0 
            #dj = degree - sum(a[1:j-1])            
            dj = degree - sum(view(a,1:j-1))            
            for i=0:a[j]-1            
                nj = n-j
                di = dj - i
                idx += binomial(nj+di, di)
            end
        end
    end
    
    end
    idx + a[n]
end

basis_index(p::Poly, degree::Int) = basis_index(vec(p.alpha), degree)

linear_index(p::Poly, d::Int) = Int[ basis_index(vec(p.alpha[i,:]), d) for i=1:p.m ]

vectorize(p::Polyopt.Poly{T}, d::Int) where {T<:Number} = sparsevec(linear_index(p, d), p.c, binomial(p.n+d,d))

function find_replacement(v)
    [x for (x,y) in filter((y->(y[2] != 0)), Tuple(enumerate(v)))]
end

function vectorize(A::AbstractArray{Poly{T}}, d::Int) where {T<:Number}
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

function momentprob(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}, peq::Array{Poly{U},1}) where {S<:Number,T<:Number,U<:Number}

    v = monomials(2*order, variables(obj.syms))
    obj.deg <= 2*order || error("obj has degree higher than 2*order")

    p = vectorize(obj, 2*order)
    mom = Array{Any}(undef, length(pineq)+1)
    mom[1] = vectorize(moment(order, obj.syms), 2*order)
    for k=1:length(pineq)
       pineq[k].deg <= 2*order || error("pineq[$(k)] has degree higher than 2*order")
       mom[k+1] = vectorize(moment(order, pineq[k].syms, pineq[k]), 2*order)
    end

    momeq = Array{Any}(undef, length(peq))
    for k=1:length(peq)
        peq[k].deg <= 2*order || error("peq[$(k)] has degree higher than 2*order")
        momeq[k] = vectorize(moment(order, peq[k].syms, peq[k]), 2*order)
    end

    MomentProb(order, v, p, mom, momeq)
end

momentprob(order::Int, obj::Poly{S}) where {S<:Number} =
    momentprob(order, obj, Poly{Int}[], Poly{Int}[])

momentprob(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}) where {S<:Number,T<:Number}  =
    momentprob(order, obj, pineq, Poly{Int}[])

function momentprob_chordal(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S},
                                   pineq::Array{Poly{T},1}, peq::Array{Poly{U},1}) where {S<:Number,T<:Number,U<:Number}

    v = monomials(2*order, variables(obj.syms))
    obj.deg <= 2*order || error("obj has degree higher than 2*order")
    
    p = vectorize(obj, 2*order)
    mom = []
    for k=1:length(cliques)
        push!(mom,vectorize(moment(order, obj.syms, cliques[k]), 2*order))
    end
    
    for k=1:length(pineq)
        pineq[k].deg <= 2*order || error("pineq[$(k)] has degree higher than 2*order")
        for j=clique_index(cliques, find_replacement(sum(pineq[k].alpha,dims=1)), obj.n)
            push!(mom,vectorize(moment(order, pineq[k].syms, pineq[k], cliques[j]), 2*order))
        end
    end

    momeq = []
    for k=1:length(peq)
        peq[k].deg <= 2*order || error("peq[$(k)] has degree higher than 2*order")
        for j=clique_index(cliques, find_replacement(sum(peq[k].alpha,dims=1)), obj.n)
            push!(momeq, vectorize(moment(order, peq[k].syms, peq[k], cliques[j]), 2*order))
        end
    end

    MomentProb(order, v,  p, mom, momeq)
end

momentprob_chordal(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S}, pineq::Array{Poly{T},1}) where {S<:Number,T<:Number} = 
    momentprob_chordal(order, cliques, obj, pineq, Poly{Int}[])

momentprob_chordal(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S}) where {S<:Number} = 
    momentprob_chordal(order, cliques, obj, Poly{Int}[], Poly{Int}[])
                                  
function correlative_sparsity(obj::Poly{S}, p::Array{Poly{T},1}) where {S<:Number,T<:Number}
    A = Matrix{Int}(LinearAlgebra.I,obj.n,obj.n)

    for k=1:obj.m
        I = find_replacement(obj.alpha[k,:])
        A[I,I] .= 1
    end

    for j=1:length(p)
        I = find_replacement(sum(p[j].alpha,dims=1))
        A[I,I] .= 1
    end
    sparse(A)
end

function momentprob_chordalembedding(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}, peq::Array{Poly{U},1}) where {S<:Number,T<:Number,U<:Number}
    C = correlative_sparsity(obj, [pineq; peq])
    cliques = chordal_embedding(C)
    momentprob_chordal(order, cliques, obj, pineq, peq)
end

momentprob_chordalembedding(order::Int, obj::Poly{S}, pineq::Array{Poly{T},1}) where {S<:Number,T<:Number} =
    momentprob_chordalembedding(order, obj, pineq, Poly{Int}[])

momentprob_chordalembedding(order::Int, obj::Poly{S}) where {S<:Number} =
    momentprob_chordalembedding(order, obj, Poly{Int}[], Poly{Int}[])

function symbol_restrict(p::Poly{T}, syms::Symbols, I::Array{Int,1}) where {T<:Number}

    if p.n == 0
        p
    else
        mask = ones(Int, p.n)
        mask[I] .= 0    
        J = find_replacement(p.alpha * mask .== 0)
        Poly{T}(syms, p.alpha[J,I], p.c[J])
    end
end

function bsosprob_chordal(degree::Int, order::Int, cliques::Array{Array{Int,1}}, 
                                    obj::Poly{S}, pineq::Array{Poly{T},1}, peq::Array{Poly{V},1})where {S<:Number,T<:Number,V<:Number}

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
        for j=clique_index(cliques, find_replacement(sum(gi.alpha,dims=1)), n)
            append!(Ji[j], i)
        end
    end

    Je = Vector{Int}[]
    for j=1:length(cliques)
        push!(Je, [])
    end

    for (i, hi) in enumerate(peq)
        for j=clique_index(cliques, find_replacement(sum(hi.alpha,dims=1)), n)
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
        #len_ab_d = length(ab_d)
        len_ab_d = length(ab_d) - (j > 1)
        
        #lbj = zeros(length(ab_d))
        lbj = zeros(len_ab_d)
        ai, aj, av = Int[], Int[], Float64[]
        lgta = binomial(obj.n+dmax,dmax)
        
        for (i, ab) in enumerate(ab_d)

            if j > 1
                if i==1 
                    continue
                else
                    i = i - 1
                end            
            end    
                            
            # for efficiency reasons we unroll   h_ab = prod(p .^ ab)
            h_ab = p[1]^ab[1]
            for l=2:length(ab)
                if ab[l] == 1
                    h_ab *= p[l]
                elseif ab[l] == 2
                    h_ab *= p[l]*p[l]
                elseif ab[l] > 2
                    h_ab *= p[l]^(ab[l])
                end
            end
            #println("h_$(ab)=$(h_ab)")
            ak = SparseVector(lgta, linear_index(h_ab, dmax), h_ab.c)

            append!(ai, i*ones(ak.nzind))
            append!(aj, ak.nzind)
            append!(av, ak.nzval)
            
            lbj[i] = (sum(ab[2*length(pj)+1:end]) > 0 ? -Inf : 0.0)
        end            
        #push!(Al, sparse(ai, aj, av, length(ab_d), binomial(length(c)+dmax,dmax)))  
        push!(Al, sparse(ai, aj, av, len_ab_d, binomial(length(c)+dmax,dmax)))  
        push!(lb, lbj)
           
        k = Array{Int,1}(undef, binomial(length(c)+dmax, dmax))
        alpha = zeros(Int, n)
        for (i,a)=enumerate(monomialpowers(length(c), dmax))            
            alpha[c] = a
            k[i] = basis_index(alpha, dmax)
            alpha[c] = 0
        end
        
        push!(El, SparseMatrixCSC{Int,Int}(m,length(k),collect(1:length(k)+1),k,ones(Int,length(k))))
                       
    end
    
    BSOSProb(degree, order, f, Al, El, As, lb)
end

bsosprob_chordal(degree::Int, order::Int, cliques::Array{Array{Int,1}}, obj::Poly{S}, pineq::Array{Poly{T},1}) where {S<:Number,T<:Number} =
    bsosprob_chordal(degree, order, cliques, obj, pineq, Poly{Int}[])

end
