

import Base: show, *, -, +, isless, ==, convert, conj, truncate, zero, one, promote_rule, transpose, adjoint, sub, round
import LinearAlgebra: dot, mul!

using SparseArrays

export variables

struct Symbols{T<:AbstractString}
    names :: Array{T,1}
end

struct Poly{T<:Number}
    n        :: Int             # number of symbols in monomials
    m        :: Int             # number of monomials
    syms     :: Symbols
    c        :: Array{T,1}      # coefficients in polynomial
    alpha    :: Array{Int, 2}   # powers of symbols for all monomials
    deg      :: Int

    function Poly{T}(syms::Symbols, alpha::Array{Int,2}, c::Array{T,1}) where {T<:Number}
        m, n = size(alpha)
        if length(c) != m error("incompatible dimensions") end
        
        deg = size(alpha, 1) > 0 ? maximum(sum(alpha,dims=2)) : 0        
        new{T}(n, m, syms, c, alpha, deg)
    end

    # we allow a constant polynomial without a proper symbol basis
    Poly{T}(c::T) where {T<:Number} = new{T}(0, 1, Symbols([""]), [c], Array{T,2}(undef, 1,0), 0)
end

#Poly{T}(syms::Symbols, alpha::Array{Int,2}, c::Array{T,1}) where {T<:Number} = Poly{T}(syms, alpha, c)
#Poly{T}(c::T) where {T<:Number} = Poly{T}(c)

convert(::Type{Poly{S}}, p::Poly{T}) where {S<:Number,T<:Number} = Poly{S}(p.syms, p.alpha, convert(Array{S}, p.c))
convert(::Type{Poly{S}}, c::T) where {S<:Number,T<:Number} = Poly{promote_type(S,T)}(p.syms, p.alpha, convert(promote_type(S,T), c))
convert(::Type{Polyopt.Poly{S}}, v::T) where {S<:Number,T<:Number} = convert(Polyopt.Poly{S},convert(S,v))

zero(::Type{Poly{T}}) where {T<:Number} = Poly{T}(zero(T))
zero(::Poly{T}) where {T<:Number} = Poly{T}(zero(T))
one(::Type{Poly{T}}) where {T<:Number} = Poly{T}(one(T))
one(::Poly{T}) where {T<:Number} = Poly{T}(one(T))

promote_rule(::Type{S}, ::Type{Poly{T}}) where {S<:Number,T<:Number} = Poly{promote_type(S,T)}
promote_rule(::Type{Poly{S}}, ::Type{Poly{T}}) where {S<:Number,T<:Number} = Poly{promote_type(S,T)}
promote_rule(::Type{Rational}, ::Type{Int}) = Rational

dot(a::Poly, b::Poly) = a*b
dot(a::Number, b::Poly) = a*b
dot(a::Poly, b::Number) = a*b

function show(io::IO, p::Poly{T}) where {T<:Number}
    if p.m == 0 print(io, zero(T)) end
    for i=1:p.m
        if i==1
            print(io, p.c[i]>=0 ? "" : "-")
        else
            print(io, p.c[i]>=0 ? "+" : "-")
        end

        ci = abs(p.c[i])
        if sum(p.alpha[i,:]) == 0
            print(io, ci)
        else
            if ci != one(ci) print(io, "$(ci)*") end
            for k=1:p.n
                if p.alpha[i,k] > 0
                    print(io, p.syms.names[k])
                    if p.alpha[i,k] > 1
                        print(io, "^$(p.alpha[i,k])")
                    end
                    if k < p.n && sum(p.alpha[i,k+1:end]) > 0
                        print(io, "*")
                    end
                end
            end
        end
    end
end

-(p::Poly{T}) where {T<:Number}= Poly{T}(p.syms, p.alpha, -p.c)
+(p::Poly{S}, a::T) where {S<:Number,T<:Number} = p + Poly{T}(p.syms, zeros(Int64, 1, p.n), T[a])
+(a::S, p::Poly{T}) where {S<:Number,T<:Number} = Poly{S}(p.syms, zeros(Int64, 1, p.n), S[a]) + p
-(p::Poly{T}, a::S) where {S<:Number,T<:Number} = p - Poly{S}(p.syms, zeros(Int64, 1, p.n), S[a])
-(a::S, p::Poly{T}) where {S<:Number,T<:Number} = Poly{S}(p.syms, zeros(Int64, 1, p.n), S[a]) - p

*(p::Poly{S}, a::T) where {S<:Number,T<:Number} = simplify(Poly{promote_type(S,T)}(p.syms, p.alpha, p.c*a))
*(a::S, p::Poly{T}) where {S<:Number,T<:Number} = *(p::Poly{T}, a::S)

function add(p1::Poly{S}, p2::Poly{T}) where {S<:Number,T<:Number}
    U = Array{promote_type(S, T),1}
    p1, p2 = promote_poly(p1, p2)
    Poly{promote_type(S,T)}(p1.syms, vcat(p1.alpha, p2.alpha), vcat(convert(U,p1.c), convert(U,p2.c)))
end

function sub(p1::Poly{S}, p2::Poly{T}) where {S<:Number,T<:Number}
    U = Array{promote_type(S, T),1}
    p1, p2 = promote_poly(p1, p2)
    Poly{promote_type(S,T)}(p1.syms, vcat(p1.alpha, p2.alpha), vcat(convert(U,p1.c), -convert(U,p2.c)))
end

+(p1::Poly{S}, p2::Poly{T}) where {S<:Number,T<:Number} = simplify( add(p1, p2) )
-(p1::Poly{S}, p2::Poly{T}) where {S<:Number,T<:Number} = simplify( sub(p1, p2) )

function *(p1::Poly{S}, p2::Poly{T}) where {S<:Number,T<:Number}
    p1, p2 = promote_poly(p1, p2)
    if p2.m == 1
        Poly{promote_type(S,T)}(p1.syms, p1.alpha .+ p2.alpha, p1.c*p2.c[1])
    elseif p1.m == 1
        Poly{promote_type(S,T)}(p1.syms, p2.alpha .+ p1.alpha, p2.c*p1.c[1])
    else
        if p1.m < p2.m
            #r = Poly(p1.syms, p2.alpha .+ p1.alpha[1:1,:], p2.c*p1.c[1])
            r = Poly{promote_type(S,T)}(p1.syms, p2.alpha .+ view(p1.alpha,1:1,:), p2.c*p1.c[1])
            for k=2:p1.m
                #r = add(r, Poly(p1.syms, p2.alpha .+ p1.alpha[k:k,:], p2.c*p1.c[k]))
                r = add(r, Poly{promote_type(S,T)}(p1.syms, p2.alpha .+ view(p1.alpha,k:k,:), p2.c*p1.c[k]))
            end
        else
            #r = Poly(p1.syms, p1.alpha .+ p2.alpha[1:1,:], p1.c*p2.c[1])
            r = Poly{promote_type(S,T)}(p1.syms, p1.alpha .+ view(p2.alpha,1:1,:), p1.c*p2.c[1])
            for k=2:p2.m
                #r = add(r, Poly(p1.syms, p1.alpha .+ p2.alpha[k:k,:], p1.c*p2.c[k]))
                r = add(r, Poly{promote_type(S,T)}(p1.syms, p1.alpha .+ view(p2.alpha,k:k,:), p1.c*p2.c[k]))
            end
        end
        simplify(r)
    end
end

function ^(p::Poly{T}, a::Int) where {T<:Number}
    if a<0 error("power must be nonnegative") end
    if a==0 return Poly{T}(p.syms, zeros(Int, 1, p.n), [one(T)]) end
    if a==1 return Poly{T}(p.syms, p.alpha, p.c) end
    if p.m == 1 return Poly{T}(p.syms, p.alpha*a, p.c.^a) end

    powers = Array{Poly{T},1}(undef, Int(1+ceil(log2(a))))
    
    @inbounds begin    
    powers[1] = p
    k = 1
    while 2^k<=a
        powers[k+1] = powers[k]*powers[k]
        k += 1
    end

    b = 2^(k-1)
    r = powers[k]
    while ( b != a)
        k -= 1
        if (b + 2^(k-1) <= a)
            b += 2^(k-1)
            r = r*powers[k]
        end            
    end
    
    end
    r
end

const MatOrVec{T} = Union{Array{T,1},Array{T,2}}

*(a::T, v::MatOrVec{Poly{S}}) where {S<:Number,T<:Number} = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))
*(v::MatOrVec{Poly{T}}, a::S) where {S<:Number,T<:Number} = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))
*(a::Poly{T}, v::MatOrVec{Poly{S}}) where {S<:Number,T<:Number} = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))
*(v::MatOrVec{Poly{S}}, a::Poly{T}) where {S<:Number,T<:Number} = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))

function ==(p1::Poly{S}, p2::Poly{T}) where {S<:Number,T<:Number} 
    r = p1-p2
    r.m == 0 || r.n == 0
end

==(p::Poly{S}, c::T) where {S<:Number,T<:Number} = p == Poly{T}(c) 
==(c::S, p::Poly{T}) where {S<:Number,T<:Number} = p == Poly{S}(c) 

conj(p::Poly{T}) where {T<:Number} = Poly{T}(p.syms, p.alpha, conj(p.c))
convert(::Type{Poly{T}}, a::T) where {T<:Number} = Poly{T}(convert(T,a))
#isconst(p::Poly{T}) where {T<:Number} = (p.m == 0 || p.n == 0)
isconst(p::Poly{T}) where {T<:Number} = p.n == 0

function mul!(alpha::Poly{S}, A::SparseMatrixCSC{T,Int}, x::Array{Poly{U},1}, beta::Poly{V}, y::Array{Poly{V},1}) where {S<:Number,T<:Number,U<:Number,V<:Number}

    @inbounds begin
    for i=1:length(y)
        if beta == zero(y[i])
            y[i] = zero(y[i])
        else
            y[i] *= beta
        end
    end

    for j=1:size(A,2)
        for k=A.colptr[j]:A.colptr[j+1]-1
            y[A.rowval[k]] += alpha*A.nzval[k]*x[j]
        end
    end
    end
    y
end

transpose(p::Poly{T}) where {T<:Number} = p
adjoint(p::Poly{T}) where {T<:Number} = p

# promote polynomial to share the same Symbol basis
function promote_poly(p1::Poly{S}, p2::Poly{T}) where {S<:Number,T<:Number}
    if (p1.syms == p2.syms) || (isconst(p1) && isconst(p2)) return (p1, p2) end
    
    if isconst(p1)  && ~isconst(p2)
        return (Poly{S}(p2.syms, zeros(Int, p1.m, p2.n), p1.c), p2)
    elseif ~isconst(p1) && isconst(p2)
        return (p1, Poly{T}(p1.syms, zeros(Int, p2.m, p1.n), p2.c))
    else
        error("incompatible monomial bases")
    end
end

function evalpoly(p::Poly{T}, x::AbstractArray{S}) where {S<:Number,T<:Number}
    val = zero(S)
    for j=1:p.m
        val += p.c[j]*prod(vec(x) .^ vec(p.alpha[j,:]))
    end
    val
end

# combine identical terms and remove zero terms
function simplify(p::Poly{T}) where {T<:Number}

    # handle simple cases explicitly to speed function up
    if p.m == 1 && p.c[1] == zero(T)
        return Poly{T}(p.syms, Array{Int,2}(undef, 0, p.n), Array{T,1}(undef, 0))
    end
    
    if p.m <= 1 
        return p
    end
    
    #y = ordermap(p)
    #perm = sortperm(y)
    
    c = 1:p.n
    rows = [ view(p.alpha,i,c) for i=1:p.m ]
    perm = sortperm(rows)
    y = p.alpha*rand(p.n)
            
    first = Array{Int,1}(undef, p.m)
    c = Array{T,1}(undef, p.m)
    
    l = 1
    @inbounds begin
    first[l] = 1
    c[l] = p.c[perm[1]]
    
    for k=2:p.m
        if y[perm[k]] == y[perm[first[l]]]
            c[l] += p.c[perm[k]]
        else
            l += 1
            first[l] = k
            c[l] = p.c[perm[k]]
        end
    end

    I = view(c,1:l) .!= 0    
    y = Poly{T}(p.syms, p.alpha[perm[view(first,1:l)[I]],:], view(c,1:l)[I])    
    end
    y
end

function truncate(p::Poly{T}, threshold=1e-10) where {T<:Number}
    lgt = 0
    labeled = falses(p.m)
    
    for k=1:p.m
        if abs(p.c[k]) < threshold
            labeled[k] = true
        else
            lgt += 1
        end
    end

    # the unlabeled terms form the reduced polynomial
    alpha = Array{Int,2}(undef, lgt, p.n)
    c = Array{T,1}(undef, lgt)
    lgt = 1
    for k=1:p.m
        if ~labeled[k]
            alpha[lgt,:] = p.alpha[k,:]
            c[lgt] = p.c[k]
            lgt += 1
        end
    end
    perm = sortperm([ vec(alpha[i,1:p.n]) for i=1:size(alpha,1)], lt=isless)
    Poly{T}(p.syms, alpha[perm,:], c[perm])
end


function round(p::Poly{T}; sigdigits=5) where {T<:Number}
    # the unlabeled terms form the reduced polynomial
    alpha = Array{Int,2}(undef, p.m, p.n)
    c = Array{T,1}(undef, p.m)
    for k=1:p.m
        alpha[k,:] = p.alpha[k,:]
        c[k] = round(p.c[k],sigdigits=sigdigits)
    end
    Poly{T}(p.syms, alpha, c)
end


# ordering for polynomials
function isless(p1::Poly{T}, p2::Poly{T}) where {T<:Number}
    isless(p1.alpha, p2.alpha)
end

# function ordermap{T<:Number}(p::Poly{T})
# 
#     @inbounds begin
# 
#     r = zeros(Float64, p.n)
#     r[p.n] = 1
#     t = p.deg
#     
#     for l=p.n-1:-1:1
#         r[l] = t+1
#         t += p.deg*r[l] 
#     end
#      
#     end
#     
#     p.alpha*r    
# end

# function order{T<:Number}(p::Poly{T})    
#     perm = sortperm(ordermap(p))    
#     Poly{T}(p.syms, p.alpha[perm,:], p.c[perm])
# end

function order(p::Poly{T}) where {T<:Number}
    c = 1:size(p.alpha,2)
    rows = [ view(p.alpha,i,c) for i=1:size(p.alpha,1) ]
    perm = sortperm(rows, order=Base.Order.Lexicographic)
    Poly{T}(p.syms, p.alpha[perm,:], p.c[perm])
end

variables(syms::Symbols) = [Poly{Int}(syms, [zeros(Int,1,k-1) 1 zeros(Int,1,length(syms.names)-k)], [1]) for k=1:length(syms.names)]
variables(syms::Vector{T}) where {T<:AbstractString}= variables(Symbols(syms))
variables(sym::T, n::Int) where {T<:AbstractString} = variables([string(sym,i) for i=1:n])