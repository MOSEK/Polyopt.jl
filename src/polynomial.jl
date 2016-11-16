import Base: show, *, -, +, isless, ==, convert, conj, truncate, zero, one, promote_rule, A_mul_B!, dot, transpose

export variables

immutable Symbols{T<:AbstractString}
    names :: Array{T,1}
end

immutable Poly{T<:Number}
    n       :: Int             # number of symbols in monomials
    m       :: Int             # number of monomials
    syms    :: Symbols
    c       :: Array{T,1}      # coefficients in polynomial
    alpha   :: Array{Int, 2}   # powers of symbols for all monomials
    deg     :: Int

    function Poly(syms::Symbols, alpha::Array{Int,2}, c::Array{T,1})
        m, n = size(alpha)
        if length(m) == 0 error("Poly must include at least one monomial term") end
        if length(c) != m error("incompatible dimensions") end
        deg = size(alpha, 1) > 0 ? maximum(sum(alpha,2)) : 0
        new(n, m, syms, c, alpha, deg)
    end

    # we allow a constant polynomial without a proper symbol basis
    Poly(c::T) = new(0, 1, Symbols([""]), [c], Array(T, 1, 0), 0)
end

Poly{T<:Number}(syms::Symbols, alpha::Array{Int,2}, c::Array{T,1}) = Poly{T}(syms, alpha, c)
Poly{T<:Number}(c::T) = Poly{T}(c)

convert{S,T}(::Type{Poly{S}}, p::Poly{T}) = Poly(p.syms, p.alpha, convert(Array{S}, p.c))
convert{S,T}(::Type{Poly{S}}, c::T) = Poly{promote_type(S,T)}(p.syms, p.alpha, convert(promote_type(S,T), c))
convert{S<:Number,T<:Number}(::Type{Polyopt.Poly{S}}, v::T) = convert(Polyopt.Poly{S},convert(S,v))

zero{T<:Number}(::Type{Poly{T}}) = Poly(zero(T))
zero{T<:Number}(::Poly{T}) = Poly(zero(T))
one{T<:Number}(::Type{Poly{T}}) = Poly(one(T))
one{T<:Number}(::Poly{T}) = Poly(one(T))

promote_rule{S<:Number,T<:Number}(::Type{S}, ::Type{Poly{T}}) = Poly{promote_type(S,T)}
promote_rule{S<:Number,T<:Number}(::Type{Poly{S}}, ::Type{Poly{T}}) = Poly{promote_type(S,T)}
promote_rule(::Type{Rational}, ::Type{Int}) = Rational

dot(a::Poly, b::Poly) = a*b
dot(a::Number, b::Poly) = a*b
dot(a::Poly, b::Number) = a*b

function show{T}(io::IO, p::Poly{T})
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

-{T<:Number}(p::Poly{T}) = Poly{T}(p.syms, p.alpha, -p.c)
+{S<:Number,T<:Number}(p::Poly{S}, a::T) = p + Poly{T}(p.syms, zeros(Int64, 1, p.n), T[a])
+{S<:Number,T<:Number}(a::S, p::Poly{T}) = Poly{S}(p.syms, zeros(Int64, 1, p.n), S[a]) + p
-{S<:Number,T<:Number}(p::Poly{T}, a::S) = p - Poly{S}(p.syms, zeros(Int64, 1, p.n), S[a])
-{S<:Number,T<:Number}(a::S, p::Poly{T}) = Poly{S}(p.syms, zeros(Int64, 1, p.n), S[a]) - p

*{S<:Number,T<:Number}(p::Poly{S}, a::T) = Poly{promote_type(S,T)}(p.syms, p.alpha, p.c*a)
*{S<:Number,T<:Number}(a::S, p::Poly{T}) = *(p::Poly{T}, a::S)

function add{S<:Number,T<:Number}(p1::Poly{S}, p2::Poly{T})
    U = Array{promote_type(S, T),1}
    p1, p2 = promote_poly(p1, p2)
    Poly{promote_type(S,T)}(p1.syms, vcat(p1.alpha, p2.alpha), vcat(convert(U,p1.c), convert(U,p2.c)))
end

function sub{S<:Number,T<:Number}(p1::Poly{S}, p2::Poly{T})
    U = Array{promote_type(S, T),1}
    p1, p2 = promote_poly(p1, p2)
    Poly{promote_type(S,T)}(p1.syms, vcat(p1.alpha, p2.alpha), vcat(convert(U,p1.c), -convert(U,p2.c)))
end

+{S<:Number,T<:Number}(p1::Poly{S}, p2::Poly{T}) = simplify( add(p1, p2) )
-{S<:Number,T<:Number}(p1::Poly{S}, p2::Poly{T}) = simplify( sub(p1, p2) )

function *{S<:Number,T<:Number}(p1::Poly{S}, p2::Poly{T})
    p1, p2 = promote_poly(p1, p2)
    if p2.m == 1
        Poly(p1.syms, p1.alpha .+ p2.alpha, p1.c*p2.c[1])
    elseif p1.m == 1
        Poly(p1.syms, p2.alpha .+ p1.alpha, p2.c*p1.c[1])
    else
        if p1.m < p2.m
            r = Poly(p1.syms, p2.alpha .+ p1.alpha[1:1,:], p2.c*p1.c[1])
            for k=2:p1.m
                r = add(r, Poly(p1.syms, p2.alpha .+ p1.alpha[k:k,:], p2.c*p1.c[k]))
            end
        else
            r = Poly(p1.syms, p1.alpha .+ p2.alpha[1:1,:], p1.c*p2.c[1])
            for k=2:p2.m
                r = add(r, Poly(p1.syms, p1.alpha .+ p2.alpha[k:k,:], p1.c*p2.c[k]))
            end
        end
        simplify(r)
    end
end

function ^{T<:Number}(p::Poly{T}, a::Int)
    if a<0 error("power must be nonnegative") end
    if a==0 return Poly{T}(p.syms, zeros(Int, 1, p.n), [one(T)]) end
    if p.m == 1 return Poly{T}(p.syms, p.alpha*a, p.c.^a) end

    r = Poly{T}(p)
    for k=2:a r *= p end
    r
end

typealias MatOrVec{T} Union{Array{T,1},Array{T,2}}

*{T<:Number,S<:Number}(a::T, v::MatOrVec{Poly{S}}) = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))
*{T<:Number,S<:Number}(v::MatOrVec{Poly{T}}, a::S) = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))
*{T<:Number,S<:Number}(a::Poly{T}, v::MatOrVec{Poly{S}}) = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))
*{T<:Number,S<:Number}(v::MatOrVec{Poly{S}}, a::Poly{T}) = reshape(Poly{promote_type(T,S)}[ a*vi for vi=v ], size(v))

=={S<:Number,T<:Number}(p1::Poly{S}, p2::Poly{T}) = (p1-p2).n == 0
=={S<:Number,T<:Number}(p::Poly{S}, c::T) = p == Poly(c) 
=={S<:Number,T<:Number}(c::S, p::Poly{T}) = p == Poly(c) 

conj{T<:Number}(p::Poly{T}) = Poly{T}(p.syms, p.alpha, conj(p.c))
convert{T<:Number}(::Type{Poly{T}}, a::T) = Poly{T}(convert(T,a))
isconst{T<:Number}(p::Poly{T}) = p.m == 1 && p.n == 0

function A_mul_B!{S<:Number,T<:Number,U<:Number,V<:Number}(alpha::Poly{S}, A::SparseMatrixCSC{T,Int}, x::Array{Poly{U},1}, beta::Poly{V}, y::Array{Poly{V},1})

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
    y
end

transpose{T<:Number}(p::Poly{T}) = p

# promote polynomial to share the same Symbol basis
function promote_poly{S<:Number,T<:Number}(p1::Poly{S}, p2::Poly{T})
    if (p1.syms == p2.syms) || (isconst(p1) && isconst(p2)) return (p1, p2) end
    if isconst(p1) && ~isconst(p2)
        return (Poly{S}(p2.syms, zeros(Int, 1, p2.n), p1.c), p2)
    elseif ~isconst(p1) && isconst(p2)
        return (p1, Poly{T}(p1.syms, zeros(Int, 1, p1.n), p2.c))
    else
        error("incompatible monomial bases")
    end
end

function evalpoly{T<:Number,S<:Number}(p::Poly{T}, x::AbstractArray{S})
    val = zero(S)
    for j=1:p.m
        val += p.c[j]*prod(vec(x) .^ vec(p.alpha[j,:]))
    end
    val
end

# combine identical terms and remove zero terms
function simplify{T<:Number}(p::Poly{T})

    lgt = 0
    g = Poly{T}(p.syms, copy(p.alpha), copy(p.c))
    labeled = falses(g.m)
    # labeled monomial terms are either zero, or have been merged with other terms.
    for k=1:g.m
        if ~labeled[k]
            for j=k+1:g.m
                if (~labeled[j] && view(g.alpha,k,:) == view(g.alpha,j,:))
                    g.c[k] += g.c[j]
                    labeled[j] = true
                end
            end
            if g.c[k] != zero(T)
                lgt += 1
            else
                labeled[k] = true
            end
        end
    end

    if lgt == 0
        return zero(Poly{T})
    end

    # the unlabeled terms form the reduced polynomial
    alpha = Array(Int, lgt, g.n)
    c = Array(T, lgt)
    lgt = 1
    for k=1:g.m
        if ~labeled[k]
            alpha[lgt,:] = view(g.alpha,k,:)
            c[lgt] = g.c[k]
            lgt += 1
        end
    end
    #perm = sortperm([ view(alpha,i,1:g.n) for i=1:size(alpha,1)], lt=lexless)
    #Poly{T}(g.syms, alpha[perm,:], c[perm])
    Poly{T}(g.syms, alpha, c)
end

function truncate{T<:Number}(p::Poly{T}, threshold=1e-10)
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
    alpha = Array(Int, lgt, p.n)
    c = Array(T, lgt)
    lgt = 1
    for k=1:p.m
        if ~labeled[k]
            alpha[lgt,:] = p.alpha[k,:]
            c[lgt] = p.c[k]
            lgt += 1
        end
    end
    perm = sortperm([ vec(alpha[i,1:p.n]) for i=1:size(alpha,1)], lt=lexless)
    Poly(p.syms, alpha[perm,:], c[perm])
end

# ordering for polynomials
function isless{T<:Number}(p1::Poly{T}, p2::Poly{T})
#     # sort polynomials by degree
#     if max(p1.m, p2.m) > 1 return p1.deg < p2.deg end
# 
#     # otherwise use graded inverse lex-order for monomials
#     grilex_isless(p1.alpha, p2.alpha)
    lexless(p1.alpha, p2.alpha)
end


# graded inverse lexicographic order (lowest total order first, then reverse lex)
# function grilex_isless(a::Array{Int}, b::Array{Int})
#     i, j = sum(a), sum(b)
#     i == j ? ~lexless(a, b) : i < j
# end

variables(syms::Symbols) = [Poly{Int}(syms, [zeros(Int,1,k-1) 1 zeros(Int,1,length(syms.names)-k)], [1]) for k=1:length(syms.names)]
variables{T<:AbstractString}(syms::Vector{T}) = variables(Symbols(syms))
variables{T<:AbstractString}(sym::T, n::Int) = variables([string(sym,i) for i=1:n])