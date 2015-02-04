using Polyopt
using Polyopt.Poly, Polyopt.chordal_embedding, Polyopt.moment
using Polyopt.Symbols, Polyopt.indexmap, Polyopt.vectorize
using Base.Test

approxzero{T<:Number}(p::Polyopt.Poly{T}, threshold=1e-6) = (norm(p.c, Inf) < threshold)

function sparsity(d, f)

    v=monomials(d, variables(f.syms))
    l=length(v)
    M=v*v'

    A = zeros(l,l)
    for j=1:l
        for i=1:l

            for k=1:f.m
                if M[i,j].alpha[1,:] == f.alpha[k,:]
                    A[i,j] = 1
                    break
                end
            end
        end
    end

    for k=1:l
        if norm(A[:,k]) > 0
            A[k,k] = 1
            for j=1:l
                for i=j+1:l
                    if v[k]^2 == v[i]*v[j]
                        A[i,j] = 1
                        A[j,i] = 1
                    end
                end
            end
        end
    end

    sparse(A)
end


function moment_new(order::Int, syms::Symbols, I::Array{Int})
    v = monomials(order, variables(syms))[I]
    v*v'
end

function momentprob_chordal_new{S}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S})

    v = monomials(2*order, variables(obj.syms))
    l = length(v)
    imap = indexmap(v)

    obj.deg <= 2*order || error("obj has degree higher than 2*order")

    p = vectorize(obj, l, imap)
    mom = Array(Any, length(cliques))

    for k=1:length(cliques)
        println("Moment-matrix for clique $(cliques[k]):\n", moment_new(order, obj.syms, cliques[k]))
        mom[k] = vectorize(moment_new(order, obj.syms, cliques[k]), l, imap)
    end

    MomentProb(order, v,  p, mom, [])
end

# The regular chordal embedding always includes the diagonal of A.
# This version filters rows and columns of A with 0 diagonal
function chordal_embedding_new{Tv<:Number,Ti<:Int}(A::SparseMatrixCSC{Tv,Ti})

    idx = find(diag(A) .> 0)
    cliques = chordal_embedding(A[idx,idx])
    Array{Int,1}[ idx[c] for c=cliques ]

end

if false
x, z = variables(["x", "z"])
#f = 2*x^4 + 2*x^3*z - x^2*z^2 + 5*z^4

f = z^2*x^10+2*z^3*x^9+z^4*x^8+2*z^2*x^8+2*z^3*x^7+z^2*x^6

prob = momentprob(f.deg >> 1, f)
X, t, y, solsta = solve_mosek(prob)

v = monomials(f.deg >> 1,[x,z])
@test approxzero( f - t -dot(v,X[1]*v) )

println("f: ", f)

println("full moment-matrix:\n", v*v')
A = sparsity(f.deg >> 1, f)
println("Sparsity of used monominals:\n", full(A))

cliques = chordal_embedding_new(A)
println("cliques for sparsity:", cliques)

prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f)

X2, t2, y2, solsta2 = solve_mosek(prob2)

end

if false
x1, x2, x3 = variables(["x1", "x2", "x3"])

f  = (x1-1)^2 + (x2-1/2)^2 + (x3-1/3)^2
g1 = 2 - x1*x2
g2 = 2 - x2*x3
prob = momentprob(1, f, [g3,g4])
X, t, y, solsta = solve_mosek(prob)
end

if false

x1, x2, x3 = variables(["x1", "x2", "x3"])

f = x1^4*x2^2*x3^2 + x1^2*x2^4*x3^2-2*x1^2*x2^2*x3^2+x3^2+x1^2*x2^2+x1^2*x2^2*x3^4

prob = momentprob(f.deg >> 1, f)
X, t, y, solsta = solve_mosek(prob)

v = monomials(f.deg >> 1,[x1,x2,x3])
@test approxzero( f - t -dot(v,X[1]*v) )

println("f: ", f)

println("full moment-matrix:\n", v*v')
A = sparsity(f.deg >> 1, f)
println("Sparsity of used monominals:\n", full(A))

cliques = chordal_embedding_new(A)
println("cliques for sparsity:", cliques)

prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f)

X2, t2, y2, solsta2 = solve_mosek(prob2)

end


# Frank's counter example
if true

x, = variables(["x"])

f = x^6 + x + 100

prob = momentprob(f.deg >> 1, f)
X, t, y, solsta = solve_mosek(prob)

v = monomials(f.deg >> 1,[x])
@test approxzero( f - t -dot(v,X[1]*v) )

println("f: ", f)

println("full moment-matrix:\n", v*v')
A = sparsity(f.deg >> 1, f)
println("Sparsity of used monominals:\n", full(A))

cliques = chordal_embedding_new(A)

println("cliques for sparsity:", cliques)

prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f)

X2, t2, y2, solsta2 = solve_mosek(prob2)

@test approxzero(f - t2  - sum([dot(v[cliques[i]],X2[i]*v[cliques[i]]) for i=1:length(cliques)]))

end


