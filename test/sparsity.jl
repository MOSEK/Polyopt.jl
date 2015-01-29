using Polyopt
using Polyopt.Poly, Polyopt.chordal_embedding, Polyopt.moment
using Polyopt.Symbols, Polyopt.indexmap, Polyopt.vectorize
using Base.Test

approxzero{T<:Number}(p::Polyopt.Poly{T}, threshold=1e-6) = (norm(p.c, Inf) < threshold)


#function sparsity{S}(f::Poly{S})
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
    # Haverly1 standard pooling instance
    x14, x24, x35, x36, x45, x46, w4 = variables(["x_{14}", "x_{24}", "x_{35}", "x_{36}", "x_{45}", "x_{46}", "w_4"])

    c14, c24, c35, c36, c45, c46 = [6, 16, 1, -9, -5, -15]
    q1, q2, q3, q5, q6 = [ 3, 1, 2, 5//2, 3//2 ]
    ub1, ub2, ub3, ub4, ub5, ub6 = [ 300, 300, 300, 300, 100, 200 ] // 300

    flow_eq    = [ x14+x24-(x45+x46) ]
    blend_eq   = [ w4*(x45+x46)-(q1*x14+q2*x24) ]
    qual_bnd   = [ q5*(x45+x36)-(w4*x45+q3*x35), q6*(x46+x36)-(w4*x46+q3*x36) ]
    cap_bnd    = [ ub1-x14, ub2-x24, ub3-(x35+x36), ub4-(x45+x46), ub5-(x35+x45), ub6-(x36+x46)]
    flow_bnd   = [ x14,x24,x35,x36,x45,x46]

    obj = dot([c14, c24, c35, c36, c45, c46], [x14, x24, x35, x36, x45, x46])
    prob = momentprob(order, obj, [qual_bnd, cap_bnd, flow_bnd], [flow_eq,blend_eq])

    X, t, y, solsta = solve_mosek(prob)

    Ao = sparsity(1, obj)
    Ai = sum([sparsity(1,v) for v=[obj, qual_bnd, cap_bnd, flow_bnd, flow_eq, blend_eq]])

    println("Sparsity of used monominals:\n", full(Ao + Ai))

    cliques = chordal_embedding(Ao + Ai)
    println("cliques for sparsity:", cliques)

    prob2 = momentprob_chordal_new(1, cliques, obj, [qual_bnd, cap_bnd, flow_bnd], [flow_eq,blend_eq])

    X2, t2, y2, solsta2 = solve_mosek(prob2)
end

if true

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


