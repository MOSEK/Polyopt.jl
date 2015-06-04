using Polyopt
using Polyopt.Poly, Polyopt.chordal_embedding, Polyopt.moment
using Polyopt.Symbols, Polyopt.indexmap, Polyopt.vectorize, Polyopt.inverse_indexmap
using Base.Test

approxzero{T<:Number}(p::Polyopt.Poly{T}, threshold=1e-6) = (norm(p.c, Inf) < threshold)

function sparsity(d, f)

    v = monomials(d + mod(d,2), variables(f.syms))
    imap  = indexmap(v)
    iimap = inverse_indexmap(v, imap)

    #n = binomial(v[end].n + v[end].deg >> 1, v[end].n)
    n = binomial(v[end].n + int(v[end].deg/2), v[end].n)
    A = zeros(Int, n, n)

    # stack of monomials still to be added
    stack = Int[ imap[f.alpha[k,:]] for k=1:f.m ]

    # indicator for monomials already added to the graph
    labeled = falses(length(v))

    # indicator for monomials in the stack
    instack = falses(length(v))
    instack[stack] = true

    while !isempty(stack)
        #println("\nNEW ITERATION")
        #println("stack: ", v[stack])

        vk = pop!(stack)
        instack[vk] = false

        #println("popping $(v[vk])")

        labeled[vk] = true
        # fill in all elements in A corresponding to vk
        for (i,j) = iimap[ vk ]

            #println("Adding $(v[vk]) as $(v[i]) * $(v[j]) to A")
            A[i,j] = A[j,i] = 1
            #println("M:\n", v[1:n]*v[1:n]')
            #println("A:\n", full(A))
            # fill the clique formed by monotome adjancency, madj_A(vk)

            aj = A[j:n, j]          # insert A[j,j]
            aj[1] = aj[i-j+1] = 1   # insert A[i,j]
            madj = j - 1 + find(aj)
            #println("madj: ", v[madj])
            #println("clique fill-in:\n", v[madj]*v[madj]')
            for k=1:length(madj)
                for l=k:length(madj)
                    uk  = v[ madj[k] ] * v[ madj[l] ]
                    iuk = imap[ uk.alpha ]
                    #println("Checking to push $(uk) to stack")
                    # Add uk to stack, if not already added
                    if !labeled[ iuk ] && !instack[ iuk ]
                        #println("pushed $(uk) to stack")
                        push!(stack, iuk)
                        instack[ iuk ] = true
                    end
                end
            end
        end

        #println("PAUSED - PRESS ENTER")
        #readline(STDIN)
    end

    sparse(A)

end



function moment_new(order::Int, syms::Symbols, I::Array{Int})
    v = monomials(order, variables(syms))[I]
    v*v'
end

function momentprob_chordal_new{S}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S})

    v = monomials(2*order, variables(obj.syms))
    imap = indexmap(v)

    obj.deg <= 2*order || error("obj has degree higher than 2*order")

    p = vectorize(obj, imap)
    mom = Array(Any, length(cliques))

    for k=1:length(cliques)
        println("Moment-matrix for clique $(cliques[k]):\n", moment_new(order, obj.syms, cliques[k]))
        mom[k] = vectorize(moment_new(order, obj.syms, cliques[k]), imap)
    end

    MomentProb(order, v,  p, mom, [])
end

# The regular chordal embedding always includes the diagonal of A.
# This version filters rows and columns of A with 0 diagonal
function chordal_embedding_new{Tv<:Number,Ti<:Int}(A::SparseMatrixCSC{Tv,Ti})

    idx = find(diag(A) .> 0)
    cliques = chordal_embedding(A[idx,idx], [1:length(idx)])
    Array{Int,1}[ idx[c] for c=cliques ]

end

if true
x, z = variables(["x", "z"])

f = 1 + x*z^4  + z^5 + x^3*z^3 + x^2*z^4 + x*z^5 + x^6 + z^6 

prob = momentprob(f.deg >> 1, f)
X, t, y, solsta = solve_mosek(prob)

v = monomials(f.deg >> 1, variables(f.syms))
@test approxzero( f - t -dot(v,X[1]*v) )

println("f: ", f)

println("full moment-matrix:\n", v*v')
A = sparsity(f.deg, f)
println("Sparsity of used monominals:\n", full(A))

cliques = chordal_embedding_new(A)
println("cliques for sparsity:", cliques)

prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f)

X2, t2, y2, solsta2 = solve_mosek(prob2)

Z = zeros(X[1])
for j=1:length(cliques)
    Z[cliques[j],cliques[j]] += X2[j]
end

@test approxzero(f - t2 - dot(v,Z*v))

end

if false

x1, x2, x3 = variables(["x1", "x2", "x3"])

f =  x1^4*x2^2*x3^2 + x1^2*x2^4*x3^2-2*x1^2*x2^2*x3^2+x3^2+x1^2*x2^2+x1^2*x2^2*x3^4

prob = momentprob(f.deg >> 1, f)
X, t, y, solsta = solve_mosek(prob)

v = monomials(f.deg >> 1, variables(f.syms))
@test approxzero( f - t -dot(v,X[1]*v) )

println("f: ", f)

println("full moment-matrix:\n", v*v')
A = sparsity(f.deg, f)
println("Sparsity of used monominals:\n", full(A))

cliques = chordal_embedding_new(A)
println("cliques for sparsity:", cliques)

prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f)

X2, t2, y2, solsta2 = solve_mosek(prob2)

Z = zeros(X[1])
for j=1:length(cliques)
    Z[cliques[j],cliques[j]] += X2[j]
end

@test approxzero(f - t2 - dot(v,Z*v))


end

if false

x1, x2, x3 = variables(["x1", "x2", "x3"])
#f = 3 - 4*(x1+x2+x3) + 6*(x1^2+x2^2+x3^2) - 4*(x1^3+x2^3+x3^3) + x1*x2*x3 + x1^4 + x2^4 + x3^4
f = - x3^2 - x1*x3^2 - x2*x3^2 + x3^3 + x1^2*x3^2 + x1*x2*x3^2 + x1*x3^3 + x2^2*x3^2 + x2*x3^3 + x3^4

prob = momentprob(f.deg >> 1, f)
X, t, y, solsta = solve_mosek(prob)

v = monomials(f.deg >> 1, variables(f.syms))
#@test approxzero( f - t -dot(v,X[1]*v) )

println("f: ", f)

println("full moment-matrix:\n", v*v')
A = sparsity(f.deg, f)

println("Sparsity of used monominals:\n", full(A))

cliques = chordal_embedding_new(A)

println("cliques for sparsity:", cliques)

prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f)

X2, t2, y2, solsta2 = solve_mosek(prob2)

Z = zeros(X[1])
for j=1:length(cliques)
    Z[cliques[j],cliques[j]] += X2[j]
end

@test approxzero(f - t2 - dot(v,Z*v))

#@test approxzero(f - t2  - sum([dot(v[cliques[i]],X2[i]*v[cliques[i]]) for i=1:length(cliques)]))

end

if false

x, = variables(["x"])
f = x^6 + x^5 + 1

prob = momentprob( int(f.deg/2), f)
X, t, y, solsta = solve_mosek(prob)
v = monomials( f.deg >> 1, variables(f.syms))
@test approxzero( f - t -dot(v,X[1]*v) )

println("f: ", f)

A = sparsity( f.deg, f)
println("Sparsity of used monominals:\n", full(A))

println("A:\n", full(A))
cliques = chordal_embedding_new(A)

println("cliques for sparsity:", cliques)

prob2 = momentprob_chordal_new( int(f.deg/2), cliques, f)

X2, t2, y2, solsta2 = solve_mosek(prob2)

@test approxzero(f - t2  - sum([dot(v[cliques[i]],X2[i]*v[cliques[i]]) for i=1:length(cliques)]))

Z = zeros(4,4)
for j=1:3
    Z[cliques[j],cliques[j]] += X2[j]
end

@test approxzero(f - t2 - dot(v,Z*v))

println("t=$(t), t2=$(t2)")
end


