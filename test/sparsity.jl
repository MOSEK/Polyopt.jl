using Polyopt
using Polyopt.Poly, Polyopt.chordal_embedding, Polyopt.moment, Polyopt.momentprob_chordal
using Polyopt.Symbols, Polyopt.indexmap, Polyopt.vectorize, Polyopt.inverse_indexmap
using Base.Test

using Base.SparseMatrix.CHOLMOD

approxzero{T<:Number}(p::Polyopt.Poly{T}, threshold=1e-6) = (norm(p.c, Inf) < threshold)

function basis_reorder{T<:Number}(d::Int, f::Polyopt.Poly{T})

    v = monomials(d + mod(d,2), variables(f.syms))
    imap  = indexmap(v)
    iimap = inverse_indexmap(v, imap)

    n = binomial(v[end].n + round(Int,v[end].deg/2), v[end].n)
    A = zeros(Int, n, n)
    # stack of monomials still to be added
    stack = Int[ imap[f.alpha[k,:]] for k=1:f.m ]

    # indicator for monomials already added to the graph
    labeled = falses(length(v))

    # indicator for monomials in the stack
    instack = falses(length(v))
    instack[stack] = true

    while !isempty(stack)

        vk = pop!(stack)
        instack[vk] = false

        labeled[vk] = true
        # fill in all elements in A corresponding to vk
        for (i,j) = iimap[ vk ]
            A[i,j] = A[j,i] = 1
        end
    end

    F  = cholfact(round(Float64, sparse(A)), shift=size(A,1))
    s = unsafe_load(F.p)
    
    Int[unsafe_load(s.Perm,i)+1 for i=1:s.n]
end

function chordal_sparsity_sos{T<:Number}(d::Int, f::Polyopt.Poly{T}, perm=[])

    v = monomials(d + mod(d,2), variables(f.syms))
    imap  = indexmap(v)
    iimap = inverse_indexmap(v, imap)

    if perm == [] perm = [1:length(v);] end
    
    iperm = zeros(perm)
    for i=1:length(perm)
        iperm[perm[i]] = i
    end

    n = binomial(v[end].n + round(Int,v[end].deg/2), v[end].n)
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

            ii, jj = max(iperm[i], iperm[j]), min(iperm[i], iperm[j])

            #println("Adding $(v[vk]) as $(v[i]) * $(v[j]) to A")
            A[ii,jj] = A[jj,ii] = 1
            #println("M:")
            #dump(v[perm[1:n]]*v[perm[1:n]]')
            #println("A:\n", full(A))
            # fill the clique formed by monotome adjancency, madj_A(vk)

            #println("ii=$(ii), jj=$(jj)")
            aj = A[jj:n, jj]
            aj[1] = 1   # insert A[j,j]
            madj = jj - 1 + find(aj)

            #println("madj: ", v[perm[madj]])
            #println("clique fill-in:\n", v[perm[madj]]*v[perm[madj]]')
            for k=1:length(madj)
                for l=k:length(madj)
                    uk  = v[ perm[madj[k]] ] * v[ perm[madj[l]] ]
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

    A = sparse(A) 
    S = CHOLMOD.Sparse(round(Float64,A + speye(A)))

    par, post, colcount = Polyopt.chm_analyze_ordering(S, Int32(1), [0:length(v);])
    ns, flag = Polyopt.pothen_sun(par+1, post+1, colcount)

    L = tril(A + speye(A))
    c = Array{Int,1}[ [L.rowval[L.colptr[i]:L.colptr[i+1]-1];] for i = find(flag .< 0) ]
    c = filter(c -> (length(c) > 1 || A[c[1],c[1]] > 0), c)

    A, c
end

function momentprob_chordal_new{S}(order::Int, cliques::Array{Array{Int,1},1}, obj::Poly{S}, perm=[])

    v = monomials(2*order, variables(obj.syms))
    imap = indexmap(v)

    obj.deg <= 2*order || error("obj has degree higher than 2*order")

    p = vectorize(obj, imap)
    mom = Array(Any, length(cliques))

    u = monomials(order, variables(obj.syms))
    if length(perm) > 0
        u = u[perm]
    end

    for k=1:length(cliques)
        uk = u[cliques[k]]
        mom[k] = vectorize(uk*uk', imap)
    end

    MomentProb(order, v,  p, mom, [])
end

if false

    x, = variables(["x"])
    f = x^6 + x^5 + 1

    prob = momentprob( int(f.deg/2), f)
    X, Z, t, y, solsta = solve_mosek(prob)
    v = monomials( f.deg >> 1, variables(f.syms))
    @test approxzero( f - t -dot(v,X[1]*v) )

    A, cliques = chordal_sparsity_sos(f.deg, f)
    prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f)
    X2, Z2, t2, y2, solsta2 = solve_mosek(prob2)

    @test approxzero(f - t2  - sum([dot(v[cliques[i]],X2[i]*v[cliques[i]]) for i=1:length(cliques)]))

    K = zeros(4,4)
    for j=1:length(cliques)
        K[cliques[j],cliques[j]] += X2[j]
    end
    @test approxzero(f - t2 - dot(v,K*v))

    println("t=$(t), t2=$(t2)")
end

# Example 3.5 "Sparse SOS Relaxations for Minimizing Functions that are Summations of Small Polynomials", Nie & Demmel, 2008.
if false
   x1, x2, x3 = variables(["x1", "x2", "x3"]);
   f1 = x1^4 + (x1*x2 - 1)^2
   f2 = x2^2*x3^2 + (x3^2-1)^2

   f  = f1 + f2

   prob = momentprob(f.deg >> 1, f)
   X, Z, t, y, solsta = solve_mosek(prob)

   # Here the CSP matrix is chordal, so Nie's method is identical to Waki's.
   #prob2 = Polyopt.momentprob_chordal(2, Array{Int,1}[ [1,2], [2,3] ], f)
   #X2, Z2, t2,ngth y2, solsta2 = solve_mosek(prob2)

   perm = basis_reorder(f.deg, f)
   perm = [1:length(perm);]
   #perm = [length(perm):-1:1;]

   A, cliques = chordal_sparsity_sos(f.deg, f, perm)
   prob3 = momentprob_chordal_new(f.deg >> 1, cliques, f, perm)
   X3, Z3, t3, y3, solsta3 = solve_mosek(prob3)
end

# Example 3.8 "Sparse SOS Relaxations for Minimizing Functions that are Summations of Small Polynomials", Nie & Demmel, 2008.
if false

   x1, x2, x3 = variables(["x1", "x2", "x3"]);
   f1 = 0.5*(x1^2 + x2^2) + 2*x1*x2
   f2 = 0.5*(x2^2 + x3^2) + 2*x2*x3
   f3 = 0.5*(x1^2 + x3^2) + 2*x1*x3
   f  = f1 + f2 + f3

   prob = momentprob(1, f)
   X, Z, t, y, solsta = solve_mosek(prob)

   prob2 = Polyopt.momentprob_chordal(1, Array{Int,1}[ [1,2], [2,3], [1,3]], f)
   X2, Z2, t2, y2, solsta2 = solve_mosek(prob2)

   A, cliques = chordal_sparsity_sos( f.deg, f)
   prob3 = momentprob_chordal_new(f.deg >> 1, cliques, f)
   X3, Z3, t3, y3, solsta3 = solve_mosek(prob3)

end

# Generalized Rosenbrock, "Sparse SOS Relaxations for Minimizing Functions that are Summations of Small Polynomials", Nie & Demmel, 2008.
if false

    n = 20
    x = variables(ASCIIString[string("x",k) for k=1:n])
    f = sum([ 100*(x[i]-x[i-1]^2)^2 + (1-x[i])^2 for i=2:n ])

    #prob = momentprob(2, f)
    #X, Z, t, y, solsta = solve_mosek(prob)

    prob2 = Polyopt.momentprob_chordal(2, Array{Int,1}[ [i,i-1] for i=2:n ], f);
    X2, Z2, t2, y2, solsta2 = solve_mosek(prob2);

    perm = basis_reorder(f.deg, f)
    A, cliques = chordal_sparsity_sos(f.deg, f, perm)
    prob3 = momentprob_chordal_new(f.deg >> 1, cliques, f, perm)
    X3, Z3, t3, y3, solsta3 = solve_mosek(prob3);

end

# Ex 3.14, "Sparse SOS Relaxations for Minimizing Functions that are Summations of Small Polynomials", Nie & Demmel, 2008.
if false

    n = 10
    x = variables(ASCIIString[string("x",k) for k=1:n])
    h = [ 2*x[1]^2 - 3*x[1] + 2*x[2] - 1, [ 2*x[i]^2 + x[i-1] - 3*x[i] + 2*x[i+1] - 1 for i=2:n-1 ], 2*x[n]^2 + x[n-1] - 3*x[n] - 1 ]
    f = sum([ hi^2 for hi=h ]) + 1e-5*dot(randn(n), x)

    I =  Array{Int,1}[Array[[1, 2]]; [[i - 1, i, i + 1] for i = 2:n-1]; Array[n-1:n]]
    prob = Polyopt.momentprob_chordal(2, I, f);
    X, Z, t, y, solsta = solve_mosek(prob2);

    perm = basis_reorder(f.deg, f)
    #perm = [1:length(perm);]
    perm = [length(perm):-1:1;]
    A, cliques = chordal_sparsity_sos(f.deg, f, perm)
    prob2 = momentprob_chordal_new(f.deg >> 1, cliques, f, perm)
    X2, Z2, t2, y2, solsta2 = solve_mosek(prob2);

    println([prob.basis[1:n+1] y[1:n+1] y2[1:n+1]])
    
    err1 = [ Polyopt.evalpoly(hi, y[2:n+1]) for hi = h ]
    err2 = [ Polyopt.evalpoly(hi, y2[2:n+1]) for hi = h ]
            
end

if true
   x1, x2, x3 = variables(["x1", "x2", "x3"]);

   f  = (x1-1)^4 + (x2-0.5)^4 + (x3-1/3)^4 + x1*x2*x3

   prob = momentprob(f.deg >> 1, f)
   X, Z, t, y, solsta = solve_mosek(prob)

   # Here the CSP matrix is chordal, so Nie's method is identical to Waki's.
   #prob2 = Polyopt.momentprob_chordal(2, Array{Int,1}[ [1,2], [2,3] ], f)
   #X2, Z2, t2,ngth y2, solsta2 = solve_mosek(prob2)

   perm = basis_reorder(f.deg, f)
   #perm = [1:length(perm);]
   #perm = [length(perm):-1:1;]

   A, cliques = chordal_sparsity_sos(f.deg, f, perm)
   prob3 = momentprob_chordal_new(f.deg >> 1, cliques, f, perm)
   X3, Z3, t3, y3, solsta3 = solve_mosek(prob3)
end
