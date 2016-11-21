using Polyopt
using Base.Test

# Simple convex quadratic problem with two variables in different cliques
let
    println("BSOS test 1")
    x = variables(["x1", "x2"])
    f = x[1] - x[2]
    g = [ 1.0-x[1]^2, 1-x[2]^2 ]

    I = Array{Int,1}[ [1], [2] ]
    prob = bsosprob_chordal(1, 1, I, f, g)
    X, t, l, y, solsta = solve_mosek(prob)
end

# P4_2 from Weisser's paper
let
    println("BSOS test 2")
    x = variables(["x1", "x2", "x3", "x4"])
    f = x[1]^2 - x[2]^2 + x[3]^2 - x[4]^2 + x[1] - x[2]
    g = [ 2*x[1]^2 + 3*x[2]^2 + 2*x[1]*x[2] + 2*x[3]^2 + 3*x[4]^2 + 2*x[3]*x[4],
          3*x[1]^2 + 2*x[2]^2 - 4*x[1]*x[2] + 3*x[3]^2 + 2*x[4]^2 - 4*x[3]*x[4],
          x[1]^2 + 6*x[2]^2 - 4*x[1]*x[2] + x[3]^2 + 6*x[4]^2 - 4*x[3]*x[4],
          x[1]^2 + 4*x[2]^2 - 3*x[1]*x[2] + x[3]^2 + 4*x[4]^2 - 3*x[3]*x[4],
          2*x[1]^2 + 5*x[2]^2 + 3*x[1]*x[2] + 2*x[3]^2 + 5*x[4]^2 + 3*x[3]*x[4],
          x[1],x[2],x[3],x[4]]
         
    I = Array{Int,1}[ [1,2,3,4] ]
    prob = bsosprob_chordal(1, 1, I, f, g)
    X, t, l, y, solsta = solve_mosek(prob)
end             

# P4_4 from Weisser's paper (we cannot reproduce bounds)
let
    println("BSOS test 3")
    x = variables(["x1", "x2", "x3", "x4"])
    f = x[1]^4 - x[2]^4 + x[3]^4 - x[4]^4
    g = [ 2*x[1]^4 + 3*x[2]^2 + 2*x[1]*x[2] + 2*x[3]^4 + 3*x[4]^2 + 2*x[3]*x[4],
          3*x[1]^2 + 2*x[2]^2 - 4*x[1]*x[2] + 3*x[3]^2 + 2*x[4]^2 - 4*x[3]*x[4],
          x[1]^2 + 6*x[2]^2 - 4*x[1]*x[2] + x[3]^2 + 6*x[4]^2 - 4*x[3]*x[4],
          x[1]^2 + 4*x[2]^4 - 3*x[1]*x[2] + x[3]^2 + 4*x[4]^4 - 3*x[3]*x[4],
          2*x[1]^2 + 5*x[2]^2 + 3*x[1]*x[2] + 2*x[3]^2 + 5*x[4]^2 + 3*x[3]*x[4],
          x[1], x[2], x[3], x[4] ]
          
    I = Array{Int,1}[ [1,2,3,4] ]
    prob = bsosprob_chordal(2, 2, I, f, g)
    X, t, l, y, solsta = solve_mosek(prob)
end

# P4_6 from Weisser's paper
let
    println("BSOS test 4")
    x = variables(["x1", "x2", "x3", "x4"])
    f = x[1]^4*x[2]^2 + x[1]^2*x[2]^4 - x[1]^2*x[2]^2 + x[3]^4*x[4]^2 + x[3]^2*x[4]^4 - x[3]^2*x[4]^2
    g = [ x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2,
          3*x[1]^2 + 2*x[2]^2 - 4*x[1]*x[2] + 3*x[3]^2 + 2*x[4]^2 - 4*x[3]*x[4],
          x[1]^2 + 6*x[2]^4 - 8*x[1]*x[2] + x[3]^2 + 6*x[4]^4 - 8*x[3]*x[4] + 2.5,
          x[1]^4 + 3*x[2]^4 + x[3]^4 + 3*x[4]^4,
          x[1]^2 + x[2]^3 + x[3]^2 + x[4]^4,
          x[1], x[2], x[3], x[4] ]
      
    I = Array{Int,1}[ [1,2,3,4] ]
    prob = bsosprob_chordal(3, 3, I, f, g)
    X, t, l, y, solsta = solve_mosek(prob)
end

# P4_8 from Weisser's paper
let
    println("BSOS test 5")
    x = variables(["x1", "x2", "x3", "x4"])
    f = x[1]^4*x[2]^2 + x[1]^2*x[2]^6 - x[1]^2*x[2]^2 + x[3]^4*x[4]^2 + x[3]^2*x[4]^6 - x[3]^2*x[4]^2
    g = [ x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2,
          3*x[1]^2 + 2*x[2]^2 - 4*x[1]*x[2] + 3*x[3]^2 + 2*x[4]^2 - 4*x[3]*x[4],
          x[1]^2 + 6*x[2]^4 - 8*x[1]*x[2] + x[3]^2 + 6*x[4]^4 - 8*x[3]*x[4] + 2.5,
          x[1]^4 + 3*x[2]^4 + x[3]^4 + 3*x[4]^4,
          x[1]^2 + x[2]^3 + x[3]^2 + x[4]^4,
          x[1], x[2], x[3], x[4] ]

    I = Array{Int,1}[ [1,2,3,4] ]
    prob = bsosprob_chordal(3, 4, I, f, g)
    X, t, l, y, solsta = solve_mosek(prob)
end

# Haverly1 from Marandi's paper 
let
    println("BSOS test 6")
    x = variables(["x1", "x2", "x3", "x4", "x5"])
    f = -200*x[2]*(15*x[1]-12) - 200*x[3]*(15*x[1]-6) + 200*x[4] - 1000*x[5]
    g = [-3//4*(x[1]-1)*(x[2]+x[3]),
          1//4*(3*x[1]-1)*(x[2]+x[3]),
          1 - 2*(x[2]+x[4]),
          1 - (x[3]+x[5]),
          1//2*(x[4]+x[2]) - 2//5*x[4] - 3//5*x[1]*x[2],
          1//2*(x[5]+x[3]) - 2//3*x[5] - x[1]*x[3],
          x[1], x[2], x[3], x[4], x[5] ]
    
    #I = Array{Int,1}[ [1,2,3,4,5] ]
    I = Polyopt.chordal_embedding(Polyopt.correlative_sparsity(f,g))

    prob = bsosprob_chordal(3, 2, I, f, g)
    X, t, l, y, solsta = solve_mosek(prob)
end

# Generalized Rosenbrock function
let
    println("BSOS test 7")
    n = 100
    x = variables("x", n)

    I = Array{Int,1}[]
    for i=1:n-1
        push!(I, [i, i+1])
    end

    f = sum([ 100*(x[i]-x[i-1]^2)^2 + (1-x[i])^2 for i=2:n ])
    g = vcat(Polyopt.Poly{Int}[x[i] for i=1:n], Polyopt.Poly{Int}[2 - sum([xk^2 for xk=x[Ik]]) for Ik in I])

    prob = bsosprob_chordal(3, 2, I, f, g);
    X, t, l, y, solsta = solve_mosek(prob);
    t
end
