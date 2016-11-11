using Polyopt
using Base.Test

approxzero{T<:Number}(p::Polyopt.Poly{T}, threshold=1e-6) = (norm(p.c, Inf) < threshold)

# example of unconstrained global optimization
let
    println("Unconstrained global optimization")
    
    x, z = variables(["x", "z"])
    f = 4.0*x^2 + x*z - 4*z^2 - 21//10*x^4 + 4*z^4 + 1//3*x^6

    # perturb problem to find global optimizer
    f = f + 1e-3*(x+z)
    prob = momentprob(3, f)
    
    X, Z, t, y, solsta = solve_mosek(prob)

    # test that (x, z) = (y[8], y[2]) is globally optimal 
    xo = Polyopt.vectorize([x,z],6)*y
    @test abs(t - Polyopt.evalpoly(f, xo)) < 1e-6
end

# gloptipoly example
let
    println("Gloptipoly example")
    
    x1, x2, x3 = variables(["x1", "x2", "x3"])

    # This is the Gloptipoly example.  
    f = -2*x1 + x2 - x3
    g = [ 24 - 20*x1 + 9*x2 - 13*x3 + 4*x1^2 - 4*x1*x2 + 4*x1*x3 + 2*x2^2 - 2*x2*x3 + 2*x3^2,
          4 - (x1 + x2 + x3),
          6 - (3*x2 + x3),
          x1, 2-x1,
          x2,
          x3, 3-x3 ]

    # perturb problem to find global optimizer        
    f = f + 1e-3*(x1+x2+x3) 
    
    prob = momentprob(4, f, g)
    X, Z, t, y, solsta = solve_mosek(prob)
        
    # test that (x1, x2, x3) extracted from y is globally optimal
    xo = Polyopt.vectorize([x1,x2,x3], 8)*y
    @test abs(t - Polyopt.evalpoly(f, xo)) < 1e-6
    @test all( [ Polyopt.evalpoly(gi, xo)  for gi=g ] .> -1e-4 )
end

# determine if f(x,z) is can be written as a SOS
let
    println("SOS example 1")

    x, z = variables(["x", "z"])
    f = 2*x^4 + 2*x^3*z - x^2*z^2 + 5*z^4

    prob = momentprob(2, f)
    X, Z, t, y, solsta = solve_mosek(prob)
    @test t > -1e-6

    v = monomials(2, [x,z])
    @test approxzero( f - dot(v,X[1]*v) )
end

# in this case f(x,z) is not SOS, but f(x,z)-t is SOS
let
    println("SOS example 2")
    
    x, z = variables(["x", "z"])

    f = (x + z + 1)^2 - 1
    prob = momentprob(1, f)

    X, Z, t, y, solsta = solve_mosek(prob)
    @test abs(t + 1) < 1e-6

    v = monomials(1, [x,z])

    # test that (f(x,z)-t) is SOS
    @test approxzero( f - t - dot(v,X[1]*v) )
end

# test duality for problem with inequality constraints
let
    println("Duality example 1")

    x, z = variables(["x", "z"])

    f = x^4 + z^2 + 1
    g = x^2 - z^2 - 1
    prob = momentprob(2, f, [g])

    X, Z, t, y, solsta = solve_mosek(prob)

    v1 = monomials(2,[x,z])
    v2 = monomials(1,[x,z])

    # test that f(x,z) - t - g(x,z)*s(x,z) is SOS,  where s(x,z) = v2'*X[2]*v2 is SOS
    @test approxzero( f - t - dot(v1,X[1]*v1) - g*dot(v2,X[2]*v2) )
end

# Test duality for problem with both inequality and equality constraints 
let
    println("Duality example 2")

    x1, x2 = variables(["x1", "x2"])
    f = x1 + x2  
    g = x1 - x2 
    h = x1^2 + x2^2 - 1 

    prob = momentprob(2, f, [g], [h])
    X, Z, t, y, solsta = solve_mosek(prob);
    
    u = monomials(2, [x1,x2])
    v = monomials(1, [x1,x2])
    
    # test optimality
    xo = Polyopt.vectorize([x1,x2], 4)*y
    @test abs(t - Polyopt.evalpoly(f, xo)) < 1e-4
    @test abs(Polyopt.evalpoly(h, xo)) < 1e-4
    @test Polyopt.evalpoly(g, xo) > -1e-4
    
    # test that f(x1,x2) - t - g(x1,x2)*s(x1,x2) - h(x1,x2)*w(x1,x2) is SOS,  where s1(x1,x2) = v'*X[2]*v is SOS, but w(x1,x2) = v'*Z[1]*z is not     
    @test approxzero(f - t - (dot(u, X[1]*u) + g*dot(v, X[2]*v) + h*dot(v, Z[1]*v)))
end    

# Example 6.23 "Sums of Squares, Moment Matrices and Polynomial Optimization", M. Laurent
let
    println("Duality example 3")
    
    x1, x2 = variables(["x1", "x2"])

    f = -x1 - x2
    g = [ 2*x1^4 - 8*x1^3 + 8*x1^2 + 2 - x2,
          4*x1^4 - 32*x1^3 + 88*x1^2 - 96*x1 + 36 - x2,
          x1, 3-x1,
          x2, 4-x2 ]

    prob = momentprob(4, f, g)

    X, Z, t, y, solsta = solve_mosek(prob)
    xo = Polyopt.vectorize([x1,x2], 8)*y
    @test abs(t - Polyopt.evalpoly(f, xo)) < 1e-6
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
end

# Example 6.25 "Sums of Squares, Moment Matrices and Polynomial Optimization", M. Laurent
let
    println("Constrained optimization, Laurent 1")
    
    x1, x2, x3 = variables(["x1", "x2", "x3"])
    
    f = x1^2*x2^2*(x1^2 + x2^2 - 3*x3^2) + x3^6 + 1e-2*(x1^6 + x2^6 + x3^6)
    g = 1 - (x1^2 + x2^2 + x3^2)

    prob = momentprob(4, f, [g])

    X, Z, t, y, solsta = solve_mosek(prob)
    xo = Polyopt.vectorize([x1,x2,x3], 8)*y
    @test abs(t - Polyopt.evalpoly(f, xo)) < 1e-6
    @test Polyopt.evalpoly(g, xo) > -1e-4 
end

# Example 6.26 "Sums of Squares, Moment Matrices and Polynomial Optimization", M. Laurent
let
    println("Constrained optimization, Laurent 2")
    
    x1, x2, x3 = variables(["x1", "x2", "x3"])
    
    f = 1 + 1e-3*(x1+x2+x3)
    h = [ 5*x1^9 - 6*x1^5*x2 + x1*x2^4 + 2*x1*x3,
          -2*x1^6*x2 + 2*x1^2*x2^3 + 2*x2*x3,
          x1^2 + x2^2 - 0.265625 ]

    prob1 = momentprob(6, f, Polyopt.Poly{Int}[], h)

    X, Z, t, y, solsta = solve_mosek(prob1)
    xo = Polyopt.vectorize([x1,x2,x3],12)*y
    @test abs(t - Polyopt.evalpoly(f, xo)) < 1e-4
    @test norm([ Polyopt.evalpoly(hi, xo) for hi=h ], Inf) < 1e-4 
end

# Example 8.8 "Sums of Squares, Moment Matrices and Polynomial Optimization", M. Laurent
let
    println("Constrained optimization, Laurent 3")
    
    x1, x2, x3 = variables(["x1", "x2", "x3"]);

    g1, g2 = x1^4 + (x1*x2-1)^2, x2^2*x3^2 + (x3^2-1)^2;
    f = g1+g2;
    g = [g1, g2, 1-x1^2, 1-x2^2, 1-x3^2];

    prob = momentprob(3, f, g);
    X, Z, t, y, solsta = solve_mosek(prob);

    probc = momentprob_chordalembedding(3, f, g);
    Xc, Zc, tc, yc, solstac = solve_mosek(probc);
    
    @test abs(t-tc) < 1e-6
end    
    
# http://gamsworld.org/global/globallib/ex2_1_1.htm
let
    println("Globallib ex2_1_1")
    
    x1, x2, x3, x4, x5 = variables(["x1", "x2", "x3", "x4", "x5"])

    f =  -0.5*(100*x1^2 + 100*x2^2 + 100*x3^2 + 100*x4^2 + 100*x5^2) + 42*x1 + 44*x2 + 45*x3 + 47*x4 + 47.5*x5
    g = [40 - (20*x1 + 12*x2 + 11*x3 + 7*x4 + 4*x5), x1, 1-x1, x2, 1-x2, x3, 1-x3, x4, 1-x4, x5, 1-x5]
    
    prob = momentprob(3, f, g)
    X, Z, t, y, solsta = solve_mosek(prob);
    xo = Polyopt.vectorize([x1,x2,x3,x4,x5],6)*y    
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
end

# http://gamsworld.org/global/globallib/ex2_1_2.htm
let
    println("Globallib ex2_1_2")
    
    x1, x2, x3, x4, x5, x6 = variables(["x1", "x2", "x3", "x4", "x5", "x6"])
    
    f = -0.5*(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5) - 10.5*x1 - 7.5*x2 - 3.5*x3 - 2.5*x4 - 1.5*x5 - 10*x6
    g = [ 6.5 - (6*x1 + 3*x2 + 3*x3 + 2*x4 + x5), 20 - (10*x1 + 10*x3 + x6), x1, 1-x1, x2, 1-x2, x3, 1-x3, x4, 1-x4, x5, 1-x5, x6 ]
    
    prob = momentprob(2, f, g)
    X, Z, t, y, solsta = solve_mosek(prob);
    xo = Polyopt.vectorize([x1,x2,x3,x4,x5,x6],4)*y
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
end

# http://gamsworld.org/global/globallib/ex2_1_4.htm
let
    println("Globallib ex2_1_4")

    x1, x2, x3, x4, x5, x6 = variables(["x1", "x2", "x3", "x4", "x5", "x6"])
    
    f = (6.5*x1 - 0.5*x1*x1) - x2 - 2*x3 - 3*x4 - 2*x5 - x6
    g = [ 16-(x1 + 2*x2 + 8*x3 + x4 + 3*x5 + 5*x6),
          -1-(-8*x1 - 4*x2 - 2*x3 + 2*x4 + 4*x5 - x6),
          24-(2*x1 + 0.5*x2 + 0.2*x3 - 3*x4 - x5 - 4*x6),
          12-(0.2*x1 + 2*x2 + 0.1*x3 - 4*x4 + 2*x5 + 2*x6),
           3-(-0.1*x1 - 0.5*x2 + 2*x3 + 5*x4 - 5*x5 + 3*x6),
           x1, 1-x1, x2, x3, x4, 1-x4, x5, 1-x5, x6, 2-x6 ]
           
    prob = momentprob(2, f, g)
    X, Z, t, y, solsta = solve_mosek(prob);
    xo = Polyopt.vectorize([x1,x2,x3,x4,x5,x6],4)*y
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
end
  
# http://gamsworld.org/global/globallib/ex3_1_2.htm
let 
    println("Globallib ex3_1_2")
    
    x1, x2, x3, x4, x5 = variables(["x1", "x2", "x3", "x4", "x5"])
 
    f = -40792.141 + 0.8356891*x1*x5 + 37.293239*x1 + 5.3578547*x3*x3 
    g = [ 6.665593-(0.0056858*x2*x5 - 0.0022053*x3*x5 + 0.0006262*x1*x4),
          85.334407-(0.0022053*x3*x5 - 0.0056858*x2*x5 - 0.0006262*x1*x4),
          29.48751-(0.0071317*x2*x5 + 0.0021813*x3*x3 + 0.0029955*x1*x2),
          -9.48751-(-0.0071317*x2*x5 - 0.0021813*x3*x3 - 0.0029955*x1*x2),
          15.599039-(0.0047026*x3*x5 + 0.0019085*x3*x4 + 0.0012547*x1*x3),
          -10.699039-((-0.0047026*x3*x5) - 0.0019085*x3*x4 - 0.0012547*x1*x3),
          x1-78, 102-x1, x2-33, 45-x2, x3-27, 45-x3, x4-27, 45-x4, x5-27, 45-x5 ]

    prob = momentprob(2, f, g)
    X, Z, t, y, solsta = solve_mosek(prob);
    xo = Polyopt.vectorize([x1,x2,x3,x4,x5],4)*y
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
end

# http://gamsworld.org/global/globallib/ex3_1_1.htm
let
    println("Globallib ex3_1_1")

    x1, x2, x3, x4, x5, x6, x7, x8 = variables(["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"])
    #f = x1 + x2 + x3
    #g = [ 1 - (0.0025*x4 + 0.0025*x6),
    #      1 - (-0.0025*x4 + 0.0025*x5 + 0.0025*x7),
    #      1 - (-0.01*x5 + 0.01*x8),
    #      83333.333 - (100*x1 - x1*x6 + 833.33252*x4),
    #      -(x2*x4 - x2*x7 - 1250*x4 + 1250*x5),
    #      -1250000 - (x3*x5 - x3*x8 - 2500*x5),
    #      x1-100,  10000-x1, 
    #      x2-1000, 10000-x2, 
    #      x3-1000, 10000-x3,
    #      x4-10,   1000-x4,                                    
    #      x5-10,   1000-x5,                                    
    #      x6-10,   1000-x6,                                    
    #      x7-10,   1000-x7,                                    
    #      x8-10,   1000-x8 ]                                    

    # rescale model to make it easier to solve
    f = x1 + x2 + x3
    g = [ 1 - (2.5*x4 + 2.5*x6),
          1 - (-2.5*x4 + 2.5*x5 + 2.5*x7),
          1 - (-10*x5 + 10*x8),
          0.08333 - (x1 - 10*x1*x6 + 0.83333*x4),
          -(x2*x4 - x2*x7 - 0.1250*x4 + 0.125*x5),
          -0.125 - (x3*x5 - x3*x8 - 0.25*x5),
          x1-0.01, 1-x1, 
          x2-0.1,  1-x2, 
          x3-0.1,  1-x3,
          x4-0.01, 1-x4,                                    
          x5-0.01, 1-x5,                                    
          x6-0.01, 1-x6,                                    
          x7-0.01, 1-x7,                                    
          x8-0.01, 1-x8 ]                                    

    probc = momentprob_chordalembedding(3, f, g)
    Xc, Zc, tc, yc, solstac = solve_mosek(probc);
    xo = Polyopt.vectorize([x1,x2,x3,x4,x5,x6,x7,x8],6)*yc
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
end

#http://gamsworld.org/global/globallib/ex5_4_2.htm
let 
    println("Globallib ex5_4_2")
    
    x1,x2,x3,x4,x5,x6,x7,x8 = variables(["x1","x2","x3","x4","x5","x6","x7","x8"])

    #f = x1 + x2 + x3 
    #g = [ 400 - (x4 + x6),
    #      300 - (-x4 + x5 + x7),
    #      100 - (-x5 + x8),
    #      83333.3333333333 - (x1 - x1*x6 + 833.333333333333*x4),
    #      -(x2*x4 - x2*x7 - 1250*x4 + 1250*x5),
    #      -1250000 - (x3*x5 - x3*x8 - 2500*x5),
    #      x1-100, 10000-x1, 
    #      x2-1000, 10000-x2, 
    #      x3-1000, 10000-x3, 
    #      x4-10, 1000-x4, 
    #      x5-10, 1000-x5, 
    #      x6-10, 1000-x6, 
    #      x7-10, 1000-x7, 
    #      x8-10, 1000-x8 ]

    # rescale model to make it easier to solve
    f = x1 + x2 + x3
    g = [ 400/1000 - (x4 + x6),
          300/1000 - (-x4 + x5 + x7),
          100/1000 - (-x5 + x8),
          .083333 - (0.01*x1 - 10*x1*x6 + 0.83333*x4),      
          -(x2*x4 - x2*x7 - 0.125*x4 + 0.125*x5),
          -0.125 - (x3*x5 - x3*x8 - 0.25*x5),
          x1-100/10000, 1-x1, 
          x2-1000/10000, 1-x2, 
          x3-1000/10000, 1-x3, 
          x4-10/1000, 1-x4, 
          x5-10/1000, 1-x5, 
          x6-10/1000, 1-x6, 
          x7-10/1000, 1-x7, 
          x8-10/1000, 1-x8 ]

    probc = momentprob_chordalembedding(3, f, g)
    Xc, Zc, tc, yc, solstac = solve_mosek(probc);
    xo = Polyopt.vectorize([x1,x2,x3,x4,x5,x6,x7,x8],6)*yc
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
end

let
    println("Chordal relaxation")
    
    x1,x2,x3 = variables(["x1","x2","x3"])

    f = x1-x2+x3 
    g = [1-x1*x2,
         1-x2*x3,
        1-x1^2,
        1-x2^2,
        1-x3^2 ]

    order = 2     
    prob = momentprob(order, f, g)
    X, Z, t, y, solsta = solve_mosek(prob);
    xo = Polyopt.vectorize([x1,x2,x3],2*order)*y
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
     
    K = Array{Int,1}[ [1,2], [2,3] ]    
    probc = momentprob_chordal(order, K, f, g) 
    Xc, Zc, tc, yc, solstac = solve_mosek(probc);
    xo = Polyopt.vectorize([x1,x2,x3],2*order)*yc
    @test all( [ Polyopt.evalpoly(gi, xo) for gi=g ] .> -1e-4 )
    
    u1 = monomials(order, [x1,x2,x3][K[1]])
    u2 = monomials(order, [x1,x2,x3][K[2]])
    v1 = monomials(order-1, [x1,x2,x3][K[1]]) 
    v2 = monomials(order-1, [x1,x2,x3][K[2]]) 
     
    r = dot(u1, Xc[1]*u1) + dot(u2, Xc[2]*u2) +  
            g[1]*dot(v1, Xc[3]*v1) +  
            g[2]*dot(v2, Xc[4]*v2) + 
            g[3]*dot(v1, Xc[5]*v1) +
            g[4]*(dot(v1, Xc[6]*v1) + dot(v2, Xc[7]*v2)) +
            g[5]*dot(v2, Xc[8]*v2)
        
    @test approxzero( f - t - r )
end
