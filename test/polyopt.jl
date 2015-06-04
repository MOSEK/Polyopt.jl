approxzero{T<:Number}(p::Polyopt.Poly{T}, threshold=1e-6) = (norm(p.c, Inf) < threshold)

# example of unconstrained global optimization
let
    x, z = variables(["x", "z"])
    f = 4*x^2 + x*z - 4*z^2 - 21//10*x^4 + 4*z^4 + 1//3*x^6
    
    # perturb problem to find global optimizer
    f = f + 1e-5*dot(randn(2),[x,z])
    
    prob = momentprob(3, f)
    X, t, y, solsta = solve_mosek(prob, 1e-10)

    # test that (x, z) = y[2:3] is globally optimal 
    @test abs(t - Polyopt.evalpoly(f, y[2:3])) < 1e-6
end

# gloptipoly example
let
    x1, x2, x3 = variables(["x1", "x2", "x3"])

    # This is the Gloptipoly example.  We normalize it to 0 <= x1 <= 1, 0 <= x3 <= 1 below for stability reasons
    #  f = -2*x1 + x2 - x3
    # g1 = 24 - 20*x1 + 9*x2 - 13*x3 + 4*x1^2 - 4*x1*x2 + 4*x1*x3 + 2*x2^2 - 2*x2*x3 + 2*x3^2
    # g2 = 4 - (x1 + x2 + x3)
    # g3 = 6 - (3*x2 + x3)
    
    f = -2*2*x1 + x2 - 3*x3

    # perturb problem to find global optimizer
    f = f + 1e-4*dot(randn(3),[x1,x2,x3]) 
    
    g1 = (24 - 20*2*x1 + 9*x2 - 13*3*x3 + 4*4*x1^2 - 4*2*x1*x2 + 4*2*3*x1*x3 + 2*x2^2 - 2*3*x2*x3 + 2*9*x3^2)
    g2 = (4 - (2*x1 + x2 + 3*x3))
    g3 = (6 - (3*x2 + 3*x3))
    
    prob = momentprob(4, f, [g1, g2, g3, 1-x1, x1, x2, 1-x3, x3])
    X, t, y, solsta = solve_mosek(prob, 1e-10)
        
    # test that (x1, x2, x3) = y[2:4] is globally optimal 
    @test abs(t - Polyopt.evalpoly(f, y[2:4])) < 1e-6
    @test all( [ Polyopt.evalpoly(g, y[2:4]) for g=[g1, g2, g3, 1-x1, x1, x2, 1-x3, x3] ] .> -1e-3 )
end

# determine if f(x,z) is can be written as a SOS
let
    x, z = variables(["x", "z"])
    f = 2*x^4 + 2*x^3*z - x^2*z^2 + 5*z^4

    prob = momentprob(2, f)
    X, t, y, solsta = solve_mosek(prob)
    @test t > -1e-6

    v = monomials(2, [x,z])
    @test approxzero( f - dot(v,X[1]*v) )
end

# in this case f(x,z) is not SOS, but f(x,z)-t is SOS
let
    x, z = variables(["x", "z"])

    f = (x + z + 1)^2 - 1
    prob = momentprob(1, f)

    X, t, y, solsta = solve_mosek(prob)
    @test abs(t + 1) < 1e-6

    v = monomials(1, [x,z])

    # test that (f(x,z)-t) is SOS
    @test approxzero( f - t - dot(v,X[1]*v) )
end

# test duality for problem with constraints
let
    x, z = variables(["x", "z"])

    f = x^4 + z^2 + 1
    g = x^2 - z^2 - 1
    prob = momentprob(2, f, [g])

    X, t, y, solsta = solve_mosek(prob)

    v1 = monomials(2,[x,z])
    v2 = monomials(1,[x,z])

    # test that f(x,z) - t - g(x,z)*s(x,z) is SOS,  where s(x,z) = v2'*X[2]*v2 is SOS
    @test approxzero( f - t - dot(v1,X[1]*v1) - g*dot(v2,X[2]*v2) )
end

