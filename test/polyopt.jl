approxzero{T<:Number}(p::Polyopt.Poly{T}, threshold=1e-6) = (norm(p.c, Inf) < threshold)

# example of unconstrained global optimization
let
    x, z = variables(["x", "z"])
    f = 4*x^2 + x*z - 4*z^2 - 21//10*x^4 + 4*z^4 + 1//3*x^6

    prob = momentprob(3, f)
    X, t, y, solsta = solve_mosek(prob)

    @test abs(t - (-1.0316)) < 1e-4
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

# Motzkin polynomial
# XXX: MOSEK cannot solve this instance
#let
#    x, z = variables(["x", "z"])
#
#    f = x^2*z^4 + x^4*z^2 - 3*x^2*z^2 + 1
#
#    prob = momentprob(3, f)
#
#    X, t, y, solsta = solve_mosek2(prob)
#end
