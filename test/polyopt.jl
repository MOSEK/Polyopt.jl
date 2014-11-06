using Polyopt
using Base.Test

# unconstrained problem
let
    x, y = variables(["x", "y"])
    obj = 4*x^2 + x*y - 4*y^2 - 21//10*x^4 + 4*y^4 + 1//3*x^6

    prob = momentprob(3, obj)
    X, y, objval, solsta = solve_mosek(prob)

    @test abs(objval - (-1.0316)) < 1e-4
end

# constrained problems
