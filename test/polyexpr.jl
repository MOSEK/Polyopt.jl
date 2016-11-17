let
    x, y = variables(["x", "y"])
    
    #@test (x*y).alpha == [1 1]
    #@test (x*y).c == [1]

    #@test (1 + 2*x*y).alpha == [0 0; 1 1]
    #@test (1 + 2*x*y).c == [1; 2]

    #@test ( (x+y)*(x-y) ).alpha == [2 0; 0 2]
    #@test ( (x+y)*(x-y) ).c == [1; -1]

    #@test ( (1+x+y)^2 * y ).deg == 3
    #@test ( (1+x+y)^2 * y ).alpha == [0 1; 1 1; 0 2; 2 1; 1 2; 0 3]
    #@test ( (1+x+y)^2 * y ).c == [1; 2; 2; 1; 2; 1]

    @test ( x*(x^2-y) == x^3 - x*y)
    @test ( (x^2-y)*x == x^3 - x*y)
    @test ( (x^2-y)*(x+y+x*y) == x^3 + x^2*y + x^3*y - x*y - y^2 - x*y^2)
    @test ( (x+y+x*y)*(x^2-y) == x^3 + x^2*y + x^3*y - x*y - y^2 - x*y^2)

    @test ( [x, y]*[x, y]' == [ x^2 x*y; x*y y^2] )
    
    @test ( [1 0; 0 2]*[x, y] == [x, 2*y])
    @test ( dot([x,y], [x,y]) == 1.0*x^2 + y^2 )

    @test ( dot([1//2*x, y, 1], [2, 1, 1]) == x + y + 1 )

    @test ( x + (Polyopt.Poly(0) + Polyopt.Poly(0)) == x )
end

let
    x, = variables(["x"])
    @test ( (1+x)^2 == 1+2*x+x^2 )
end

let
    x, y, z = variables(["x", "y", "z"])
    @test Polyopt.monomials(3,[x,y,z]) == sort([1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2,x^3,x^2*y,x^2*z,x*y^2,x*y*z,x*z^2,y^3,y^2*z,y*z^2,z^3])
end
