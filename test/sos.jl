using Base.Test

let

    using Polyopt
    using Sos

    eps = 10e-5;

    function check_decomp(f,gramMatrix,basis)
        errPoly = f-(basis'*gramMatrix*basis)[1];
        norm(errPoly.c) < eps
    end

    x = variables(["x","y","z"])

    f = (1.0*x[1])^4+(x[1])^2 + 2.0*x[1]^2   
    (isSos,gramMatrix,basis) = issos(f);
    @test check_decomp(f,gramMatrix,basis);
    @test isSos

    f = x[1]^4+x[2]^4+x[3]^4-4*x[1]*x[2]*x[3]+x[1]+x[2]+x[3]+3
    (isSos,gramMatrix,basis) = issos(f,eps);
    @test check_decomp(f,gramMatrix,basis);
    @test isSos

    f = (1.0*x[1])^4+(x[1])^2 + 2.0*x[1]^2-10  
    (isSos,gramMatrix,basis) = issos(f);
    @test !isSos

end






