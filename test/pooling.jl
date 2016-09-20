let
    order = 2

    # Haverly1 standard pooling instance
    x14, x24, x35, x36, x45, x46, w4 = variables(["x_{14}"; "x_{24}"; "x_{35}"; "x_{36}"; "x_{45}"; "x_{46}"; "w_4"])

    c14, c24, c35, c36, c45, c46 = [6; 16; 1; -9; -5; -15]
    q1, q2, q3, q5, q6 = [ 3; 1; 2; 5//2; 3//2 ]
    ub1, ub2, ub3, ub4, ub5, ub6 = [ 300; 300; 300; 300; 100; 200 ] // 300 

    flow_eq    = [ x14+x24-(x45+x46) ]
    blend_eq   = [ w4*(x45+x46)-(q1*x14+q2*x24) ]
    qual_bnd   = [ q5*(x45+x36)-(w4*x45+q3*x35); q6*(x46+x36)-(w4*x46+q3*x36) ]
    cap_bnd    = [ ub1-x14; ub2-x24; ub3-(x35+x36); ub4-(x45+x46); ub5-(x35+x45); ub6-(x36+x46) ]
    flow_bnd   = [ x14; x24; x35; x36; x45; x46 ]

    obj = dot([c14; c24; c35; c36; c45; c46], [x14; x24; x35; x36; x45; x46])
    prob = momentprob(order, obj, [qual_bnd; cap_bnd; flow_bnd], [flow_eq; blend_eq])

    X, Z, t, y, solsta = solve_mosek(prob)
    xo = Polyopt.vectorize([x14,x24,x35,x36,x45,x46,w4],2*order)*y

    @test norm(xo-[0; 1/3; 0; 1/3; 0; 1/3; 1]) < 1e-4
    
    probc = momentprob_chordalembedding(order, obj, [qual_bnd; cap_bnd; flow_bnd], [flow_eq; blend_eq])
    Xc, Zc, tc, yc, solstac = solve_mosek(probc)
    xo = Polyopt.vectorize([x14,x24,x35,x36,x45,x46,w4],2*order)*yc

    @test norm(xo-[0; 1/3; 0; 1/3; 0; 1/3; 1]) < 1e-4    
end


