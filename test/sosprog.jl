using Base.Test
using Polyopt
using Sos

function SolveAndVerify(prog,eps,fminsosRef)

    (optval,cert,s,lambda) = solve(prog)
    cnstTerm = evalpoly(prog.obj,vec(zeros(1,prog.obj.n)));
    
    fminsos =  cnstTerm + optval
    @test norm(fminsosRef-fminsos) < eps
    
    certVerify = prog.obj-fminsos
    
    (isSos,temp,temp) = issos(cert,eps)
    @test isSos
   
    for i=1:length(prog.pineq) 
        certVerify = certVerify - s[i]*prog.pineq[i].p
        (isSos,temp,temp) = issos(s[i],eps)
        @test isSos
    end
  
    (isSos,temp,temp) = issos(cert,eps)
    @test isSos

    for i=1:length(prog.peq) 
        certVerify = certVerify - lambda[i]*prog.peq[i].p
    end

    certErr = certVerify-cert
    @test norm(certErr.c) < eps
    
end

let
    
    #######################################
    eps = 10e-5;
    x = variables(["x","y","z"])
    f = x[1]*x[1]+1
    multDegree = 2;
  
    pineq = [Constraint(x[1]-1,multDegree)]; 
    peq = [Constraint(x[1]-1,multDegree)]; 
    optval = 2
    prog = sosprog(f,pineq,peq);
    
    SolveAndVerify(prog,eps,optval)
    #######################################
    
    #######################################
    eps = 10e-5;
    x = variables(["x","y","z"])
    f = x[1]*x[1]+1
    multDegree = 2;
  
    pineq = [Constraint(x[1]-1,multDegree)]; 
    peq = Array(Constraint,0);
    optval = 2
    prog = sosprog(f,pineq,peq);
    
    SolveAndVerify(prog,eps,optval)
    #######################################
   
    #######################################
    eps = 10e-5;
    x = variables(["x","y","z"])
    f = x[1]*x[1]+1
    multDegree = 2;
  
    peq = [Constraint(x[1]-1,multDegree)]; 
    pineq = Array(Constraint,0);
    optval = 2
    prog = sosprog(f,pineq,peq);
    
    SolveAndVerify(prog,eps,optval)
    #######################################
    
    #######################################
    
    multDegree = 3;
    pineq = [Constraint(1-x[1]^2-x[2]^2,multDegree)];
    peq = Array(Constraint,0);
    optval = 9
    
    prog = sosprog(10-x[1],pineq,peq)
    SolveAndVerify(prog,eps,optval)
   
    ######################################
    
    multDegree = 2;
    peq = [Constraint(1-x[1]^2-x[2]^2,multDegree)];
    pineq = Array(Constraint,0);
    optval = -1/sqrt(2)*2
    
    prog = sosprog(x[1]+x[2],pineq,peq)
    SolveAndVerify(prog,eps,optval)
   
 
    #######################################
      
end






