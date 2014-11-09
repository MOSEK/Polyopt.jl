using Sos
using Polyopt



x = variables(["x","y"])
f = (1.0*x[1])^4+(x[1])^2 + 2.0*x[1]^2   
issos(f,10e-5)

f = (1.0*x[1])^4+(x[1])^2 + 2.0*x[1]^2-10  
issos(f,10e-5)

x = variables(["x","y","z"])

f = x[1]^4+x[2]^4+x[3]^4-4*x[1]*x[2]*x[3]+x[1]+x[2]+x[3]+3
issos(f,10e-5)
