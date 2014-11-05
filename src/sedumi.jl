using Polyopt, MATLAB

function write(matfile::ASCIIString, prob::MomentProb)

    numcon = length(prob.obj) - 1
    barvardim = Int[ sqrt(size(prob.mom[k],1)) for k=1:length(prob.mom) ]

    eqdim = Int[ sqrt(size(prob.eq[k],1)) for k=1:length(prob.eq) ]
    eqidx = Array(Int, length(prob.eq)+1)
    eqidx[1] = 0
    for k=1:length(prob.eq)
        eqidx[k+1] = eqidx[k] + eqdim[k]*(eqdim[k]+1)>>1
    end
    numvar = eqidx[ end ]

    b = float64(prob.obj[2:end])

    # off-diagonal coefficients are scaled by 2
    At = vcat([offdiag_scale(mk) for mk = prob.mom]...)

    # add free variables from equality constraints
    subj = Array(Int,0)
    val  = Array(Float64,0)
    for j=1:length(prob.eq)
        k = prob.eq[j].colptr[1]:prob.eq[j].colptr[2]-1
        append!(subj, Solver.trilind( prob.eq[j].rowval[k], eqdim[j] ) + eqidx[j])
        append!(val, -float64(prob.eq[j].nzval[k]))
    end

    c = [sparsevec(subj, val, numvar), At[:,1]]

    subi = Array(Int,0)
    subj = Array(Int,0)
    val  = Array(Float64,0)
    for j=1:length(prob.eq)
        for i=1:numcon
            k = prob.eq[j].colptr[i+1]:prob.eq[j].colptr[i+2]-1
            append!(subi, i*ones(Int, length(k)))
            append!(subj, Solver.trilind( prob.eq[j].rowval[k], eqdim[j] ) + eqidx[j])
            append!(val, float64(prob.eq[j].nzval[k]))
        end
    end
    At = [sparse(subj,subi,val,numvar,numcon), At[:,2:end]]
    At = (At')'

    MATLAB.write_matfile(matfile;
                         c  = c,
                         At = At,
                         b  = b,
                         K  = Dict( "f" => float(numvar), "s" => float(barvardim')) )
end

function offdiag_scale(A::AbstractArray)
    n = int(sqrt(size(A,1)))
    spdiagm(vec(2 - eye(n)))*A
end

