using Mosek

# The 'prob' structure specifies the moment problem
#
# minimize    prob.obj'*y
# subject to  prob.mom[k]*y is PSD ,  k=1,...,length(prob.mom)
#             y[1] = 1
#
# We formulate the dual problem for MOSEK
#
# maximize     t
# subject to   sum_j dot(prob.mom[j][:,1], Xj) = prob.obj[1] - t
#              sum_j dot(prob.mom[j][:,i], Xj) = prob.obj[i], i=2,...,length(prob.obj)
#              Xj is PSD, j=1,...,length(prob.mom)
#
function solve_mosek(prob::MomentProb)

    printstream(msg::String) = print(msg)

    # Create a task object and attach log stream printer
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,printstream)

    # The momemt problem is in dual form, so we dualize it
    numcon = length(prob.obj)
    numbarvar = length(prob.mom)
    barvardim = Int[ sqrt(size(prob.mom[k],1)) for k=1:numbarvar ]

    eqdim = Int[ sqrt(size(prob.eq[k],1)) for k=1:length(prob.eq) ]
    eqidx = Array(Int, length(prob.eq)+1)
    eqidx[1] = 0
    for k=1:length(prob.eq)
        eqidx[k+1] = eqidx[k] + eqdim[k]*(eqdim[k]+1)>>1
    end
    numvar = 1 + eqidx[ end ]

    # Append 'numcon' empty constraints.
    appendcons(task, Int32(numcon))

    # add free variables from equality constraints
    appendvars(task, Int32(numvar))

    putvarboundslice(task, 1, numvar+1,
                     round(Int32, [ MSK_BK_FR  for i in 1:numvar ]),
                     [ -Inf       for i in 1:numvar ],
                     [ +Inf       for i in 1:numvar ])

    putcj(task, 1, 1.0)
    for j=1:length(prob.eq)
        k = prob.eq[j].colptr[1]:prob.eq[j].colptr[2]-1
        subj = trilind( prob.eq[j].rowval[k], eqdim[j] ) + eqidx[j] + 1
        putclist(task, subj, -map(Float64,prob.eq[j].nzval[k]))
    end

    putaij(task, 1, 1, 1.0)

    for j=1:length(prob.eq)
        for i=1:numcon
            k = prob.eq[j].colptr[i]:prob.eq[j].colptr[i+1]-1
            subj = trilind( prob.eq[j].rowval[k], eqdim[j] ) + eqidx[j] + 1
            putaijlist(task, i*ones(Int, length(subj)), subj, map(Float64,prob.eq[j].nzval[k]))
        end
    end

    # Append matrix variables of sizes in 'BARVARDIM'.
    appendbarvars(task, round(Int32,barvardim))

    bkc = Int32[ MSK_BK_FX for k=1:numcon ]
    blc = map(Float64,prob.obj)
    buc = map(Float64,prob.obj)

    # Set the bounds on constraints.
    putconboundslice(task, 1, numcon+1, bkc, blc, buc)

    # Add constraints
    for i=1:numcon
        for j=1:numbarvar
            nj = Int64(barvardim[j])
            k = prob.mom[j].colptr[i]:prob.mom[j].colptr[i+1]-1
            subk, subl = ind2sub( (nj, nj), prob.mom[j].rowval[k] )
            aij = appendsparsesymmat(task, Int32(nj), round(Int32,subk), round(Int32,subl), map(Float64,prob.mom[j].nzval[k]))
            putbaraij(task, Int32(i), Int32(j), [aij], [1.0])
        end
    end

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    # Write .task file
    writetask(task, "polyopt.task")

    # Solve the problem and print summary
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    # Get status information about the solution
    solsta = getsolsta(task,MSK_SOL_ITR)

    if solsta == MSK_SOL_STA_OPTIMAL
        X = [ symm(getbarxj(task, MSK_SOL_ITR, j), int(sqrt(size(prob.mom[j],1)))) for j=1:length(prob.mom) ]
        t = getxxslice(task, MSK_SOL_ITR, 1, 2)[1]
        y = gety(task, MSK_SOL_ITR)
        return (X, t, y, "Optimal")
    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
        X = [ symm(getbarxj(task, MSK_SOL_ITR, j), int(sqrt(size(prob.mom[j],1)))) for j=1:length(prob.mom) ]
        t = getxxslice(task, MSK_SOL_ITR, 1, 2)[1]
        y = gety(task, MSK_SOL_ITR)
        return (X, t, y, "Near optimal")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        X = [ symm(getbarxj(task, MSK_SOL_ITR, j), int(sqrt(size(prob.mom[j],1)))) for j=1:length(prob.mom) ]
        t = getxxslice(task, MSK_SOL_ITR, 1, 2)[1]
        y = gety(task, MSK_SOL_ITR)
        return (X, t, y, "Unknown")
    else
        error("Other solution status")
    end
end

trilind(k::Vector{Int}, n::Int) = Int[i + (j-1)*(n-1) - (j-1)*(j-2)>>1 for (i,j) = zip(ind2sub((n, n), k)...) ]

function symm{T<:Number}(x::Array{T,1}, n::Int)
    X = zeros(n,n)
    k = 0
    for j=1:n
        X[j:n,j] = x[k + (1:n-j+1)]
        k += n-j+1
    end

    full(Symmetric(X,:L))
end

