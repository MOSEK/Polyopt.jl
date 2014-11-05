using Mosek

# The 'prob' structure specifies the moment problem
#
# minimize    prob.obj'*y
# subject to  prob.mom[k]*y is PSD ,  k=1,...,length(prob.mom)
#             y[1] = 1
#
# We formulate the dual problem for MOSEK
#
# maximize    -sum_j dot(prob.mom[j][:,1], Xj)
# subject to   sum_j dot(prob.mom[j][:,i], Xj) = prob.obj[i],  i=2,...,length(prob.obj)
#              Xj is PSD, j=1,...,length(prob.mom)
function solve_mosek(prob::MomentProb, scaling=false)

    printstream(msg::String) = print(msg)

    # Create a task object and attach log stream printer
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,printstream)

    # The momemt problem is in dual form, so we dualize it
    numcon = length(prob.obj) - 1
    numbarvar = length(prob.mom)
    barvardim = Int[ sqrt(size(prob.mom[k],1)) for k=1:numbarvar ]

    eqdim = Int[ sqrt(size(prob.eq[k],1)) for k=1:length(prob.eq) ]
    eqidx = Array(Int, length(prob.eq)+1)
    eqidx[1] = 0
    for k=1:length(prob.eq)
        eqidx[k+1] = eqidx[k] + eqdim[k]*(eqdim[k]+1)>>1
    end
    numvar = eqidx[ end ]

    # Append 'numcon' empty constraints.
    appendcons(task, int32(numcon))

    # find norm of each row in (A, b)
    norm_rows = ones(numcon)
    if scaling
        for i=1:numcon
            norm_rows[i] = max(norm_rows[i], abs(prob.obj[i+1]))
            for j=1:length(prob.mom)
                norm_rows[i] = max(norm_rows[i], norm(prob.mom[j].nzval[prob.mom[j].colptr[i+1]:prob.mom[j].colptr[i+2]-1], Inf))
            end

            for j=1:length(prob.eq)
                norm_rows[i] = max(norm_rows[i], norm(prob.eq[j].nzval[prob.eq[j].colptr[i+1]:prob.eq[j].colptr[i+2]-1], Inf))
            end
        end
    end

    # add free variables from equality constraints
    if numvar > 0
        appendvars(task, int32(numvar))

        putvarboundslice(task, 1, numvar+1,
                         [ MSK_BK_FR::Int32 for i in 1:numvar ],
                         [ -Inf             for i in 1:numvar ],
                         [ +Inf             for i in 1:numvar ])

        for j=1:length(prob.eq)
            k = prob.eq[j].colptr[1]:prob.eq[j].colptr[2]-1
            subj = trilind( prob.eq[j].rowval[k], eqdim[j] ) + eqidx[j]
            putclist(task, subj, -float64(prob.eq[j].nzval[k])/norm_rows[subj])
        end

        for j=1:length(prob.eq)
            for i=1:numcon
                k = prob.eq[j].colptr[i+1]:prob.eq[j].colptr[i+2]-1
                subj = trilind( prob.eq[j].rowval[k], eqdim[j] ) + eqidx[j]
                putaijlist(task, i*ones(Int, length(subj)), subj, float64(prob.eq[j].nzval[k])/norm_rows[i])
            end
        end
    end

    # Append matrix variables of sizes in 'BARVARDIM'.
    appendbarvars(task, int32(barvardim))

    bkc = Int32[ MSK_BK_FX for k=1:numcon ]
    blc = float64(prob.obj[2:end])./norm_rows
    buc = float64(prob.obj[2:end])./norm_rows

    # Set the bounds on constraints.
    putconboundslice(task, 1, numcon+1, bkc, blc, buc)

    # Add objective
    for j=1:numbarvar
        nj = int64(barvardim[j])
        k = prob.mom[j].colptr[1]:prob.mom[j].colptr[2]-1
        subk, subl = ind2sub( (nj, nj), prob.mom[j].rowval[k] )
        cj = appendsparsesymmat(task,
                                int32(nj),
                                int32(subk),
                                int32(subl),
                                float64(-prob.mom[j].nzval[k]))
        putbarcj(task, j, [cj], [1.0])
    end

    # Add constraints
    for i=1:numcon
        for j=1:numbarvar
            nj = int64(barvardim[j])
            k = prob.mom[j].colptr[i+1]:prob.mom[j].colptr[i+2]-1
            subk, subl = ind2sub( (nj, nj), prob.mom[j].rowval[k] )
            aij = appendsparsesymmat(task, int32(nj), int32(subk), int32(subl), float64(prob.mom[j].nzval[k])/norm_rows[i])
            putbaraij(task, int32(i), int32(j), [aij], [1.0])
        end
    end

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    # Write .task file
    #writetask(task, "polyopt.task")

    # Solve the problem and print summary
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    # Get status information about the solution
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    #return getbarxj(task, MSK_SOL_ITR, 1)

    if solsta == MSK_SOL_STA_OPTIMAL
        return ([1, gety(task, MSK_SOL_ITR)], getprimalobj(task, MSK_SOL_ITR), "Optimal")
    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
        return ([1, gety(task, MSK_SOL_ITR)], getprimalobj(task, MSK_SOL_ITR), "Near optimal")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
        error("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
        return ([1, gety(task, MSK_SOL_ITR)],  getprimalobj(task, MSK_SOL_ITR), "Unknown")
    else
        error("Other solution status")
    end

end

trilind(k::Vector{Int}, n::Int) = Int[i + (j-1)*(n-1) - (j-1)*(j-2)>>1 for (i,j) = zip(ind2sub((n, n), k)...) ]
