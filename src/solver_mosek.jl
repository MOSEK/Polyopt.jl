using Mosek

# The 'prob' structure specifies the moment problem
#
# minimize    prob.obj'*y
# subject to  prob.mom[k]*y is PSD,  k=1,...,length(prob.mom)
#             prob.eq[k]*y = 0,      k=1,...,length(prob.eq) 
#             y[1] = 1
#
# We formulate the dual problem for MOSEK
#
# maximize     t
# subject to   sum_j dot(prob.mom[j][:,1], Xj) + sum_k dot(prb.eq[k][:,1], Zk)) = prob.obj[1] - t
#              sum_j dot(prob.mom[j][:,i], Xj) + sum_k dot(prb.eq[k][:,1], Zk)) = prob.obj[i], i=2,...,length(prob.obj)
#              Xj is PSD, j=1,...,length(prob.mom)
#              Zk is symmetric but free,  k=1,...,length(prob.eq)   
#
function solve_mosek(prob::MomentProb, tolrelgap=1e-10; showlog=true)

    printstream(msg::String) = print(msg)

    # Create a task object and attach log stream printer
    task = maketask()
    if showlog  putstreamfunc(task,MSK_STREAM_LOG,printstream)  end

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
    
    # add free variables from equality constraints
    appendvars(task, Int32(numvar))

    putvarboundslice(task, 1, numvar+1,
                     round(Int32, [ MSK_BK_FR  for i in 1:numvar ]),
                     [ -Inf       for i in 1:numvar ],
                     [ +Inf       for i in 1:numvar ])

    putcj(task, 1, 1.0)

    # Append matrix variables of sizes in 'BARVARDIM'.
    appendbarvars(task, round(Int32,barvardim))

    # Add constraints
    numconst = 1
    appendcons(task, 1)
    putaij(task, 1, 1, 1.0)
    
    for i=1:numcon

        added_const = ( i == 1);
        
        for j=1:numbarvar
            nj = Int64(barvardim[j])            
            k1, k2 = prob.mom[j].colptr[i], prob.mom[j].colptr[i+1]-1
            if k2 >= k1
            
                if !added_const
                    appendcons(task, 1)
                    added_const = true
                end                
                subk, subl = ind2sub( (nj, nj), prob.mom[j].rowval[k1:k2] )
                aij = appendsparsesymmat(task, Int32(nj), round(Int32,subk), round(Int32,subl), map(Float64,prob.mom[j].nzval[k1:k2]))
                putbaraij(task, Int32(numconst), Int32(j), [aij], [1.0])
            end
        end
        
        for j=1:length(prob.eq)
            k1, k2 = prob.eq[j].colptr[i], prob.eq[j].colptr[i+1]-1
            if k2 >= k1
                if !added_const
                    appendcons(task, 1)
                    added_const = true
                end                
                
                subj = trilind( prob.eq[j].rowval[k1:k2], eqdim[j] ) + eqidx[j] + 1
                putaijlist(task, numconst*ones(Int, length(subj)), subj, map(Float64,prob.eq[j].nzval[k1:k2]))
            end
        end
        
        if added_const
            putconbound(task, numconst, MSK_BK_FX, prob.obj[i], prob.obj[i])
            numconst = numconst + 1
        end
    end

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    putparam(task, "MSK_IPAR_NUM_THREADS", "8")
    putparam(task, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", string(tolrelgap))
    
    # Write .task file
    writetask(task, "polyopt.task")

    # Solve the problem and print summary
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    # Get status information about the solution
    solsta = getsolsta(task,MSK_SOL_ITR)

    X = [ symm(getbarxj(task, MSK_SOL_ITR, j), Int(sqrt(size(prob.mom[j],1)))) for j=1:length(prob.mom) ]
    t = getxxslice(task, MSK_SOL_ITR, 1, 2)[1]
    Z = [ symm(getxxslice(task, MSK_SOL_ITR, 2+eqidx[k], 2+eqidx[k+1]), eqdim[k]) for k=1:length(prob.eq) ]
    Z = [ 0.5*(Zk + diagm(diag(Zk))) for Zk = Z ]

    y = gety(task, MSK_SOL_ITR)
    if solsta == MSK_SOL_STA_OPTIMAL
        return (X, Z, t, y, "Optimal")
    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
        return (X, Z, t, y, "Near optimal")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        return (X, Z, t, y, "Dual infeasibility")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        return (X, Z, t, y, "Primal infeasibility")
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
        return (X, Z, t, y, "Near dual infeasibility")
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
        return (X, Z, t, y, "Near primal infeasibility")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        return (X, Z, t, y, "Unknown")
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

