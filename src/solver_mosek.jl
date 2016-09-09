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

    printstream(msg::AbstractString) = print(msg)

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
    
    I = Array(Int,0)
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
		        trilidx = subk .>= subl            
                aij = appendsparsesymmat(task, nj, subk[trilidx], subl[trilidx], prob.mom[j].nzval[k1:k2][trilidx])
                putbaraij(task, numconst, j, [aij], [1.0])
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
                putaijlist(task, numconst*ones(Int, length(subj)), subj, prob.eq[j].nzval[k1:k2])
            end
        end
        
        if added_const
            push!(I, i)
            putconbound(task, numconst, MSK_BK_FX, prob.obj[i], prob.obj[i])
            numconst = numconst + 1
        end
    end

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    putparam(task, "MSK_IPAR_NUM_THREADS", "4")
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
    if length(y) < length(prob.obj) y = sparsevec(I, y) end
    
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

#
# min  t.
# st.  fj - Al[j]'*lj + As[j]'*Xj = 0,  j=1,...,size(Al,1)
#      sum Ej * fj = f - t*e1
#      lj >= 0,   Xj >= 0
function solve_mosek(prob::BSOSProb, tolrelgap=1e-10; showlog=true)
    
    printstream(msg::AbstractString) = print(msg)

    # Create a task object and attach log stream printer
    task = maketask()
    if showlog  putstreamfunc(task,MSK_STREAM_LOG,printstream)  end
    
    f, Al, El, As = prob.obj, prob.Al, prob.El, prob.As
    
    m = length(Al)
    diml = [ size(a,1) for a in Al ]
    dimf = [ size(a,2) for a in Al ]    
    dimX = Int[ round(Int, sqrt(size(as,1))) for as in As ]
    
    # variables are indexed as (l1,...,lm,f1,...,fm,t)
    idx_l = zeros(Int, length(diml)+1)
    idx_l[1] = 1
    for j=1:length(diml)
        idx_l[j+1] = idx_l[j] + diml[j]
    end
    
    idx_f = zeros(Int, length(dimf)+1)
    idx_f[1] = idx_l[end]
    for j=1:length(dimf)
        idx_f[j+1] = idx_f[j] + dimf[j]
    end    
    
    numvar = sum(diml) + sum(dimf) + 1
    numcon = sum(dimf) + length(f)
    numbarvar = m
    
    #println("diml:", diml)
    #println("dimf:", dimf)
    #println("dimX:", dimX)
    #println("NUMVAR=", numvar)
    #println("NUMBARVAR=", numbarvar)
    #println("DIMBARVAR=", dimX)
    #println("NUMCON=", numcon)
        
    appendvars(task, numvar)

    # lj >= 0
    for j=1:sum(diml)
        putvarbound(task, j, MSK_BK_LO, 0.0, Inf)
        #println("VARBOUND($(j)): LOWER 0.0")
    end 
    
    # fj and t free
    for j=sum(diml)+(1:sum(dimf)+1)
        putvarbound(task, j, MSK_BK_FR, -Inf, Inf)
        #println("VARBOUND($(j)): FREE")
    end 

    putcj(task, numvar, 1.0)

    # Append matrix variables
    appendbarvars(task, dimX)

    # Add constraints
    appendcons(task, numcon)

    for i=1:sum(dimf)
        putconbound(task, i, MSK_BK_FX, 0.0, 0.0)
        #println("CONBOUND A($(i)): FX 0.0")
    end
    
    for i=1:length(f)
        putconbound(task, sum(dimf)+i, MSK_BK_FX, f[i], f[i])
        #println("CONBOUND B($(sum(dimf)+i)): FX $(f[i])")
    end

    idx_const = 1
    for j=1:m
        #println("Block $(j)")
        for i=1:size(Al[j],2)
            k1, k2 = Al[j].colptr[i], Al[j].colptr[i+1]-1        
            #if k1>=k2 && i>size(As[j],2)
            #    println("XXX:  Constraint can be simplified!");
            #end
            
            subj = [idx_l[j]-1 + Al[j].rowval[k1:k2]; idx_f[j]-1 + i] 
            val  = [Al[j].nzval[k1:k2]; -1.0]                
            putarow(task, idx_const, subj, val)    
            #println("CONSTRAINT A($(idx_const)): $(subj), $(val)")
            idx_const += 1            
        end
    end

    Et = vcat(El...)
    for i=1:size(Et,2)
        k1, k2 = Et.colptr[i], Et.colptr[i+1]-1
        #if k1>=k2
        #    println("XXX:  Constraint can be simplified!");
        #end
        subj = idx_f[1]-1 + Et.rowval[k1:k2]
        val  = Et.nzval[k1:k2]
        if i==1
            push!(subj, numvar)
            push!(val, 1.0)
        end
        #println("CONSTRAINT B($(idx_const)): $(subj), $(val)")
        putarow(task, idx_const, subj, val)    
        idx_const += 1
    end
              
    idx_const = 1
    for j=1:m
        nj = dimX[j]            
        for i=1:size(As[j],2)
            k1, k2 = As[j].colptr[i], As[j].colptr[i+1]-1
            subk, subl = ind2sub( (nj, nj), As[j].rowval[k1:k2] )
            I = subk .>= subl            
            aij = appendsparsesymmat(task, nj, subk[I], subl[I], As[j].nzval[k1:k2][I])
            #println("BARAIJ($(idx_const+i-1),$(j)): $(subk[I]), $(subl[I]), $(As[j].nzval[k1:k2][I])")
            putbaraij(task, idx_const+i-1, j, [aij], [1.0])
        end
    idx_const += dimf[j]
    end
    
    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    putparam(task, "MSK_IPAR_NUM_THREADS", "4")
    putparam(task, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", string(tolrelgap))
    
    # Write .task file
    writetask(task, "polyopt.task")

    # Solve the problem and print summary
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    # Get status information about the solution
    solsta = getsolsta(task,MSK_SOL_ITR)

    X = [ symm(getbarxj(task, MSK_SOL_ITR, j), Int(sqrt(size(As[j],1)))) for j=1:m ]    
    l = [ getxxslice(task, MSK_SOL_ITR, idx_l[j], idx_l[j+1]) for j=1:m ]        
    f = [ getxxslice(task, MSK_SOL_ITR, idx_f[j], idx_f[j+1]) for j=1:m ]         
    t = getxxslice(task, MSK_SOL_ITR, numvar, numvar+1)[1]
    y = gety(task, MSK_SOL_ITR)
    X, t, l, f, y    
end

function symm{T<:Number}(x::Array{T,1}, n::Int)
    X = zeros(n,n)
    k = 0
    for j=1:n
        X[j:n,j] = x[k + (1:n-j+1)]
        k += n-j+1
    end

    full(Symmetric(X,:L))
end

