using Mosek
using LinearAlgebra


function ind2sub_replacement(a, b)
    indices = CartesianIndices(a)[b]
    subk = Array{Int32}(undef, length(indices))
    subl = Array{Int32}(undef, length(indices))
    for (i, ind) in enumerate(indices)
        subk[i], subl[i] = Tuple(ind)
    end
    return (subk, subl)
end

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
function solve_mosek(prob::MomentProb; tolrelgap=1e-10, showlog=true)

    printstream(msg::AbstractString) = print(msg)

    # Create a task object and attach log stream printer
    task = maketask()
    if showlog  putstreamfunc(task,MSK_STREAM_LOG,printstream)  end

    # The momemt problem is in dual form, so we dualize it
    numcon = length(prob.obj)
    numbarvar = length(prob.mom)
    barvardim = Int[ sqrt(size(prob.mom[k],1)) for k=1:numbarvar ]
    
    eqdim = Int[ sqrt(size(prob.eq[k],1)) for k=1:length(prob.eq) ]
    eqidx = Array{Int}(undef, length(prob.eq)+1)
    eqidx[1] = 0
    for k=1:length(prob.eq)
        eqidx[k+1] = eqidx[k] + eqdim[k]*(eqdim[k]+1)/2
    end
    numvar = 1 + eqidx[ end ]
    
    # add free variables from equality constraints
    appendvars(task, Int32(numvar))

    putvarboundslice(task, 1, numvar+1,
                     [ MSK_BK_FR  for i in 1:numvar ],
                     [ -Inf       for i in 1:numvar ],
                     [ +Inf       for i in 1:numvar ])

    putcj(task, 1, 1.0)

    # Append matrix variables of sizes in 'BARVARDIM'.
    appendbarvars(task, barvardim)

    # Add constraints
    numconst = 1
    appendcons(task, 1)
    putaij(task, 1, 1, 1.0)
    
    I = Array{Int}(undef,0)
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
                    
                subk, subl = ind2sub_replacement((nj,nj), prob.mom[j].rowval[k1:k2])
                    
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

                subk, subl = ind2sub_replacement((eqdim[j], eqdim[j]), prob.eq[j].rowval[k1:k2])
                trilidx = subk .>= subl
                subj = trilind( prob.eq[j].rowval[k1:k2][trilidx], eqdim[j] ) .+ (eqidx[j] + 1)

                putaijlist(task, numconst*ones(Int, length(subj)), subj, prob.eq[j].nzval[k1:k2][trilidx])
#                 subj = trilind( prob.eq[j].rowval[k1:k2], eqdim[j] ) + eqidx[j] + 1
#                 putaijlist(task, numconst*ones(Int, length(subj)), subj, prob.eq[j].nzval[k1:k2])
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
#    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
#        return (X, Z, t, y, "Near optimal")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        return (X, Z, t, y, "Dual infeasibility")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        return (X, Z, t, y, "Primal infeasibility")
#    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
#        return (X, Z, t, y, "Near dual infeasibility")
#    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
#        return (X, Z, t, y, "Near primal infeasibility")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        return (X, Z, t, y, "Unknown")
    else
        error("Other solution status")
    end
end

trilind(k::Vector{Int}, n::Int) = Int[i + (j-1)*(n-1) - (j-1)*(j-2)/2 for (i,j) = zip(ind2sub_replacement((n, n),k)...) ]

function solve_mosek_no_elim(prob::BSOSProb; tolrelgap=1e-10, showlog=true)
    
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
    
    l = 1
    for k=1:length(El)
        l += nnz(El[k])
    end    
    const_idx = Array{Int,1}(undef, l)
    const_idx[1] = 1
    l = 1
    for k=1:length(El)
        const_idx[l + (1:nnz(El[k]))] = El[k].rowval
        l += nnz(El[k])
    end        
    const_idx = unique(const_idx)
#     
#     const_idx = Int[ 1 ]
#     for i=1:length(Al)
#         push!(const_idx, prob.El[i].rowval...)
#     end
#     const_idx = unique(const_idx)
    const_map = Dict(zip(const_idx, collect(1:length(const_idx))))

    numvar = sum(diml) + sum(dimf) + 1
    numcon = sum(dimf) + length(const_idx)
    numbarvar = m
            
    appendvars(task, numvar)

    # bounds on lj 
    offs = 0
    for l=1:length(diml)
        for j=1:diml[l]
            if prob.lb[l][j] == -Inf
                putvarbound(task, offs+j, MSK_BK_FR, -Inf, Inf)            
                #println("VARBOUND($(offs+j)): FREE")
            else
                putvarbound(task, offs+j, MSK_BK_LO, 0.0, Inf)
                #println("VARBOUND($(offs+j)): LOWER 0.0")
            end
        end
        offs += diml[l]
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
     
    idx_const = 1
    
    for j=1:m
        #println("Block $(j)")
        for i=1:dimf[j]
            k1, k2 = Al[j].colptr[i], Al[j].colptr[i+1]-1
            subj = [idx_l[j]-1 + Al[j].rowval[k1:k2]; idx_f[j]-1 + i]            
            val  = [Al[j].nzval[k1:k2]; -1.0]                
            putarow(task, idx_const, subj, val)    
            #println("CONSTRAINT A($(idx_const)): $(subj), $(val)")
            
            k1, k2 = As[j].colptr[i], As[j].colptr[i+1]-1
            if k2>=k1
                subk, subl = ind2sub_replacement( (dimX[j], dimX[j]), As[j].rowval[k1:k2] )
                I = subk .>= subl            
                aij = appendsparsesymmat(task, dimX[j], subk[I], subl[I], As[j].nzval[k1:k2][I])
                #println("BARAIJ($(idx_const+i-1),$(j)): $(subk[I]), $(subl[I]), $(As[j].nzval[k1:k2][I])")
                putbaraij(task, idx_const, j, [aij], [1.0])
            end
                        
            putconbound(task, idx_const, MSK_BK_FX, 0.0, 0.0)                        
            idx_const += 1            
        end
        
    end
    
    putaij(task, idx_const, numvar, 1.0)
    for k=1:length(El)
        for (j,i) in enumerate(prob.El[k].rowval)
            #println("putaij: $(idx_const-1+const_map[i]), $(idx_f[k]-1+j), $(prob.El[k].nzval[j])")
            putaij(task, idx_const-1+const_map[i], idx_f[k]-1+j, prob.El[k].nzval[j])
        end            
    end
    
    for j=idx_const:numcon
        putconbound(task, j, MSK_BK_FX, 0.0, 0.0)
    end
    
    for (j,i) in enumerate(f.nzind)
        #println("putconbound: $(idx_const-1+const_map[i]), $(f.nzval[j])")        
        putconbound(task, idx_const-1+const_map[i], MSK_BK_FX, f.nzval[j], f.nzval[j])    
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
    if solsta == MSK_SOL_STA_OPTIMAL
        return (X, t, l, f, y, "Optimal")
#    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
#        return (X, t, l, f, y, "Near optimal")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        return (X, t, l, f, y, "Dual infeasibility")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        return (X, t, l, f, y, "Primal infeasibility")
#    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
#        return (X, t, l, f, y, "Near dual infeasibility")
#    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
#        return (X, t, l, f, y, "Near primal infeasibility")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        return (X, t, l, f, y, "Unknown")
    else
        error("Other solution status")
    end
end


function solve_mosek(prob::BSOSProb; tolrelgap=1e-10, showlog=true)
 
    function getrow(Ak, Ek, j)    
        sub, val = Int[], Float64[]
        i = findfirst(Ek.rowval, j)
        if i>0
            k1, k2 = Ak.colptr[i], Ak.colptr[i+1]-1
            sub = view(Ak.rowval, k1:k2)
            val = view(Ak.nzval, k1:k2)
        else
            sub, val = Int[], Float64[]
        end
        sub, val                
    end
    
    function constraint_size(A, E, j)
        sz = (j == 1 ? 1 : 0)
        for k=1:length(A)
            sz += length(getrow(A[k], E[k], j)[1])            
        end
        sz
    end

    printstream(msg::AbstractString) = print(msg)

    # Create a task object and attach log stream printer
    task = maketask()
    if showlog  putstreamfunc(task,MSK_STREAM_LOG,printstream)  end
    
    Al, El, As = prob.Al, prob.El, prob.As
    
    m = length(Al)
    diml = [ size(a,1) for a in Al ]
    dimX = Int[ round(Int, sqrt(size(as,1))) for as in As ]

    # variables are indexed as (l1,...,lm,t)
    idx_l = zeros(Int, length(diml)+1)
    idx_l[1] = 1
    for j=1:length(diml)
        idx_l[j+1] = idx_l[j] + diml[j]
    end
        
    l = 1
    for k=1:length(El)
        l += nnz(El[k])
    end    
    const_idx = Array{Int,1}(undef, l)
    const_idx[1] = 1
    l = 1
    for k=1:length(El)
        const_idx[l + (1:nnz(El[k]))] = El[k].rowval
        l += nnz(El[k])
    end    
    
    const_idx = sort(unique(const_idx))
    
    nonzero_markers = falses(length(const_idx))
    for (i,k)=enumerate(const_idx)

        if (nonzero_markers[i]) continue end        
    
        if i==1
            nonzero_markers[i] = true
            continue
        end
        
        for j=1:m
            subj, valj = getrow(Al[j], El[j], k)                
            if length(subj) > 0
                nonzero_markers[i] = true
                break
            end
        end
        
        if (nonzero_markers[i]) continue end
            
        for j=1:m
            subj, valj = getrow(As[j], El[j], k)
            if length(subj) >0
                nonzero_markers[i] = true
                break            
            end
        end 
    end   
    const_idx = const_idx[ nonzero_markers ]

    numvar = sum(diml) + 1
    numcon = length(const_idx)
    numbarvar = m
            
    appendvars(task, numvar)

    # bounds on lj 
    offs = 0
    for l=1:m
        for j=1:diml[l]
            if prob.lb[l][j] == -Inf
                putvarbound(task, offs+j, MSK_BK_FR, -Inf, Inf)
                #println("VARBOUND($(offs+j)): FREE")            
            else
                putvarbound(task, offs+j, MSK_BK_LO, 0.0, Inf)
                #println("VARBOUND($(offs+j)): LOWER 0.0")
            end
        end
        offs += diml[l]
    end
    
    # t free
    putvarbound(task, numvar, MSK_BK_FR, -Inf, Inf)
    
    putcj(task, numvar, 1.0)

    # Append matrix variables
    appendbarvars(task, dimX)

    # Add constraints
    appendcons(task, numcon)
     
    f = prob.obj[const_idx]
       
    for (i,k)=enumerate(const_idx)    
        offs = 0
        for j=1:m
            subj, valj = getrow(Al[j], El[j], k)                
            if length(subj) > 0
                for r=1:length(subj)
                    putaij(task, i, subj[r]+offs, valj[r])
                end
            end
            offs += size(Al[j],1)
        end
        if i==1
            putaij(task, i, numvar, 1.0)
        end
        putconbound(task, i, MSK_BK_FX, f[i], f[i])                        
       
        for j=1:m
            subj, valj = getrow(As[j], El[j], k)
            if length(subj) >0
                subk, subl = ind2sub_replacement( (dimX[j], dimX[j]), subj )
                I = subk .>= subl            
                aij = appendsparsesymmat(task, dimX[j], subk[I], subl[I], valj[I])
                putbaraij(task, i, j, [aij], [1.0])
            end
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

    X = [ symm(getbarxj(task, MSK_SOL_ITR, j), Int(sqrt(size(As[j],1)))) for j=1:m ]    
    l = [ getxxslice(task, MSK_SOL_ITR, idx_l[j], idx_l[j+1]) for j=1:m ]        
    t = getxxslice(task, MSK_SOL_ITR, numvar, numvar+1)[1]
    y = sparsevec(const_idx, gety(task, MSK_SOL_ITR), length(prob.obj))
    if solsta == MSK_SOL_STA_OPTIMAL
        return (X, t, l, y, "Optimal")
    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
        return (X, t, l,  y, "Near optimal")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        return (X, t, l,  y, "Dual infeasibility")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        return (X, t, l,  y, "Primal infeasibility")
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
        return (X, t, l,  y, "Near dual infeasibility")
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
        return (X, t, l,  y, "Near primal infeasibility")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        return (X, t, l,  y, "Unknown")
    else
        error("Other solution status")
    end   
end


function symm(x::Array{T,1}, n::Int) where {T<:Number}
    X = zeros(n,n)
    k = 0
    for j=1:n
        X[j:n,j] = x[k .+ (1:n-j+1)]
        k += n-j+1
    end

    Matrix(Symmetric(X,:L))
end

