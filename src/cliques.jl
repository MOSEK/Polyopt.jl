using Base.SparseArrays.CHOLMOD

function chm_analyze_ordering{Tv<:CHOLMOD.VTypes}(A::CHOLMOD.Sparse{Tv}, ordering::Cint, Perm::Vector{Int})
    s = unsafe_load(A.p)
    Parent   = Array(Int,s.nrow)
    Post     = Array(Int,s.nrow)
    ColCount = Array(Int,s.nrow)
    First    = Array(Int,s.nrow)
    Level    = Array(Int,s.nrow)

    ccall((@CHOLMOD.cholmod_name("analyze_ordering",Int),:libcholmod), Cint,
          (Ptr{CHOLMOD.C_Sparse{Tv}}, Cint, Ptr{Int}, Ptr{Int}, Csize_t,
           Ptr{Int}, Ptr{Int}, Ptr{Int},
           Ptr{Int}, Ptr{Int}, Ptr{UInt8}),
           A.p, ordering, Perm, C_NULL, zero(Csize_t), Parent, Post, ColCount, First, Level, CHOLMOD.common())
    Parent, Post, ColCount
end

function pothen_sun{T<:Integer}(par::Vector{T}, post::Vector{T}, colcount::Vector{T})
    n = length(par)
    flag = -ones(T,n)
    ns = n
    for j=post
        mdeg = colcount[j] - 1
        if par[j] != 0 && mdeg == colcount[par[j]] && flag[par[j]] == -1
            # par[j] not assigned to supernode
            ns -= 1
            if flag[j] < 0    # j is a repr. vertex
                flag[par[j]] = j
                flag[j] -= 1
            else              # j is not a repr. vertex
                flag[par[j]] = flag[j]
                flag[flag[j]] -= 1
            end
        end
    end
    ns, flag
end

function chordal_embedding{Tv<:Number,Ti<:Int}(A::SparseMatrixCSC{Tv,Ti}, Perm::Vector{Ti})
    m, n = size(A)
    S = CHOLMOD.Sparse(round(Float64,A))
    
    F = cholfact(round(Float64,A), shift=n, perm=Perm)    
    s = unsafe_load(F.p)
    p = round(Int64,[ unsafe_load(s.Perm, i) for i=1:s.n])

    par, post, colcount = chm_analyze_ordering(S, s.ordering, p)
    ns, flag = pothen_sun(par+1, post+1, colcount)

    # extract cliques
    L = sparse(CHOLMOD.Sparse(F))
    Array{Int,1}[ sort(1+p[L.rowval[L.colptr[i]:L.colptr[i+1]-1]]) for i = find(flag .< 0) ]
end

# If no permutation is specified, use the CHOLMOD ordering
chordal_embedding{Tv<:Number,Ti<:Int}(A::SparseMatrixCSC{Tv,Ti}) = chordal_embedding(A, Array(Int,0))

function clique_index(I::Array{Array{Int,1},1}, S::Array{Int,1}, n::Int)

    J = Int[]
    mask = zeros(Int, n)
    s = sparsevec(S, 1, n)    
    for k=1:length(I)
        if length(S) <= length(I[k]) 
            mask[I[k]] = 1    
            if dot(mask, s) == length(S)
                push!(J, k)
            end
            mask[I[k]] = 0
        end
    end
    J
end
