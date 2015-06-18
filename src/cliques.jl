using Base.SparseMatrix.CHOLMOD

for Ti in CHOLMOD.IndexTypes
    @eval begin
        function chm_analyze_ordering{Tv<:CHOLMOD.VTypes}(A::CHOLMOD.Sparse{Tv,$Ti},
                                                          ordering::Cint,
                                                          Perm::Vector{$Ti})
            s = unsafe_load(A.p)
            Parent   = Array($Ti,s.nrow)
            Post     = Array($Ti,s.nrow)
            ColCount = Array($Ti,s.nrow)
            First    = Array($Ti,s.nrow)
            Level    = Array($Ti,s.nrow)

            ccall((@CHOLMOD.cholmod_name("analyze_ordering",$Ti),:libcholmod), Cint,
                  (Ptr{CHOLMOD.C_Sparse{Tv,$Ti}}, Cint, Ptr{$Ti}, Ptr{$Ti}, Csize_t,
                   Ptr{$Ti}, Ptr{$Ti}, Ptr{$Ti},
                   Ptr{$Ti}, Ptr{$Ti}, Ptr{Uint8}),
                   A.p, ordering, Perm, C_NULL, zero(Csize_t), Parent, Post, ColCount, First, Level, CHOLMOD.common($Ti))
            Parent, Post, ColCount
        end

        function chm_analyze_p{Tv<:CHOLMOD.VTypes}(A::CHOLMOD.Sparse{Tv,$Ti}, Perm::Vector{$Ti})
            F = Factor(ccall((@CHOLMOD.cholmod_name("analyze_p",$Ti),:libcholmod), Ptr{CHOLMOD.C_Factor{Tv,$Ti}},
                             (Ptr{CHOLMOD.C_Sparse{Tv,$Ti}}, Ptr{$Ti}, Ptr{$Ti}, Csize_t, Ptr{Uint8}),
                              A.p, Perm, C_NULL, zero(Csize_t), CHOLMOD.common($Ti)))
            F
        end

        function chm_factorize_p!{Tv<:CHOLMOD.VTypes}(A::CHOLMOD.Sparse{Tv,$Ti}, beta::Tv, L::Factor{Tv,$Ti})
             ccall((@CHOLMOD.cholmod_name("factorize_p",$Ti),:libcholmod), Cint,
                  (Ptr{CHOLMOD.C_Sparse{Tv,$Ti}}, Ptr{Tv}, Ptr{$Ti}, Csize_t,
                   Ptr{CHOLMOD.C_Factor{Tv,$Ti}}, Ptr{Uint8}),
                   A.p, &beta, C_NULL, zero(Csize_t), L.p, CHOLMOD.common($Ti))
        end

    end
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
    
    Perm = [1:n;]
    F = chm_analyze_p(S, length(Perm)>0 ? Perm-1 : Perm)
    chm_factorize_p!(S, 1.0, F)
    
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

function clique_index(I::Array{Array{Int,1},1}, S::Array{Int,1})

    s = IntSet(S)
    for k=1:length(I)
        if s <= IntSet(I[k]) return k end
    end
    throw(ArgumentError)
end

