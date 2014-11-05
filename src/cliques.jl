using Base.LinAlg.CHOLMOD

for Ti in (:Int32,:Int64)
    @eval begin
        function chm_analyze_ordering{Tv<:CHOLMOD.CHMVTypes}(A::CHOLMOD.CholmodSparse{Tv,$Ti},
                                                             ordering::Cint,
                                                             Perm::Vector{$Ti})
            Parent   = Array($Ti,A.c.m)
            Post     = Array($Ti,A.c.m)
            ColCount = Array($Ti,A.c.m)
            First    = Array($Ti,A.c.m)
            Level    = Array($Ti,A.c.m)

            ccall((@CHOLMOD.chm_nm("analyze_ordering",$Ti),:libcholmod), Cint,
                  (Ptr{CHOLMOD.c_CholmodSparse{Tv,$Ti}}, Cint, Ptr{$Ti}, Ptr{$Ti}, Csize_t,
                   Ptr{$Ti}, Ptr{$Ti}, Ptr{$Ti},
                   Ptr{$Ti}, Ptr{$Ti}, Ptr{Uint8}),
                   &A.c, ordering, Perm, C_NULL, zero(Csize_t), Parent, Post, ColCount, First, Level, CHOLMOD.cmn($Ti))
            Parent, Post, ColCount
        end

        function chm_analyze_p{Tv<:CHOLMOD.CHMVTypes}(A::CholmodSparse{Tv,$Ti}, Perm::Vector{$Ti})
            F = ccall((@CHOLMOD.chm_nm("analyze_p",$Ti),:libcholmod), Ptr{CHOLMOD.c_CholmodFactor{Tv,$Ti}},
                      (Ptr{CHOLMOD.c_CholmodSparse{Tv,$Ti}}, Ptr{$Ti}, Ptr{$Ti}, Csize_t, Ptr{Uint8}),
                      &A.c, Perm, C_NULL, zero(Csize_t), CHOLMOD.cmn($Ti))
            CholmodFactor(F)
        end

        function chm_factorize_p!{Tv<:CHOLMOD.CHMVTypes}(A::CholmodSparse{Tv,$Ti}, beta::Tv, L::CholmodFactor{Tv,$Ti})
             ccall((@CHOLMOD.chm_nm("factorize_p",$Ti),:libcholmod), Cint,
                  (Ptr{CHOLMOD.c_CholmodSparse{Tv,$Ti}}, Ptr{Tv}, Ptr{$Ti}, Csize_t,
                   Ptr{CHOLMOD.c_CholmodFactor{Tv,$Ti}}, Ptr{Uint8}),
                   &A.c, &beta, C_NULL, zero(Csize_t), &L.c, CHOLMOD.cmn($Ti))
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
    S = CholmodSparse(SparseMatrixCSC(m,n,A.colptr,A.rowval,zeros(length(A.nzval))))
    F = chm_analyze_p(S, length(Perm)>0 ? Perm-1 : Perm)
    chm_factorize_p!(S, 1.0, F)

    par, post, colcount = chm_analyze_ordering(S, F.c.ordering, F.Perm)
    ns, flag = pothen_sun(par+1, post+1, colcount)

    # extract cliques
    L = sparse(F)
    [ sort(1+F.Perm[L.rowval[L.colptr[i]:L.colptr[i+1]-1]]) for i = find(flag .< 0) ]
end

# If no permutation is specified, use the CHOLMOD ordering
chordal_embedding{Tv<:Number,Ti<:Int}(A::SparseMatrixCSC{Tv,Ti}) = chordal_embedding(A, Array(Int,0))

function clique_index(I::Array{Array{Int,1},1}, S::Array{Int,1})
    s = IntSet(S...)
    for k=1:length(I)
        if s <= IntSet(I[k]...) return k end
    end
    throw(ArgumentError)
end

