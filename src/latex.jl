texstring{T<:Number}(a::T) = string(a)
texstring(a::Float64) = string(@sprintf("%1.2f",a))
texstring(a::Rational) = ( a.den == one(a) ? string(a.num) : string("\\frac{", a.num, "}{", a.den, "}") )

# convert a polynomial to tex format, e.g.,  2 x_1 x_2 + x_3
function latex(p::Poly)
    if p.m == 0 return("0") end

    s = ""
    for i=1:p.m
        if i==1
            s *= (p.c[i]>=0 ? "" : "- ")
        else
            s *= (p.c[i]>=0 ? "+ " : "- ")
        end

        ci = abs(p.c[i])
        if sum(p.alpha[i,:]) == 0
            s *= texstring(ci)
        else
            if ci != one(ci) s *= texstring(ci) * " " end
            for k=1:p.n
                if p.alpha[i,k] > 0
                    s *= p.syms.names[k]
                    if p.alpha[i,k] > 1
                        s *= "^{$(p.alpha[i,k])} "
                    else
                        s *= " "
                    end
                end
            end
        end
    end
    s
end

# convert a polynomial to vectorized tex format, e.g.,  2 y_{110} + y_{001}
function latex(p::Poly, imap::Dict{Array{Int,2},Int}, texsymbols::Array{String,1})
    if p.m == 0 return("0") end

    s = ""
    for i=1:p.m
        if i==1
            s *= (p.c[i]>=0 ? "" : "- ")
        else
            s *= (p.c[i]>=0 ? "+ " : "- ")
        end

        ci = abs(p.c[i])
        if sum(p.alpha[i,:]) == 0
            s *= texstring(ci) * " "
        else
            if ci != one(ci) s *= texstring(ci) * " " end
            s *= texsymbols[imap[p.alpha[i,:]]] * " "
        end
    end
    s
end

function latex{T}(M::Array{Poly{T}}, imap::Dict{Array{Int,2},Int}, texsymbols::Array{String,1}, matrixlimit = 10)
    n = size(M,1)
    s = "\\left[\n\\begin{array}{*{$(min(matrixlimit,n))}c}\n"
    if matrixlimit < n
        for i=1:matrixlimit-2
            s *= join([ [ latex(m, imap, texsymbols) for m=M[i,1:matrixlimit-2] ],
                            ["\\dots", latex(M[i,n], imap, texsymbols)] ], " & ") * "\\\\\n"
        end
        s *= join([ ["\\vdots" for k=1:matrixlimit-2], "\\ddots", "\\vdots\\\\\n"], " & ")

        s *= join([ [ latex(m, imap, texsymbols) for m=M[n,1:matrixlimit-2] ],
                        ["\\dots", latex(M[n,n], imap, texsymbols)] ], " & ") * "\\\\\n"
    else
        for i=1:n
            s *= join([ latex(m, imap, texsymbols) for m=M[i,:]], " & ") * "\\\\\n"
        end
    end
    string(s, "\\end{array}\n\\right]\n")
end

function latex_dual(prob::MomentProb, imap::Dict{Array{Int,2},Int}, texsymbols::Array{String,1}, matrixlimit = 10)

    s = string("Dual problem:\n\\[\n\\begin{array}{ll}\n\\text{minimize} &",
               latex(dot(prob.obj, prob.basis), imap, texsymbols), "\\\\\n")
    s *= "\\text{subject to} "

    for k=1:length(prob.mom)
        nk = int(sqrt(size(prob.mom[k],1)))
        Mk = reshape(prob.basis'*prob.mom[k]', (nk, nk))
        Mk = Mk + Mk' - diagm(diag(Mk))

        if nk > 1
            s = string(s, " & ", latex(Mk, imap, texsymbols, matrixlimit), "\\succeq 0\\\\\n\n")
        else
            s = string(s, " & ", latex(Mk[1], imap, texsymbols), " \\geq 0\\\\\n")
        end
    end

    for k=1:length(prob.eq)
        nk = int(sqrt(size(prob.eq[k],1)))
        Mk = reshape(prob.basis'*prob.eq[k]', (nk, nk))
        Mk = Mk + Mk' - diagm(diag(Mk))

        if nk > 1
            s = string(s, " & ", latex(Mk, imap, texsymbols, matrixlimit), " = 0\\\\\n\n")
        else
            s = string(s, " & ", latex(Mk[1], imap, texsymbols), " = 0\\\\\n")
        end
    end

    string(s, "\\end{array}\n\\]\n")
end

function latex_dual(prob::MomentProb, vectorizednames::Bool, matrixlimit::Int)

    imap = Polyopt.indexmap(prob.basis)
    syms = Array(String, length(imap))

    if vectorizednames == true
        for (key,k) = imap
            syms[k] = "y_{" * join([string(i) for i=key]) * "}"
        end
    else
        for k=1:length(imap)
            syms[k] = latex(prob.basis[k])
        end
    end

    latex_dual(prob, imap, syms, matrixlimit)
end

function latex_probstats(prob::MomentProb)
    sdpblksize  = [ int(sqrt(size(m,1))) for m=prob.mom ]
    numsdpblk   = length(sdpblksize)
    maxblksize  = maximum(sdpblksize)
    avgblksize  = int(mean(sdpblksize))

    eqsize      = [ int(sqrt(size(m,1))) for m=prob.eq ]
    numfree     = sum(Int[ n*(n+1) >> 1 for n=eqsize ])
    numvar      = sum([ n*(n+1) >> 1 for n=sdpblksize ]) + numfree
    numcon      = length(prob.basis)-1

    maxrhs  = maximum(prob.obj)
    maxcoef = zeros(length(prob.mom) + length(prob.eq))
    for k=1:length(prob.mom)
        maxcoef[k] = norm(prob.mom[k].nzval, Inf)
    end

    for k=1:length(prob.eq)
        maxcoef[k+numsdpblk] = norm(prob.eq[k].nzval, Inf)
    end

    s = "\\begin{quote}\n\\begin{tabular}{ll}\n"
    s = string(s,"Problem statistics\\\\\n\\hline\n")
    s = string(s,"constraints & $(numcon) \\\\\n")
    s = string(s,"sdp blocks & $(numsdpblk) \\\\\n")
    s = string(s,"largest block size & $(maxblksize) \\\\\n")
    s = string(s,"avg. block size & $(avgblksize) \\\\\n")
    s = string(s,"free variables & $(numfree) \\\\\n")
    s = string(s,"total number of scalar vars& $(numvar) \\\\\n")
    s = string(s, @sprintf("largest matrix coeff & %1.1e \\\\", maximum(maxcoef)))
    s = string(s, @sprintf("largest rhs & %1.1e\\\\\n", maxrhs))
    s = string(s, "\\hline\n")
    string(s,"\\end{tabular}\n\\end{quote}\n\n")
end

function latex_header()
    s = "\\documentclass{standalone}\n"
    s = string(s,"\\usepackage{varwidth}\n")
    s = string(s,"\\usepackage[paperwidth=\\maxdimen,paperheight=\\maxdimen]{geometry}\n")
    s = string(s,"\\usepackage{amsmath}\n")
    s = string(s,"\\begin{document}\n")
    string(s,"\\begin{varwidth}{\\linewidth}\n")
end

function latex_footer()
    s = "\\end{varwidth}\n"
    string(s, "\\end{document}")
end

function latex(texfile::String, prob::MomentProb, header="", body="", vectorizednames=false, matrixlimit=10)

    f = open(texfile,"w")

    tex = string(latex_header(),
                 header,
                 latex_probstats(prob),
                 body,
                 latex_dual(prob, vectorizednames, matrixlimit),
                 latex_footer())

    print(f, tex)
    close(f)
end
