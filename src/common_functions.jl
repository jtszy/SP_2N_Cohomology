function quotient_homomorphism(free_group, sp_n, gen_dict)
    function F(i, source, target)
        if source([i]) == one(free_group)
            return one(sp_n)
        end
        for gen in keys(gen_dict)
            if source([i]) == gen
                return word(gen_dict[gen])
            elseif source([i]) == gen^(-1)
                return word(gen_dict[gen]^(-1))
            end
        end
    end
    Groups.Homomorphism(F, free_group, sp_n)
end

function com(x,y)
    x*y*x^(-1)*y^(-1)
end

function minimalistic_support(relations, quotient_hom)
    jacobian_matrix = LowCohomologySOS.jacobian_matrix(relations)

    RF_G = parent(first(jacobian_matrix))
    Sp_N = quotient_hom.target

    min_support = union!(
        [one(Sp_N)],
        gens(Sp_N),
        inv.(gens(Sp_N))
    )

    for entry in jacobian_matrix
        for i in SparseArrays.nonzeroinds(entry.coeffs)
            push!(min_support, quotient_hom(RF_G.basis[i]))
        end
    end

    min_support = unique!([min_support; inv.(min_support)])
    return min_support
end

function support_jacobian(relations, quotient_hom)
    Sp_N = quotient_hom.target
    F_G = quotient_hom.source

    support_jacobian = [one(Sp_N)]
    for r in relations
        current_factor = one(Sp_N)
        for i in 1:length(word(r))
            g = quotient_hom(F_G(word(r)[i:i]))
            push!(support_jacobian,current_factor*g)
            current_factor *= g
        end
    end
    support_jacobian = unique(support_jacobian)

    return support_jacobian
end

function symplectic_min_supports(
    quotient_hom,
    S;
    rels = "all"
)   
    Sp_N = quotient_hom.target
    F_Sp_N_Steinberg = quotient_hom.source

    N = div(size(MatrixGroups.matrix_repr(first(S)))[1],2)

    Steinberg_relations = relations_St(F_Sp_N_Steinberg, S, N, sq_adj_ = rels)

    for r in Steinberg_relations
        print(r)
        @assert quotient_hom(r) == one(Sp_N)
    end 

    sup_jacobian = support_jacobian(vcat(Steinberg_relations, S), quotient_hom)

    min_support = minimalistic_support(Steinberg_relations, quotient_hom)

    return sup_jacobian, min_support

end

function sq_adj_op(
    Δ₁⁻,
    S # gens of Sp_N without z
)
    RG = parent(first(Δ₁⁻))
    Sp2N = parent(first(RG.basis))
    N = Int8(sqrt(length(gens(Sp2N))/2))

    sq_pairs = []
    adj_pairs = []
    op_pairs = []

    M = length(S)
    A = alphabet(Sp2N)
    
    function index_pair(s::Int)
        if s <= M
            s_i, s_j = mod(A[word(S[s])[1]].i, N), mod(A[word(S[s])[1]].j, N)
        elseif s <= Int(3*M/2)
            s_i, s_j = Int(floor((s - 2*N^2 + 2*N) / (N-1))) + 1, mod(s - 1 - 2*N^2 + 2*N, N-1) + 1 + 
                Int(Int(floor((s - 2*N^2 + 2*N) / (N-1))) + 1 <= mod(s - 1 - 2*N^2 + 2*N, N-1))
        else
            s_i, s_j = Int(floor((s - 3*N^2 + 3*N) / (N-1))) + 1, mod(s - 1 - 3*N^2 + 3*N, N-1) + 1 +
                Int(Int(floor((s - 3*N^2 + 3*N) / (N-1))) + 1 <= mod(s - 1 - 3*N^2 + 3*N, N-1) )
        end
        return s_i, s_j
    end
    # needs to be checked

    for s in 1:2*M
        for t in 1:2*M
            s_i, s_j = index_pair(s)
            t_i, t_j = index_pair(t)
            if sort([s_i,s_j]) == sort([t_i, t_j])
                push!(sq_pairs,(s,t))
            elseif length(intersect!([s_i,s_j],[t_i,t_j])) == 1
                push!(adj_pairs,(s,t))
            else
                push!(op_pairs,(s,t))
            end
        end
    end

    sq = [(i,j) in sq_pairs ? Δ₁⁻[i,j] : zero(RG) for i in 1:2*M, j in 1:2*M]
    adj = [(i,j) in adj_pairs ? Δ₁⁻[i,j] : zero(RG) for i in 1:2*M, j in 1:2*M]
    op = [(i,j) in op_pairs ? Δ₁⁻[i,j] : zero(RG) for i in 1:2*M, j in 1:2*M]

    @assert sq+adj+op == Δ₁⁻

    return sq, adj, op
end