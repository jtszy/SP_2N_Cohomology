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

function mono_sq_adj_op(
    Δ₁⁻,
    S # generating set indexing Δ₁⁻
)
    RG = parent(first(Δ₁⁻))
    Sp2N = parent(first(RG.basis))
    N = Int8(sqrt(length(gens(Sp2N))/2))
    mono_pairs = []
    sq_pairs = []
    adj_pairs = []
    op_pairs = []
    A = alphabet(Sp2N)
    for s in eachindex(S)
        for t in eachindex(S)
            s_i, s_j = mod(A[word(S[s])[1]].i,N), mod(A[word(S[s])[1]].j,N)
            t_i, t_j = mod(A[word(S[t])[1]].i,N), mod(A[word(S[t])[1]].j,N)
            if length(union!([s_i,s_j],[t_i,t_j])) == 1
                push!(mono_pairs,(s,t))
            elseif length(union!([s_i,s_j],[t_i,t_j])) == 2
                push!(sq_pairs,(s,t))
            elseif length(union!([s_i,s_j],[t_i,t_j])) == 3
                push!(adj_pairs,(s,t))
            else
                push!(op_pairs,(s,t))
            end
        end
    end
    mono = [(i,j) in mono_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    sq = [(i,j) in sq_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    adj = [(i,j) in adj_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    op = [(i,j) in op_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]

    @assert mono+sq+adj+op == Δ₁⁻

    return mono, sq, adj, op
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
    Steinberg_relations,
    quotient_hom,
    S
)   
    sup_jacobian = SP_4_Cohomology.support_jacobian(vcat(Steinberg_relations, S), quotient_hom)
    min_support = SP_4_Cohomology.minimalistic_support(Steinberg_relations, quotient_hom)

    return sup_jacobian, min_support
end
