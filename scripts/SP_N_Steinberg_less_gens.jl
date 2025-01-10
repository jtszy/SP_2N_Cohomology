using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS ÷ 2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS ÷ 2)

using Groups
# using IntervalArithmetic
using JuMP
using LowCohomologySOS
using PermutationGroups
using SCS
using Serialization
using SP_4_Cohomology
using SparseArrays
using StarAlgebras
using SymbolicWedderburn

n = 3

Sp_N = MatrixGroups.SymplecticGroup{2*n}(Int8)

function extended_f_sp_2n(n::Integer)
    range_as_list = [i for i in 1:n]
    ordered_pairs = [(i,j) for i ∈ 1:n for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    unordered_pairs = [(i,j) for i ∈ 1:n for j ∈ deleteat!(copy(range_as_list), findall(j->j<=i,copy(range_as_list)))]
    
    x_z_gens = vcat(
        [(:x,i,j) for (i,j) in ordered_pairs],
        [(:z,i,j) for (i,j) in unordered_pairs],
        [(:zt,i,j) for (i,j) in unordered_pairs],
        [(:X,i,j) for (i,j) in ordered_pairs],
        [(:Z,i,j) for (i,j) in unordered_pairs],
        [(:Zt,i,j) for (i,j) in unordered_pairs],
    )

    gens_number = div(length(x_z_gens),2)
    alphabet = Alphabet(x_z_gens, vcat([i+gens_number for i in 1:gens_number],[i for i in 1:gens_number]))

    return FreeGroup(alphabet)
end

Extended_f_sp_2n = extended_f_sp_2n(n)

S = gens(Sp_N)
gen_dict_Sp_2n = Dict(LowCohomologySOS.determine_letter(S[i]) => i for i in eachindex(S))

function quotient_hom_gens(i, F, G)
    if i <= length(gens(F))
        (type,i,j) = LowCohomologySOS.determine_letter(gens(F,i))

        if type == :x
            new_type, new_i, new_j = :A, i, j
        elseif type == :z
            new_type, new_i, new_j = :B, i, i + n
        else # type :zt
            new_type, new_i, new_j = :B, i + n, i
        end

        sp_2n_i_index = gen_dict_Sp_2n[MatrixGroups.ElementarySymplectic{2*n}(new_type,new_i,new_j)]

        return Groups.word_type(G)([sp_2n_i_index]) # seems to be the same as "word([sp_2n_i_index])"
    else
        return Groups.word_type(G)(word(G(quotient_hom_gens(i-length(gens(F)),F,G))^(-1)))
    end
end

quotient_hom = let source = Extended_f_sp_2n, target = Sp_N
    Groups.Homomorphism((i, F, G) -> quotient_hom_gens(i, F,G ), source, target)
end

com(a,b) = SP_4_Cohomology.com(a,b)

function relations_St(
    F_G::Groups.FreeGroup,
    N::Integer;
    sq_adj = "all"
)
    gen_dict = Dict(LowCohomologySOS.determine_letter(gens(F_G, i)) => gens(F_G, i) for i in eachindex(gens(F_G)))

    range_as_list = [i for i in 1:N]
    unordered_pairs = [(i,j) for i ∈ 1:n for j ∈ deleteat!(copy(range_as_list), findall(j->j<=i,copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ range_as_list
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    unordered_triples_ij = [(i,j,k) for i ∈ range_as_list
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j<=i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    unordered_triples_ik = [(i,j,k) for i ∈ range_as_list
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k<=i,copy(range_as_list))) 
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
    unordered_triples_jk = [(i,j,k) for j ∈ range_as_list
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k<=j,copy(range_as_list))) 
                    for i ∈ deleteat!(copy(range_as_list), findall(i->k∈[j,k],copy(range_as_list)))]
    
    x(i,j) = gen_dict[(:x,i,j)]
    z(i,j) = gen_dict[(:z,i,j)]
    zt(i,j) = gen_dict[(:zt,i,j)]

    relations_sq = vcat(
        [z(i,j)*zt(i,j)^(-1)*z(i,j)*z(i,j)*zt(i,j)^(-1)*z(i,j)*z(i,j)*zt(i,j)^(-1)*z(i,j)*z(i,j)*zt(i,j)^(-1)*z(i,j) for (i,j) in unordered_pairs]
    )
    
    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j), x(i,k)) for (i,j,k) in triples],
        [com(x(i,j), x(k,j)) for (i,j,k) in triples],

        [com(z(i,j), x(j,k)) for (i,j,k) in unordered_triples_ij],
        [com(zt(i,j), x(j,k)) for (i,j,k) in unordered_triples_ij],

        [com(x(i,j), z(i,k)) for (i,j,k) in unordered_triples_ik],
        [com(z(i,k), x(j,k)) for (i,j,k) in unordered_triples_ik],
        [com(zt(i,k), x(j,k)) for (i,j,k) in unordered_triples_ik],

        [com(x(i,j), zt(j,k)) for (i,j,k) in unordered_triples_jk]
    )

    if sq_adj == "sq" 
        return relations_sq
    elseif sq_adj == "adj"
        return relations_adj
    elseif sq_adj == "all"
        return vcat(relations_sq, relations_adj)
    end
end

function symplectic_min_supports(
    Steinberg_relations,
    quotient_hom,
    S
)   
    Sp_N = quotient_hom.target

    for r in Steinberg_relations
        @assert quotient_hom(r) == one(Sp_N)
    end 

    sup_jacobian = SP_4_Cohomology.support_jacobian(vcat(Steinberg_relations, S), quotient_hom)
    min_support = SP_4_Cohomology.minimalistic_support(Steinberg_relations, quotient_hom)

    return sup_jacobian, min_support
end

Steinberg_relations = relations_St(Extended_f_sp_2n, n, sq_adj = "all")
support_jacobian, min_support = symplectic_min_supports(Steinberg_relations, quotient_hom, S)

Δ₁, I_N, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom, Steinberg_relations, support_jacobian);

function sq_adj_op(
    Δ₁⁻
)
    sq_pairs = []
    adj_pairs = []
    op_pairs = []
    A = alphabet(Extended_f_sp_2n)
    S = gens(Extended_f_sp_2n)
    RG = parent(first(Δ₁⁻))
    for s in eachindex(S)
        for t in eachindex(S)
            s_i, s_j = A[word(S[s])[1]][2], A[word(S[s])[1]][3]
            t_i, t_j = A[word(S[t])[1]][2], A[word(S[t])[1]][3]
            if length(intersect!([s_i,s_j],[t_i,t_j])) == 2
                push!(sq_pairs,(s,t))
            elseif length(intersect!([s_i,s_j],[t_i,t_j])) == 1
                push!(adj_pairs,(s,t))
            else
                push!(op_pairs,(s,t))
            end
        end
    end
    sq = [(i,j) in sq_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    adj = [(i,j) in adj_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    op = [(i,j) in op_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]

    @assert sq+adj+op == Δ₁⁻

    return sq, adj, op
end

Δm_sq, Δm_adj, Δm_op  = sq_adj_op(Δ₁⁻)
Δ = Δm_sq + Δm_adj + Δm_op + Δ₁⁺

RG = LowCohomologySOS.group_ring(Sp_N, min_support, star_multiplication = true)

Δ = LowCohomologySOS.embed.(identity, Δ, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N, Ref(RG))

@time begin
    sos_problem = LowCohomologySOS.sos_problem(
        Δ,
        I_N
    )
end

# Find a numerical spectral gap
JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-6, max_iters = 500))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem)

result_bool, _ = LowCohomologySOS.certify_sos_decomposition(Δ, I_N, λ, Q, min_support)
