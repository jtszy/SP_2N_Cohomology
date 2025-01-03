using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS ÷ 2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS ÷ 2)

using Groups
using IntervalArithmetic
using JuMP
using LowCohomologySOS
using PermutationGroups
using SCS
using Serialization
using SP_4_Cohomology
using SparseArrays
using SymbolicWedderburn

n = 3

Sp_N = MatrixGroups.SymplecticGroup{2*n}(Int8)
F_Sp_N_Steinberg = FreeGroup(alphabet(Sp_N))

function extended_f_sp_2n(n::Integer)
    range_as_list = [i for i in 1:n]
    ordered_pairs = [(i,j) for i ∈ 1:n for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    unordered_pairs = [(i,j) for i ∈ 1:n for j ∈ deleteat!(copy(range_as_list), findall(j->j<i,copy(range_as_list)))]
    
    x_y_z_gens = vcat(
        [(:x,i,j) for (i,j) in ordered_pairs],
        [(:y,i,j) for (i,j) in unordered_pairs],
        [(:yt,i,j) for (i,j) in unordered_pairs],
        [(:z,i,j) for (i,j) in unordered_pairs],
        [(:zt,i,j) for (i,j) in unordered_pairs],
        [(:X,i,j) for (i,j) in ordered_pairs],
        [(:Y,i,j) for (i,j) in unordered_pairs],
        [(:Yt,i,j) for (i,j) in unordered_pairs],
        [(:Z,i,j) for (i,j) in unordered_pairs],
        [(:Zt,i,j) for (i,j) in unordered_pairs],
    )

    gens_number = div(length(x_y_z_gens),2)
    alphabet = Alphabet(x_y_z_gens, vcat([i+gens_number for i in 1:gens_number],[i for i in 1:gens_number]))

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
        elseif type == :y
            new_type, new_i, new_j = :B, max(i,j), min(i,j) + n
        elseif type == :yt
            new_type, new_i, new_j = :B, min(i,j) + n, max(i,j)
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

# below code to change

for i in eachindex(S)
    @assert quotient_hom_Steinberg(gens(F_Sp_N_Steinberg,i)) == S[i]
    @assert quotient_hom_Steinberg(gens(F_Sp_N_Steinberg,i)^(-1)) == S[i]^(-1)
end

Solution = Dict()
Solution["lambda_list"] = []
Solution["Q_list"] = []
Solution["result_list"] = []


# support_jacobian, min_support = SP_4_Cohomology.symplectic_min_supports(quotient_hom_Steinberg, S; rels = "adj")

# Steinberg_relations = SP_4_Cohomology.relations_St(F_Sp_N_Steinberg, S, N; sq_adj_ = "adj")

support_jacobian, min_support = SP_4_Cohomology.symplectic_min_supports(quotient_hom_Steinberg, S)

Steinberg_relations = SP_4_Cohomology.relations_St(F_Sp_N_Steinberg, S, N)

for r in Steinberg_relations
    @assert quotient_hom_Steinberg(r) == one(Sp_N)
end

Δ₁, I_N, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom_Steinberg, Steinberg_relations, support_jacobian);

Δm_mono, Δm_sq, Δm_adj_mi, Δm_adj_db, Δm_op  = SP_4_Cohomology.mono_sq_adj_op(Δ₁⁻, S)

I_N_mono, I_N_sq = SP_4_Cohomology.mono_sq_adj_op(I_N, S)

# Δ = Δm_adj_db + Δ₁⁺
Δ = Δm_sq+ Δm_adj_mi+ Δm_adj_db+ Δm_op + Δ₁⁺
Δ = Δ₁

RG = LowCohomologySOS.group_ring(Sp_N, min_support, star_multiplication = true)

Δ = LowCohomologySOS.embed.(identity, Δ, Ref(RG))
# I_N = LowCohomologySOS.embed.(identity, I_N_sq, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N, Ref(RG))

constraints_basis, psd_basis, Σ, action = SP_4_Cohomology.wedderburn_data(RG.basis, min_support, S);

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Δ, σ, action.alphabet_perm, S) == Δ
    @assert LowCohomologySOS.act_on_matrix(I_N, σ, action.alphabet_perm, S) == I_N
end

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    sos_problem, P = LowCohomologySOS.sos_problem(
        Δ,
        I_N,
        w_dec_matrix
        # 0.7 / 0.05
    )
end

# Find a numerical spectral gap
JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-6, max_iters = 500))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)

# Q, λ = deserialize("./Steinberg_Solution_Sp_6.sjl")
# Q = Q[2]
# λ = λ[2]

result_bool, _ = LowCohomologySOS.certify_sos_decomposition(Δ, I_N, λ, Q, min_support)


Solution["lambda_list"] = push!(Solution["lambda_list"], λ)
Solution["Q_list"] = push!(Solution["Q_list"], Q)
Solution["result_list"] = push!(Solution["result_list"], result_bool)



# serialize("./Steinberg_Solution_Sp_6_comm_relations.sjl", Solution)