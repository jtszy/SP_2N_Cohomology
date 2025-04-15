using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS ÷ 2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS ÷ 2)

using Groups
using JuMP
using LowCohomologySOS
using PermutationGroups
using SCS
using Serialization
using SP_4_Cohomology
using SparseArrays
using SymbolicWedderburn


# Load the paremeters
N - parse(Integer, ARGS[1])
precomputed = parse(Bool, ARGS[2])
sq_adj_all = parse(String, ARGS[3])
# N = 3
# precomputed = false
# quotient_flag = true

# Define Sp₂ₙ(Z) and the quotient homomorphism on it from the free group on Sp₂ₙ(Z)'s gens
Sp_2N = MatrixGroups.SymplecticGroup{2*N}(Int8)
F_Sp_2N_Steinberg = FreeGroup(alphabet(Sp_2N))
S = gens(Sp_2N)
quotient_hom_Steinberg = let source = F_Sp_2N_Steinberg, target = Sp_2N
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end
for i in eachindex(S)
    @assert quotient_hom_Steinberg(gens(F_Sp_2N_Steinberg,i)) == S[i]
    @assert quotient_hom_Steinberg(gens(F_Sp_2N_Steinberg,i)^(-1)) == S[i]^(-1)
end

# Compute the relations as words in the free group on gens of Sp₂ₙ(Z)
quotient_flag = (sq_adj_all == "all" ? false : true)
Steinberg_relations = SP_4_Cohomology.relations_St(
    F_Sp_2N_Steinberg, 
    S, 
    N,
    quotient_flag = quotient_flag
)

# Compute the Laplacian of interest
support_jacobian, min_support = SP_4_Cohomology.symplectic_min_supports(
    Steinberg_relations, 
    quotient_hom_Steinberg, 
    S
)
Δ₁, I_N_whole, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(
    quotient_hom_Steinberg, 
    Steinberg_relations, 
    support_jacobian
);
Δm_mono, Δm_sq, Δm_adj, Δm_op  = SP_4_Cohomology.mono_sq_adj_op(Δ₁⁻, S)
if sq_adj_all == "adj"
    laplacian = Δm_adj + Δ₁⁺
    descr = "_H_"*string(N)*"_adj.sjl"
elseif sq_adj_all == "sq"
    laplacian = Δm_sq + Δ₁⁺
    descr = "_H_"*string(N)*"_sq.sjl"
else
    laplacian = Δ₁
    descr = "_Sp_"*string(2*N)*"_delta_1.sjl"
end
RG = LowCohomologySOS.group_ring(Sp_2N, min_support, star_multiplication = true)
laplacian = LowCohomologySOS.embed.(identity, laplacian, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N_whole, Ref(RG))

# Either compute or load the numerical solution of the corresponding SDP
if !precomputed
    # Wedderburn symmetrization
    constraints_basis, psd_basis, Σ, action = SP_4_Cohomology.wedderburn_data(
        RG.basis, 
        min_support, 
        S
    )
    for σ in Σ
        @assert LowCohomologySOS.act_on_matrix(laplacian, σ, action.alphabet_perm, S) == laplacian
        @assert LowCohomologySOS.act_on_matrix(I_N, σ, action.alphabet_perm, S) == I_N
    end
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(
        Float64, 
        Σ, 
        action, 
        constraints_basis, 
        psd_basis
    )

    # Define and solve the corresponding SDP problem
    if N == 2
        upper_bound_ = Inf
    elseif sq_adj_all != "adj"
        upper_bound_ = 1.0
    else
        upper_bound_ = 0.25
    end
    sos_problem, P = LowCohomologySOS.sos_problem(
        laplacian,
        I_N,
        w_dec_matrix,
        upper_bound_
    )
    JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-6, max_iters = 200_000))
    JuMP.optimize!(sos_problem)
    λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)

    # Serialize the solution
    solution = Dict("lambda" => λ, "Q" => Q)
    serialize("./scripts/Steinberg_Solution"*descr, solution)
else
    # Load the precomputed numerical solution
    solution = deserialize("./scripts/Steinberg_Solution"*descr)
    λ, Q = solution["lambda"], solution["Q"]
end

# Turn into a rigorous proof - certify the numerical estimate
certified_flag, certified_interval = SP_4_Cohomology.certify_sos_decomposition(laplacian, I_N, λ, Q, min_support)
if certified_flag == true
    @info "Spectral gap fits within the interval:" certified_interval
end
