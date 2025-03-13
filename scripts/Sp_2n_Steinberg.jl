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


# sq_adj = parse(String, ARGS[1])
# precomputed = parse(String, ARGS[2])
N = 3
sq_adj = "adj"
precomputed = false

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
Steinberg_relations = SP_4_Cohomology.relations_St(
    F_Sp_2N_Steinberg, 
    S, 
    N
)
# for r in Steinberg_relations
#     @assert quotient_hom_Steinberg(r) == one(Sp_2N)
# end

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
laplacian = (sq_adj == "adj" ? Δm_adj + Δ₁⁺ : Δm_sq + Δ₁⁺)
RG = LowCohomologySOS.group_ring(Sp_2N, min_support, star_multiplication = true)
laplacian = LowCohomologySOS.embed.(identity, laplacian, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N_whole, Ref(RG))

# Either compute or load the numerical solution of the corresponding SDP
descr = (sq_adj == "adj" ? "_Sp_6_adj.sjl" : "_Sp_6_sq.sjl")
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
    sos_problem, P = LowCohomologySOS.sos_problem(
        laplacian,
        I_N,
        w_dec_matrix,
        0.003
    )
    JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-6, max_iters = 200_000))
    JuMP.optimize!(sos_problem)
    λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)

    # Serialize the solution
    solution = Dict("lambda" => λ, "Q" => Q)
    serialize("./scripts/Steinberg_Solution"*descr, solution)
else
    # Load the precomputed numerical solution
    # solution = deserialize("./scripts/Steinberg_Solution"*descr)
    # λ, Q = solution["lambda"], solution["Q"]
end

# Turn into a rigorous proof - certify the numerical estimate
SP_4_Cohomology.certify_sos_decomposition(laplacian, I_N, λ, Q, min_support)


# function set_optimal_start_values(model::Model)
#     # Store a mapping of the variable primal solution
#     variable_primal = Dict(x => value(x) for x in all_variables(model))
#     # In the following, we loop through every constraint and store a mapping
#     # from the constraint index to a tuple containing the primal and dual
#     # solutions.
#     constraint_solution = Dict()
#     for (F, S) in list_of_constraint_types(model)
#         # We add a try-catch here because some constraint types might not
#         # support getting the primal or dual solution.
#         try
#             for ci in all_constraints(model, F, S)
#                 constraint_solution[ci] = (value(ci), dual(ci))
#             end
#         catch
#             @info("Something went wrong getting $F-in-$S. Skipping")
#         end
#     end
#     # Now we can loop through our cached solutions and set the starting values.
#     for (x, primal_start) in variable_primal
#         set_start_value(x, primal_start)
#     end
#     for (ci, (primal_start, dual_start)) in constraint_solution
#         set_start_value(ci, primal_start)
#         set_dual_start_value(ci, dual_start)
#     end
#     return
# end

# set_optimal_start_values(sos_problem)
# JuMP.optimize!(sos_problem)