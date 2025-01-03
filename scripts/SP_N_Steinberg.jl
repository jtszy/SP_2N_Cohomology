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

N = 3

Sp_N = MatrixGroups.SymplecticGroup{2*N}(Int8)

S = gens(Sp_N)

F_Sp_N_Steinberg = FreeGroup(alphabet(Sp_N))

F_z = FreeGroup(2 * N*(N-1))

G = Groups.Constructions.DirectProduct(F_Sp_N_Steinberg,F_z)

g = vcat(gens(G)[1 : N^2 - N], gens(G)[N^2+1 : Int((3 *N^2 - N) / 2)], gens(G)[Int((3 *N^2 + N) / 2) + 1: 4*N^2 - 2*N])
#removing z_i from original generators of Sp_N

g = vcat(g, inv.(g))

K = 4 * N^2 - 4 * N # number of new generators

F = FreeGroup(Alphabet(g, [(K+1):(2K); 1:K]))

S_without_z = vcat(S[1 : N^2 - N], S[N^2+1 : Int((3 *N^2 - N) / 2)], S[Int((3 *N^2 + N) / 2) + 1: 2*N^2])

function index_projection(i)
    if i <= N^2 - N # x variables
        return i
    elseif i <= Int(3 * (N^2 - N) / 2) # yt variables
        return i + N
    elseif i <= 2*N^2 - 2*N # y variables
        return i + 2*N
    elseif i <= 3*N^2 - 3*N #zt variables
        ind =  Int(floor((i - 2*N^2 + 2*N - 1) / (N-1))) + 1
        return N^2 - N + ind
    elseif i <= 4*N^2 - 4*N #z variables
        ind =  Int(floor((i - 2*N^2 + 2*N - 1) / (N-1))) + 1 - N
        return Int((3*N^2 - N)/ 2) + ind
    else #inverses of generators
        return index_projection(i - 4*N^2 + 4*N) + 2*N^2 
    end
end

quotient_hom_Steinberg = let source = F, target = Sp_N
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([index_projection(i)]), source, target)
end

# for i in eachindex(S)
#     @assert quotient_hom_Steinberg(gens(F_Sp_N_Steinberg,i)) == S[i]
#     @assert quotient_hom_Steinberg(gens(F_Sp_N_Steinberg,i)^(-1)) == S[i]^(-1)
# end


support_jacobian, min_support = SP_4_Cohomology.symplectic_min_supports(quotient_hom_Steinberg, S_without_z; rels = "adj")

Steinberg_relations = SP_4_Cohomology.relations_St(F, S_without_z, N; sq_adj_ = "adj")

for r in Steinberg_relations
    @assert quotient_hom_Steinberg(r) == one(Sp_N)
end

Δ₁, I_N, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom_Steinberg, Steinberg_relations, support_jacobian);

Δm_sq, Δm_adj, Δm_op = SP_4_Cohomology.sq_adj_op(Δ₁⁻, S_without_z)

# I_N_mono, I_N_sq = SP_4_Cohomology.mono_sq_adj_op(I_N, S)

# Δ₁ = Δm_ + Δ₁⁺
# I_N = I_N_

RG = LowCohomologySOS.group_ring(Sp_N, min_support, star_multiplication = true)

Δ₁ = LowCohomologySOS.embed.(identity, Δ₁, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N, Ref(RG))

constraints_basis, psd_basis, Σ, action = SP_4_Cohomology.wedderburn_data(RG.basis, min_support, S); 
# probably here should be sth different than S

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Δ₁, σ, action.alphabet_perm, S) == Δ₁
    @assert LowCohomologySOS.act_on_matrix(I_N, σ, action.alphabet_perm, S) == I_N
end

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    sos_problem, P = LowCohomologySOS.sos_problem(
        Δ₁,
        I_N,
        w_dec_matrix
        # 0.7 / 0.05
    )
end

# Find a numerical spectral gap
JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-6, max_iters = 1000))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)

# Q, λ = deserialize("./Steinberg_Solution_Sp_6.sjl")
# Q = Q[2]
# λ = λ[2]

result_bool, _ = LowCohomologySOS.certify_sos_decomposition(Δ, I_N, λ, Q, min_support)


# serialize("./Steinberg_Solution_Sp_6_comm_relations.sjl", Solution)