using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using SP_4_Cohomology
# using StarAlgebras
using Groups


# TODO: adjust the whole code below to sp2ns

const N = 3
const M = 4

i = SP_4_Cohomology.sp2n_sp2m_embedding(2*N, 2*M)

Sp2N = i.source
Sp2M = i.target

const half_radius = 2

Sp2N_S = gens(Sp2N);
Sp2N_S_inv = let s = Sp2N_S
    [s; inv.(s)]
end;
Sp2M_S = gens(Sp2M);
Sp2M_S_inv = let s = Sp2M_S
    [s; inv.(s)]
end;

Sp2N_half_basis, Sp2N_sizes = Groups.wlmetric_ball(Sp2N_S_inv, radius = half_radius);
Sp2M_half_basis, Sp2M_sizes = Groups.wlmetric_ball(Sp2M_S_inv, radius = half_radius);

Sp2N_Δ₁, Sp2N_Iₙ, Sp2N_Δ₁⁺, Sp2N_Δ₁⁻ = LowCohomologySOS.laplacians(Sp2N, Sp2N_half_basis, Sp2N_S, sq_adj_ = "adj");
Sp2M_Δ₁, Sp2M_Iₙ, Sp2M_Δ₁⁺, Sp2M_Δ₁⁻ = LowCohomologySOS.laplacians(Sp2M, Sp2M_half_basis, Sp2M_S, sq_adj_ = "adj");
Sp2N_sq, Sp2N_adj, Sp2N_op = LowCohomologySOS.sq_adj_op(Sp2N_Δ₁⁻, Sp2N_S)
Sp2M_sq, Sp2M_adj, Sp2M_op = LowCohomologySOS.sq_adj_op(Sp2M_Δ₁⁻, Sp2M_S)

RG_prime = parent(first(Sp2M_Δ₁⁺))

Δ₁⁺_emb = LowCohomologySOS.embed_matrix(Sp2N_Δ₁⁺, i, RG_prime)
adj_emb = LowCohomologySOS.embed_matrix(Sp2N_adj, i, RG_prime)

@assert parent(first(Δ₁⁺_emb)) == parent(first(adj_emb)) == parent(first(Sp2M_Δ₁⁺)) == parent(first(Sp2M_adj))

using PermutationGroups

Δ₁⁺_emb_symmetrized = let
    Σ = PermutationGroups.SymmetricGroup(3)
    LowCohomologySOS.weyl_symmetrize_matrix(Δ₁⁺_emb, Σ, LowCohomologySOS._conj, Sp2M_S)
end

adj_emb_symmetrized = let
    Σ = PermutationGroups.SymmetricGroup(3)
    LowCohomologySOS.weyl_symmetrize_matrix(adj_emb, Σ, LowCohomologySOS._conj, Sp2M_S)
end

# TODO change scalars below:
24*Sp2M_Δ₁⁺-Δ₁⁺_emb_symmetrized # it looks like symmetrization works for upper Laplacians!
24*Sp2M_adj-adj_emb_symmetrized # Adj symmetrizes as well with the same pace!!