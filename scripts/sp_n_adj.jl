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
using SymbolicWedderburn

N = 3

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

com(a,b) = SP_4_Cohomology.com(a,b)

function relations_St(
    F_G::Groups.FreeGroup,
    S, # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    N::Integer;
    sq_adj_op = "all"
)
    gen_dict = Dict(LowCohomologySOS.determine_letter(S[i]) => gens(F_G, i) for i in eachindex(S))

    range_as_list = [i for i in 1:N]
    pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ range_as_list
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    
    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,max(i,j),min(i,j) + N)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,min(i,j) + N,max(i,j))]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i,i + N)]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i + N,i)]
    
    relations_mono = vcat(
        # [(z(i)*zt(i)^(-1)*z(i))^4 for i in 1:N],
        [com(z(i),zt(i)) for i in 1:N]
    )

    relations_sq = vcat(
        [com(x(i,j),y(i,j))*z(i)^(-2) for (i,j) in pairs],
        [com(x(i,j),yt(i,j))*zt(j)^2 for (i,j) in pairs],
        [com(x(i,j),z(j))*(y(i,j)*z(i))^(-1) for (i,j) in pairs],
        [com(z(i),y(i,j)) for (i,j) in pairs],
        [com(x(i,j),zt(i))*yt(i,j)*zt(j)^(-1) for (i,j) in pairs],
        [com(zt(j),yt(i,j)^(-1)) for (i,j) in pairs],
        [com(y(i,j),zt(i))*z(j)*x(j,i)^(-1) for (i,j) in pairs],
        [com(x(j,i),z(j)^(-1)) for (i,j) in pairs],
        [com(yt(i,j), z(i))*zt(j)*x(i,j) for (i,j) in pairs],
        [com(x(i,j),zt(j)) for (i,j) in pairs],

        # [com(x(i,j),x(j,i)) for (i,j) in pairs],
        # [com(y(i,j),yt(i,j)) for (i,j) in pairs],
    )

    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples]                
    )
    
    pair_commutators = []
    for i in 1:N
        for j in 1:N
            if i != j
                push!(pair_commutators, com(x(i,j), z(i)))
                push!(pair_commutators, com(x(i,j), zt(j)))

                push!(pair_commutators, com(y(i,j), z(j)))
                push!(pair_commutators, com(yt(i,j), zt(j)))

                push!(pair_commutators, com(y(i,j), z(i)))
                push!(pair_commutators, com(yt(i,j), zt(i)))
            end
        end
    end

    triplet_commutators = []
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if i != j && j != k && i != k
                    push!(triplet_commutators, com(x(i,j), x(i,k)))
                    push!(triplet_commutators, com(x(i,j), x(k,j)))

                    push!(triplet_commutators, com(x(i,j), y(i,k)))

                    push!(triplet_commutators, com(x(i,j), yt(k,j)))

                    push!(triplet_commutators, com(y(i,j), y(j,k)))
                    push!(triplet_commutators, com(y(i,j), y(i,k)))
                    push!(triplet_commutators, com(y(i,j), y(k,j)))

                    push!(triplet_commutators, com(yt(i,j), yt(j,k)))
                    push!(triplet_commutators, com(yt(i,j), yt(j,k)))
                    push!(triplet_commutators, com(yt(i,j), yt(j,k)))

                    push!(triplet_commutators, com(z(i), y(j,k)))
                    push!(triplet_commutators, com(z(i), yt(j,k)))
                    push!(triplet_commutators, com(z(i), x(j,k)))

                    push!(triplet_commutators, com(zt(i), y(j,k)))
                    push!(triplet_commutators, com(zt(i), yt(j,k)))
                    push!(triplet_commutators, com(zt(i), x(j,k)))
                end
            end
        end
    end

    quad_commutators = []
    for i in 1:N
        for j in 1:N
            for k in 1:N
                for l in 1:N
                    if i != j && j != k && k != l && i != k && i != l && j != l
                        push!(quad_commutators, com(x(i,j), x(k,l)))

                        push!(quad_commutators, com(x(i,j), y(k,l)))
                        push!(quad_commutators, com(x(i,j), yt(k,l)))

                        push!(quad_commutators, com(y(i,j), y(k,l)))

                        push!(quad_commutators, com(y(i,j), yt(k,l)))
                        
                        push!(quad_commutators, com(yt(i,j), yt(k,l)))
                    end
                end
            end
        end
    end

    if sq_adj_op == "sq" 
        return vcat(relations_sq, pair_commutators)
    elseif sq_adj_op == "adj"
        return vcat(relations_adj, triplet_commutators)
    elseif sq_adj_op == "op"
        return vcat(quad_commutators)
        # return vcat(relations_sq, relations_adj, pair_commutators, triplet_commutators, quad_commutators)
        # return vcat(relations_adj, pair_commutators, triplet_commutators, quad_commutators)
    else
        return vcat(relations_mono, relations_sq, relations_adj, pair_commutators, triplet_commutators, quad_commutators)
        # return vcat(relations_mono, relations_adj, pair_commutators, triplet_commutators, quad_commutators)
    end
end

function symplectic_min_supports(
    Steinberg_relations,
    quotient_hom,
    S
)   
    Sp_N = quotient_hom.target

    sup_jacobian = SP_4_Cohomology.support_jacobian(vcat(Steinberg_relations, S), quotient_hom)
    min_support = SP_4_Cohomology.minimalistic_support(Steinberg_relations, quotient_hom)

    return sup_jacobian, min_support
end

Steinberg_relations = relations_St(F_Sp_2N_Steinberg, S, N, sq_adj_op = "all")

# for r in Steinberg_relations
#     @assert quotient_hom_Steinberg(r) == one(Sp_2N)
# end

support_jacobian, min_support = symplectic_min_supports(Steinberg_relations, quotient_hom_Steinberg, S)

min_support

# gens_sp4_inv = [S; inv.(S)]
# Ball_4, sizes = Groups.wlmetric_ball(gens_sp4_inv, radius=2)

# min_support = collect(union(Set(min_support),Set(Ball_4)))

Δ₁, I_N_whole, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom_Steinberg, Steinberg_relations, support_jacobian);

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
            if sort([s_i,s_j]) == sort([t_i, t_j])
                if s_i == s_j
                    push!(mono_pairs,(s,t))
                else
                    push!(sq_pairs,(s,t))
                end
            elseif length(intersect!([s_i,s_j],[t_i,t_j])) == 1
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

Δm_mono, Δm_sq, Δm_adj, Δm_op  = mono_sq_adj_op(Δ₁⁻,S)
I_mono, I_sq, I_adj, I_op  = mono_sq_adj_op(I_N_whole,S)

# Δ = Δm_adj+Δm_op + 10*Δ₁⁺
# Δ = Δm_sq + Δm_adj + Δm_op + Δ₁⁺
Δ = Δm_adj + Δm_op + Δ₁⁺

RG = LowCohomologySOS.group_ring(Sp_2N, min_support, star_multiplication = true)

Δ = LowCohomologySOS.embed.(identity, Δ, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N_whole, Ref(RG))

constraints_basis, psd_basis, Σ, action = SP_4_Cohomology.wedderburn_data(RG.basis, min_support, S);

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
    @info "SDP problem definition"
    sos_problem, P = LowCohomologySOS.sos_problem(
        Δ,
        I_N,
        w_dec_matrix,
        0.003
    )
end

# Find a numerical spectral gap
JuMP.set_optimizer(sos_problem, SP_4_Cohomology.scs_opt(eps = 1e-6, max_iters = 13000))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
LowCohomologySOS.certify_sos_decomposition(Δ, I_N, λ, Q, min_support)

solution = Dict("lambda" => λ, "Q" => Q)
serialize("./scripts/Steinberg_Solution_Sp_6_adj.sjl", solution)