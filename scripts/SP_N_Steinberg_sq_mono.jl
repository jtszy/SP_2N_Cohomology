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

n = 2

Sp_N = MatrixGroups.SymplecticGroup{2*n}(Int8)

function extended_f_sp_2n(n::Integer)

    x_z_gens = vcat(
        [(:x,i) for i in 1:(n-1)],
        [(:xt,i) for i in 1:(n-1)],
        [(:z,i) for i in 1:n],
        [(:zt,i) for i in 1:n],
        [(:X,i) for i in 1:(n-1)],
        [(:Xt,i) for i in 1:(n-1)],
        [(:Z,i) for i in 1:n],
        [(:Zt,i) for i in 1:n]
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
        (type,i) = LowCohomologySOS.determine_letter(gens(F,i))

        if type == :x
            new_type, new_i, new_j = :A, 1, i+1
        elseif type == :xt
            new_type, new_i, new_j = :A, i+1, 1
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
    N::Integer
)
    gen_dict = Dict(LowCohomologySOS.determine_letter(gens(F_G, i)) => gens(F_G, i) for i in eachindex(gens(F_G)))

    x(i,j) = let
        if i == 1
            result = gen_dict[(:x,j-1)]
        elseif j == 1
            result = gen_dict[(:xt,i-1)]
        else
            result = com(x(i,1),x(1,j))
        end
        
        result
    end
    z(i) = gen_dict[(:z,i)]
    zt(i) = gen_dict[(:zt,i)]
    y(i,j) = com(x(i,j),z(j))*z(i)^(-1)
    yt(i,j) = (com(x(i,j),zt(i))*zt(j)^(-1))^(-1)
    
    range_as_list = [i for i in 1:N]
    pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ range_as_list
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]

    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples]                
    )

    relations_sq = vcat(
        [(z(i)*zt(i)^(-1)*z(i))^4 for i in 1:N],

        [com(x(i,j),y(i,j))*z(i)^(-2) for (i,j) in pairs],
        [com(x(i,j),yt(i,j))*zt(j)^2 for (i,j) in pairs],
        [com(x(i,j),z(j))*(y(i,j)*z(i))^(-1) for (i,j) in pairs],
        [com(x(i,j),z(j))*(z(i)*y(i,j))^(-1) for (i,j) in pairs],
        [com(z(i),y(i,j)) for (i,j) in pairs],
        [com(x(i,j),zt(i))*yt(i,j)*zt(j)^(-1) for (i,j) in pairs],
        [com(x(i,j),zt(i))*zt(j)^(-1)*yt(i,j) for (i,j) in pairs],
        [com(zt(j),yt(i,j)^(-1)) for (i,j) in pairs],
        [com(y(i,j),zt(i))*z(j)*x(j,i)^(-1) for (i,j) in pairs],
        [com(y(i,j),zt(i))*x(j,i)^(-1)*z(j) for (i,j) in pairs],
        [com(x(j,i),z(j)^(-1)) for (i,j) in pairs],
        [com(yt(i,j), z(i))*zt(j)*x(i,j) for (i,j) in pairs],
        [com(yt(i,j), z(i))*x(i,j)*zt(j) for (i,j) in pairs],
        [com(x(i,j),zt(j)) for (i,j) in pairs]
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

    commutators = vcat(pair_commutators, triplet_commutators, quad_commutators)

    return vcat(relations_sq, relations_adj, commutators)
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

Steinberg_relations = relations_St(Extended_f_sp_2n, n)
support_jacobian, min_support = symplectic_min_supports(Steinberg_relations, quotient_hom, S)

Δ₁, I_N, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom, Steinberg_relations, support_jacobian);

function mono_sq(
    Δ₁⁻
)
    mono_pairs = []
    sq_pairs = []

    A = alphabet(Extended_f_sp_2n)
    S = gens(Extended_f_sp_2n)
    RG = parent(first(Δ₁⁻))
    for s in eachindex(S)
        for t in eachindex(S)
            s_i = A[word(S[s])[1]][2]
            t_i = A[word(S[t])[1]][2]
            if s_i == t_i
                push!(mono_pairs,(s,t))
            else
                push!(sq_pairs,(s,t))
            end
        end
    end
    mono = [(i,j) in mono_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    sq = [(i,j) in sq_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]

    @assert mono+sq == Δ₁⁻

    return mono, sq
end

Δm_mono, Δm_sq = mono_sq(Δ₁⁻)
Δ = Δm_sq + Δ₁⁺

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
