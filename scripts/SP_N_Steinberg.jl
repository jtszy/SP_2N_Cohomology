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

n = 3

Sp_N = MatrixGroups.SymplecticGroup{2*n}(Int8)

function extended_f_sp_2n(n::Integer)
    range_as_list = [i for i in 1:n]
    ordered_pairs = [(i,j) for i ∈ 1:n for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    unordered_pairs = [(i,j) for i ∈ 1:n for j ∈ deleteat!(copy(range_as_list), findall(j->j<i,copy(range_as_list)))]
    
    x_y_z_gens = vcat(
        [(:x,i,j) for (i,j) in ordered_pairs],
        [(:y,i,j) for (i,j) in unordered_pairs],
        [(:yt,i,j) for (i,j) in unordered_pairs],
        [(:z,i,j) for (i,j) in ordered_pairs],
        [(:zt,i,j) for (i,j) in ordered_pairs],
        [(:X,i,j) for (i,j) in ordered_pairs],
        [(:Y,i,j) for (i,j) in unordered_pairs],
        [(:Yt,i,j) for (i,j) in unordered_pairs],
        [(:Z,i,j) for (i,j) in ordered_pairs],
        [(:Zt,i,j) for (i,j) in ordered_pairs],
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

com(a,b) = SP_4_Cohomology.com(a,b)

function relations_St(
    F_G::Groups.FreeGroup,
    N::Integer;
    sq_adj_op = "all"
)
    gen_dict = Dict(LowCohomologySOS.determine_letter(gens(F_G, i)) => gens(F_G, i) for i in eachindex(gens(F_G)))

    range_as_list = [i for i in 1:N]
    pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ range_as_list
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    quadruples = [(i,j,k,l) for i ∈ range_as_list
                        for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                        for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))
                        for l ∈ deleteat!(copy(range_as_list), findall(l->l∈[i,j,k],copy(range_as_list)))]
    
    x(i,j) = gen_dict[(:x,i,j)]
    y(i,j) = gen_dict[(:y,min(i,j),max(i,j))]
    yt(i,j) = gen_dict[(:yt,min(i,j),max(i,j))]
    z(i,j) = gen_dict[(:z,i,j)]
    zt(i,j) = gen_dict[(:zt,i,j)]

    relations_sq = vcat(
        [z(i,j)*zt(i,j)^(-1)*z(i,j)*z(i,j)*zt(i,j)^(-1)*z(i,j)*z(i,j)*zt(i,j)^(-1)*z(i,j)*z(i,j)*zt(i,j)^(-1)*z(i,j) for (i,j) in pairs],

        # [com(x(i,j),y(i,j))*z(i,j)^(-2) for (i,j) in pairs],
        # [com(x(i,j),yt(i,j))*zt(j,i)^2 for (i,j) in pairs],
        # [com(x(i,j),z(j,i))*(y(i,j)*z(i,j))^(-1) for (i,j) in pairs],
        # [com(x(i,j),z(j,i))*(z(i,j)*y(i,j))^(-1) for (i,j) in pairs],
        # [com(z(i,j),y(i,j)) for (i,j) in pairs],
        # [com(x(i,j),zt(i,j))*yt(i,j)*zt(j,i)^(-1) for (i,j) in pairs],
        # [com(x(i,j),zt(i,j))*zt(j,i)^(-1)*yt(i,j) for (i,j) in pairs],
        # [com(zt(j,i),yt(i,j)^(-1)) for (i,j) in pairs],
        # [com(y(i,j),zt(i,j))*z(j,i)*x(j,i)^(-1) for (i,j) in pairs],
        # [com(y(i,j),zt(i,j))*x(j,i)^(-1)*z(j,i) for (i,j) in pairs],
        # [com(x(j,i),z(j,i)^(-1)) for (i,j) in pairs],
        # [com(yt(i,j), z(i,j))*zt(j,i)*x(i,j) for (i,j) in pairs],
        # [com(yt(i,j), z(i,j))*x(i,j)*zt(j,i) for (i,j) in pairs],
        # [com(x(i,j),zt(j,i)) for (i,j) in pairs],

        # [com(x(i,j), z(i,j)) for (i,j) in pairs],
        # [com(z(i,j), x(i,j)) for (i,j) in pairs],
        # [com(x(i,j), zt(j,i)) for (i,j) in pairs],
        # [com(zt(j,i), x(i,j)) for (i,j) in pairs],
        # [com(y(i,j), z(j,i)) for (i,j) in pairs],
        # [com(z(j,i), y(i,j)) for (i,j) in pairs],
        # [com(yt(i,j), zt(j,i)) for (i,j) in pairs],
        # [com(zt(j,i), yt(i,j)) for (i,j) in pairs],
        # [com(y(i,j), z(i,j)) for (i,j) in pairs],
        # [com(z(i,j), y(i,j)) for (i,j) in pairs],
        # [com(yt(i,j), zt(i,j)) for (i,j) in pairs],
        # [com(zt(i,j), yt(i,j)) for (i,j) in pairs],
    )
    
    relations_adj = vcat(
        # [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        # [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        # [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples],
        
        # [com(x(i,j),y(i,j))*z(i,k)^(-2) for (i,j,k) in triples],
        # [com(x(i,j),yt(i,j))*zt(j,k)^2 for (i,j,k) in triples],
        # [com(x(i,j),z(j,k))*(y(i,j)*z(i,k))^(-1) for (i,j,k) in triples],
        # [com(x(i,j),z(j,k))*(z(i,k)*y(i,j))^(-1) for (i,j,k) in triples],
        [com(z(i,k),y(i,j)) for (i,j,k) in triples],
        # [com(x(i,j),zt(i,k))*yt(i,j)*zt(j,k)^(-1) for (i,j,k) in triples],
        # [com(x(i,j),zt(i,k))*zt(j,k)^(-1)*yt(i,j) for (i,j,k) in triples],
        [com(zt(j,k),yt(i,j)^(-1)) for (i,j,k) in triples],
        # [com(y(i,j),zt(i,k))*z(j,k)*x(j,i)^(-1) for (i,j,k) in triples],
        # [com(y(i,j),zt(i,k))*x(j,i)^(-1)*z(j,k) for (i,j,k) in triples],
        [com(x(j,i),z(j,k)^(-1)) for (i,j,k) in triples],
        # [com(yt(i,j), z(i,k))*zt(j,k)*x(i,j) for (i,j,k) in triples],
        # [com(yt(i,j), z(i,k))*x(i,j)*zt(j,k) for (i,j,k) in triples],

        [com(x(i,j), z(i,k)) for (i,j,k) in triples],
        [com(z(i,k), x(i,j)) for (i,j,k) in triples],
        [com(x(i,j), zt(j,k)) for (i,j,k) in triples],
        [com(zt(j,k), x(i,j)) for (i,j,k) in triples],
        [com(y(i,j), z(j,k)) for (i,j,k) in triples],
        [com(z(j,k), y(i,j)) for (i,j,k) in triples],
        [com(yt(i,j), zt(j,k)) for (i,j,k) in triples],
        [com(zt(j,k), yt(i,j)) for (i,j,k) in triples],
        [com(y(i,j), z(i,k)) for (i,j,k) in triples],
        [com(z(i,k), y(i,j)) for (i,j,k) in triples],
        [com(yt(i,j), zt(i,k)) for (i,j,k) in triples],
        [com(zt(i,k), yt(i,j)) for (i,j,k) in triples],

        [com(x(i,j), x(i,k)) for (i,j,k) in triples],
        [com(x(i,j), x(k,j)) for (i,j,k) in triples],
        [com(x(i,j), y(i,k)) for (i,j,k) in triples],
        [com(y(i,k), x(i,j)) for (i,j,k) in triples],
        [com(x(i,j), yt(k,j)) for (i,j,k) in triples],
        [com(yt(k,j), x(i,j)) for (i,j,k) in triples],
        [com(y(i,j), y(j,k)) for (i,j,k) in triples],
        [com(y(i,j), y(i,k)) for (i,j,k) in triples],
        [com(y(i,j), y(k,j)) for (i,j,k) in triples],
        [com(yt(i,j), yt(j,k)) for (i,j,k) in triples],
        [com(yt(i,j), yt(i,k)) for (i,j,k) in triples],
        [com(yt(i,j), yt(k,j)) for (i,j,k) in triples],

        [com(z(i,j), y(j,k)) for (i,j,k) in triples],
        [com(z(i,k), y(j,k)) for (i,j,k) in triples],
        [com(y(j,k), z(i,j)) for (i,j,k) in triples],
        [com(y(j,k), z(i,k)) for (i,j,k) in triples],
        [com(z(i,j), yt(j,k)) for (i,j,k) in triples],
        [com(z(i,k), yt(j,k)) for (i,j,k) in triples],
        [com(yt(j,k), z(i,j)) for (i,j,k) in triples],
        [com(yt(j,k), z(i,k)) for (i,j,k) in triples],
        [com(z(i,j), x(j,k)) for (i,j,k) in triples],
        [com(z(i,k), x(j,k)) for (i,j,k) in triples],
        [com(x(j,k), z(i,j)) for (i,j,k) in triples],
        [com(x(j,k), z(i,k)) for (i,j,k) in triples],
        [com(zt(i,j), y(j,k)) for (i,j,k) in triples],
        [com(zt(i,k), y(j,k)) for (i,j,k) in triples],
        [com(y(j,k), zt(i,j)) for (i,j,k) in triples],
        [com(y(j,k), zt(i,k)) for (i,j,k) in triples],
        [com(zt(i,j), yt(j,k)) for (i,j,k) in triples],
        [com(zt(i,k), yt(j,k)) for (i,j,k) in triples],
        [com(yt(j,k), zt(i,j)) for (i,j,k) in triples],
        [com(yt(j,k), zt(i,k)) for (i,j,k) in triples],
        [com(zt(i,j), x(j,k)) for (i,j,k) in triples],
        [com(zt(i,k), x(j,k)) for (i,j,k) in triples],
        [com(x(j,k), zt(i,j)) for (i,j,k) in triples],
        [com(x(j,k), zt(i,k)) for (i,j,k) in triples],

        [com(z(i,j), y(j,k)) for (i,j,k) in triples],
        [com(z(i,k), y(j,k)) for (i,j,k) in triples],
        [com(y(j,k), z(i,j)) for (i,j,k) in triples],
        [com(y(j,k), z(i,k)) for (i,j,k) in triples],
        [com(z(i,j), yt(j,k)) for (i,j,k) in triples],
        [com(z(i,k), yt(j,k)) for (i,j,k) in triples],
        [com(yt(j,k), z(i,j)) for (i,j,k) in triples],
        [com(yt(j,k), z(i,k)) for (i,j,k) in triples],
        [com(z(i,j), x(j,k)) for (i,j,k) in triples],
        [com(z(i,k), x(j,k)) for (i,j,k) in triples],
        [com(x(j,k), z(i,j)) for (i,j,k) in triples],
        [com(x(j,k), z(i,k)) for (i,j,k) in triples],
        [com(zt(i,j), y(j,k)) for (i,j,k) in triples],
        [com(zt(i,k), y(j,k)) for (i,j,k) in triples],
        [com(y(j,k), zt(i,j)) for (i,j,k) in triples],
        [com(y(j,k), zt(i,k)) for (i,j,k) in triples],
        [com(zt(i,j), yt(j,k)) for (i,j,k) in triples],
        [com(zt(i,k), yt(j,k)) for (i,j,k) in triples],
        [com(yt(j,k), zt(i,j)) for (i,j,k) in triples],
        [com(yt(j,k), zt(i,k)) for (i,j,k) in triples],
        [com(zt(i,j), x(j,k)) for (i,j,k) in triples],
        [com(zt(i,k), x(j,k)) for (i,j,k) in triples],
        [com(x(j,k), zt(i,j)) for (i,j,k) in triples],
        [com(x(j,k), zt(i,k)) for (i,j,k) in triples]
    )

    relations_op = vcat(
        [com(z(i,l), y(j,k)) for (i,j,k,l) in quadruples],
        [com(y(j,k), z(i,l)) for (i,j,k,l) in quadruples],
        [com(z(i,l), yt(j,k)) for (i,j,k,l) in quadruples],
        [com(yt(j,k), z(i,l)) for (i,j,k,l) in quadruples],
        [com(z(i,l), x(j,k)) for (i,j,k,l) in quadruples],
        [com(x(j,k), z(i,l)) for (i,j,k,l) in quadruples],
        [com(zt(i,l), y(j,k)) for (i,j,k,l) in quadruples],
        [com(y(j,k), zt(i,l)) for (i,j,k,l) in quadruples],
        [com(zt(i,l), yt(j,k)) for (i,j,k,l) in quadruples],
        [com(yt(j,k), zt(i,l)) for (i,j,k,l) in quadruples],
        [com(zt(i,l), x(j,k)) for (i,j,k,l) in quadruples],
        [com(x(j,k), zt(i,l)) for (i,j,k,l) in quadruples],

        [com(x(i,j), x(k,l)) for (i,j,k,l) in quadruples],
        [com(x(i,j), y(k,l)) for (i,j,k,l) in quadruples],
        [com(y(k,l), x(i,j)) for (i,j,k,l) in quadruples],
        [com(x(i,j), yt(k,l)) for (i,j,k,l) in quadruples],
        [com(yt(k,l), x(i,j)) for (i,j,k,l) in quadruples],
        [com(y(i,j), y(k,l)) for (i,j,k,l) in quadruples],
        [com(y(i,j), yt(k,l)) for (i,j,k,l) in quadruples],
        [com(yt(k,l), y(i,j)) for (i,j,k,l) in quadruples],
        [com(yt(i,j), yt(k,l)) for (i,j,k,l) in quadruples],

        [com(z(i,l), y(j,k)) for (i,j,k,l) in quadruples],
        [com(y(j,k), z(i,l)) for (i,j,k,l) in quadruples],
        [com(z(i,l), yt(j,k)) for (i,j,k,l) in quadruples],
        [com(yt(j,k), z(i,l)) for (i,j,k,l) in quadruples],
        [com(z(i,l), x(j,k)) for (i,j,k,l) in quadruples],
        [com(x(j,k), z(i,l)) for (i,j,k,l) in quadruples],
        [com(zt(i,l), y(j,k)) for (i,j,k,l) in quadruples],
        [com(y(j,k), zt(i,l)) for (i,j,k,l) in quadruples],
        [com(zt(i,l), yt(j,k)) for (i,j,k,l) in quadruples],
        [com(yt(j,k), zt(i,l)) for (i,j,k,l) in quadruples],
        [com(zt(i,l), x(j,k)) for (i,j,k,l) in quadruples],
        [com(x(j,k), zt(i,l)) for (i,j,k,l) in quadruples]
    )

    if sq_adj_op == "sq" 
        return relations_sq
    elseif sq_adj_op == "adj"
        return relations_adj
    elseif sq_adj_op == "op"
        return relations_op
    elseif sq_adj_op == "all"
        return vcat(relations_sq, relations_adj, relations_op)
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

Steinberg_relations = relations_St(Extended_f_sp_2n, n, sq_adj_op = "adj")
# support_jacobian, min_support = symplectic_min_supports(Steinberg_relations, quotient_hom, S)

# Steinberg_relations = relations_St(Extended_f_sp_2n, n)
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

# Δ = Δm_adj + Δ₁⁺
Δ = Δm_sq + Δm_adj + Δm_op + Δ₁⁺
# Δ = Δ₁

RG = LowCohomologySOS.group_ring(Sp_N, min_support, star_multiplication = true)
RG.mstructure
Δ = LowCohomologySOS.embed.(identity, Δ, Ref(RG))
I_N = LowCohomologySOS.embed.(identity, I_N, Ref(RG))

# Wedderburn has to be customized for zij's ########################################################
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
###########################################################################################

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
# λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)

result_bool, _ = LowCohomologySOS.certify_sos_decomposition(Δ, I_N, λ, Q, min_support)
