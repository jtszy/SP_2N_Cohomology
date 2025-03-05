function relations_St(
    F_G::Groups.FreeGroup,
    S,
    N::Integer
)
    gen_dict = Dict(LowCohomologySOS.determine_letter(S[i]) => gens(F_G, i) for i in eachindex(S))
    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,max(i,j),min(i,j) + N)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,min(i,j) + N,max(i,j))]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i,i + N)]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i + N,i)]

    range_as_list = [i for i in 1:N]
    pairs = [
        (i,j) for i ∈ range_as_list 
              for j ∈ deleteat!(copy(range_as_list),findall(j->j==i,copy(range_as_list)))
    ]
    triples = [
        (i,j,k) for i ∈ range_as_list
                for j ∈ deleteat!(copy(range_as_list),findall(j->j==i,copy(range_as_list))) 
                for k ∈ deleteat!(copy(range_as_list),findall(k->k∈[i,j],copy(range_as_list)))
    ]
    
    relations_mono = vcat(
        [(z(i)*zt(i)^(-1)*z(i))^4 for i in 1:N] # Matsumoto's relation
    )
    relations_sq = vcat(
        [com(x(i,j),y(i,j))*z(i)^(-2) for (i,j) in pairs],
        [com(x(i,j),yt(i,j))*zt(j)^2 for (i,j) in pairs],
        [com(x(i,j),z(j))*(z(i)*y(i,j))^(-1) for (i,j) in pairs],
        [com(x(i,j),z(j))*(y(i,j)*z(i))^(-1) for (i,j) in pairs],
        [com(z(i),y(i,j)) for (i,j) in pairs],
        [com(x(i,j),zt(i))*yt(i,j)*zt(j)^(-1) for (i,j) in pairs],
        [com(x(i,j),zt(i))*zt(j)^(-1)*yt(i,j) for (i,j) in pairs],
        [com(zt(j),yt(i,j)^(-1)) for (i,j) in pairs],
        [com(y(i,j),zt(i))*z(j)*x(j,i)^(-1) for (i,j) in pairs],
        [com(y(i,j),zt(i))*x(j,i)^(-1)*z(j) for (i,j) in pairs],
        [com(x(j,i),z(j)^(-1)) for (i,j) in pairs],
        [com(yt(i,j),z(i))*zt(j)*x(i,j) for (i,j) in pairs],
        [com(yt(i,j),z(i))*x(i,j)*zt(j) for (i,j) in pairs],
        [com(x(i,j),zt(j)) for (i,j) in pairs]
    )
    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples],
    )
    
    lift = Dict(S[i] => gens(F_G, i) for i in eachindex(S))
    commutator_relations = []
    for i in eachindex(S)
        for j in eachindex(S)
            if i != j
                s, t = S[i], S[j]
                if com(s,t) == one(parent(first(S)))
                    push!(commutator_relations, com(lift[s], lift[t]))
                end
            end
        end
    end

    if N == 2 
        return vcat(relations_mono, relations_sq, commutator_relations)
    else
        # Matsumoto's relation yields too big support
        return vcat(relations_sq, relations_adj, commutator_relations)
    end
end
