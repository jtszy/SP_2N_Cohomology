function relations_St(
    F_G::Groups.FreeGroup,
    S, # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    N::Integer,
    W = [];
    sq_adj_ = "all"
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

    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples]                
    )

    relations_list = [
        [com(x(i,j),y(i,j))*z(i)^(-2) for (i,j) in pairs],
        [com(x(i,j),yt(i,j))*zt(j)^2 for (i,j) in pairs],
        reduce(vcat, [[com(x(i,j),z(j))*(y(i,j)*z(i))^(-1), com(x(i,j),z(j))*(z(i)*y(i,j))^(-1), com(z(i),y(i,j))] for (i,j) in pairs]),
        reduce(vcat, [[com(x(i,j),zt(i))*yt(i,j)*zt(j)^(-1), com(x(i,j),zt(i))*zt(j)^(-1)*yt(i,j), com(zt(j),yt(i,j)^(-1))] for (i,j) in pairs]),
        reduce(vcat, [[com(y(i,j),zt(i))*z(j)*x(j,i)^(-1), com(y(i,j),zt(i))*x(j,i)^(-1)*z(j), com(x(j,i),z(j)^(-1))] for (i,j) in pairs]),
        reduce(vcat, [[com(yt(i,j), z(i))*zt(j)*x(i,j), com(yt(i,j), z(i))*x(i,j)*zt(j), com(x(i,j),zt(j))] for (i,j) in pairs])
    ]
    
    relations_sq = vcat([relations_list[k] for k in 1:length(relations_list) if !(k in W)]...)

    if sq_adj_ == "sq" 
        return relations_sq
    elseif sq_adj_ == "adj"
        return relations_adj
    elseif sq_adj_ == "all"
        return vcat(relations_sq, relations_adj)
    end
end
