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

    relations_mono = vcat(
        [z(i)*zt(i)^(-1)*z(i)*z(i)*zt(i)^(-1)*z(i)*z(i)*zt(i)^(-1)*z(i)*z(i)*zt(i)^(-1)*z(i) for i in 1:N]
    )

    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples]                
    )

    # relations_list = [
    #     [com(x(i,j),y(i,j))*z(i)^(-2) for (i,j) in pairs],
    #     [com(x(i,j),yt(i,j))*zt(j)^2 for (i,j) in pairs],
    #     reduce(vcat, [[com(x(i,j),z(j))*(y(i,j)*z(i))^(-1), com(x(i,j),z(j))*(z(i)*y(i,j))^(-1), com(z(i),y(i,j))] for (i,j) in pairs]),
    #     reduce(vcat, [[com(x(i,j),zt(i))*yt(i,j)*zt(j)^(-1), com(x(i,j),zt(i))*zt(j)^(-1)*yt(i,j), com(zt(j),yt(i,j)^(-1))] for (i,j) in pairs]),
    #     reduce(vcat, [[com(y(i,j),zt(i))*z(j)*x(j,i)^(-1), com(y(i,j),zt(i))*x(j,i)^(-1)*z(j), com(x(j,i),z(j)^(-1))] for (i,j) in pairs]),
    #     reduce(vcat, [[com(yt(i,j), z(i))*zt(j)*x(i,j), com(yt(i,j), z(i))*x(i,j)*zt(j), com(x(i,j),zt(j))] for (i,j) in pairs])
    # ]
    
    # relations_sq = vcat([relations_list[k] for k in 1:length(relations_list) if !(k in W)]...)

        relations_sq = [
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
        ]
    relations_sq = vcat([relations_sq[k] for k in 1:length(relations_sq) if !(k in W)]...)

    sq_commutators = []
    for i in 1:N
        for j in 1:N
            if i != j
                push!(sq_commutators, com(x(i,j), z(i)))
                push!(sq_commutators, com(z(i), x(i,j)))

                push!(sq_commutators, com(x(i,j), zt(j)))
                push!(sq_commutators, com(zt(j), x(i,j)))


                push!(sq_commutators, com(y(i,j), z(j)))
                push!(sq_commutators, com(z(j), y(i,j)))

                push!(sq_commutators, com(yt(i,j), zt(j)))
                push!(sq_commutators, com(zt(j), yt(i,j)))


                push!(sq_commutators, com(y(i,j), z(i)))
                push!(sq_commutators, com(z(i), y(i,j)))

                push!(sq_commutators, com(yt(i,j), zt(i)))
                push!(sq_commutators, com(zt(i), yt(i,j)))
            end
        end
    end

    adj_commutators = []
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if i != j && j != k && i != k
                    push!(adj_commutators, com(x(i,j), x(i,k)))
                    push!(adj_commutators, com(x(i,j), x(k,j)))

                    push!(adj_commutators, com(x(i,j), y(i,k)))
                    push!(adj_commutators, com(y(i,k), x(i,j)))

                    push!(adj_commutators, com(x(i,j), yt(k,j)))
                    push!(adj_commutators, com(yt(k,j), x(i,j)))

                    push!(adj_commutators, com(y(i,j), y(j,k)))
                    push!(adj_commutators, com(y(i,j), y(i,k)))
                    push!(adj_commutators, com(y(i,j), y(k,j)))

                    push!(adj_commutators, com(yt(i,j), yt(j,k)))
                    push!(adj_commutators, com(yt(i,j), yt(j,k)))
                    push!(adj_commutators, com(yt(i,j), yt(j,k)))

                    push!(adj_commutators, com(z(i), y(j,k)))
                    push!(adj_commutators, com(y(j,k), z(i)))
                    push!(adj_commutators, com(z(i), yt(j,k)))
                    push!(adj_commutators, com(yt(j,k), z(i)))
                    push!(adj_commutators, com(z(i), x(j,k)))
                    push!(adj_commutators, com(x(j,k), z(i)))

                    push!(adj_commutators, com(zt(i), y(j,k)))
                    push!(adj_commutators, com(y(j,k), zt(i)))
                    push!(adj_commutators, com(zt(i), yt(j,k)))
                    push!(adj_commutators, com(yt(j,k), zt(i)))
                    push!(adj_commutators, com(zt(i), x(j,k)))
                    push!(adj_commutators, com(x(j,k), zt(i)))

                end
            end
        end
    end

    op_commutators = []
    for i in 1:N
        for j in 1:N
            for k in 1:N
                for l in 1:N
                    if i != j && j != k && k != l && i != k && i != l && j != l
                        push!(op_commutators, com(x(i,j), x(k,l)))

                        push!(op_commutators, com(x(i,j), y(k,l)))
                        push!(op_commutators, com(y(k,l), x(i,j)))

                        push!(op_commutators, com(x(i,j), yt(k,l)))
                        push!(op_commutators, com(yt(k,l), x(i,j)))

                        push!(op_commutators, com(y(i,j), y(k,l)))

                        push!(op_commutators, com(y(i,j), yt(k,l)))
                        push!(op_commutators, com(yt(k,l), y(i,j)))

                        push!(op_commutators, com(yt(i,j), yt(k,l)))
                    end
                end
            end
        end
    end

    z_commutators = sq_commutators
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if i != j && j != k && i != k
                    push!(z_commutators, com(z(i), y(j,k)))
                    push!(z_commutators, com(y(j,k), z(i)))
                    push!(z_commutators, com(z(i), yt(j,k)))
                    push!(z_commutators, com(yt(j,k), z(i)))
                    push!(z_commutators, com(z(i), x(j,k)))
                    push!(z_commutators, com(x(j,k), z(i)))

                    push!(z_commutators, com(zt(i), y(j,k)))
                    push!(z_commutators, com(y(j,k), zt(i)))
                    push!(z_commutators, com(zt(i), yt(j,k)))
                    push!(z_commutators, com(yt(j,k), zt(i)))
                    push!(z_commutators, com(zt(i), x(j,k)))
                    push!(z_commutators, com(x(j,k), zt(i)))

                end
            end
        end
    end

    relations_sq = vcat(relations_sq, sq_commutators)
    relations_adj = vcat(relations_adj, adj_commutators)
    relations_op = vcat(op_commutators)

    if sq_adj_ == "mono" 
        return relations_mono
    elseif sq_adj_ == "sq" 
        return relations_sq
    elseif sq_adj_ == "adj"
        return relations_adj
    elseif sq_adj_ == "op"
        return relations_op
    elseif sq_adj_ == "op_comm"
        return vcat(relations_op, z_commutators, relations_mono)
    elseif sq_adj_ == "adj_comm"
        return vcat(relations_adj, relations_sq, relations_op)
    elseif sq_adj_ == "all"
        return vcat(relations_mono, relations_sq, relations_adj, relations_op)
    end
end