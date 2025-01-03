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
    zt(i,j) = gens(F_G)[2*N^2 - 2*N + (i-1)*(N-1) + j - Int(i<j)] # adjustment of indexing 
    z(i,j) = gens(F_G)[3*N^2 - 3*N + (i-1)*(N-1) + j - Int(i<j)]

    # relations_mono = vcat(
    #     [z(i)*zt(i)^(-1)*z(i)*z(i)*zt(i)^(-1)*z(i)*z(i)*zt(i)^(-1)*z(i)*z(i)*zt(i)^(-1)*z(i) for i in 1:N]
    # )
    # ignored for now

    relations_adj = vcat(
        [com(x(i,j),x(j,k))*x(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),y(j,k))*y(i,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),yt(i,k))*yt(j,k) for (i,j,k) in triples], #original adj relations
        
        [com(x(i,j),y(i,j))*z(i,k)^(-2) for (i,j,k) in triples],
        [com(x(i,j),yt(i,j))*zt(j,k)^2 for (i,j,k) in triples],
        [com(x(i,j),z(j,k))*(y(i,j)*z(i,j))^(-1) for (i,j,k) in triples],
        [com(x(i,j),z(j,k))*(z(i,j)*y(i,j))^(-1) for (i,j,k) in triples],
        [com(x(i,j),z(j,k))*(y(i,j)*z(i,k))^(-1) for (i,j,k) in triples],
        [com(x(i,j),z(j,k))*(z(i,k)*y(i,j))^(-1) for (i,j,k) in triples],
        [com(z(i,k),y(i,j)) for (i,j,k) in triples],
        [com(x(i,j),zt(i,k))*yt(i,j)*zt(j,i)^(-1) for (i,j,k) in triples],
        [com(x(i,j),zt(i,k))*zt(j,i)^(-1)*yt(i,j) for (i,j,k) in triples],
        [com(x(i,j),zt(i,k))*yt(i,j)*zt(j,k)^(-1) for (i,j,k) in triples],
        [com(x(i,j),zt(i,k))*zt(j,k)^(-1)*yt(i,j) for (i,j,k) in triples],
        [com(zt(j,k),yt(i,j)^(-1)) for (i,j,k) in triples],
        [com(y(i,j),zt(i,k))*z(j,i)*x(j,i)^(-1) for (i,j,k) in triples],
        [com(y(i,j),zt(i,k))*x(j,i)^(-1)*z(j,i) for (i,j,k) in triples],
        [com(y(i,j),zt(i,k))*z(j,k)*x(j,i)^(-1) for (i,j,k) in triples],
        [com(y(i,j),zt(i,k))*x(j,i)^(-1)*z(j,k) for (i,j,k) in triples],
        [com(x(j,i),z(j,k)^(-1)) for (i,j,k) in triples],
        [com(yt(i,j), z(i,k))*zt(j,i)*x(i,j) for (i,j,k) in triples],
        [com(yt(i,j), z(i,k))*x(i,j)*zt(j,i) for (i,j,k) in triples],
        [com(yt(i,j), z(i,k))*zt(j,k)*x(i,j) for (i,j,k) in triples],
        [com(yt(i,j), z(i,k))*x(i,j)*zt(j,k) for (i,j,k) in triples],
        [com(x(i,j),zt(j,k)) for (i,j,k) in triples] # new adj relations given by double indexed z's
    )
    relations_sq = vcat(
            [com(x(i,j),y(i,j))*z(i,j)^(-2) for (i,j) in pairs],
            [com(x(i,j),yt(i,j))*zt(j,i)^2 for (i,j) in pairs],
            [com(x(i,j),z(j,i))*(y(i,j)*z(i,j))^(-1) for (i,j) in pairs],
            [com(x(i,j),z(j,i))*(z(i,j)*y(i,j))^(-1) for (i,j) in pairs],
            [com(z(i,j),y(i,j)) for (i,j) in pairs],
            [com(x(i,j),zt(i,j))*yt(i,j)*zt(j,i)^(-1) for (i,j) in pairs],
            [com(x(i,j),zt(i,j))*zt(j,i)^(-1)*yt(i,j) for (i,j) in pairs],
            [com(zt(j,i),yt(i,j)^(-1)) for (i,j) in pairs],
            [com(y(i,j),zt(i,j))*z(j,i)*x(j,i)^(-1) for (i,j) in pairs],
            [com(y(i,j),zt(i,j))*x(j,i)^(-1)*z(j,i) for (i,j) in pairs],
            [com(x(j,i),z(j,i)^(-1)) for (i,j) in pairs],
            [com(yt(i,j), z(i,j))*zt(j,i)*x(i,j) for (i,j) in pairs],
            [com(yt(i,j), z(i,j))*x(i,j)*zt(j,i) for (i,j) in pairs],
            [com(x(i,j),zt(j,i)) for (i,j) in pairs]
    )

    sq_commutators = []
    for i in 1:N
        for j in 1:N
            if i != j
                push!(sq_commutators, com(x(i,j), z(i,j)))
                push!(sq_commutators, com(z(i,j), x(i,j)))

                push!(sq_commutators, com(x(i,j), zt(j,i)))
                push!(sq_commutators, com(zt(j,i), x(i,j)))

                push!(sq_commutators, com(y(i,j), z(j,i)))
                push!(sq_commutators, com(z(j,i), y(i,j)))

                push!(sq_commutators, com(yt(i,j), zt(j,i)))
                push!(sq_commutators, com(zt(j,i), yt(i,j)))

                push!(sq_commutators, com(y(i,j), z(i,j)))
                push!(sq_commutators, com(z(i,j), y(i,j)))

                push!(sq_commutators, com(yt(i,j), zt(i,j)))
                push!(sq_commutators, com(zt(i,j), yt(i,j))) 
            end
        end
    end

    adj_commutators = []
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if i != j && j != k && i != k
                    push!(sq_commutators, com(x(i,j), z(i,k)))
                    push!(sq_commutators, com(z(i,k), x(i,j)))
    
                    push!(sq_commutators, com(x(i,j), zt(j,k)))
                    push!(sq_commutators, com(zt(j,k), x(i,j)))
    
                    push!(sq_commutators, com(y(i,j), z(j,k)))
                    push!(sq_commutators, com(z(j,k), y(i,j)))
    
                    push!(sq_commutators, com(yt(i,j), zt(j,k)))
                    push!(sq_commutators, com(zt(j,k), yt(i,j))) 

                    push!(sq_commutators, com(y(i,j), z(i,k)))
                    push!(sq_commutators, com(z(i,k), y(i,j)))
    
                    push!(sq_commutators, com(yt(i,j), zt(i,k)))
                    push!(sq_commutators, com(zt(i,k), yt(i,j))) # new adj relations from sq
    
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

                    push!(adj_commutators, com(z(i,j), y(j,k)))
                    push!(adj_commutators, com(y(j,k), z(i,j)))
                    push!(adj_commutators, com(z(i,j), yt(j,k)))
                    push!(adj_commutators, com(yt(j,k), z(i,j)))
                    push!(adj_commutators, com(z(i,j), x(j,k)))
                    push!(adj_commutators, com(x(j,k), z(i,j)))

                    push!(adj_commutators, com(zt(i,j), y(j,k)))
                    push!(adj_commutators, com(y(j,k), zt(i,j)))
                    push!(adj_commutators, com(zt(i,j), yt(j,k)))
                    push!(adj_commutators, com(yt(j,k), zt(i,j)))
                    push!(adj_commutators, com(zt(i,j), x(j,k)))
                    push!(adj_commutators, com(x(j,k), zt(i,j)))

                    push!(adj_commutators, com(z(i,k), y(j,k)))
                    push!(adj_commutators, com(y(j,k), z(i,k)))
                    push!(adj_commutators, com(z(i,k), yt(j,k)))
                    push!(adj_commutators, com(yt(j,k), z(i,k)))
                    push!(adj_commutators, com(z(i,k), x(j,k)))
                    push!(adj_commutators, com(x(j,k), z(i,k)))

                    push!(adj_commutators, com(zt(i,k), y(j,k)))
                    push!(adj_commutators, com(y(j,k), zt(i,k)))
                    push!(adj_commutators, com(zt(i,k), yt(j,k)))
                    push!(adj_commutators, com(yt(j,k), zt(i,k)))
                    push!(adj_commutators, com(zt(i,k), x(j,k)))
                    push!(adj_commutators, com(x(j,k), zt(i,k))) # original / adjusted adj relations

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

                        push!(op_commutators, com(yt(i,j), yt(k,l))) # original op relations

                        push!(adj_commutators, com(z(i,l), y(j,k)))
                        push!(adj_commutators, com(y(j,k), z(i,l)))
                        push!(adj_commutators, com(z(i,l), yt(j,k)))
                        push!(adj_commutators, com(yt(j,k), z(i,l)))
                        push!(adj_commutators, com(z(i,l), x(j,k)))
                        push!(adj_commutators, com(x(j,k), z(i,l)))
    
                        push!(adj_commutators, com(zt(i,l), y(j,k)))
                        push!(adj_commutators, com(y(j,k), zt(i,l)))
                        push!(adj_commutators, com(zt(i,l), yt(j,k)))
                        push!(adj_commutators, com(yt(j,k), zt(i,l)))
                        push!(adj_commutators, com(zt(i,l), x(j,k)))
                        push!(adj_commutators, com(x(j,k), zt(i,l))) # new op relations from adj

                        
                    end
                end
            end
        end
    end

    relations_sq = vcat(relations_sq, sq_commutators)
    relations_adj = vcat(relations_adj, adj_commutators)
    relations_op = vcat(op_commutators)

    if sq_adj_ == "sq" 
        return relations_sq
    elseif sq_adj_ == "adj"
        return relations_adj
    elseif sq_adj_ == "op"
        return relations_op
    elseif sq_adj_ == "all"
        return vcat(relations_sq, relations_adj, relations_op)
    end
end