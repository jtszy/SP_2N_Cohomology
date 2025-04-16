function conjugation(
    l, # l is an element of Sp_2n
    σ::PermutationGroups.AbstractPerm
)
    l_matrix = MatrixGroups.matrix_repr(l)
    N = size(l_matrix)[1]
    n = div(N,2)

    for i in 1:N
        for j in 1:N
            if l_matrix[i,j] < 0 && (i <= n || j <= n)
                return elementary_conjugation(l^(-1), σ)^(-1)
            end
        end
    end

    return elementary_conjugation(l,σ)
    
end

function elementary_conjugation(
    l, # l is an element of Sp_2n
    σ::PermutationGroups.AbstractPerm 
)
    l_matrix = MatrixGroups.matrix_repr(l)
    N = size(l_matrix)[1]
    n = div(N,2)

    type, l_i, l_j, l_transpose = elementary_element_info(l_matrix)

    if type == :A
        return Groups.MatrixGroups.ElementarySymplectic{N}(type, l_i^σ, l_j^σ)
    end

    if l_transpose == 0
        i = max(l_i^σ, l_j^σ)
        j = min(l_i^σ, l_j^σ) + n
        return Groups.MatrixGroups.ElementarySymplectic{N}(type, i, j)
    
    else
        i = min(l_i^σ, l_j^σ) + n
        j = max(l_i^σ, l_j^σ)
        return Groups.MatrixGroups.ElementarySymplectic{N}(type, i, j)
    end

end

function elementary_element_info(l_matrix)
    N = size(l_matrix)[1]
    n = div(N,2)
    type = :B
    l_transpose = 1
    l_i, l_j = -1, -1

    for i in 1:N
        for j in 1:N

            if l_matrix[i,j] < 0
                type = :A
            end

            if l_matrix[i,j] == 1 && i != j
                l_i = i % n
                l_j = j % n

                l_i = (l_i == 0) ? n : l_i
                l_j = (l_j == 0) ? n : l_j

                if i < j
                    l_transpose = 0
                end
            end
        end
    end

    return (type, l_i, l_j, l_transpose)
end

function wedderburn_data(basis, half_basis, S)
    @time begin
        N = div(size(first(S))[1], 2)

        Σ = PermutationGroups.SymmetricGroup(N)

        actions = LowCohomologySOS.WedderburnActions(
            alphabet(parent(first(S))), 
            Σ, 
            SP_2N_Cohomology.conjugation, 
            S, 
            basis
        )
        constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(
            basis, 
            half_basis, 
            S
        )
    end
    
    return constraints_basis, psd_basis, Σ, actions
end