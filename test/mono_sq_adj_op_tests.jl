@testset "mono_sq_adj_op" begin
    N = 4
    Sp_2N = MatrixGroups.SymplecticGroup{2*N}(Int8)
    F_Sp_2N_Steinberg = FreeGroup(alphabet(Sp_2N))
    S = gens(Sp_2N)

    Steinberg_relations = SP_4_Cohomology.relations_St(
        F_Sp_2N_Steinberg, 
        S, 
        N
    )
    quotient_hom_Steinberg = let source = F_Sp_2N_Steinberg, target = Sp_2N
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end
    support_jacobian, min_support = SP_4_Cohomology.symplectic_min_supports(
        Steinberg_relations, 
        quotient_hom_Steinberg, 
        S
    )
    Δ₁, I_N_whole, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(
        quotient_hom_Steinberg, 
        Steinberg_relations, 
        support_jacobian
    );

    gen_dict = Dict(LowCohomologySOS.determine_letter(S[i]) => i for i in eachindex(S))
    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,max(i,j),min(i,j) + N)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,min(i,j) + N,max(i,j))]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i,i + N)]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{2*N}(:B,i + N,i)]

    RG = LowCohomologySOS.group_ring(Sp_2N, min_support, star_multiplication = false)
    Δ₁⁻ = LowCohomologySOS.embed.(identity, Δ₁⁻, Ref(RG))
    gen_dict_2 = Dict(LowCohomologySOS.determine_letter(s) => RG(s) for s in S)
    X(i,j) = gen_dict_2[MatrixGroups.ElementarySymplectic{2*N}(:A,i,j)]
    Y(i,j) = gen_dict_2[MatrixGroups.ElementarySymplectic{2*N}(:B,max(i,j),min(i,j) + N)]
    Yt(i,j) = gen_dict_2[MatrixGroups.ElementarySymplectic{2*N}(:B,min(i,j) + N,max(i,j))]
    Z(i) = gen_dict_2[MatrixGroups.ElementarySymplectic{2*N}(:B,i,i + N)]
    Zt(i) = gen_dict_2[MatrixGroups.ElementarySymplectic{2*N}(:B,i + N,i)]

   mono, sq, adj, op = SP_4_Cohomology.mono_sq_adj_op(Δ₁⁻, S)

   for i in 1:20
        i1, j1 = rand(1:N), rand(1:N)
        while i1 == j1
            i1, j1 = rand(1:N), rand(1:N)
        end
        i2, j2 = rand(1:N), rand(1:N)
        while i2 == j2
            i2, j2 = rand(1:N), rand(1:N)
        end
        if length(union!([i1,j1],[i2,j2])) == 2
            @test sq[x(i1,j1),x(i2,j2)] == (one(RG)-X(i1,j1))*(one(RG)-X(i2,j2))'
            @test mono[x(i1,j1),x(i2,j2)] == adj[x(i1,j1),x(i2,j2)] == op[x(i1,j1),x(i2,j2)] == zero(RG)
        elseif length(union!([i1,j1],[i2,j2])) == 3
            @test adj[x(i1,j1),x(i2,j2)] == (one(RG)-X(i1,j1))*(one(RG)-X(i2,j2))'
            @test mono[x(i1,j1),x(i2,j2)] == sq[x(i1,j1),x(i2,j2)] == op[x(i1,j1),x(i2,j2)] == zero(RG)
        else
            length(union!([i1,j1],[i2,j2])) == 4
            @test op[x(i1,j1),x(i2,j2)] == (one(RG)-X(i1,j1))*(one(RG)-X(i2,j2))'
            @test mono[x(i1,j1),x(i2,j2)] == sq[x(i1,j1),x(i2,j2)] == adj[x(i1,j1),x(i2,j2)] == zero(RG)
        end
        @test sq[z(i1),y(i1,j1)] == (one(RG)-Z(i1))*(one(RG)-Y(i1,j1))'
        @test sq[zt(i1),yt(i1,j1)] == (one(RG)-Zt(i1))*(one(RG)-Yt(i1,j1))'
        @test mono[z(i1),zt(i1)] == (one(RG)-Z(i1))*(one(RG)-Zt(i1))'
   end
end