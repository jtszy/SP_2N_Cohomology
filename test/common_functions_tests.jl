@testset "com" begin
    F_2 = Groups.FreeGroup(2)
    a, b = Groups.gens(F_2)
    Z2 = Groups.FPGroup(F_2, [a * b => b * a])
    A, B = Groups.gens(Z2)

    @test SP_2N_Cohomology.com(A, B) == one(Z2)

    sl_3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
    e12, e13, e21, e23, e31, e32 = Groups.gens(sl_3)

    @test SP_2N_Cohomology.com(e12, e13) == one(sl_3)
    @test SP_2N_Cohomology.com(e12, e32) == one(sl_3)
    @test SP_2N_Cohomology.com(e13, e23) == one(sl_3)
    @test SP_2N_Cohomology.com(e21, e23) == one(sl_3)
    @test SP_2N_Cohomology.com(e21, e31) == one(sl_3)
    @test SP_2N_Cohomology.com(e31, e32) == one(sl_3)
    @test SP_2N_Cohomology.com(e12, e23) == e13
    @test SP_2N_Cohomology.com(e13, e32) == e12
    @test SP_2N_Cohomology.com(e21, e13) == e23
    @test SP_2N_Cohomology.com(e23, e31) == e21
    @test SP_2N_Cohomology.com(e31, e12) == e32
    @test SP_2N_Cohomology.com(e32, e21) == e31
end

@testset "quotient_homomorphism" begin
    F_2 = Groups.FreeGroup(2)
    a, b = Groups.gens(F_2)
    Z2 = Groups.FPGroup(F_2, [a * b => b * a])
    A, B = Groups.gens(Z2)

    test_group = Dict()
    test_group[a] = A
    test_group[b] = B

    hom_Z2 = SP_2N_Cohomology.quotient_homomorphism(F_2, Z2, test_group)

    @test hom_Z2(a) == A
    @test hom_Z2(b) == B
    @test hom_Z2(a * b) == A * B
    @test hom_Z2(b * a) == A * B
end

@testset "symplectic_min_supports" begin
    Sp_4 = MatrixGroups.SymplecticGroup{4}(Int8)
    F_Sp_4_Steinberg = FreeGroup(alphabet(Sp_4))

    S = gens(Sp_4)
    gen_dict = Dict(LowCohomologySOS.determine_letter(S[i]) => gens(F_Sp_4_Steinberg, i) for i in eachindex(S))
    x(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:A,i,j)]
    y(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,max(i,j),min(i,j) + 2)]
    yt(i,j) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,min(i,j) + 2,max(i,j))]
    z(i) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,i,i + 2)]
    zt(i) = gen_dict[MatrixGroups.ElementarySymplectic{4}(:B,i + 2,i)]

    selected_relations = [
        z(1)*y(1,2)*z(1)^(-1)*y(1,2)^(-1),
        zt(2)*yt(1,2)^(-1)*zt(2)^(-1)*yt(1,2),
        x(2,1)*zt(1)*x(2,1)^(-1)*zt(1)^(-1)
    ]
    quotient_hom_Steinberg = let source = F_Sp_4_Steinberg, target = Sp_4
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    support_jacobian, min_support = SP_2N_Cohomology.symplectic_min_supports(selected_relations, quotient_hom_Steinberg, S)

    gen_dict_2 = Dict(LowCohomologySOS.determine_letter(s) => s for s in S)
    X(i,j) = gen_dict_2[MatrixGroups.ElementarySymplectic{4}(:A,i,j)]
    Y(i,j) = gen_dict_2[MatrixGroups.ElementarySymplectic{4}(:B,max(i,j),min(i,j) + 2)]
    Yt(i,j) = gen_dict_2[MatrixGroups.ElementarySymplectic{4}(:B,min(i,j) + 2,max(i,j))]
    Z(i) = gen_dict_2[MatrixGroups.ElementarySymplectic{4}(:B,i,i + 2)]
    Zt(i) = gen_dict_2[MatrixGroups.ElementarySymplectic{4}(:B,i + 2,i)]
    I = one(Sp_4)

    test_min_support = vcat(
        [
            I,
            Y(1,2),
            Z(1)
        ],
        [
            I,
            Zt(2),
            Zt(2)*Yt(1,2)^(-1),
            Yt(1,2)^(-1)
        ],
        [
            I,
            X(2,1),
            Zt(1)
        ],
        S
    )
    test_min_support = unique!([test_min_support;inv.(test_min_support)])

    @info length(min_support), length(test_min_support)

    @test Set(min_support) == Set(test_min_support)
end
