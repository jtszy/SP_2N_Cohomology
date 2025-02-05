using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS ÷ 2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS ÷ 2)

using Groups
using IntervalArithmetic
using JuMP
using LowCohomologySOS
using SCS
using SP_4_Cohomology
using SparseArrays

sp4 = MatrixGroups.SymplecticGroup{4}(Int8)

F_sp4_Behr = FreeGroup(6)
x_a, x_b, x_ab, x_2ab, w_a, w_b = gens(F_sp4_Behr)

Behr_matrix = Dict(
    x_a => [
        1 1 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 -1 1
    ],
    x_b => [
        1 0 0 0;
        0 1 0 1;
        0 0 1 0;
        0 0 0 1
    ],
    x_ab => [
        1 0 0 1;
        0 1 1 0;
        0 0 1 0;
        0 0 0 1
    ],
    x_2ab => [
        1 0 1 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1
    ],
    w_a => [
        0 -1 0 0;
        1 0 0 0;
        0 0 0 -1;
        0 0 1 0
    ],
    w_b => [
        1 0 0 0;
        0 0 0 -1;
        0 0 1 0;
        0 1 0 0
    ]
)

gens_sp4_inv = [gens(sp4); inv.(gens(sp4))]

Ball_4, sizes = Groups.wlmetric_ball(gens_sp4_inv, radius=4)

Behr_group = Dict()
for g in Ball_4
    for gen in keys(Behr_matrix)
        if MatrixGroups.matrix_repr(g) == Behr_matrix[gen]
            Behr_group[gen] = g
        end
    end
end

quotient_hom_Behr = SP_4_Cohomology.quotient_homomorphism(F_sp4_Behr, sp4, Behr_group)

Behr_relations = [
    x_a * x_b * x_a^(-1) * x_b^(-1) * x_2ab^(-1) * x_ab^(-1),
    x_a * x_ab * x_a^(-1) * x_ab^(-1) * x_2ab^(-2),
    x_a * x_2ab * x_a^(-1) * x_2ab^(-1),
    x_b * x_ab * x_b^(-1) * x_ab^(-1),
    x_b * x_2ab * x_b^(-1) * x_2ab^(-1),
    x_ab * x_2ab * x_ab^(-1) * x_2ab^(-1),
    w_a * w_b^2 * w_a * w_b^(-2),
    w_b * w_a^2 * w_b^(-1) * w_a^(-2),
    (w_a * w_b)^2 * (w_b * w_a)^(-2),
    w_b^4,
    w_a * x_b * w_a^(-1) * x_2ab^(-1),
    w_a * x_2ab * w_a^(-1) * x_b^(-1),
    w_a * x_ab * w_a^(-1) * x_ab,
    w_b * x_a * w_b^(-1) * x_ab^(-1),
    w_b * x_ab * w_b^(-1) * x_a,
    w_b * x_2ab * w_b^(-1) * x_2ab^(-1),
    w_a * x_a * w_a^(-1) * x_a * w_a * x_a,
    w_b * x_b * w_b^(-1) * x_b * w_b * x_b,
]

for r in Behr_relations
    @assert quotient_hom_Behr(r) == one(sp4)
end

jacobian_matrix_Behr = LowCohomologySOS.jacobian_matrix(Behr_relations)

min_support = SP_4_Cohomology.minimalistic_support(Behr_relations, quotient_hom_Behr)

Δ₁, I_ = LowCohomologySOS.spectral_gap_elements(
    quotient_hom_Behr,
    Behr_relations,
    min_support,
)

sos_problem_Behr = LowCohomologySOS.sos_problem(
    Δ₁,
    I_,
)

JuMP.set_optimizer(sos_problem_Behr, SP_4_Cohomology.scs_opt(eps=1e-7, max_iters=15000))

JuMP.optimize!(sos_problem_Behr)

λ, Q = LowCohomologySOS.get_solution(sos_problem_Behr)

result_bool, result = LowCohomologySOS.certify_sos_decomposition(
    Δ₁,
    I_,
    λ,
    Q,
    min_support,
)