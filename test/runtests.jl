using Test

include("../src/mol_dyn_md.jl")
using .MolDyn

function setup_test_environment()
    r_ab_eq_hcl = 1.57e-10
    num_steps = 1000
    qs = zeros(Float64, num_steps, 2, 3)
    qs[1, 2, :] = [r_ab_eq_hcl, 0.0, 0.0]
    return Dict(:qs => qs, :r_ab_eq_hcl => r_ab_eq_hcl)
end

@testset "Example Test Set" begin
    env = setup_test_environment()
    result = u_stretch(env[:qs][1,:,:], env[:qs][2,:,:], 1.0, env[:r_ab_eq_hcl])
    @test result ≈ 0.0 atol = 0.1
end

