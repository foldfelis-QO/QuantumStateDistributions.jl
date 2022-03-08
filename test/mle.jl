@testset "mle" begin
    # gen data
    ρ = SqueezedState(0.8, π/4, Matrix, dim=100)
    d = GaussianStateBHD(ρ)
    data = rand(d, 40960)

    # mle
    ppit = PositivePreservingIterator(data, 50, dim=35)
    run!(ppit)
end
