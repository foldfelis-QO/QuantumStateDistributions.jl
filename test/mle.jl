@testset "mle" begin
    # gen data
    ρ = SqueezedState(0.8, π/4, Matrix, dim=100)
    d = GaussianStateBHD(ρ)
    n = 8192
    data = rand(d, n)

    # mle
    t = 50
    ppit = PositivePreservingIterator(data, t, dim=35)
    run!(ppit)
end
