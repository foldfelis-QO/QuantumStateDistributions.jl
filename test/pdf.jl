@testset "qpdf" begin
    r, θ, dim = 0.8, 2π, 100
    ρ = SqueezedState(r, θ, Matrix, dim=dim)
    d = GaussianStateBHD(ρ)

    θs = LinRange(0, 2π, 101)
    xs = LinRange(-3, 3, 101)

    @test size(qpdf(d, θs, xs)) == (101, 101)
    @test length(qpdf(d, π/2, 1)) == 1
    @test size(qpdf(Float32, d, θs, xs)) == (101, 101)
    @test length(qpdf(Float32, d ,π/3, 2)) == 1
    @test size(qpdf(d, LinRange(0, 2π, 10), LinRange(-3, 3, 10))) == (10, 10)
end
