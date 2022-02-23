@testset "gaussian state bhd" begin
    dim = 100
    ρ = SqueezedState(0.8, π, Matrix, dim=dim)
    d = GaussianStateBHD(ρ)

    θ, x = rand(d)
    @test 0 ≤ θ < 2π

    points = rand(d, 4096)
    @test all(0 ≤ θ < 2π for θ in points[1, :])
end

@testset "operator" begin
    T = Float64
    dim = 100

    ρ = Matrix{T}(I, dim, dim)
    c = QSB.Creation(T, dim=dim)
    a = QSB.Annihilation(T, dim=dim)

    @test QSD.mean_of_create(ρ) ≈ tr(c * ρ)
    @test QSD.mean_of_create²(ρ) ≈ tr(c^2 * ρ)
    @test QSD.mean_of_annihilate(ρ) ≈ tr(a * ρ)
    @test QSD.mean_of_annihilate²(ρ) ≈ tr(a^2 * ρ)
    @test QSD.mean_of_create_annihilate(ρ) ≈ tr(c*a * ρ)

    θ = π/8
    @test QSD.mean_of_π̂ₓ²(ρ, θ) ≈ (
        exp(-2im*θ) * tr(a^2 * ρ) +
        exp(2im*θ) * tr(c^2 * ρ) +
        1 + 2tr(c*a * ρ)
    ) / 4
    @test QSD.mean_of_π̂ₓ(ρ, θ) ≈ (
        exp(-im*θ) * tr(a * ρ) +
        exp(im*θ) * tr(c * ρ)
    ) / 2

    θs = LinRange(0, 2π, 4096)
    @test QSD.mean_of_π̂ₓ²(ρ, θs) ≈ (
        exp.(-2im*θs) .* tr(a^2 * ρ) .+
        exp.(2im*θs) .* tr(c^2 * ρ) .+
        1 .+ 2tr(c*a * ρ)
    ) ./ 4
    @test QSD.mean_of_π̂ₓ(ρ, θs) ≈ (
        exp.(-im*θs) .* tr(a * ρ) .+
        exp.(im*θs) .* tr(c * ρ)
    ) ./ 2
end
