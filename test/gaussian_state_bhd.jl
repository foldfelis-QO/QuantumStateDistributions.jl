@testset "operator" begin
    T = Float64
    dim = 100

    ρ = Matrix{T}(I, dim, dim)
    c = QSB.Creation(T, dim)
    a = QSB.Annihilation(T, dim)

    @test QSD.mean_of_create(ρ) ≈ tr(c * ρ)
    @test QSD.mean_of_create²(ρ) ≈ tr(c^2 * ρ)
    @test QSD.mean_of_annihilate(ρ) ≈ tr(a * ρ)
    @test QSD.mean_of_annihilate²(ρ) ≈ tr(a^2 * ρ)
    @test QSD.mean_of_create_annihilate(ρ) ≈ tr(c*a * ρ)

    θs = LinRange(0, 2π, 100)
    @test QSD.mean_of_π̂ₓ²(θs, ρ) ≈ (
        exp.(-2im*θs) .* tr(a^2 * ρ) .+
        exp.(2im*θs) .* tr(c^2 * ρ) .+
        1 .+ 2tr(c*a * ρ)
    ) ./ 4
    @test QSD.mean_of_π̂ₓ(θs, ρ) ≈ (
        exp.(-im*θs) .* tr(a * ρ) .+
        exp.(im*θs) .* tr(c * ρ)
    ) ./ 2
end
