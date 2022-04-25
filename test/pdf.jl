using LsqFit

@testset "qpdf for GaussianStateBHD" begin
    r, θ, dim = 0.8, 2π, 100
    ρ = SqueezedState(r, θ, Matrix, dim=dim)
    d = GaussianStateBHD(ρ)

    m = 10
    n = 1000
    θs = LinRange(0, 2π, m)
    xs = LinRange(-3, 3, n)

    g_pdf = Vector{Float64}(undef, n)

    for i in 1:n
        g_pdf[i] = qpdf(d, π/3, xs[i])
    end

    g_pdfs = qpdf(d, θs, xs)

    moments = [QuantumStateDistributions.mean(d, π/3), QuantumStateDistributions.std(d, π/3)]
    @. model(x, p) = 1 / (p[2] * √(2π)) * exp(-(x - p[1])^2 / (2 * p[2]^2))
    p0 = [0.5, 0.5]

    g_fit = curve_fit(model, xs, g_pdf, p0)
    @test coef(g_fit) ≈ moments

    for i in 1:m
        g_fits = curve_fit(model, xs, g_pdfs[i, :], p0)
        moments = [QuantumStateDistributions.mean(d, θs[i]), QuantumStateDistributions.std(d, θs[i])]

        @test coef(g_fits) ≈ moments
    end
end

@testset "qpdf for QuantumStateBHD" begin
    r, θ, dim = 0.8, 2π, 100
    ρ = SqueezedState(r, θ, Matrix, dim=dim)
    q = QuantumStateBHD(ρ)
    d = GaussianStateBHD(ρ)

    m = 10
    n = 1000
    θs = LinRange(0, 2π, m)
    xs = LinRange(-3, 3, n)

    q_pdf = Vector{Float64}(undef, n)

    for i in 1:n
        q_pdf[i] = qpdf(q, π/3, xs[i])
    end

    q_pdfs = qpdf(q, θs, xs)

    moments = [QuantumStateDistributions.mean(d, π/3), QuantumStateDistributions.std(d, π/3)]
    @. model(x, p) = 1 / (p[2] * √(2π)) * exp(-(x - p[1])^2 / (2 * p[2]^2))
    p0 = [0.5, 0.5]

    q_fit = curve_fit(model, xs, q_pdf, p0)
    @test coef(q_fit) ≈ moments

    for i in 1:m
        q_fits = curve_fit(model, xs, q_pdfs[i, :], p0)
        moments = [QuantumStateDistributions.mean(d, θs[i]), QuantumStateDistributions.std(d, θs[i])]

        @test coef(q_fits) ≈ moments
    end
end

@testset "pdf for density matrix" begin
    r, θ, dim = 0.8, 2π, 100
    ρ = SqueezedState(r, θ, Matrix, dim=dim)
    d = GaussianStateBHD(ρ)

    m = 10
    n = 1000
    θs = LinRange(0, 2π, m)
    xs = LinRange(-3, 3, n)

    ρ_pdf = Vector{Float64}(undef, n)

    for i in 1:n
        ρ_pdf[i] = qpdf(ρ, π/3, xs[i])
    end

    ρ_pdfs = qpdf(ρ, θs, xs)

    moments = [QuantumStateDistributions.mean(d, π/3), QuantumStateDistributions.std(d, π/3)]
    @. model(x, p) = 1 / (p[2] * √(2π)) * exp(-(x - p[1])^2 / (2 * p[2]^2))
    p0 = [0.5, 0.5]

    ρ_fit = curve_fit(model, xs, ρ_pdf, p0)
    @test coef(ρ_fit) ≈ moments

    for i in 1:m
        ρ_fits = curve_fit(model, xs, ρ_pdfs[i, :], p0)
        moments = [QuantumStateDistributions.mean(d, θs[i]), QuantumStateDistributions.std(d, θs[i])]

        @test coef(ρ_fits) ≈ moments
    end
end
