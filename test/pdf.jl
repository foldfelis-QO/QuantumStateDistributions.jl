using LsqFit

@testset "qpdf" begin
    r, θ, dim = 0.8, 2π, 100
    ρ = SqueezedState(r, θ, Matrix, dim=dim)
    d = GaussianStateBHD(ρ)
    xs = LinRange(-3, 3, 1000)
    θs = LinRange(0, 2π, 10)
    g_pdf = zeros(1000)
    g_pdfs = zeros(10, 1000)

    for i in 1:1000
        g_pdf[i] = qpdf(d, π/3, xs[i])
    end

    g_pdfs = qpdf(d, θs, xs)

    param = [QuantumStateDistributions.mean(d, π/3), QuantumStateDistributions.std(d, π/3)]
    @. model(x, p) = 1 / (p[2] * √(2π)) * exp(-(x - p[1])^2 / (2 * p[2]^2))
    p0 = [0.5, 0.5]

    g_fit = curve_fit(model, xs, g_pdf, p0)
    @test coef(g_fit) ≈ param

    for i in 1:10
        g_fits = curve_fit(model, xs, g_pdfs[i, :], p0)
        params = [QuantumStateDistributions.mean(d, θs[i]), QuantumStateDistributions.std(d, θs[i])]

        @test coef(g_fits) ≈ params
    end

    q = QuantumStateBHD(ρ)
    q_pdf = zeros(1000)
    q_pdfs = zeros(10, 1000)

    for i in 1:1000
        q_pdf[i] = qpdf(q, π/3, xs[i])
    end

    q_pdfs = qpdf(q, θs, xs)

    @test q_pdf ≈ g_pdf
    @test q_pdfs ≈ g_pdfs
end
