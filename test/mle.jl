# @testset "mle" begin

# end

using QuantumStateBase, QuantumStateDistributions, Plots

function gen_data()
    ρ = SqueezedState(0.8, π/4, Matrix, dim=100)
    d = GaussianStateBHD(ρ)

    return rand(d, 40960)
end

function main()
    data = gen_data()
    ppit = PositivePreservingIterator(data, 35, 50)
    run!(ppit)

    ρ = ppit.ρ

    # w = wigner(ρ, LinRange(-3, 3, 101), LinRange(-3, 3, 101))

    # lim = maximum(abs.(w.𝐰_surface))
    # heatmap(w.𝐰_surface, color=:coolwarm, clim=(-lim, lim))

    return ρ
end
