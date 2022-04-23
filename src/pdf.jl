export
    qpdf,
    qpdf!

"""
    qpdf([T=Float64], d::QuantumStateBHD, θ::Real, x::Real)
    qpdf([T=Float64], d::QuantumStateBHD, θs::AbstractRange, xs::AbstractRange)
    qpdf([T=Float64], d::GaussianStateBHD, θ::Real, x::Real)
    qpdf([T=Float64], d::GaussianStateBHD, θs::AbstractRange, xs::AbstractRange)

Quadrature prabability at point (θ, x) or points (θs, xs)
"""
qpdf(d, θ, x) = qpdf(Float64, d, θ, x)

function qpdf(T::Type{<:Real}, d::GaussianStateBHD, θ::Real, x::Real)
    μ = QuantumStateDistributions.mean(d, θ)
    σ = QuantumStateDistributions.std(d, θ)

    return T(pdf(Normal(μ, σ), x))
end

function qpdf(T::Type{<:Real}, d::GaussianStateBHD, θs::AbstractRange, xs::AbstractRange)
    m, n = length(θs), length(xs)

    μs = Vector{T}(undef, m)
    σs = Vector{T}(undef, m)
    gaussians = Vector{Normal}(undef, m)
    𝐩 = Matrix{T}(undef, m, n)

    return qpdf!(μs, σs, gaussians, 𝐩, d, θs, xs)
end

function qpdf!(μs, σs, gaussians, 𝐩, d::GaussianStateBHD, θs::AbstractRange, xs::AbstractRange)
    μs .= QuantumStateDistributions.mean(d, θs)
    σs .= QuantumStateDistributions.std(d, θs)

    gaussians .= Normal.(μs, σs)

    for i in 1:length(θs)
        for (j, x) in enumerate(xs)
		    𝐩[i, j] = pdf(gaussians[i], x)
        end
    end

    return 𝐩
end

qpdf(T::Type{<:Real}, d::QuantumStateBHD, θ, x) = qpdf(T, d.ρ, θ, x)

"""
    qpdf([T=Float64], ρ::AbstractArray, θ::Real, x::Real)
    qpdf([T=Float64], ρ::AbstractArray, θs::AbstractRange, xs::AbstractRange)

Quadrature prabability in intensity-to-measurement-phase quadrature coordinate.

``p(\\rho, \\theta, x) = tr(\\hat{\\Pi}(\\theta, x) \\rho)``
"""
qpdf(ρ::AbstractArray, θ, x) = qpdf(Float64, ρ, θ, x)

function qpdf(T::Type{<:Real}, ρ::AbstractArray, θ::Real, x::Real)
    dim = size(ρ, 1)
    𝛑̂_res = Matrix{Complex{T}}(undef, dim, dim)

    return qpdf!(𝛑̂_res, ρ, θ, x)
end

function qpdf!(𝛑̂_res::AbstractMatrix, ρ::AbstractArray, θ::Real, x::Real)
    dim = size(ρ, 1)

    return real_tr_mul(𝛑̂!(𝛑̂_res, θ, x, dim=dim), ρ)
end

function qpdf(T::Type{<:Real}, ρ::AbstractArray, θs::AbstractRange, xs::AbstractRange)
    dim = size(ρ, 1)
    𝛑̂_res_vec = [Matrix{Complex{T}}(undef, dim, dim) for _ in 1:Threads.nthreads()]
    𝐩 = Matrix{T}(undef, length(θs), length(xs))

    return qpdf!(𝛑̂_res_vec, 𝐩, ρ, θs, xs)
end

function qpdf!(
    𝛑̂_res_vec::AbstractVector{Matrix{Complex{T}}}, 𝐩::Matrix{T},
    ρ::AbstractArray, θs::AbstractRange, xs::AbstractRange
) where {T}
    @sync for (j, x) in enumerate(xs)
        for (i, θ) in enumerate(θs)
            Threads.@spawn 𝐩[i, j] = qpdf!(𝛑̂_res_vec[Threads.threadid()], ρ, θ, x)
        end
    end

    return 𝐩
end

##### for arb. state in intensity-to-measurement-phase quadrature coordinate #####

"""
    ψₙ(n::Integer, θ::Real, x::Real)

Eigenstate of BHD measurement operator.

``\\psi_n(\\theta, x) = \\langle n | \\theta, x \\rangle``
"""
function ψₙ(n::Integer, θ::Real, x::Real)
    # |θ, x⟩ = ∑ₙ |n⟩ ⟨n|θ, x⟩ = ∑ₙ ψₙ(θ, x) |n⟩
    # ⟨n|θ, x⟩ = ψₙ(θ, x) = exp(im n θ) (2/π)^(1/4) exp(-x^2) Hₙ(√2 x)/√(2^n n!)

    return (2/π)^(1/4) * exp(im*n*θ - x^2) * hermiteh(n, sqrt(2)x) / sqrt(2^n * factorial(n))
end

"""
    𝛑̂(θ::Real, x::Real; dim::Integer)

BHD measurement operator.

``\\hat{\\Pi}_{m, n}(\\theta, x) = \\langle m | \\hat{\\Pi}(\\theta, x) | n \\rangle = \\langle m | \\theta, x \\rangle \\langle \\theta, x | n \\rangle``
"""
𝛑̂(θ::Real, x::Real; dim) = 𝛑̂(ComplexF64, θ, x, dim=dim)

function 𝛑̂(T::Type{<:Complex}, θ::Real, x::Real; dim)
    result = Matrix{T}(undef, dim, dim)

    return 𝛑̂!(result, θ, x, dim=dim)
end

function 𝛑̂!(result::AbstractMatrix{<:Complex}, θ::Real, x::Real; dim)
    view(result, :, 1) .= ψₙ.(big.(0:dim-1), θ, x)
    result .= view(result, :, 1) * view(result, :, 1)'

    return result
end

# #########
# # utils #
# #########

real_tr_mul(𝐚, 𝐛) = sum(real(𝐚[i, :]' * 𝐛[:, i]) for i in 1:size(𝐚, 1))
