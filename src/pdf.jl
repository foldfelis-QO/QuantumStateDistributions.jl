export
    qpdf,
    qpdf!

"""
    q_pdf([T=Float64], ρ::AbstractArray, θ::Real, x::Real)

Quadrature prabability at point (θ, x)

    q_pdf([T=Float64], ρ::AbstractArray, θs::AbstractRange, xs::AbstractRange)

Quadrature prabability at points (θs, xs)
"""
qpdf(ρ, θ, x) = qpdf(Float64, ρ, θ, x)

function qpdf(T::Type{<:Real}, ρ, θ::Real, x::Real)
    dim = size(ρ, 1)
    𝛑̂_res = Matrix{Complex{T}}(undef, dim, dim)

    return qpdf!(𝛑̂_res, ρ, θ, x)
end

function qpdf!(𝛑̂_res::AbstractMatrix, ρ::AbstractArray, θ::Real, x::Real)
    dim = size(ρ, 1)

    return real_tr_mul(𝛑̂!(𝛑̂_res, θ, x, dim=dim), ρ)
end

function qpdf(T::Type{<:Real}, ρ, θs::AbstractRange, xs::AbstractRange)
    dim = size(ρ, 1)
    𝛑̂_res_vec = [Matrix{Complex{T}}(undef, dim, dim) for _ in 1:Threads.nthreads()]
    𝐩 = Matrix{T}(undef, length(θs), length(xs))

    return qpdf!(𝛑̂_res_vec, 𝐩, ρ, θs, xs)
end

function qpdf!(𝛑̂_res_vec::AbstractVector{Matrix{Complex{T}}}, 𝐩::Matrix{T}, ρ::AbstractArray, θs::AbstractRange, xs::AbstractRange) where {T}
    @sync for (j, x) in enumerate(xs)
        for (i, θ) in enumerate(θs)
            Threads.@spawn 𝐩[i, j] = qpdf!(𝛑̂_res_vec[Threads.threadid()], ρ, θ, x)
        end
    end

    return 𝐩
end

##### for arb. state in intensity-to-measurement-phase quadrature coordinate #####

# |θ, x⟩ = ∑ₙ |n⟩ ⟨n|θ, x⟩ = ∑ₙ ψₙ(θ, x) |n⟩
# ⟨n|θ, x⟩ = ψₙ(θ, x) = exp(im n θ) (2/π)^(1/4) exp(-x^2) Hₙ(√2 x)/√(2^n n!)
function ψₙ(n::Integer, θ::Real, x::Real)
    return (2/π)^(1/4) * exp(im*n*θ - x^2) * hermiteh(n, sqrt(2)x) / sqrt(2^n * factorial(n))
end

function 𝛑̂!(result::AbstractMatrix{<:Complex}, θ::Real, x::Real; dim)
    view(result, :, 1) .= ψₙ.(big.(0:dim-1), θ, x)
    result .= view(result, :, 1) * view(result, :, 1)'

    return result
end

function 𝛑̂(T::Type{<:Complex}, θ::Real, x::Real; dim)
    result = Matrix{T}(undef, dim, dim)

    return 𝛑̂!(result, θ, x, dim=dim)
end

𝛑̂(θ::Real, x::Real; dim) = 𝛑̂(ComplexF64, θ, x, dim=dim)

# #########
# # utils #
# #########

real_tr_mul(𝐚, 𝐛) = sum(real(𝐚[i, :]' * 𝐛[:, i]) for i in 1:size(𝐚, 1))
