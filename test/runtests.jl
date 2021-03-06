using QuantumStateDistributions
using QuantumStateBase
using LinearAlgebra
using Test

const QSB = QuantumStateBase
const QSD = QuantumStateDistributions

@testset "QuantumStateDistributions.jl" begin
    include("pdf.jl")
    include("mle.jl")
    include("gaussian_state_bhd.jl")
end
