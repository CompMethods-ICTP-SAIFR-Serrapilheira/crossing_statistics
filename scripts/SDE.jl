# Final project for the 2022 QBio Computational Methods course
# ---------------------------------
# In this module, named SDE.jl, we implement integration of trajectories on two stochastic differential equations:
# The Ornstein-Uhlenbeck (OU) process, and Brownian motion.

module SDE

export IsoOU2D, Brownian2D, StochasticProcess, updateSP, trajectorySP

include("PRNG.jl")
using .PRNG

abstract type StochasticProcess end

mutable struct IsoOU2D <: StochasticProcess
    # Represents the state of an isotropic two-dimensional OU process

    x::Tuple{Float64, Float64} # current position
    const λ::Tuple{Float64, Float64} # attraction center
    const τ::Float64 # return time-scale
    const g::Float64 # noise amplitude
end

mutable struct Brownian2D <: StochasticProcess
    # Represents the state of a two-dimensional Brownian motion

    x::Tuple{Float64, Float64} # current position
    const g::Float64 # noise amplitude
end

function updateSP(rng::LCG, X::StochasticProcess, h::Float64)::Nothing

    if typeof(X) == IsoOU2D
        u = rand2vecn(rng)

        X.x = (exp(-h/X.τ) .* X.x) .+ ((1-exp(-h/X.τ)) .* X.λ) .+ (√(X.g*X.τ/2 * (1 - exp(-2h/X.τ))) .* u)
    elseif typeof(X) == Brownian2D
        u = rand2vecn(rng)

        X.x = X.x .+ (√(X.g*h) .* u)
    end

    nothing
end

function trajectorySP(rng::LCG, X::StochasticProcess, h::Float64, maxT::Float64)::Vector{Tuple{Float64,Float64,Float64}}
    N = floor(Int, maxT/h)
    out = Array{Tuple{Float64,Float64,Float64}}(undef, N)
    for i in 0:N
        out[i] = (i*h, X.x[1],X.x[2])
        updateSP(rng, X, h)
    end

    out
end
end