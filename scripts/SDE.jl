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
    # Evolves a stochastic process forward a time-step
    # rng : LCG random number generator
    # X : stochastic process to evolve
    # h : time step

    if typeof(X) == IsoOU2D # Ornstein-Uhlenbeck case
        u = rand2vecn(rng)

        X.x = (exp(-h/X.τ) .* X.x) .+ ((1-exp(-h/X.τ)) .* X.λ) .+ (√(X.g*X.τ/2 * (1 - exp(-2h/X.τ))) .* u)
    elseif typeof(X) == Brownian2D # Brownian motion case
        u = rand2vecn(rng)

        X.x = X.x .+ (√(X.g*h) .* u)
    end

    nothing
end

function trajectorySP(rng::LCG, X::StochasticProcess, h::Float64, maxT::Float64)::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    # Generates a trajectory of an stochastic process
    # rng : LCG random number generator
    # X : stochastic process to evolve
    # h : time step size
    # maxT : total time
    # output : A triplet of arrays ([x1,...],[y1,..],[t1,...]) representing the trajetory in time
    N = floor(Int, maxT/h) # number of time steps
    outX = Array{Float64}(undef, N) # initialize output arrays
    outY = Array{Float64}(undef, N)
    outT = Array{Float64}(undef, N) 
    for i in 1:N
        outT[i] = (i-1)*h
        outX[i] = X.x[1]
        outY[i] = X.x[2] # record time and position
        updateSP(rng, X, h) # evolve the process by a time step
    end

    (outX, outY, outT)
end
end