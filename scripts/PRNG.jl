# Final project for the 2022 QBio Computational Methods course
# ---------------------------------
# In this module, named PRNG.jl, we implement some very basic pseudo-random number generation.
# In particular, we're aiming at implementing the generation of a 2D standard normal distirbution.
# 
# We take the simplest route: Implement a U[0,1] distribution with a linear congruential generator.
# The 2D standard normal distribution is then obtained from a U[0,1] by the Box-Muller transformation.

module PRNG

export LCG, rand, rand2vecn

mutable struct LCG # Linear Congruential Generator
    # The linear congruential generator (LCG) is a simple pseudo-random number generator.
    # It works by producing a sequence of integers x(i+1) = a*x(i) + c mod M
    # Here x(i), a, c, M are integers. The initialization "seed" is x(0).
    # In a computer, the modulus operation is automatically done by integer overflow, so that
    # instead of choosing an M, we will work with the maximum representable value of an 64-bit integer.
    
    # Since x(i) is a number between 0 and M - 1, if x(i) is uniform in this integer range
    # we would find, for large M, for x(i)/M to be uniformly distributed in [0,1].
    # It is known that the LCG has a series of flaws (e.g. Marsaglia's theorem)
    # but for a long time it was the standard PRNG for most applications.
    
    x::Int64 # current state
    const a::Int64
    const c::Int64

    LCG(seed::Int64) = new(seed, 2862933555777941757, 3037000493) # Initializing to some "good" default values
end

function rand(rng::LCG, n::Int64 = 1)
    # Producing (sequences of n) U[0,1] pseudo-random samples from an LCG.   
    # rng : The linear congruential generator to use.
    # n : Size of the sample (if n > 1 returns an array of independent U[0,1] samples)
    # output: if n == 1 u::Float64 which is U[0,1]
    #         if n > 1, an Array{Float64, n} of independent U[0,1] samples.
    # ------------------------
    # The procedure is: we update the current state x(i+1) = a*x(i) + c mod M
    # At first glance, if would suffice to take u = x(i+1)/M.
    # Note however, due to integer overflow, x(i+1) can be negative. If that's the case,
    # we instead use x(i+1) = x(i) + M to generate the float u = x(i+1)/M.
    if n == 1
        rng.x = (rng.a*rng.x + rng.c)
        u = rng.x * 2.0^(-64) # normalize integer to the interval [-1, 1]
        u < 0 ? u + 1 : u # return u in [0,1]
    else
        nn = Array{Float64}(undef, n) # if n > 1 we allocate an n-element array
        for i in 1:n # same as above for each sample
            rng.x = (rng.a*rng.x + rng.c)
            u = rng.x * 2.0^(-64)
            nn[i] = u < 0 ? u + 1 : u
        end 
        nn
    end
end

function rand2vecn(rng::LCG, μ::Tuple{Float64,Float64} = (0.0,0.0), σ::Float64 = 1.0)::Tuple{Float64, Float64}
    # Produces a vector of 2 independent N(0,1) samples from an LCG.
    # rng : The linear congruential generator to use.
    # output: a pair (x,y) of iid N(0,1) variables.      
    # ------------------------
    # The producedure is to utilize the Box-Muller transform
    # It boils down to exploiting the following fact: if x, y are iid N(0,1) distributed
    # then the variables exp(-r^2/2) and θ/(2π) are iid U[0,1] distributed,
    # where r = ||(x,y)|| and θ is the angle between (x,y) and the y = 0 axis.

    u1, u2 = rand(rng, 2)
    r = sqrt(-2log(u1))
    θ = 2π*u2

    σ .* (r*cos(θ), r*sin(θ)) .+ μ
end
end