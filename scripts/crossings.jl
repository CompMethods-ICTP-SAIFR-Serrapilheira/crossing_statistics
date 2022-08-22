# Final project for the 2022 QBio Computational Methods course
# ---------------------------------
# In this script, we generate encounter statistics between two stationary OU processes

include("SDE.jl")
using .SDE

function collisionCheck(x,y,ϵ)
	# Hard sphere collision test.
	# x, y : Numerical 2-tuples representing spheres on the plane
	# ϵ : Numerical value representing the diameter of either "sphere"
	# output: true if the two spheres intersect

	(x[1]-y[1])^2 + (x[2]-y[2])^2 < ϵ^2
end

function crossingEvent(rng::SDE.LCG, X::StochasticProcess, Y::StochasticProcess, ϵ::Float64, h::Float64, cutOffT::Float64)
	# Calcualates the first-crossing event at a distance ϵ for two stochastic processes X, Y
	# rng : Pseudo-random number generator (here an LCG)
	# X, Y: StochasticProcess objects, with their initial positions already initialized (X.x, Y.x)
	# ϵ : Collision-distance
	# h : Time-step
	# cutOffT : Maximum time to run for.
	# output: (position, time) of the first-crossing event if it happens, otherwise (0, -1).

	maxN = floor(Int, cutOffT/h) # maximum number of timesteps

	for i in 0:maxN
		if collisionCheck(X.x, Y.x, ϵ)
			return (X.x, i * h) # return the event data if a collision has been detected
		end

		# evolve the system forward a time step
		updateSP(rng, X, h)
		updateSP(rng, Y, h)
	end

	((0,0), -1)
end

rng = SDE.LCG(9999); # initilizing the RNG

# We set the parameters for the simulation (identical processes)
λ = (0.0,0.0) 
g = 1.0
τ = 1.0

X = IsoOU2D(λ, λ, g, τ)
Y = IsoOU2D(λ, λ, g, τ)

# Simulation parameters
realization_N = 10000 # number of realizations of the process to calculate

# We estimate time step h so that step size is of order ϵq in the OU process
timeStep(q, ϵ, g, τ) = - τ * log(1 - q^2*ϵ^2/(2g*τ))

q = 1
ϵ = 0.05
h = timeStep(q, ϵ, 1, 1)

maxTime = 100.

# Opening the output file
outName = string("./data/crossings.dat")
outputFile = open(outName, "w")

write(outputFile, string("x y t\n"))

# Simulating N realizations of the stochastic process
for i in 1:realization_N

	# Positions are initialized at the stationary distribution
	X.x = SDE.rand2vecn(rng, λ, g*τ/2) 
	Y.x = SDE.rand2vecn(rng, λ, g*τ/2)

	# Record and write the crossing data
	(x, y), t = crossingEvent(rng, X, Y, ϵ, h, maxTime)

	write(outputFile, string(x, " ", y, " ", t, "\n"))
end

close(outputFile)