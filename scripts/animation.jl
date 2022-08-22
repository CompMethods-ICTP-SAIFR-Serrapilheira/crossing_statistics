# Final project for the 2022 QBio Computational Methods course
# ---------------------------------
# In this script, we generate an animation of the time evolution of an OU process

include("SDE.jl")
using Plots
using .SDE

rng = SDE.LCG(10)

x = IsoOU2D((0,0),(0,0),1,1) # initiate an OU process at the origin

trajX, trajY, trajT = trajectorySP(rng, x, 0.01, 10.0) # sample trajectory

anim = @animate for i in 1:length(trajX) # produce a plot for each time step
	plot(trajX[1:i],trajY[1:i], xlims = (-1,1), ylims = (-1,1), labels = ["x(t)"])
end

gif(anim, "./figs/trajectory_fps15.gif", fps = 15) # output .gif animation