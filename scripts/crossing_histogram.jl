using DelimitedFiles
using Plots

allDat = readdlm("./data/crossings.dat") # read the output of crossings.jl

encounterTimes = view(allDat, 2:floor(Int, length(allDat)/3), 3) # extract the encounter times

hist = histogram(encounterTimes, labels = "First-encounter time") # plot histogram

png(hist, "./figs/first_encounters.pdf")