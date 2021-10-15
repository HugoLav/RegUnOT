# Generate the figures

include("computations_prox.jl")
include("dynamical_RUOT.jl")

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------


# Uncomment below depending if you want OT, ROT, UOT, RUOT

# For OT
# nu = 0.0
# isGrowth = false
# nameFile = "store_OT.txt"

# For ROT
# nu = 0.015
# isGrowth = false
# nameFile = "store_ROT.txt"

# For UOT
# nu = 0.0
# isGrowth = true
# nameFile = "store_UOT.txt"

# For RUOT
nu = 0.015
isGrowth = true
nameFile = "store_RUOT.txt"



# Parameters for the growth penalization
growthInvIntensity = 1e2

# Size of the grids
nTime = 24
nSpace = 100


# Parameters Douglas Rachford
nIter = 10000
alpha = 1.0
gamma = 1e-1


# ------------------------------------------------------------
# Auxilliary functions for the boundary conditions and Psi
# ------------------------------------------------------------

function gaussian(x,m,sigma)
    return exp( - (x-m)^2 / (2 * sigma^2) )
end

# Definition of the function Psi*, together with its derivatives

function psiStar(s)
    return growthInvIntensity * (cosh(s) - 1.)
end

function dPsiStar(s)
    return growthInvIntensity*sinh(s)
end

function d2PsiStar(s)
    return growthInvIntensity*cosh(s)
end

function proxPsiStar(s,gamma)
    return proxCosh(s, gamma*growthInvIntensity)
end


# ------------------------------------------------------------
# Define the boundary conditions
# ------------------------------------------------------------

gridSpace = LinRange(0,1-1/nSpace,nSpace)
gridSpaceAvg = LinRange(1/(2*nSpace),1-1/(2*nSpace),nSpace)

# Gaussian with different weights
rho0 = 2*gaussian.( gridSpace, 0.2, 0.03 ) + gaussian.( gridSpace, 0.8, 0.03 )
rho1 = gaussian.( gridSpace, 0.4, 0.03 ) + 2*gaussian.( gridSpace, 0.6, 0.03 )

# Normalize
# Note that a priori it is not necessary for UOT and RUOT
rho0 /= sum(rho0)
rho1 /= sum(rho1)

# ------------------------------------------------------------
# Compute the solution
# ------------------------------------------------------------


println("start")

rho, momentum, growth, simP = dynamical_RUOT(nTime,nSpace,rho0,rho1;
    nu = nu,
    isGrowth = isGrowth,
    psiStar=psiStar, dPsiStar = dPsiStar, d2PsiStar = d2PsiStar, proxPsiStar=proxPsiStar,
    nIter=nIter,
    alpha = alpha,
    gamma = gamma,
    verbose = true)

# Plot the error made

println("Final error update Z")
println(simP["errorZ"])
println("Final error update W")
println(simP["errorW"])
println("Discrepancy continuity equation")
println(simP["errorCE"])

# ------------------------------------------------------------
# Do some plotting
# ------------------------------------------------------------

# Plot the trajectories at some instant in time

p = plot(rho0, title =  "rho0" )
display(p)

p = plot(rho1, title =  "rho1" )
display(p)


# for i = 1:nTime+1
for i = 1:nTime
    p = plot(rho[i,:], title =  string("time: ",(i-1)/nTime) )
    display(p)
end

# ------------------------------------------------------------
# Store the result in a file which is readible by pgfplots
# ------------------------------------------------------------

# Write the initial and final conditions on a file

fileObj = open("temporal_boundary.txt", "w")

write(fileObj, "grid rho0 rho1")
write(fileObj, "\n")

# Then write the values
for j = 1:nSpace
    write(fileObj, string( gridSpace[j], " "))
    write(fileObj, string(rho0[j]*nSpace, " "))
    write(fileObj, string(rho1[j]*nSpace, " "))
    write(fileObj, "\n")
end

close(fileObj)

# Then write the interpolation

fileObj = open(nameFile, "w")

# First line: time

write(fileObj, "grid " )
for i = 1:nTime
    write(fileObj, string("rho", i, " " ))
end
write(fileObj, "\n")

# Other lines: value of the function in space
for j = 1:nSpace
    write(fileObj, string( gridSpaceAvg[j], " "))
    for i = 1:nTime
        write(fileObj, string( rho[i,j]*nSpace, " "))
    end
    write(fileObj, "\n")
end

close(fileObj)
