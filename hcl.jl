###########################################################
# MODULES                                                 #
###########################################################

include("src/mol_dyn_md.jl")
using .MolDyn

###########################################################
# CONSTANTS                                               #
###########################################################

kg_per_amu = 1.661e-27
num_steps = 1000

###########################################################
# ARRAYS HOLDING ATOM INFORMATION                         #
###########################################################

# Positions (qs), velocities (vs), and accelerations (accels) arrays:
# First axis is timestep
# Second axis are atoms
# Third axis is x,y,z (meters for position, m/s for velocities)

vs = zeros(Float64, num_steps, 2, 3)
qs = zeros(Float64, num_steps, 2, 3)
accels = zeros(Float64, num_steps, 2, 3)

# Masses: The masses of each atom (kg)
ms = zeros(Float64, num_steps, 2, 1)

###########################################################
# INITIALIZE SIMULATION                                   #
###########################################################

# Equilibrium bond length for HCl
r_ab_eq_hcl = 1.57e-10

# Assume Cl is at 0,0,0 and H lies along the x-axis

# HCl equilibrium bond length
qs[1, 2, :] = [r_ab_eq_hcl, 0.0, 0.0]

# Masses, Cl first then H
ms[1] = 35 * kg_per_amu
ms[2] = 1 * kg_per_amu

###########################################################
# PROPAGATE ONE TIME STEP                                 #
###########################################################

dt = 1e-15

qs[2,1,:] = qs[1,1,:] + vs[1,1,:].*dt + accels[1,1,:].*dt^2
accels[2,1,:] = -one_bond_stretch_gradient(qs[1,1,:], qs[1,2,:], 1.0, r_ab_eq_hcl) / ms[1]
vs[2,1,:] = vs[1,1,:] + (accels[1,1,:]+accels[2,1,:]).*dt.*0.5

qs[2,2,:] = qs[1,2,:] + vs[1,2,:].*dt + accels[1,2,:].*dt^2
accels[2,2,:] = -one_bond_stretch_gradient(qs[1,2,:], qs[1,1,:], 1.0, r_ab_eq_hcl) / ms[2]
vs[2,2,:] = vs[1,2,:] + (accels[1,2,:]+accels[2,2,:]).*dt.*0.5

###########################################################
# PRINT RESULT                                            #
###########################################################

println("Cl start $(qs[1,1,:]), Cl end $(qs[2,1,:])")
println("H start $(qs[1,2,:]), H end $(qs[2,2,:])")
