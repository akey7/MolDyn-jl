###########################################################
# MODULES                                                 #
###########################################################

include("lib/mol_dyn_md.jl")
using .MolDyn

###########################################################
# CONSTANTS                                               #
###########################################################

kg_per_amu = 1.661e-27
num_steps = 1000

###########################################################
# ARRAYS HOLDING ATOM INFORMATION                         #
###########################################################

# Positions (qs) and velocities (vs) arrays:
# First axis is timestep
# Second axis are atoms
# Third axis is x,y,z (meters for position, m/s for velocities)

vs = zeros(Float64, num_steps, 2, 3)
qs = zeros(Float64, num_steps, 2, 3)

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
# INITIALIZE SIMULATION                                   #
###########################################################

println(r_ab(qs[1,:,:], qs[2,:,:]))
println(u_stretch(qs[1,:,:], qs[2,:,:], 1.0, r_ab_eq_hcl))
