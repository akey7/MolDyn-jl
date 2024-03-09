###########################################################
# MODULES                                                 #
###########################################################

using Plots

include("src/mol_dyn_md.jl")
using .MolDyn

###########################################################
# CONSTANTS                                               #
###########################################################

kg_per_amu = 1.661e-27
num_steps = 1000

###########################################################
# ARRAYS HOLDING ATOM AND BOND INFORMATION                #
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
qs[1, 2, :] = [r_ab_eq_hcl*0.9, 0.0, 0.0]

# Masses, Cl first then H
ms[1] = 35 * kg_per_amu
ms[2] = 1 * kg_per_amu

# Reduced masses are identical
ms[1] = ((35*1)/(35+1)) * kg_per_amu
ms[2] = ((35*1)/(35+1)) * kg_per_amu

# 1-2 Bonds
# Rows are bonds, columns are atoms participating in bond
# Note: This is specifying edges on a graph, so 1-2 also has 2-1

one_two_bonds = [1 2; 2 1]

# 1-2 Bonds, stretch constants
# Note: There is the same constant for each direction of the bond
# HCl bond constant 516 N/m according to Atkins and de Paula, pg. 454

one_two_bonds_kab = [1.0 1.0]

# 1-2 Bonds, equilibrium distances
# Note: There is a distance for each direction of the bond

one_two_bonds_req = [r_ab_eq_hcl r_ab_eq_hcl]

###########################################################
# VELOCITY VERLET                                         #
###########################################################

dt = 1e-15

stretch_velocity_verlet(qs, vs, accels, one_two_bonds, one_two_bonds_kab, one_two_bonds_req, ms, dt, num_steps)

###########################################################
# PLOT H X-AXIS TRAJECTORY                                #
###########################################################

display(plot(eachindex(qs[:, 2, 1]), qs[:, 2, 1]))
println("When done looking at the plot, press enter to exit.")
readline()
