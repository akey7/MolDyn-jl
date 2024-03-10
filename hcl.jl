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
num_steps = 25000

###########################################################
# ARRAYS HOLDING ATOM AND BOND INFORMATION                #
###########################################################

# Positions (xs), velocities (vs), and accelerations (accels) arrays:
# First axis is timestep
# Second axis are atoms
# Third axis is x,y,z (meters for position, m/s for velocities)

vs = zeros(Float64, num_steps, 2, 3)
xs = zeros(Float64, num_steps, 2, 3)
accels = zeros(Float64, num_steps, 2, 3)

# Masses: The masses of each atom (kg)
ms = zeros(Float64, 2)

# Total kinetic energies: The total kinetic energy of the system
# at each timestep.

tkes = zeros(Float64, num_steps)

###########################################################
# INITIALIZE SIMULATION                                   #
###########################################################

# Equilibrium bond length for HCl
r_ab_eq_hcl = 1.57e-10

# Assume Cl is at 0,0,0 and H lies along the x-axis

# HCl equilibrium bond length
xs[1, 2, :] = [r_ab_eq_hcl*0.999, 0.0, 0.0]

# # Masses, Cl first then H
# ms[1] = 35 * kg_per_amu
# ms[2] = 1 * kg_per_amu

# Make each mass the reduced mass of 1H35Cl
ms[1] = (35*1) / (35+1) * kg_per_amu
ms[2] = (35*1) / (35+1) * kg_per_amu

# 1-2 Bonds
# Rows are bonds, columns are atoms participating in bond
# Note: This is specifying edges on a graph, so 1-2 also has 2-1

one_two_bonds = [1 2; 2 1]

# 1-2 Bonds, stretch constants
# Note: There is the same constant for each direction of the bond
# HCl bond constant 516 N/m according to Atkins and de Paula, pg. 454

one_two_bonds_kab = [516.0 516.0]

# 1-2 Bonds, equilibrium distances
# Note: There is a distance for each direction of the bond

one_two_bonds_req = [r_ab_eq_hcl r_ab_eq_hcl]

###########################################################
# VELOCITY VERLET                                         #
###########################################################

dt = 1e-18

stretch_velocity_verlet(xs, vs, accels, tkes, one_two_bonds, one_two_bonds_kab, one_two_bonds_req, ms, dt, num_steps)

###########################################################
# PLOT H X-AXIS TRAJECTORY                                #
###########################################################

display(plot(eachindex(xs[:, 2, 1])/1000, xs[:, 2, 1], xlabel="Time (fs)", ylabel="H x position (m)"))
println("When done looking at the trajectory plot, press enter to continue.")
readline()

display(plot(eachindex(xs[:, 2, 1])/1000, tkes, xlabel="Time (fs)", ylabel="Kinetic Energy (J)"))
println("When done looking the kinetic energies plot, press enter to exit.")
readline()
