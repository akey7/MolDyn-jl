###########################################################
# MODULES                                                 #
###########################################################

include("src/mol_dyn_md.jl")
using .MolDyn

###########################################################
# CONSTANTS                                               #
###########################################################

kg_per_amu = 1.661e-27
num_steps = 100

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

# 1-2 Bonds
# Rows are bonds, columns are atoms participating in bond
# Note: This is specifying edges on a graph, so 1-2 also has 2-1
one_two_bonds = [1 2; 2 1]

# 1-2 Bonds, stretch constants
# Note: There is a constant for each direction of the bond
one_two_bonds_kab = [1.0 1.0]

# 1-2 Bonds, equilibrium distances
# Note: There is a distance for each direction of the bond
one_two_bonds_req = [1.57e-10 1.57e-10]

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
# VELOCITY VERLET                                         #
###########################################################

dt = 1e-15

for time_i in 2:num_steps
    for bond_i in [1 2]
        k_ab = one_two_bonds_kab[bond_i]
        r_eq = one_two_bonds_req[bond_i]
        a_i = one_two_bonds[bond_i, 1]
        b_i = one_two_bonds[bond_i, 2]
        qs[time_i, a_i, :] = qs[time_i-1, a_i, :] + vs[time_i-1, a_i, :] .* dt + accels[time_i-1, a_i, :].*dt^2
        accels[time_i, a_i, :] = -one_bond_stretch_gradient(qs[time_i-1, a_i, :], qs[time_i-1, b_i, :], k_ab, r_eq) / ms[a_i]
        vs[time_i, a_i, :] = vs[time_i-1, a_i, :] + (accels[time_i-1, a_i, :] + accels[time_i, a_i, :]) .* dt .* 0.5
    end
end

###########################################################
# PRINT RESULT                                            #
###########################################################

println("Cl start $(qs[1,1,:]), Cl end $(qs[100,1,:])")
println("H start $(qs[1,2,:]), H end $(qs[100,2,:])")
