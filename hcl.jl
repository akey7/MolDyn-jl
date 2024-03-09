###########################################################
# CONSTANTS                                               #
###########################################################

kg_per_amu = 1.661e-27

###########################################################
# ARRAYS HOLDING ATOM INFORMATION                         #
###########################################################

# Positions: Each row is an atom with columns for x, y, z (meters)
qs = zeros(Float64, 2, 3)

# Velocities: Each row is an atom with columns for x, y, z (meters)
vs = zeros(Float64, 2, 3)

# Masses: The masses of each atom (kg)
ms = zeros(Float64, 2, 1)

###########################################################
# INITIALIZE SIMULATION                                   #
###########################################################

# Assume Cl is at 0,0,0 and H lies along the x-axis

# HCl equilibrium bond length
qs[2, :] = [1.57e-10, 0.0, 0.0]

# Masses, Cl first then H
ms[1] = 35 * kg_per_amu
ms[2] = 1 * kg_per_amu

###########################################################
# INITIALIZE SIMULATION                                   #
###########################################################

println(qs)
