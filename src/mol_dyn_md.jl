module MolDyn

# Export everything: Minimially it will all need testing.
export r_ab, u_stretch, one_bond_stretch_gradient

# Distance from atom a to b
function r_ab(a, b)
    sqrt(sum((a-b).^2))
end

# Stretch energy
function u_stretch(a::Matrix{Float64}, b::Matrix{Float64}, k_ab::Float64, r_ab_eq::Float64) 
    0.5*k_ab*(r_ab(a,b)-r_ab_eq)^2
end

function one_bond_stretch_gradient(a, b, k_ab, r_ab_eq)
    du_drab = k_ab*(r_ab(a, b)-r_ab_eq)
    drab_dxa = (a[1]-b[1])/r_ab(a, b)
    drab_dya = (a[2]-b[2])/r_ab(a, b)
    drab_dza = (a[3]-b[3])/r_ab(a, b)

    [drab_dxa, drab_dya, drab_dza] .* du_drab
end

end
