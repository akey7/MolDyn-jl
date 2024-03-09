module MolDyn

export r_ab, u_stretch

# Distance from atom a to b
function r_ab(a::Matrix{Float64}, b::Matrix{Float64})
    sqrt(sum((a-b).^2))
end

# Stretch energy
function u_stretch(a::Matrix{Float64}, b::Matrix{Float64}, k_ab::Float64, r_ab_eq::Float64) 
    0.5*k_ab*(r_ab(a,b)-r_ab_eq)^2
end

end
