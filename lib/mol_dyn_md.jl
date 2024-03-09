module MolDyn

# Stretch energy
function stretch(a::Vector{Float64}, b::Vector{Float64}, k_ab::Float64, r_ab_eq::Float64)
    r_ab = sqrt(sum((a-b).^2))
    0.5*k_ab*(r_ab-r_ab_eq)^2
end

# First derivative of stretch energy
function stretch_prime(a::Vector{Float64}, b::Vector{Float64}, k_ab::Float64, r_ab_eq::Float64)
end

end
