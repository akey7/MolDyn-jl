module MolDyn

# Distance from atom a to b
function r_ab(a::Vector{Float64}, b::Vector{Float64})
    sqrt(sum((a-b).^2))
end

# Stretch energy
function stretch(a::Vector{Float64}, b::Vector{Float64}, k_ab::Float64, r_ab_eq::Float64) 
    0.5*k_ab*(r_ab(a,b)-r_ab_eq)^2
end

end
