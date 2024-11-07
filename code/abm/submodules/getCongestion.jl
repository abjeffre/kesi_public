# Calculate the congestion at a specific group `gid` based on the provided location data,
# effort data, and congestion alpha.
# Note that if no agents harvest in the patch the function returns nan and thus is converted to zero. 

# # Arguments
# - `loc::Vector{Int}`: A vector containing group locations.
# - effort::Matrix{Float64}`: A matrix containing effort data, where each row corresponds
#   to an agents effort and the second column represents the effort values for harvesting the resource.
# - `gid::Int`: The group ID for the agents which congestion is to be calculated.
# - `ngroups::Int`: The total number of groups.
# - `congestion_alpha::Float64`: A scaling factor for congestion calculation.

# # Returns
# - `Float64`: The calculated congestion cost for each agent.

# # Example
# loc = [1, 2, 1, 3, 2, 3]
# effort = [0.8 0.2; 0.6 0.4; 0.9 0.1; 0.7 0.3; 0.5 0.5; 0.8 0.2]
# gid = [1, 1, 2, 2, 3, 3]
# ngroups = 3
# congestion_alpha = 0.5

# congestion_value = getCongestion(loc, effort, gid, ngroups, congestion_alpha)
# println("Congestion for group $gid: $congestion_value")

function getCongestion(loc, effort, gid, ngroups, congestion_alpha)
    cong=(congestion_alpha .*([sum(loc .== i) for i in 1:ngroups].*
    [mean(effort[loc .==i]) for i in 1:ngroups]))[gid]
    return(ifelse.(isnan.(cong), 0, cong))
end


