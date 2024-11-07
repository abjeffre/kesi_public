#################################
######### UPDATING ##############

# Case where each agent has their own wage rate. 
function GetWages(effort::DataFrame, wages::Vector{Vector{Float64}}, agents::DataFrame)
    n, L = size(effort)
    wage_results = zeros(n, L)
    for j in 1:L
        wage_results[:, j] = wages[j][agents.gid] .* effort[:, j].^wage_elasticities[j]
    end
    return wage_results
end

# The case where each good has its own wage rate. 
function GetWages(effort::DataFrame, wages::Vector{Float64}, wage_elasticities)
    L = size(effort, 2)
    output = zeros(size(effort))
    for j in 1:L
        output[:, j] = wages[j] .* effort[:, j].^wage_elasticities[j]
    end
    return output
end

# This is for labor market dynamics

function GetWages(effort::Matrix{Float64}, wages::Vector{Float64}, gid::Vector{Int}, labor_market::Vector{Float64}, ngroups::Int)
    L = size(effort, 2)
    output = zeros(size(effort))
    
    for j in 1:L
        group_counts = [sum(gid .== i) for i in 1:ngroups]
        group_means = [mean(effort[gid .== i, j]) for i in 1:ngroups]
        adjusted_wages = wages[j] .- (labor_market .* (group_counts .* group_means))[gid]
        output[:, j] = adjusted_wages .* effort[:, j].^wage_elasticities[j]
    end
    
    return output
end


# This is for environmental damage

function GetWages(effort::Matrix{Float64}, wages::Vector{Float64}, gid::Vector{Int}, ngroups::Int, K::Vector{Float64}, kmax::Vector{Float64}, eco_sys_dam::Vector{Float64})
    L = size(effort, 2)
    output = zeros(size(effort))
    
    for j in 1:L
        # Compute the ecosystem damage adjustment for each group
        eco_adjustment = [eco_sys_dam[j] * (1 - K[i] / kmax[i]) for i in 1:ngroups]
        
        # Compute the adjusted wages for each group
        adjusted_wages = (wages[j] .- eco_adjustment)[gid]
        
        # Compute the output for the j-th column
        output[:, j] = adjusted_wages .* effort[:, j]
    end
    
    return output
end

function GetWages(effort, wages::Float64, wage_elasticities) wages.*effort.^wage_elasticities end