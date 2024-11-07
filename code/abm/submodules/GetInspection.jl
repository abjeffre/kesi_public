#  function GetInspection(harvest, x, loc, gid, policy, monitor_tech, def, type)
#   caught = asInt.(zeros(length(x)))
#   for i = 1:length(x)
#     if type == "nonlocal"
#       prob_caught = (sum(x[gid .== loc[i]])/def[gid[i]]).^monitor_tech
#       prob_caught =prob_caught > 1 ? 1 : prob_caught
#       if loc[i] != gid[i]
#                 caught[i] = rbinom(1, 1, prob_caught)[1]
#       end
#     end
#     if type == "local"
#        prob_caught = (sum(x[gid .== loc[i]])/def[gid[i]]).^monitor_tech
#        prob_caught = prob_caught > 1 ? 1 : prob_caught
#        if loc[i] == gid[i]
#          if harvest[i] >= policy[gid[i]] caught[i] = rbinom(1, 1, prob_caught)[1]end
#        end
#      end
#     end
#     return(caught)
#   end

  ######################################################
  ################## TESTING OPTIMIZATION ##############


# function GetInspection(harvest, x, loc, gid, policy, monitor_tech, def, type)
#   n = length(x)
#   caught = zeros(Int, n)
#   prob_caught = zeros(Float64, n)
  
#   # Precompute probabilities for each location and group
#   for i in 1:n
#       prob_caught[i] = (sum(x[gid .== loc[i]]) / def[gid[i]]) ^ monitor_tech
#   end
#   prob_caught .= min.(prob_caught, 1.0)  # Clamp values to 1

#   if type == "nonlocal"
#       for i in 1:n
#           if loc[i] != gid[i]
#               caught[i] = rand() < prob_caught[i] ? 1 : 0
#           end
#       end
#   elseif type == "local"
#       for i in 1:n
#           if loc[i] == gid[i] && harvest[i] >= policy[gid[i]]
#               caught[i] = rand() < prob_caught[i] ? 1 : 0
#           end
#       end
#   end

#   return caught
# end


using DataFrames
using Random

function precompute_prob_caught(x, gid, loc, def, monitor_tech)
    groups = unique(gid)
    prob_caught_dict = Dict{Int, Float64}()
    
    for group in groups
        prob_caught_dict[group] = (sum(x[gid .== group]) / def[group]) ^ monitor_tech
    end
    
    # Clamp probabilities to a maximum of 1
    for group in keys(prob_caught_dict)
        prob_caught_dict[group] = min(prob_caught_dict[group], 1.0)
    end
    
    return prob_caught_dict
end

function GetInspection(harvest, x, loc, gid, policy, monitor_tech, def, type)
    n = length(x)
    caught = zeros(Int, n)

    prob_caught_dict = precompute_prob_caught(x, gid, loc, def, monitor_tech)
    
    if type == "nonlocal"
        for i in 1:n
            if loc[i] != gid[i]
                prob_caught = prob_caught_dict[gid[i]]
                caught[i] = rand() < prob_caught ? 1 : 0
            end
        end
    elseif type == "local"
        for i in 1:n
            if loc[i] == gid[i] && harvest[i] >= policy[gid[i]]
                prob_caught = prob_caught_dict[gid[i]]
                caught[i] = rand() < prob_caught ? 1 : 0
            end
        end
    end
    
    return caught
end



