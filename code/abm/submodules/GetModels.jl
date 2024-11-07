# # The Function first Samples a number of groups for each individual and stores that to be called later if needed
# # The function then loops through indiviudals. It creates a reference list, other_agents_in_group that has a list of all individuals who are not the focal individual in that persons group.
# # It then checks the Group learning stratgey. 
# # If it is randome then it simply samples either all the agents in a an individuals 'model_group' or thier home group
# # 


# function GetModels(agents, ngroups, gmean, nmodels, out, learn_type, glearn_strat)
#       n = size(agents)[1]
#       models = zeros(n)
#       model_groups = sample(1:ngroups,gmean, n)
#       for i = 1:n
#         other_agents_in_group = agents.id[(agents.gid .∈ Ref(agents.gid[i])) .== (agents.id .!= i)]
#         if glearn_strat == false
#           model_list = ifelse(out[i]==1,
#             sample(agents.id[agents.gid.==model_groups[i]], nmodels, replace = false),
#             sample(other_agents_in_group, nmodels, replace = false))
#           else
#             model_list = ifelse(out[i]==1,
#               sample(agents.id[1:end .!= agents.id[i]], nmodels, replace = false),
#               sample(other_agents_in_group, nmodels, replace = false))
#           end

#         if learn_type == "wealth"
#           candidate_models = model_list[findmax(agents.payoff[model_list])[2]]
#           models[i] = ifelse(agents.payoff[candidate_models] > agents.payoff[i], candidate_models, i) end
#         if learn_type == "income"
#             candidate_models = model_list[findmax(agents.payoff_round[model_list])[2]]
#             models[i] = ifelse(agents.payoff_round[candidate_models] > agents.payoff_round[i], candidate_models, i) end
#       end # End finding models
#       models=convert.(Int64, models)
# end


# MINE IS FASTER

# # #####################################################
# # ################ TESTING OPTIMIZATION ###############


# function GetModels(agents::DataFrame, ngroups::Int, gmean, nmodels::Int, out::Vector{Int}, learn_type::String, glearn_strat::Bool)::Vector{Int}
#     n = size(agents, 1)
#     models = zeros(Int, n)
#     model_groups = rand(1:ngroups, n)  # Sample group for each individual

#     for i in 1:n
#         other_agents_in_group = agents.id[(agents.gid .== agents.gid[i]) .& (agents.id .!= i)]
        
#         model_list = if glearn_strat == false
#             if out[i] == 1
#                 sample(agents.id[agents.gid .== model_groups[i]], nmodels; replace=false)
#             else
#                 sample(other_agents_in_group, nmodels; replace=false)
#             end
#         else
#             if out[i] == 1
#                 sample(agents.id[1:end .!= agents.id[i]], nmodels; replace=false)
#             else
#                 sample(other_agents_in_group, nmodels; replace=false)
#             end
#         end

#         if learn_type == "wealth"
#             candidate_models = model_list[argmax(agents.payoff[model_list])]
#             models[i] = if agents.payoff[candidate_models] > agents.payoff[i]
#                 candidate_models
#             else
#                 i
#             end
#         elseif learn_type == "income"
#             candidate_models = model_list[argmax(agents.payoff_round[model_list])]
#             models[i] = if agents.payoff_round[candidate_models] > agents.payoff_round[i]
#                 candidate_models
#             else
#                 i
#             end
#         end
#     end
#     return models
# end



# function GetModels(agents::DataFrame, ngroups::Int, gmean, nmodels::Int, out::Vector{Int}, learn_type::String, glearn_strat::Bool)::Vector{Int}
#     n = size(agents, 1)
#     models = Vector{Int}(undef, n)
#     model_groups = rand(1:ngroups, n)  # Sample group for each individual

#     # Precompute indices for each group to avoid recomputation in the loop
#     group_indices = Dict{Int, Vector{Int}}()
#     for gid in unique(agents.gid)
#         group_indices[gid] = findall(x -> x == gid, agents.gid)
#     end

#     for i in 1:n
#         gid = agents.gid[i]
#         agent_id = agents.id[i]
#         other_agents_in_group = filter(x -> x != agent_id, agents.id[group_indices[gid]])

#         if out[i] == 1
#             model_list = if glearn_strat
#                 sample(setdiff(agents.id, [agent_id]), nmodels; replace=false)
#             else
#                 sample(agents.id[agents.gid .== model_groups[i]], nmodels; replace=false)
#             end
#         else
#             model_list = sample(other_agents_in_group, nmodels; replace=false)
#         end

#         if learn_type == "wealth"
#             candidate_models = model_list[argmax(agents.payoff[model_list])]
#             models[i] = if agents.payoff[candidate_models] > agents.payoff[i]
#                 candidate_models
#             else
#                 i
#             end
#         elseif learn_type == "income"
#             candidate_models = model_list[argmax(agents.payoff_round[model_list])]
#             models[i] = if agents.payoff_round[candidate_models] > agents.payoff_round[i]
#                 candidate_models
#             else
#                 i
#             end
#         end
#     end
#     return models
# end



function GetModels(agents::DataFrame, ngroups::Int, gmean, nmodels::Int, out::Vector{Int}, learn_type::String, glearn_strat::Bool)::Vector{Int}
  n = size(agents, 1)
  models = Vector{Int}(undef, n)
  model_groups = rand(1:ngroups, n)  # Sample group for each individual

  # Precompute indices and other_agents_in_group for each group
  group_indices = Dict{Int, Vector{Int}}()
  other_agents_in_group = Dict{Int, Vector{Int}}()
  for gid in unique(agents.gid)
      group_indices[gid] = findall(x -> x == gid, agents.gid)
      other_agents_in_group[gid] = filter(x -> x != gid, agents.id[group_indices[gid]])
  end

  for i in 1:n
      gid = agents.gid[i]
      agent_id = agents.id[i]
      group_id = model_groups[i]

      if out[i] == 1
          model_list = if glearn_strat
              sample(setdiff(agents.id, [agent_id]), nmodels; replace=false)
          else
              sample(agents.id[group_indices[group_id]], nmodels; replace=false)
          end
      else
          model_list = sample(other_agents_in_group[gid], nmodels; replace=false)
      end

      if learn_type == "wealth"
          candidate_models = model_list[argmax(agents.payoff[model_list])]
          models[i] = if agents.payoff[candidate_models] > agents.payoff[i]
              candidate_models
          else
              i
          end
      elseif learn_type == "income"
          candidate_models = model_list[argmax(agents.payoff_round[model_list])]
          models[i] = if agents.payoff_round[candidate_models] > agents.payoff_round[i]
              candidate_models
          else
              i
          end
      end
  end
  return models
end
