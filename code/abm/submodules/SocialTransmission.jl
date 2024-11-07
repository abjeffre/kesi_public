# # function SocialTransmission(trait, models, fidelity, types)
# #         x=copy(trait)
# #         n=size(x)[1]
# #         k=size(x)[2]
# #         trans_error = rand(n,k) .< fidelity
# #         self = models .== collect(1:n)
# #         #trans_error[self,:] .=0
# #         if types == "Dirichlet"
# #                 temp=(Matrix(x).*5).+.01
# #                 n_error  = zeros(n,k)
# #                 for i in 1:n
# #                         n_error[i,:] = rand(Dirichlet(temp[i,:]), 1)'
# #                 end
# #                 x =  ifelse.(trans_error[:,1], n_error[models,:], x[models,:])
# #         else
# #                 for i in 1:k
# #                         if types[i] == "binary" x[:,i]=ifelse.(trans_error[:,i], ifelse.(x[models,i].==1,0,1), x[models,i])  end
# #                         if types[i] == "prob"
# #                                 n_error = rand(Normal(), n)
# #                                 proposal = inv_logit.(logit.(x[models,i]) .+ n_error)
# #                                 proposal = ifelse.(proposal .> .999, .9, proposal)                                
# #                                 new_value = ifelse.(proposal .< 0.001, .1, proposal)
# #                                 x[:,i] =  ifelse.(trans_error[:,i], new_value, x[models,i])
# #                         end
# #                         if types[i] == "positivecont"
# #                                  n_error = rand(Normal(0, .2), n)
# #                                  proposed = x[models,i] .+ n_error
# #                                  x[:,i] =  ifelse.(trans_error[:,i], ifelse.(proposed .<= 0, .1, proposed), x[models,i])
# #                                  #new = rand(Uniform(0, 6), n)
# #                                  #x[:,i] =  ifelse.(trans_error[:,i], n_error, x[models,i])
# #                         end
# #                 end
# #         end
# #         return(x)
# # end




# # function SocialTransmission(trait, models, fidelity, types)
# #     x = copy(Matrix(trait))  # Convert DataFrame to a matrix for computation
# #     n, k = size(x)
# #     trans_error = rand(n, k) .< fidelity
# #     col_names = names(trait)  # Save the column names to restore later

# #     if types == "Dirichlet"
# #         temp = (x .* 5) .+ 0.01
# #         n_error = zeros(n, k)
# #         perturbation = rand(Normal(0, 0.1), n, k)
# #         temp += perturbation
# #         for i in 1:n
# #         n_error[i, :] = softmax(temp[i, :])
# #         end
# #         for i in 1:n
# #             if trans_error[i, 1]
# #                 x[i, :] .= n_error[models[i], :]
# #             else
# #                 x[i, :] .= x[models[i], :]
# #             end
# #         end
# #     else
# #         for i in 1:k
# #             if types[i] == "binary"
# #                 x[:, i] .= ifelse.(trans_error[:, i], 1 .- x[models, i], x[models, i])
# #             elseif types[i] == "prob"
# #                 n_error = rand(Normal(), n)
# #                 proposal = clamp.(inv_logit.(logit.(x[models, i]) .+ n_error), 0.001, 0.999)
# #                 x[:, i] .= ifelse.(trans_error[:, i], proposal, x[models, i])
# #             elseif types[i] == "positivecont"
# #                 n_error = rand(Normal(0, 0.2), n)
# #                 proposed = x[models, i] .+ n_error
# #                 x[:, i] .= ifelse.(trans_error[:, i], ifelse.(proposed .<= 0, 0.1, proposed), x[models, i])
# #             end
# #         end
# #     end

# #     return DataFrame(x, col_names)  # Convert matrix back to DataFrame with original column names
# # end



# using Distributions
# using StatsFuns: softmax
# using DataFrames

# function SocialTransmission(trait::DataFrame, models, fidelity, types)
#     x = copy(Matrix(trait))  # Convert DataFrame to a matrix for computation
#     n, k = size(x)
#     trans_error = rand(n, k) .< fidelity
#     col_names = names(trait)  # Save the column names to restore later

#     if types == "Dirichlet"
#         temp = (x .* 5) .+ 0.01
#         to_modify = findall(trans_error[:, 1])  # Pre-calculate indices to modify
        
#         for i in to_modify
#             perturbation = rand(Normal(0, 0.2), k)
#             temp[models[i], :] += perturbation
#             x[i, :] .= softmax(temp[models[i], :])
#         end
#     else
#         for i in 1:k
#             if types[i] == "binary"
#                 x[:, i] .= ifelse.(trans_error[:, i], 1 .- x[models, i], x[models, i])
#             elseif types[i] == "prob"
#                 n_error = rand(Normal(), n)
#                 proposal = clamp.(inv_logit.(logit.(x[models, i]) .+ n_error), 0.001, 0.999)
#                 x[:, i] .= ifelse.(trans_error[:, i], proposal, x[models, i])
#             elseif types[i] == "positivecont"
#                 n_error = rand(Normal(0, 0.2), n)
#                 proposed = x[models, i] .+ n_error
#                 x[:, i] .= ifelse.(trans_error[:, i], ifelse.(proposed .<= 0, 0.1, proposed), x[models, i])
#             end
#         end
#     end

#     return DataFrame(x, col_names)  # Convert matrix back to DataFrame with original column names
# end



# function SocialTransmission2(trait, models, fidelity, types)
#     x = copy(Matrix(trait))  # Convert DataFrame to a matrix for computation
#     n, k = size(x)
#     trans_error = rand(n, k) .< fidelity

#     if types == "Dirichlet"
#         temp = (x .* 5) .+ 0.01
#         to_modify = findall(trans_error[:, 1])  # Pre-calculate indices to modify
        
#         for i in to_modify
#             perturbation = rand(Normal(0, 0.2), k)
#             temp[models[i], :] += perturbation
#             x[i, :] .= softmax(temp[models[i], :])
#         end
#     end
    
#     return x
# end




# function SocialTransmissionGroup(trait, models, fidelity, types, agents, ngroups, out)
#     x = copy(Matrix(trait))  # Convert to matrix for computation
#     n, k = size(x)
#     trans_error = rand(n, k) .< fidelity

#     if types == "Dirichlet"
#         temp = (x .* 5) .+ 0.01
#         to_modify = findall(trans_error[:, 1])  # Pre-calculate indices to modify
        
#         for i in to_modify
#             perturbation = rand(Normal(0, 0.2), k)
#             temp[models[i], :] += perturbation
#             x[i, :] .= softmax(temp[models[i], :])
#         end
#     else
#         for i in 1:k
#             if types[i] == "binary"
#                 x[:, i] .= ifelse.(trans_error[:, i], 1 .- x[models, i], x[models, i])
#             elseif types[i] == "prob"
#                 n_error = rand(Normal(), n)
#                 proposal = clamp.(inv_logit.(logit.(x[models, i]) .+ n_error), 0.001, 0.999)
#                 x[:, i] .= ifelse.(trans_error[:, i], proposal, x[models, i])
#             elseif types[i] == "positivecont"
#                 n_error = rand(Normal(0, 0.2), n)
#                 proposed = x[models, i] .+ n_error
#                 x[:, i] .= ifelse.(trans_error[:, i], ifelse.(proposed .<= 0, 0.1, proposed), x[models, i])
#             elseif types[i] == "grouppositivecont"
#                 medians = reportMedian(x[:, i], agents.gid, ngroups)
#                 x[:, i][Bool.(out)] .= medians[agents.gid][Bool.(out)]
#                 n_error = rand(Normal(0, 0.2), n)
#                 proposed = x[models, i] .+ n_error
#                 x[:, i] .= ifelse.(trans_error[:, i], ifelse.(proposed .<= 0, 0.1, proposed), x[models, i])
#             end
#         end
#     end

#     return x
# end





# # function SocialTransmission(trait::DataFrame, models, fidelity, types)
# #     x = copy(Matrix(trait))  # Convert DataFrame to a matrix for computation
# #     n, k = size(x)
# #     trans_error = rand(n, k) .< fidelity
# #     col_names = names(trait)  # Save the column names to restore later

# #     if types == "Dirichlet"
# #         temp = (x .* 5) .+ 0.01
# #         n_error = [rand(Dirichlet(temp[i, :])) for i in 1:n]
# #         for i in 1:n
# #             if trans_error[i, 1]
# #                 x[i, :] .= n_error[models[i]]
# #             else
# #                 x[i, :] .= x[models[i], :]
# #             end
# #         end
# #     else
# #         for i in 1:k
# #             if types[i] == "binary"
# #                 x[:, i] .= ifelse.(trans_error[:, i], 1 .- x[models, i], x[models, i])
# #             elseif types[i] == "prob"
# #                 n_error = rand(Normal(), n)
# #                 proposal = clamp.(inv_logit.(logit.(x[models, i]) .+ n_error), 0.001, 0.999)
# #                 x[:, i] .= ifelse.(trans_error[:, i], proposal, x[models, i])
# #             elseif types[i] == "positivecont"
# #                 n_error = rand(Normal(0, 0.2), n)
# #                 proposed = x[models, i] .+ n_error
# #                 x[:, i] .= ifelse.(trans_error[:, i], ifelse.(proposed .<= 0, 0.1, proposed), x[models, i])
# #             end
# #         end
# #     end

# #     return DataFrame(x, col_names)  # Convert matrix back to DataFrame with original column names
# # end


function SocialTransmission(trait, models, fidelity, types)
    x=copy(trait)
    n=size(x)[1]
    k=size(x)[2]
    trans_error = rand(n,k) .< fidelity
    self = models .== collect(1:n)
    #trans_error[self,:] .=0
    if types == "Dirichlet"
            temp=(Matrix(x).*5).+.01
            n_error  = zeros(n,k)
            for i in 1:n
                    n_error[i,:] = rand(Dirichlet(temp[i,:]), 1)'
            end
            x =  ifelse.(trans_error[:,1], n_error[models,:], x[models,:])
    else
            for i in 1:k
                    if types[i] == "binary" x[:,i]=ifelse.(trans_error[:,i], ifelse.(x[models,i].==1,0,1), x[models,i])  end
                    if types[i] == "prob"
                            n_error = rand(Normal(), n)
                            proposal = inv_logit.(logit.(x[models,i]) .+ n_error)
                            proposal = ifelse.(proposal .> .999, .9, proposal)                                
                            new_value = ifelse.(proposal .< 0.001, .1, proposal)
                            x[:,i] =  ifelse.(trans_error[:,i], new_value, x[models,i])
                    end
                    if types[i] == "positivecont"
                             n_error = rand(Normal(0, .2), n)
                             proposed = x[models,i] .+ n_error
                             x[:,i] =  ifelse.(trans_error[:,i], ifelse.(proposed .<= 0, .1, proposed), x[models,i])
                             #new = rand(Uniform(0, 6), n)
                             #x[:,i] =  ifelse.(trans_error[:,i], n_error, x[models,i])
                    end
            end
    end
    return(x)
end
