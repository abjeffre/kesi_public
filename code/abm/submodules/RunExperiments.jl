# This function runs the experiments settings the values of particular traits to fixed values
# Note tha the begin_leakage_experiment variable controls the timing of the leakage experiment - its default is to never set it so you should pay attention
# Note that this special experimental leakage group also allows for some groups to be set diffrent from other groups
function RunExperiment(;
 experiment = experiment,
 experiment_group = experiment_group,
 traits = traits,
 effort = effort, agents = agents,
 t= t,
 begin_leakage_experiment = begin_leakage_experiment,
 special_experiment_leak = special_experiment_leak,
 special_leakage_group = special_leakage_group,
 experiment_effort =  experiment_effort,
 experiment_limit =  experiment_limit,
 experiment_punish1 =  experiment_punish1,
 experiment_punish2 =  experiment_punish2,
 experiment_leak =  experiment_leak,
 special_experiment_effort = special_experiment_effort,
 special_effort_group = special_effort_group,
 ngoods = ngoods, 
 )
 if experiment
    for i = 1:length(experiment_group)
        if experiment_leak != 0
            if t >= begin_leakage_experiment 
                    traits.leakage_type[agents.gid .== experiment_group[i]] = rand(Binomial(1, experiment_leak), sum(agents.gid.==experiment_group[i]))
            end
        end
        if experiment_punish1 != 0
            traits.punish_type[agents.gid .==experiment_group[i]] .= experiment_punish1
        end
        if experiment_punish2 != 0
            traits.punish_type2[agents.gid .==experiment_group[i]] .= experiment_punish2
        end
        if experiment_limit !=0
            traits.harv_limit[agents.gid .==experiment_group[i]] .=experiment_limit
        end
        if experiment_effort != 0
            for j in 1:ngoods
                #println(experiment_group)
                #println(j)
                #println(experiment_effort)
                effort[agents.gid.==experiment_group[i], j] .= experiment_effort[j]
            end
        end

        if special_leakage_group != nothing
            for j in 1:length(special_leakage_group)
                # println(special_leakage_group)
                traits.leakage_type[agents.gid .== special_leakage_group[j]] = rbinom(sum(agents.gid.==special_leakage_group[j]),1, special_experiment_leak)
            end
        end
        if special_effort_group != nothing
            for j in 1:length(special_effort_group)
                # println(special_effort_group)
                for j in 1:ngooods
                    effort[agents.gid .== special_leakage_group[j], j] .= special_experiment_effort[j] 
                end
            end
        end
        
    end
 end
 return effort, traits
end