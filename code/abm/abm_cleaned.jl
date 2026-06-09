using DataFrames
using Statistics
using Distributions
using Random
using Distributions
using StatsBase
# include("cpr/code/abm/submodules/SplitGroupsTest.jl")
# include("cpr/code/abm/submodules/SocialTransmission.jl")
# include("cpr/code/abm/submodules/ResourceDynamics.jl")
# include("cpr/code/abm/submodules/MutateAgents.jl")
# include("cpr/code/abm/submodules/MakeBabies.jhisl")
# include("cpr/code/abm/submodules/KillAgents.jl")
# include("cpr/code/abm/submodules/GetSeizedPay.jl")
# include("cpr/code/abm/submodules/GetPollution.jl")
# include("cpr/code/abm/submodules/GetPolicy.jl")
# include("cpr/code/abm/submodules/GetPatch.jl")
# include("cpr/code/abm/submodules/GetModels.jl")
# include("cpr/code/abm/submodules/GetInspection.jl")
# include("cpr/code/abm/submodules/GetIndvHarvest.jl")
# include("cpr/code/abm/submodules/GetHarvest.jl")
# include("cpr/code/abm/submodules/GetGroupHarvest.jl")
# include("cpr/code/abm/submodules/GetGroupSeized.jl")
# include("cpr/code/abm/submodules/GetFinesPay.jl")
# include("cpr/code/abm/submodules/GetEcoSysServ.jl")
# include("cpr/code/abm/submodules/TransferWealth.jl")
# include("cpr/code/abm/submodules/GetModelsn1.jl")
# include("cpr/code/abm/submodules/RWLearn.jl")
# include("functions/utility.jl")

function cpr_abm(
  ;nsim = 1,                       # Number of simulations per call
  nrounds = 500,                  # Number of rounds per generation
  n = 300,                        # Size of the population
  ngroups = 2,                    # Number of Groups in the population
  lattice = (1, 2),               # This controls the dimensions of the lattice that the world exists on
  mortality_rate = 0.03,          # The number of deaths per 100 people
  mutation = 0.01,                # Rate of mutation on traits
  wages = .1,                     # Wage rate in other sectors - opportunity costs
  max_forest = 350000,            # Average max stock
  var_forest = 0,                 # Controls athe heterogeneity in forest size across diffrent groups
  degrade = 1.0,                    # This measures how degradable a resource is(when invasion the resource declines linearly with size and as it increase it degrades more quickly, if negative it decreases the rate of degredation), degradable resource means that as the resouce declines in size beyond its max more additional labor is required to harvest the same amount
  regrow = .01,                   # The regrowth rate
  volatility = Normal(0, 0),      # The volatility of the resource each round - set as variance on a normal
  pollution = false,              # Pollution provides a public cost based on 
  pol_slope = .1,                 # As the slope increases the rate at which pollution increases as the resource declines increase
  pol_C = .1,                     # As the constant increases the total amount of polution increases
  ecosys = false,                 # Eco-system services provide a public good to all members of a group proportionalm to size of resource
  eco_slope = 1,                  # As the slope increases the resource will continue to produce ecosystem servies
  eco_C = .01,                    # As the constant increases the total net benifit of the ecosystem services increases
  tech = .00001,                  # Used for scaling Cobb Douglas production function
  labor = .7,                     # The elasticity of labor on harvesting production
  price = 1.0,                      # This sets the price of the resource on the market
  ngoods = 2,                     # Specifiies the number of sectors
  necessity = 0,                  # This sets the minimum amount of the good the household requires
  monitor_tech = 1,           # This controls the efficacy of monitnoring, higher values increase the detection rate -  to understand the functio check plot(curve(pbeta(i, 1, x), 0, 5), where i is the proportion of monitors in a pop
  defensibility = 1,              # This sets the maximum possible insepction rate if all indiviudals participate in guarding it.
  def_perc = true,                # This sets the maximum possible insepction rate if all indiviudals participate in guarding it.
  punish_cost = 0.1,              # This is the cost that must be paid for individuals <0 to monitor their forests - For the default conditions this is about 10 percent of mean payoffs
  fine = 0.0,                     # This controls the size of the fine issued when caught, note that in a real world situation this could be recouped by the injured parties but it is not
  limit_seed = [0, 5],              # This is the average harvest limit. If a person is a punisher it controls the max effort a person is allowed to allocate
  harvest_var = 1.5,              # Harvest limit group offset 
  harvest_var_ind = .5,           # Harvest limit individual offset
  pun2_on = true,                 # Turns punishment on or off. 
  pun1_on = true,                 # Turns group borders on or lff
  seized_on = true,               # Turns seizures on or nff
  fines_evolve = false,            # If false fines stay fixed at initial value
  fines1_on = false,              # Turns on fines for local agents
  fines2_on = false,              # Turns on fines for non locals
  fine_start = 1,                 # Determine mean fine value for all populations at the beginiing SET TO NOHTING TO TURN OFF
  fine_var = .2,                  # Determines the group offset for fines at the begining
  distance_adj =0.9,              # This affects the proboabilty of sampling a more close group.
  travel_cost = .01,              # This basically controls the travel time for individuals who travel to neightboring communities to steal from Note that higher values mean less leakage
  groups_sampled = 1,             # When leakers sample candidate wards to visit they draw this number of groups to compare forest sizes
  social_learning = true,         # Toggels whether Presitge Biased Tranmission is on
  nmodels = 3,                    # The number of models sampled to copy from in social learning
  fidelity = 0.02,                # This is the fidelity of social transmission
  learn_type = "income",          # Two Options - "wealth" and "income" indiviudals can choose to copy wealth or income if copy wealth they copy total overall payoffs, if copy income they copy payoffs from the previous round
  outgroup = 0.01,                # This is the probability that the individual samples from the whole population and not just his group when updating0...
  baseline = .01,                 # Baseline payoff to be added each round -
  leak = true,                    # This controls whether individuals are able to move into neightboring territory to harvest
  verbose = false,                # verbose reporting for debugging
  seed = 1984,
  og_on = false,                  # THe evolution of listening to outgroup members.
  experiment_leak = false,        # THIS SETS THE VALUE OF THE OTHER GROUP LEAKAGE and Punish
  experiment_punish1 = false,     # THIS SETS THE VALUE OF THE OTHER GROUPS PUNISHMENT
  experiment_punish2 = false,     # THIS SETS THE VALUE OF THE OTHER GROUPS PUNISHMENT
  experiment_limit = false,       # THIS SETS THE VALUE OF THE OTHER GROUPS LIMIT
  experiment_effort = false,      # THIS SETS THE VALUE OF THE OTHER GROUPS LIMIT
  experiment_price = false,      # THIS SETS THE VALUE OF THE OTHER GROUPS LIMIT
  experiment_group = 1,           # Determines the set of groups which the experiment will be run on
  experiment_stock = false,
  cmls = false,                   # Determines whether cmls will operate
  invasion = false,                   # Checks invasion criteria by setting all trait start values to near invasion
  glearn_strat = false,           # options: "wealth", "income", "env"
  split_method = "random",        # How groups split is CLMS is on
  kmax_data = nothing,            # This takes data for k
  back_leak = false,              # Determines whether or not individuals can back_invade
  fines_on = false,               # Turns fines on or off!
  inspect_timing = nothing,       # options: nothing, "before", "after", if nothing then it randomizes to half and half
  inher = false,                  # Turns wealth inheretence on or off
  tech_data = nothing,            # The modle can recieve data specifiying the technological capacity of the system over time
  harvest_type = "individual",    # Individuals can pool labor before harvests 
  policy_weight = "equal",        # Typically takes wealth as a weight, but can take any-data that can be used to rank people.
  rec_history =  false,           # You can record the history of wealth but it is costly. 
  resource_zero = false,          # Sets resource to invasion to observe regrowth dynamics
  harvest_zero = false,           # Automatically sets harvest to invasion to observe simple regrowth dynamics
  wealth_degrade = nothing,       # When wealth is passed on it degrades by some percent
  slearn_freq = 1,                # Not opperational - it defines the frequency of social learning 
  reset_stock = nothing,           # Is the year in which resources are reset to max
  socialLearnYear = nothing,      # Which years individuals are allowed to socially learn in  - specify as a vector of dates
  αlearn = 1,                     # Controls learning rate
  indvLearn = false,              # Controls whehter individual learning is turned on
  full_save = false,              # Saves everything if true
  compress_data = true,           # Compresses data to Float64 if true
  control_learning = false,       # Agents can learn from experimental groups if true
  learn_group_policy = false,     # Agents learn the policy of targeted out group not trait of inidivudal
  bsm = "individual",              # Defines how gains from seizures are split options = "Collective" or "individual"
  genetic_evolution = true,        # defines whether or not genetic evolution opperates. 
  special_leakage_group = nothing, # If the leakage experiment group cannot be the same as the other groups. 
  begin_leakage_experiment = 1, 
  population_growth = false,
  pgrowth_data = nothing,
  unitTest = false,
  special_experiment_leak = nothing,
  special_experiment_effort = nothing,
  α = 1,
  population_data = nothing,
  new_leakage_experiment = nothing,
  kseed = nothing,
  limit_seed_override = nothing,
  effort_seed = nothing,
  density_alpha = 0,
  wage_data = nothing,
  congestion = 0, 
  labor_market = false, # This can be false or a scalar between zero and
  set_stock = nothing,
  seizure_pun2_correlation = nothing
) 
  # Make sure all potential parameters are converted into floats for multiple dispatch
  wages=Float64.(wages)
  tech=Float64.(tech)
  labor=Float64.(labor)
  degrade=Float64.(degrade)
  price=Float64.(price)
  
  # Establish group population sizes
  if population_data != nothing
      population_data = convert.(Int64, population_data)
      n = convert(Int64, sum(population_data))
      gs_init = convert.(Int64, population_data)
  else
    gs_init = fill(convert(Int64, n/ngroups), ngroups)
  end
  pgrowth_data == nothing ? pgrowth_data = fill(0, nrounds) : nothing # If no populationg growth data is supplied then simply fill in zeros
  price_copy = copy(price) # Record inital price
  nstart = copy(n) # Record inital n
  history=getHistoryBook(n = n, nrounds = nrounds, ngroups = ngroups,
  nsim = nsim, full_save = full_save, pgrowth_data = pgrowth_data,
  rec_history = rec_history, population_growth = population_growth) # Make history books
    
  ##############################
  ####Start Simulations #######
  for sim = 1:nsim
    seed = seed + sim # Set seed so that it is diffrent on each run for each simulation run
    n = nstart        # Reset n to nstart
    #############################
    #### EXPERIMENT SETUP #######
    temp=[experiment_effort, experiment_limit, experiment_leak, experiment_punish1,
    experiment_punish2, experiment_price]
    experiment, exclude_patches, experiment_learning_groups, learnfromcontrol, normal_groups  = setupExperiments(temp = temp, special_leakage_group = special_leakage_group,
    experiment_group = experiment_group, control_learning = control_learning, ngroups = ngroups)
    #############################
    ##### Create the world ######
    world =  A = reshape(collect(1:ngroups), lattice[1], lattice[2])
    distances = distance(world)
    if kmax_data !== nothing     #make the forests and DIVIDE them amongst the groups
      kmax = kmax_data
      K = copy(kmax)
    else
      kmax = rand(Normal(max_forest, var_forest), ngroups)/ngroups
      kmax[kmax .< 0] .=max_forest/ngroups
      K = copy(kmax)
    end
    resource_zero == true ?  K = rand(ngroups).+1 : nothing

    ################################
    ### Give birth to humanity #####
    if population_data == nothing    # Generate group ID
      tempID=repeat(1:ngroups, gs_init[1])
    else
      tempID=reduce(vcat, [fill(i, population_data[i]) for i in 1:length(population_data)])
    end
    agents = DataFrame( # Create Agents
      id = collect(1:n),
      gid = tempID,
      payoff = zeros(n),
      payoff_round = zeros(n),
      age = sample(1:2, n, replace=true) #notice a small variation in age is necessary for our motrality risk measure
      )
    effort, traits, traitTypes, traitTypesGroup = SeedTraits(n=n, ngroups=ngroups, agents = agents,  ngoods=ngoods,
                                                            gs_init = gs_init, limit_seed = limit_seed, harvest_var = harvest_var,
                                                            harvest_var_ind = harvest_var_ind, seed = seed, fine_start = fine_start,
                                                             learn_group_policy = learn_group_policy, invasion = invasion, fine_var = fine_var,
                                                              limit_seed_override = limit_seed_override, effort_seed = effort_seed)
   
    groups = DataFrame(gid = collect(1:ngroups), # Create Group tracking info
                        fine1 = zeros(ngroups),
                        fine2 = zeros(ngroups),
                        limit = zeros(ngroups),
                        group_status = ones(ngroups),
                        def = zeros(ngroups))                           
    children = Vector{Int}[]  # Set up Array To contain Family history
    for i in 1:n push!(children, Vector{Int}[]) end
    # Build in heterogentiy - fed in through data - 
   # if length(labor) == 1 labor = fill(labor, ngroups) end
    # if length(tech) == 1 tech = fill(tech, ngroups) end
    # if length(degrade) == 1 degrade = fill(degrade, ngroups) end
    #if length(wages) == 1 wages = fill(wages, ngroups) end
   # if length(labor) == ngroups labor= labor[agents.gid] end
    #if length(wages) == ngroups wages= wages[agents.gid] end   
    #set defensibility
    def = zeros(ngroups)
    for i in 1:ngroups
      if def_perc==true
        groups.def[i] = gs_init[i]/defensibility
      else
        groups.def[i] =defensibility
      end
    end
    # Implement K seed
    if kseed !==nothing K = copy(kseed) end
    
    if verbose ==true println(string("Sim: ", sim, ", Initiation: COMPLETED ")) end

    #for testing purposes
    caught2 = zeros(n)

    ############################################
    ############# Begin years ##################
    year = 1
    for t ∈ 1:nrounds
      agents.payoff_round = zeros(n)
      ########Impose Restrictions #######
      if leak == false  traits.leakage_type = zeros(n) end
      if pun1_on == false  traits.punish_type = zeros(n) end
      if pun2_on == false  traits.punish_type2 = zeros(n) end
      if fines1_on == false traits.fines1 = zeros(n)end
      if fines2_on == false traits.fines2 = zeros(n)end
      if og_on == false traits.og_type = zeros(n)end
      if tech_data  !== nothing tech = tech_data[year] end

      ############################
      ### RUN EXPERIMENT #########
       if t > 1 
          if seizure_pun2_correlation !== nothing
            last_round_caught2 =round.(report(caught2, agents.gid, ngroups), digits =2)
            experiment_punish2 =last_round_caught2[1]*seizure_pun2_correlation # Note that this .2 is a scalar 
          end
        end

      effort, traits = RunExperiment(;experiment = experiment, experiment_group = experiment_group, traits = traits,
                                      effort = effort, agents = agents, t= t, begin_leakage_experiment = begin_leakage_experiment,
                                      special_experiment_leak = special_experiment_leak, special_leakage_group = special_leakage_group,
                                      experiment_effort =  experiment_effort,
                                      experiment_limit =  experiment_limit,
                                      experiment_punish1 =  experiment_punish1,
                                      experiment_punish2 =  experiment_punish2,
                                      experiment_leak =  experiment_leak,
                                      experiment_price = experiment_price,
                                      special_experiment_effort = special_experiment_effort
                                      )
      if verbose== true print(string("Experiment: COMPLETED"))end      
      ####################
      #### Politics #######
      groups.limit=GetPolicy(traits.harv_limit, policy_weight, agents.payoff, ngroups, agents.gid, groups.group_status, t)
      if fines_evolve == true
        groups.fine1=GetPolicy(traits.fines1, policy_weight, agents.payoff, ngroups, agents.gid, groups.group_status, t)
        groups.fine2=GetPolicy(traits.fines2, policy_weight, agents.payoff, ngroups, agents.gid, groups.group_status, t)
      else
        groups.fine1 = ones(ngroups).*fine 
        groups.fine2 = ones(ngroups).*fine 
      end     
      ########################
      #### Patch Selection ###
      loc=GetPatch(ngroups, agents.gid, groups.group_status, distances, distance_adj,
                traits.leakage_type, K, groups_sampled, experiment, exclude_patches, back_leak)
      TC=travel_cost.*traits.leakage_type
      # Determine the timing of inspections
      if inspect_timing == nothing
        catch_before = sample([true, false]) # Randomize Whether people this round are caught before or after harvests
      else
        if inspect_timing == "before" catch_before = true end
        if inspect_timing == "after" catch_before = false end
      end
      ###########################
      #### Harvesting ###########
      if catch_before == true
        temp_hg = zeros(n)
        caught1=GetInspection(temp_hg, traits.punish_type, loc, agents.gid, groups.limit, monitor_tech, groups.def, "nonlocal")
        temp_effort = effort[:,2] .* (1 .-caught1)
        if harvest_type == "collective"
          GH=GetGroupHarvest(temp_effort, loc, K, kmax, tech, labor, degrade, ngroups)
          HG=GetIndvHarvest(GH, temp_effort, loc, necessity, ngroups)
        else
          HG=GetHarvest(temp_effort, loc, K, kmax, tech, labor, degrade, necessity, ngroups, agents)
          GH=reportSum(HG, loc, ngroups)
        end
      else
        if harvest_type == "collective"
          GH=GetGroupHarvest(effort[:,2], loc, K, kmax, tech, labor, degrade, ngroups)
          HG=GetIndvHarvest(GH, effort[:,2], loc, necessity, ngroups)
        else
          HG=GetHarvest(effort[:,2], loc, K, kmax, tech, labor, degrade, necessity, ngroups, agents)
          GH=reportSum(HG, loc, ngroups)
        end
      end
      ########################################
      #### Inspection, Seizure and Fines ######
      if catch_before == false  
        caught1=GetInspection(HG, traits.punish_type, loc, agents.gid, groups.limit, monitor_tech, groups.def, "nonlocal") 
      end
      caught2=GetInspection(HG, traits.punish_type2, loc, agents.gid, groups.limit, monitor_tech, groups.def, "local")
      pun1_on ? nothing : caught1 .= 0
      pun2_on ? nothing : caught2 .= 0
      caught_sum = ifelse.((caught1 + caught2) .> 0, 1, 0)
      if catch_before == true  
        seized1=GetGroupSeized(caught1, caught1, loc, ngroups) 
      else
        seized1=GetGroupSeized(HG, caught1, loc, ngroups) #REVERT
      end
      seized2=GetGroupSeized(HG, caught2, loc, ngroups)
      if bsm == "individual"
        SP1=GetSeizedPay(seized1, traits.punish_type, agents.gid, ngroups)
        SP2=GetSeizedPay(seized2, traits.punish_type2, agents.gid, ngroups)
        FP1=GetFinesPay(SP1, groups.fine1, agents.gid, ngroups)
        FP2=GetFinesPay(SP2, groups.fine2, agents.gid, ngroups)
        
      end

      if bsm == "collective"
        group_size = convert(Int64, n/ngroups)
        SP1 = seized1./group_size
        SP2 = seized2./group_size
        FP1 = SP1.*groups.fine1
        FP2 = SP2.*groups.fine2
        SP1 = SP1[agents.gid]
        SP2 = SP2[agents.gid]
        FP1 = FP1[agents.gid]
        FP2 = FP2[agents.gid]
      end
      MC1 = punish_cost*traits.punish_type
      MC2 = punish_cost*traits.punish_type2
      if seized_on == false 
        SP2 = SP1 = zeros(n)
       end
      if fines_on == false FP2 = FP1 = zeros(n) end
      if catch_before == true SP1 .=0 end
      ##################################
      ######## ECOSYSTEM SERVICES ######
      ecosys ? ECO =  GetEcoSysServ(ngroups, eco_slope, eco_C, K, kmax) :  ECO = zeros(n)
      pollution ? POL =  GetPollution(effort[:,2], loc, ngroups, pol_slope, pol_C, K, kmax) : POL = zeros(n)     
      
      ##################################
      ########## CONGESTION ############

      cong=GetCongestion(effort[:,2], loc, ngroups, congestion)
      
      ##################################
      ###### PAYOFFS ###################
      #Wage Labor Market
      if wage_data !== nothing wages = wage_data[t] end
      if labor_market == false 
        WL = GetWages(effort[:,1], wages) #tech # NOTE THIS BROKEN FOR WHEN WE HAVE DIFFRENT WAGES
      else
        WL = GetWages(effort[:,1], wages, n, ngroups, labor_market) #tech
      end

      agents.payoff_round = 
      (HG .*(1 .- caught_sum).*price + # Payoff from harvesting
      SP1.*price + # Payoff from Seizures Access Rights
      SP2.*price + # Payoff from Seizures USe Rights
      FP1.*price + # Payoff from Fines Access Right
      FP2.*price # Payoff from Fines Use Rights
      ).^α + # Scaled by some diminishing marginal returns on harvesting
       WL - # Plus wages
       MC1 - # Minus cost of investment in access rights
       MC2 - # Minus cost of investment in use rights
       TC- # Minus Travel costs
       POL[agents.gid] + # Minus Polution costs
       ECO[agents.gid] - # Plus Ecosystem Services
       cong[loc].*effort[:,2] - # Minus congestion
       ifelse(catch_before == true, (caught1).*groups.fine1[loc], HG .*(caught1).*groups.fine1[loc]) - # Minuse fines if must be paid
       HG .*(caught2).*groups.fine2[agents.gid] # Minus fines from ingroup
      
      # Store Payoffs
      agents.payoff += agents.payoff_round .+ baseline
      agents.payoff_round[isnan.(agents.payoff_round)] .=0
      agents.payoff[isnan.(agents.payoff)] .=0
      if verbose == true
        println("harvest, ", sum(isnan.(HG.*(1 .- caught_sum))))
        println("sg_og, ", sum(isnan.(SP1)))
        println("sg_ig, ", sum(isnan.(SP2)))
        println("fp_ig, ", sum(isnan.(FP1)))
        println("fp_og, ", sum(isnan.(FP1)))
        println("travel_cost:", sum(isnan.(TC)))
        println("wl, ", sum(isnan.(WL)))
        println("mc1,", sum(isnan.(MC1)))
        println("mc2,", sum(isnan.(MC2)))
        println("pol,", sum(isnan.(POL)))
        println("eco,", sum(isnan.(ECO)))
      end
      ################################################
      ########### ECOSYSTEM DYNAMICS #################
      K=ResourceDynamics(GH, K, kmax, regrow, volatility, ngroups, harvest_zero)
      if new_leakage_experiment!==nothing
        K = K .- new_leakage_experiment
        K=ifelse.(K .<=0, 1, K)
      end
      if reset_stock !== nothing
         if reset_stock[t]
            K = kmax 
         end
      end
      
      if set_stock !== nothing
        K = kmax.*set_stock
      end

      if experiment_stock !== false
        K.=experiment_stock
      end
      #################################################
      ############# RECORD HISTORY ####################
      history =RecordHistory(history = history, K = K, kmax = kmax, ngroups = ngroups, effort = effort,
      traits = traits, agents = agents, seized2 = seized2, SP1 = SP1, SP2 = SP2,
      rec_history = rec_history, full_save = full_save, HG = HG, loc = loc, caught2 = caught2,
       n = n, sim=sim, year = year)

      ###############################################
      ##############UPDATE MEMORY ###################
      # mem[:effort] = effort
      # mem[:payoff_round] = agents.payoff_round
      #################################################
      ########### SOCIAL LEARNING #####################


    if social_learning == true
      #Set Up group rankings
      if glearn_strat == "income" gmean=AnalyticWeights(report(agents.payoff_round, agents.gid, ngroups)) end
      if glearn_strat == "wealth" gmean=AnalyticWeights(report(agents.payoff, agents.gid, ngroups)) end
      if glearn_strat == "env" gmean=AnalyticWeights(K./kmax) end
      if glearn_strat == "random" gmean=AnalyticWeights(ones(ngroups))end
      if glearn_strat == false gmean=AnalyticWeights(ones(ngroups))end
      if experiment== true
        if control_learning== false
          if special_leakage_group == nothing
              if length(experiment_group)==1 gmean[experiment_group]=0 end  #this ensures all learning happens not from the experimental group.
              if length(experiment_group)>1 gmean[experiment_group].=0 end
          else
              if length(special_leakage_group)==1 gmean[special_leakage_group]=0 end  #this ensures all learning happens not from the experimental group.
              if length(special_leakage_group)>1 gmean[special_leakage_group].=0 end
          end
        end
      end

      out = zeros(n)
      #Determine if individuals learn from outgroup or ingroup.
      if og_on == false
        out = rbinom(n, 1, outgroup)
      else
        for i in 1:n
          out[i] = rbinom(1,1, traits.og_type[i])[1]
        end
      end
      # Determine whether social learning opperates this year
      sly =
        if socialLearnYear == nothing
            sly = true
        else
            if length(socialLearnYear) > 1
                sly =ifelse(year ∈ socialLearnYear, true, false)
            else
                sly= ifelse(year < socialLearnYear, false, true)
            end
        end

      if sly == true
        if n == ngroups # This for when the number of individuals equals the number of groups. 
          if year > 1 
            #models=GetModelsn1(agents, ngroups, gmean, nmodels, out, learn_type, glearn_strat)
            models=GetModelsn1(agents, history, learn_type, year)
            #println("Before: ", sum(effort[:,2][models] .>= effort[:,2][agents.id]))
            if learn_group_policy
              traits=SocialTransmissionGroup(traits, models, fidelity, traitTypesGroup, agents, ngroups, out)
            else                          
              traits=SocialTransmission(traits, models, fidelity, traitTypes)
            end                    
            # effortT = copy(effort)
            effort2 = zeros(n,2)
            effort2[:,2] = history[:effort][year-1,:,1]
            effort2[:,1] = 1 .-effort2[:,2]
            effort=SocialTransmission2(effort, effort2, models, fidelity, "Dirichlet")
            #   println("Mutation: ", sum(effort[:,2] .> effortT[:,2][models])) 
          end 
        else # This is for normal social learning
          models=GetModels(agents, ngroups, gmean, nmodels, out, learn_type, glearn_strat)
          if learn_group_policy
            traits=SocialTransmissionGroup(traits, models, fidelity, traitTypesGroup, agents, ngroups, out)
          else
            traits=SocialTransmission(traits, models, fidelity, traitTypes)
          end                    
          effort=SocialTransmission(effort, models, fidelity, "Dirichlet")
          if full_save == true history[:models][year,:,sim] .= models end 
        end # End ngroup = n check
      end#end social learning year 
    end # End Social Learning

    
      #########################################
      ############# GENERAL UPDATING ##########
      agents.age .+= 1

      ####################################
      ###### Evolutionary dynamics #######
      # println(agents.payoff)
      # First we check to see if genetic evolution is turned on
      # Second we set all negative or payoffs = to 0.0001 to aid in sampling
      # If genetic evolution is on we check to see if agents are learning from from the experimental groups.
      # First then, if learning from the control group is allowed then we first deal with agents in the experimental group
      # All agents within the experimental group only reporduce amongst themselves.
      # Agents outside of the experimental 
      if genetic_evolution
        agents.payoff_round[agents.payoff_round.<=0] .=0
        agents.payoff[agents.payoff.<=0] .=0
        sample_payoff = ifelse.(agents.payoff .!=0, agents.payoff, 0.0001)
        if learnfromcontrol==false
          # Learning from controls ensure that groups learn from the whole population
          # If learn from controls is false, then groups only update from non-experimental groups
          # If no experiment is being run then learnfromcontrol should equal true  to ensure complete updating is on 
          # Experimental Group
          pop = agents.id[agents.gid .∈  [experiment_learning_groups]]
          died =  KillAgents(pop, agents.id, agents.age, mortality_rate, sample_payoff)
          babies = MakeBabies(pop, agents.id, sample_payoff, died)
          traits[died,:]=MutateAgents(traits[babies, :], mutation, traitTypes)
          effort[died, :]=MutateAgents(effort[babies, :], mutation, "Dirichlet")
          if inher == true
            GetChildList(babies, died, children)
            agents.payoff = DistributWealth(died, agents.payoff, children)
            for i in 1:length(died) children[died[i]] = Vector{Int64}[] end #remove children
          end
          agents.payoff[died] .= 0
          agents.age[died]  .= 0
          #Non-experimental Group
          pop = agents.id[agents.gid .∉  [experiment_learning_groups]]
          died =  KillAgents(pop, agents.id, agents.age, mortality_rate, sample_payoff)
          babies = MakeBabies(pop, agents.id, sample_payoff, died)
          traits[died,:]=MutateAgents(traits[babies, :], mutation, traitTypes)
          effort[died, :]=MutateAgents(effort[babies, :], mutation, "Dirichlet")
          if inher == true
            GetChildList(babies, died, children)
            agents.payoff = DistributWealth(died, agents.payoff, children)
            for i in 1:length(died) children[died[i]] = Vector{Int64}[] end #remove children
          end
          agents.payoff[died] .= 0
          agents.age[died]  .= 0
          if cmls == true  agents.gid[died] = agents.gid[babies] end
        else # If learn from control is TRUE 
          died =  KillAgents(agents.id, agents.id, agents.age, mortality_rate, sample_payoff)
          babies = MakeBabies(agents.id, agents.id, sample_payoff, died)
          if inher == true
            GetChildList(babies, died, children)
            agents.payoff = DistributWealth(died, agents.payoff, children)
            for i in 1:length(died) children[died[i]] = Vector{Int64}[] end #remove children
          end
          traits[died,:]=MutateAgents(traits[babies, :], mutation, traitTypes)
          effort[died, :]=MutateAgents(effort[babies, :], mutation, "Dirichlet")
          agents.payoff[died] .= 0
          agents.age[died]  .= 0
        end
      end

      ###########################################
      ############# POPULATION GROWTH ###########
      # A variable pgrowth_data defines the number of individauls to add to each group at each timestep
      # The population is sampled as above for parents for the total number of new children
      # These children are then mutated and are assigned a random group id, but such that each group gains the same number of new agents
      # N is now updated and gs_init must be updated as well. 
      # There are a set of conditions that must be fufilled.
      # Control learning on 
      if population_growth == true
        if learnfromcontrol == true
            new_agents = convert(Int, pgrowth_data[t]*ngroups)
            #println(new_agents)
            babies = MakeBabies(agents.id, agents.id, sample_payoff, died, new_agents)
            # Add new rows to traits dataframes
            append!(traits, MutateAgents(traits[babies, :], mutation, traitTypes))
            append!(effort, MutateAgents(effort[babies, :], mutation, "Dirichlet"))
            new_ids = collect((n+1):(n+new_agents))
            # Update agents dataframe
            temp = DataFrame(
                id = new_ids,
                gid = repeat(1:ngroups, convert(Int, pgrowth_data[t])),
                payoff = zeros(new_agents),
                payoff_round = zeros(new_agents),
                age = sample(1:2, new_agents, replace=true) #notice a small variation in age is necessary for our motrality risk measure
                )
            append!(agents, temp)  
            if inher == true
                # Ensure wealth updating
                # First grow the child list
                # Add in the new children to the child list
                # Note that no resources need to be distributed to them. 
                new_child_list = Vector{Int}[]
                for i in 1:new_agents push!(new_child_list, Vector{Int}[]) end
                children=[children; new_child_list]
                GetChildList(babies, new_ids, childrentemp) # Update the child list
            end
            #update n
            n = size(agents)[1]    
        else # IF learn from control is false
            # FOR EXPERIMENTAL GROUPS 
            new_agents_experiment = convert(Int, pgrowth_data[t]*length(experiment_group))
            pop = agents.id[agents.gid .∈  [experiment_learning_groups]]
            babies = MakeBabies(pop, agents.id, sample_payoff, died, new_agents_experiment)
            # Add new rows to traits dataframes
            append!(traits, MutateAgents(traits[babies, :], mutation, traitTypes))
            append!(effort, MutateAgents(effort[babies, :], mutation, "Dirichlet"))
            new_ids = collect((n+1):(n+new_agents_experiment))
            # Update agents dataframe
            temp_experiment = DataFrame(
                id = new_ids,
                gid = repeat(collect(experiment_group), convert(Int, pgrowth_data[t]*length(experiment_group))),
                payoff = zeros(new_agents_experiment),
                payoff_round = zeros(new_agents_experiment),
                age = sample(1:2, new_agents_experiment, replace=true) #notice a small variation in age is necessary for our motrality risk measure
                )
              if inher == true
                # Ensure wealth updating
                # First grow the child list
                # Add in the new children to the child list
                # Note that no resources need to be distributed to them. 
                new_child_list = Vector{Int}[]
                for i in 1:new_agents_experiment push!(new_child_list, Vector{Int}[]) end
                children=[children; new_child_list]
                GetChildList(babies, new_ids, childrentemp) # Update the child list
            end
            #update n
            # For non-Experimental groups do the exact same.
            new_agents_nonexp = convert(Int,pgrowth_data[t]*(ngroups-length(experiment_group)))
            pop = agents.id[agents.gid .∉  [experiment_learning_groups]]
            babies = MakeBabies(pop, agents.id, sample_payoff, died, new_agents_nonexp)
            # Add new rows to traits dataframes
            append!(traits, MutateAgents(traits[babies, :], mutation, traitTypes))
            append!(effort, MutateAgents(effort[babies, :], mutation, "Dirichlet"))
            new_ids = collect((n+1+new_agents_experiment):(n+new_agents_experiment+new_agents_nonexp))
            # Update agents dataframe
            temp_nonexperiment = DataFrame(
                id = new_ids,
                gid = repeat(normal_groups, convert(Int, pgrowth_data[t]*length(normal_groups))),
                payoff = zeros(new_agents_nonexp),
                payoff_round = zeros(new_agents_nonexp),
                age = sample(1:2, new_agents_nonexp, replace=true) #notice a small variation in age is necessary for our motrality risk measure
                )
            if inher == true
                # Ensure wealth updating
                # First grow the child list
                # Add in the new children to the child list
                # Note that no resources need to be distributed to them. 
                new_child_list = Vector{Int}[]
                for i in 1:new_agents_nonexp push!(new_child_list, Vector{Int}[]) end
                children=[children; new_child_list]
                GetChildList(babies, new_ids, childrentemp) # Update the child list
            end # Inheretence
            #update n
            append!(agents, temp_experiment)  
            append!(agents, temp_nonexperiment)   
            n = size(agents)[1]
        end # End ControlLearning Group    
        if unitTest == true 
            println("All IDs Unique:", length(unique(agents.id)) == n)
            println("All Groups are appropiate size:", all(values(countmap(agents.gid)).== n/ngroups))
            println("n =", n)
        end
    end # End Population Growth
        

      ############################################
      ############ WEALTH DYNAMICS ###############
  
      wealth_degrade !== nothing ? agents.payoff = agents.payoff .* wealth_degrade : nothing

      #
      year += 1
    end #End the year
  end #End Sims
  if compress_data == true
    ky=keys(history)
    for key in ky
      history[key]=convert.(Float16, history[key]) 
    end#loop
  end#Compress data
  return(history)
end#End Function
