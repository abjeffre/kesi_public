# Function to generate a squiggly wave
@everywhere function getWeather(t, nperiods, amplitude, frequency_mod, amplitude_mod)
    # Base sine wave with the desired nperiods
    base_wave = amplitude .* sin.(2 .* π .* t ./ nperiods)
    
    # Modifying the base wave with additional sine and cosine waves
    mod_wave = amplitude_mod .* sin.(2 .* π .* frequency_mod .* t ./ nperiods) .+
               amplitude_mod .* cos.(2 .* π .* frequency_mod .* t ./ nperiods .+ π/4) .+
               amplitude_mod .* sin.(2 .* π .* frequency_mod .* t ./ nperiods .+ π/2)
    
    # Combining the base wave and the modulating wave
    squiggly = base_wave .+ mod_wave
    
    
    return squiggly[1:(end-1)]
end

@everywhere function clean_path(file_path) replace(file_path, " " => "_") end

@everywhere function getWeatherSim(
    gamma11 = 0.01, 
    wage_noise = .01
    
)

    # FIXED parameters for Sim
    years = 2000
    nperiods = 26
    M = 2 # Types of weather variables
    weather_noise = 1
    nrounds = years * nperiods
    ngoods = 3

    
    # Generate weather parameters
    amplitude = [0, -1]
    frequency_mod = [-2, 1]
    amplitude_mod = [2.1, -1.5]
    
    # Define wages 
    μ = [.8, .8]  # mean wage
    γ = [0.01 0.02; -0.02 0.036]
    γ[1,1] = gamma11

    # Time vector
    t = 0:1:(years*nperiods)  # Adjust the range and step size as needed

    weather = fill(1.0, nrounds, M)
    for i in 1:M
        weather[:,i] = getWeather(t, nperiods, amplitude[i], frequency_mod[i], amplitude_mod[i]) + rand(Normal(0, weather_noise), nrounds)
    end

    # Make weather data
    wage_data = fill(0.0, nrounds, ngoods - 1)
    for i in 1:(ngoods - 1)
        wage_data[:,i] = μ[i] .+ sum(γ[i,:] .* weather', dims = 1)[1,:] + rand(Normal(0, wage_noise), nrounds)
    end

     plot(plot(wage_data[1:26,:]), plot(weather[1:26,:]))

    # Do Simulation
    observed = cpr_abm(nrounds = nrounds,
        ngoods = ngoods,
        nperiods = nperiods,
        wage_data = wage_data,
        leak = false,
        set_stock = .3, 
        pun1_on = false,
        experiment_group = 1,
        experiment_effort = [.3 ,.3 ,.4],
        price = 1.5,
        labor = 1,
        pun2_on = false,
        wage_elasticities = fill(1, ngoods - 1))

    # Plot outputs
    # wages = plot(observed[:wages][(26*1502):(26*1503), :, 1], alpha = 1, title = "wage income")
    # labor = plot(observed[:effort][:, :, 1, 1], alpha = 1, title = "effort")
    # stock = plot(observed[:stock][:, :, 1], title = "stock")
    # wages2 = plot(wage_data[1:26, :], title = "wage rate")
    # weather2 = plot(weather[1:26, :], title = "weather")
    # plot(labor, wages, wages2, weather2)

    # wages = plot(observed[:wages][(26*1502):(26*1503), :, 1], alpha = 1, title = "wage income")
    # plot!(observed[:harvest][(26*1502):(26*1503), 1, 1] .* 2, alpha = 1, title = "wage income")

    # Save data to CSVs
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/effort_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(observed[:effort][:,:,1,1], :auto))
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/weather_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(weather, :auto))
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(observed[:wages][:,:,1], :auto))

    # Earnings plus harvest
    observed[:wages][:,ngoods,1,1] = observed[:harvest][:,1,1] .* 3
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_full_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(observed[:wages][:,:,1,1], :auto))
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/waves_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame([amplitude frequency_mod amplitude_mod], :auto))
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/gamma_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(γ, :auto))

    ##########################################
    ########## Make counterfactual ###########

    weather[:,1] = weather[:,1] .+ 3

    # Make counterfactual weather data
    wage_data = fill(0.0, nrounds, ngoods - 1)
    for i in 1:(ngoods - 1)
        wage_data[:,i] = μ[i] .+ sum(γ[i,:] .* weather', dims = 1)[1,:] + rand(Normal(0, wage_noise), nrounds)
    end

    counter_factual = cpr_abm(nrounds = nrounds,
        ngoods = ngoods,
        nperiods = nperiods,
        wage_data = wage_data,
        leak = false,
        set_stock = .3,
        experiment_group = 1,
        experiment_effort = [.3 ,.3 ,.4], 
        pun1_on = false,
        price = 1.5,
        labor = 1,
        pun2_on = false,
        wage_elasticities = fill(1, ngoods - 1))

    # Plot outputs
    # wages = plot(counter_factual[:wages][:,:,1], alpha = .2, title = "wage income")
    # labor = plot(counter_factual[:effort][:,:,1,1], alpha = .2, title = "effort")
    # stock = plot(counter_factual[:stock][:,:,1], title = "stock")
    # wages2 = plot(wage_data[1:26, :], title = "wage rate")
    # weather2 = plot(weather[1:26, :], title = "weather")
    # plot(labor, wages, wages2, weather2)

    # Save counterfactual data to CSVs
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/effort_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(counter_factual[:effort][:,:,1,1], :auto))
    counter_factual[:wages][:,ngoods,1,1] = counter_factual[:harvest][:,1,1] .* 2
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_full_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(counter_factual[:wages][:,:,1,1], :auto))
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/weather_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(weather, :auto))
    CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(counter_factual[:wages][:,:,1], :auto))

end


γ₁_sweep = collect(.01:.02:.2)
wage_var_sweep  = collect(.01:.02:.2)
S=expand_grid(γ₁_sweep, wage_var_sweep)
CSV.write(clean_path("data/sweeps/weather_on_earnings/sweep_list.csv"), DataFrame(S, :auto))
# addprocs(10)
pmap(getWeatherSim, S[:,1], S[:,2])


#########################################
########### GRAVE YARD ###################



# function getWeatherSim(
#     years = 2000,
#     nperiods =26,
#     M = 2, # Types of weather variabales. 
#     # M2 = 2, # Types of unobserved periodic forcers
#     ngoods = 3,
#     wage_noise = .01,
#     weather_noise = .01,
#     gamma11 = 0.01
#     )
    
#     nrounds = years*nperiods
#     ##################################################
#     ################ FOR LOOP ########################

#     # Strategy change the size of the effect of Gamma_mk = gamma_11
#     # Change the size of the noise term on the relationship between weather and wages!
#     # Generate the weather
#     #amplitude = round.(rand(Normal(0, 1), M))
#     #frequency_mod = round.(rand(Normal(0, 2), M))  # Number of modulations per nperiods
#     #amplitude_mod = rand(Normal(0, 2), M)  # Amplitude of the modulations
#     amplitude = [0, -1]
#     frequency_mod = [-2, 1]
#     amplitude_mod = [2.1, -1.5]
#     # Generate unobserved forcers
#     # amplitude2 = round.(rand(Normal(0, 1), M2))
#     # frequency_mod2 = round.(rand(Normal(0, 2), M2))  # Number of modulations per nperiods
#     # amplitude_mod2 = rand(Normal(0, 2), M2)  # Amplitude of the modulations

#     # Define wages 
#     μ = [.3, .3] # mean wage
#     #γ = rand(Uniform(-.02, .02),ngoods-1, M) # thing that links weather to wages
#     γ = [0.01 .02; -0.02 .036]
#     # Time vector
#     t = 0:1:(years*nperiods)  # Adjust the range and step size as needed
#     γ[1,1] = gamma11

#     weather = fill(1.0, nrounds, M)
#     for i in 1:M
#         weather[:,i] = getWeather(t, nperiods, amplitude[i], frequency_mod[i], amplitude_mod[i]) + rand(Normal(0, weather_noise), nrounds)
#     end


#     # unobserved = fill(1.0, nrounds, M)
#     # for i in 1:M
#     #     unobserved[:,i] = getWeather(t, nperiods, amplitude2[i], frequency_mod2[i], amplitude_mod2[i]) + rand(Normal(0, w_noise), nrounds)
#     # end
#     # γ2 = rand(Uniform(-.01, .01),ngoods-1, M) # thing that links weather to wages


#     # Make weather data
#     wage_data = fill(0.0, nrounds, ngoods-1)
#     for i in 1:(ngoods -1)
#         wage_data[:,i] = μ[i] .+ sum(γ[i,:].*weather', dims = 1)[1,:] + rand(Normal(0, wage_noise), nrounds)
#         # + sum(γ2[i,:].*unobserved', dims = 1)[1,:]
#     end

#     plot(plot(wage_data[1:26,:]), plot(weather[1:26,:]))


#     # Do Simulation
#     observed =cpr_abm(nrounds = nrounds,
#     ngoods = ngoods,
#     nperiods = nperiods,
#     wage_data = wage_data,
#     leak = false,
#     set_stock  = .07, 
#     pun1_on = false,
#     price = 2.,
#     labor = 1,
#     pun2_on = false,
#     wage_elasticities = fill(1, ngoods-1))

#     26*1500
#     # Plot outputs
#     wages=plot(observed[:wages][(26*1502):(26*1503),:,1], alpha = 1, title = "wage income")
#     labor = plot(observed[:effort][:,:,1,1], alpha = .2, title = "effort")
#     stock = plot(observed[:stock][:,:,1], title = "stock")
#     wages2 = plot(wage_data[1:26,:], title = "wage rate")
#     weather2 = plot(weather[1:26,:], title = "weather" )
#     plot(labor, wages, wages2, weather2)

#     wages=plot(observed[:wages][(26*1502):(26*1503),:,1], alpha = 1, title = "wage income")
#     plot!(observed[:harvest][(26*1502):(26*1503),1,1].*2, alpha = 1, title = "wage income")
#     # Make data 
#     # So we are going to save, weather and wages in CSVs so that we can extract it
#     # in R. 

#     # plot(observed[:effort][collect(20:26:nrounds),2:3,1,1], ylim = (0, 1))
#     using CSV

#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/effort_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(observed[:effort][:,:,1,1], :auto))
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/weather_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(weather, :auto))
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(observed[:wages][:,:,1], :auto))
#     # Earnings plus harvest
#     observed[:wages][:,ngoods,1,1] = observed[:harvest][:,1,1].*2
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_full_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(observed[:wages][:,:,1,1], :auto))
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/waves_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame([amplitude frequency_mod amplitude_mod], :auto))
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/gamma_observed_noise_$wage_noise gamma $gamma11.csv"), DataFrame(γ, :auto))

#     ##########################################
#     ########## Make counterfactual ###########


#     weather[:,1]=weather[:,1].+3

#     # Make weather data
#     wage_data = fill(0.0, nrounds, ngoods-1)
#     for i in 1:(ngoods -1)
#         wage_data[:,i] = μ[i] .+ sum(γ[i,:].*weather', dims = 1)[1,:] + rand(Normal, (0, wage_noise), nrounds)
#         # sum(γ2[i,:].*unobserved', dims = 1)[1,:]
#     end

#     counter_factual=cpr_abm(nrounds = nrounds,
#     ngoods = ngoods,
#     nperiods = nperiods,
#     wage_data = wage_data,
#     leak = false,
#     set_stock  = .07, 
#     pun1_on = false,
#     price = 2.,
#     labor = 1,
#     pun2_on = false,
#     wage_elasticities = fill(1, ngoods-1))

#     # Plot outputs

#     wages=plot(counter_factual[:wages][:,:,1], alpha = .2, title = "wage income")
#     labor = plot(counter_factual[:effort][:,:,1,1], alpha = .2, title = "effort")
#     stock = plot(counter_factual[:stock][:,:,1], title = "stock")
#     wages2 = plot(wage_data[1:26,:], title = "wage rate")
#     weather2 = plot(weather[1:26,:], title = "weather" )
#     plot(labor, wages, wages2, weather2)

#     # Counterfactual
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/effort_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(counter_factual[:effort][:,:,1,1], :auto))
#     # Earnings plus harvest
#     counter_factual[:wages][:,ngoods,1,1] = counter_factual[:harvest][:,1,1].*2
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_full_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(counter_factual[:wages][:,:,1,1], :auto))
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/weather_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(weather, :auto))
#     CSV.write(clean_path("data/sweeps/weather_on_earnings/abm/earnings_counterfactual_noise_$wage_noise gamma $gamma11.csv"), DataFrame(counter_factual[:wages][:,:,1], :auto))

# end

