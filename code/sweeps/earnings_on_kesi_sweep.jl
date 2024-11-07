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

@everywhere function getEarningsonKesi(
    pr = 0.01, 
    ins = .01
    
)

    # FIXED parameters for Sim
    years = 2000
    nperiods = 26
    M = 2 # Types of weather variables
    weather_noise = .01
    nrounds = years * nperiods
    ngoods = 3

    # Generate weather parameters
    amplitude = [0, -1]
    frequency_mod = [-2, 1]
    amplitude_mod = [2.1, -1.5]
    
    # Define wages 
    μ = [.8, .8]  # mean wage
    γ = [0.01 0.02; -0.02 0.036]

    # Time vector
    t = 0:1:(years*nperiods)  # Adjust the range and step size as needed

    weather = fill(1.0, nrounds, M)
    for i in 1:M
        weather[:,i] = getWeather(t, nperiods, amplitude[i], frequency_mod[i], amplitude_mod[i]) + rand(Normal(0, weather_noise), nrounds)
    end

    # Make weather data
    wage_data = fill(0.0, nrounds, ngoods - 1)
    for i in 1:(ngoods - 1)
        wage_data[:,i] = μ[i] .+ sum(γ[i,:] .* weather', dims = 1)[1,:] + rand(Normal(0, .03), nrounds)
    end

    # plot(plot(wage_data[1:26,:]), plot(weather[1:26,:]))

    # Do Simulation
    observed = cpr_abm(nrounds = nrounds,
        ngoods = ngoods,
        nperiods = nperiods,
        wage_data = wage_data,
        experiment_limit = .3,
        experiment_punish2 = ins,
        leak = false,
        set_stock = .25,
        pun1_on = false,
        seized_on = false,
        price = pr, 
        labor = 1,
        pun2_on = true,
        wage_elasticities = fill(1, ngoods - 1))
    # Plot outputs
    # wages = plot(observed[:wages][(26*1502):(26*1503), :,1, 1], alpha = 1, title = "wage income")
    # labor = plot(observed[:effort][:, :, 1, 1], alpha = .2, title = "effort")
    # stock = plot(observed[:stock][:, :, 1], title = "stock")
    # wages2 = plot(wage_data[1:26, :], title = "wage rate")
    # weather2 = plot(weather[1:26, :], title = "weather")
    # plot(labor, wages, wages2, weather2)

    # wages = plot(observed[:wages][(26*1502):(26*1503), :, 1], alpha = 1, title = "wage income")
    # plot(observed[:harvest][(26*1502):(26*1503), 1, 1] .* 2, alpha = 1, title = "wage income")
    # Save data to CSVs
    CSV.write(clean_path("data/sweeps/earnings_on_kesi/abm/effort_price_$pr inspect $ins.csv"), DataFrame(observed[:effort][:,:,1,1], :auto))
    CSV.write(clean_path("data/sweeps/earnings_on_kesi/abm/weather_price_$pr inspect $ins.csv"), DataFrame(weather, :auto))
    CSV.write(clean_path("data/sweeps/earnings_on_kesi/abm/earnings_price_$pr inspect $ins.csv"), DataFrame(observed[:wages][:,:,1], :auto))
    observed[:punish2][:,:,1,1]

    # Earnings plus harvest
    observed[:wages][:,ngoods,1,1] = observed[:harvest][:,1,1] .* pr
    CSV.write(clean_path("data/sweeps/earnings_on_kesi/abm/earnings_full_price_$pr inspect $ins.csv"), DataFrame(observed[:wages][:,:,1,1], :auto))
    CSV.write(clean_path("data/sweeps/earnings_on_kesi/abm/waves_price_$pr inspect $ins.csv"), DataFrame([amplitude frequency_mod amplitude_mod], :auto))
    CSV.write(clean_path("data/sweeps/earnings_on_kesi/abm/gamma_price_$pr inspect $ins.csv"), DataFrame(γ, :auto))
    CSV.write(clean_path("data/sweeps/earnings_on_kesi/abm/kesi_price_$pr inspect $ins.csv"), DataFrame(ins=Int64.(round.(observed[:caught2][:,1,1].*150))))


end


inspect = collect(.01 : .1 : 1)
price = collect(1.1:.12:2.2)

S=expand_grid(price, inspect)
CSV.write(clean_path("data/sweeps/earnings_on_kesi/sweep_list.csv"), DataFrame(S, :auto))
# addprocs(10)
pmap(getEarningsonKesi, S[:,1], S[:,2])

