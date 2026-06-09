###################################################### 
########## Set up Seasonal Variation in Wages ########

nrounds = 1500
reset_stock=fill(true, nrounds)
wage_scalar =  0.1
wave_frequency = .03
amplitude =5
wage1 =  wage_scalar .* ((sin.(wave_frequency .*collect(1:nrounds)) .+ 1)*amplitude/2)
plot(wage1)

a = cpr_abm(tech =.0001,
 leak = false,
 experiment_limit = .2,
 regrow = 0.007,
 experiment_punish2 = .2,
 outgroup = 0.001,
 wage_data = wage1,
 price = .25, 
 nrounds = nrounds,
 max_forest = 20000,
 #reset_stock = reset_stock,
 nsim = 1)

 high=plot(mean(a[:caught2][:,1,:], dims = 2), label = "Y = Reported Illegal Activity", xlab = "Time", c = "#6C0006")
 plot!(wage1, label = "W=Opportunity Cost", c = "black")
 plot!(mean(a[:stock][:,1,:], dims = 2), label = "B=Resource Stock", c = "#006c66")



b = cpr_abm(tech =.0001,
 leak = false,
 experiment_limit = .2,
 regrow = 0.007,
 experiment_punish2 = .2,
 outgroup = 0.001,
 wages = .25,
 price = .25, 
 nrounds = nrounds,
 max_forest = 20000,
 #reset_stock = reset_stock,
 nsim = 1)

low=plot(mean(b[:caught2][:,1,:], dims = 2), label = "", xlab = "Time", c = "#6C0006")
 hline!((.25, .25), label = "", c = "black")
 plot!(mean(b[:stock][:,1,:], dims = 2), label = "", c= "#006c66")


plot(low, high, size = (900, 400), bottom_margin = 20px, grid = false, lw =2)
savefig("kesi\\figures\\simulation_predictions.pdf")