
files<-dir("data/sweeps/weather_on_earnings/stan")
sweep_list <- as.data.frame(read_csv("data/sweeps/weather_on_earnings/sweep_list.csv"))

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


M <- rep(0, nrow(sweep_list))

for(i in 1:nrow(sweep_list)){
  print(i)
  # First we get whether or not there has been illegal activity
  base_name <- paste0("noise_",sweep_list[i,1], "_gamma_",sweep_list[i,2],".csv")
  # df<-read.csv(paste0("cpr/data/kesi_sweep/", base_name))
  par_recover<- read.csv(paste0("data/sweeps/weather_on_earnings/stan/y_counterfactual_1_observed_years_10_", base_name))[,2]
  par_known<- colMeans(read.csv(paste0("data/sweeps/weather_on_earnings/stan/y_estimated_1_observed_years_10_", base_name)))[2:261]
  
  M[i] = mean((par_recover - par_known)/par_known)
  #M[i] = abs(mean(par_recover - par_known)*(sd(par_known)/mean(par_known)))
}

dat = data.frame(x = sweep_list[,1], y =sweep_list[,2], z = M[])

# Define the specific color and a contrasting color
base_color <- "#006c66"
contrasting_color <- "#DD1970"

# Create a color gradient function
color_gradient <- colorRampPalette(c(base_color, contrasting_color))

# Generate a sequence of colors in the gradient
gradient_colors <- color_gradient(100)  # Generate 100 colors in the gradient

par_cont_heatmap = ggplot(data = dat) +
  geom_raster(aes(x = x, y = y, fill=z)) +
  theme_classic() + 
  scale_fill_gradient(low = base_color, high = contrasting_color, name = addline_format(c("Predicted Error (%)")), limits =c(min(dat$z), max(dat$z))) +
  xlab("Noise (Episilon)") +
  ylab("Weather on Earnings (Gamma)") +
  ggtitle("Error in Prediction") 


ggsave("figures/prediction_error.pdf", par_cont_heatmap, width = 5, height = 4)
