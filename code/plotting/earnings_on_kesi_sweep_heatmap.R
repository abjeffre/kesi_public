
files<-dir("data/sweeps/earnings_on_kesi/stan")
sweep_list <- as.data.frame(read_csv("data/sweeps/earnings_on_kesi/sweep_list.csv"))

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


M <- rep(0, nrow(sweep_list))

for(i in 1:nrow(sweep_list)){
    # First we get whether or not there has been illegal activity
    base_name <- paste0("price_",sweep_list[i,1], "_inspect_",sweep_list[i,2],".csv")
    # df<-read.csv(paste0("cpr/data/kesi_sweep/", base_name))
    par_recover<- read.csv(paste0("data/sweeps/earnings_on_kesi/stan/good2_imputed_years2_", base_name))[,2]
    par_known<- read.csv(paste0("data/sweeps/earnings_on_kesi/stan/good2_observed_years2_", base_name))[,2]
    
    M[i] = abs(mean(par_recover - par_known)*(sd(par_recover - par_known)/mean(par_known)))
    #M[i] = abs(mean(par_recover - par_known)*(sd(par_known)/mean(par_known)))
}

dat = data.frame(x = sweep_list[,1], y =sweep_list[,2], z = log(M[]))

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
  scale_fill_gradient(low = base_color, high = contrasting_color, name = addline_format(c("Log mean contast")), limits =c(min(dat$z), max(dat$z))) +
  xlab("Resource Benifit (b)") +
  ylab("Inspection Probability (I)") +
  ggtitle("Posterior Constrast") 



######################################
########## Known Parameter ###########


M1 <- rep(0, nrow(sweep_list))

for(i in 1:nrow(sweep_list)){
  par_recover <- c() 
  par_known <- c() 
  for(j in 1:1){
    # First we get whether or not there has been illegal activity
    base_name <- paste0("price_",sweep_list[i,1], "_inspect_",sweep_list[i,2],".csv")
    # df<-read.csv(paste0("cpr/data/kesi_sweep/", base_name))
    par_recover<- c(par_recover, read.csv(paste0("data/sweeps/earnings_on_kesi/stan/good2_imputed_years2_", base_name))[,2])
    par_known<- c(par_known, read.csv(paste0("data/sweeps/earnings_on_kesi/stan/good2_observed_years2_", base_name))[,2])
  }
  M1[i] = mean(par_known)
}

dat = data.frame(x = sweep_list[,1], y =sweep_list[,2], z = M1[])

par_est_heatmap= ggplot(data = dat) +
  geom_raster(aes(x = x, y = y, fill=z)) + 
  theme_classic() + 
  scale_fill_gradient(low = contrasting_color, high = base_color, name = addline_format(c("Parameter Estimate")), limits =c(-100, 7)) +
  xlab("Resource Benifit (b)") +
  ylab("Inspection Probability (I)") +
  ggtitle("Parameter Estimate")


# Combine the plots side by side using cowplot
combined_plot <- plot_grid(par_cont_heatmap, par_est_heatmap, labels = "AUTO")

# Save the combined plot as a PDF
ggsave("figures/combined_plot.pdf", combined_plot, width = 10, height = 4)

# To view the combined plot in the RStudio viewer
print(combined_plot)

