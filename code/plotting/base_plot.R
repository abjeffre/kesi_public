####################################################################
################## MAKE THE BASE PLOT ##############################

gdp_impute<-colMeans(post2$gdp_merged)
seq = seq(50,150, length.out = 100)
sim<- function(seq){
  rnbinom(10000, mu = with(post2, 
                           log1p(exp(a +
                                       mu_gdp*log1p(seq) +
                                       rowMedians(kDEnv[,2,]) +
                                       rowMedians(kDEnv[,1,]) +
                                       rowMedians(kp[,]) +
                                       rowMedians(kX[,]) +
                                       rowMedians(br_k[,]) +
                                       btimber*mean(data$timber_prices)))), post2$omega)
}

plot(data$kesi ~ gdp_impute[1:length(data$kesi)] , xlab = "Average household earnings (Two week period) ", ylab = "Observed Cases",  pch = 16,
     col =col.alpha("#606161", .5), xlim = c(50, 145), cex.lab = font_size, cex = 1.2, cex.axis = 1.2)
ests <- sapply(seq, sim)
rethinking::shade(apply(ests, 2, PI, .68), seq, col = col.alpha("#7fc5bd", .2))
lines(colMeans(ests) ~ (seq), lwd = 4, col  = "#006c66")

lines(c(0, 160), c(colMeans(ests)[1], colMeans(ests)[1]), lwd = 2, col  = "#6C0006", lty = 2)
lines(c(0, 160), c(colMeans(ests)[100], colMeans(ests)[100]), lwd = 2, col  = "#6C0006", lty = 2)
#points(data$kesi ~ (colMeans(post$gdp_merged)), pch = 16,  col = col.alpha("#606161", .1))
text(x = 142, y = c(colMeans(ests)[1])+.4, col  = "#6C0006", cex = 1.3, # Coordinates
     label = paste0("~ ", round(colMeans(ests)[1],2)))
text(x = 142, y = c(colMeans(ests)[100])-.4,  col  = "#6C0006", cex = 1.3, # Coordinates
     label = paste0("~ ", round(colMeans(ests)[100],2)))

text(x = 128, y = (colMeans(ests)[1] +colMeans(ests)[100])/2 ,  col  = "#6C0006", # Coordinates
     label = paste0("~ ", round((1-(colMeans(ests)[100]/colMeans(ests)[1]))*100,0) ,"% change"), cex = 1.5)




Arrows(142,colMeans(ests)[1]-.3,142,colMeans(ests)[100]+.3, col  = "#6C0006",  arr.width=0.25, arr.type="triangle", arr.length =.25)
Arrows(142,colMeans(ests)[100]+.3,142,colMeans(ests)[1]-.3, col  = "#6C0006", arr.width=0.25, arr.type="triangle", arr.length =.25)

#dev.off()
#pdf(file ="figures/catapiller.pdf", width = 7, height = 6)

df <- NULL
for(i in 1:data$K){
  df<-cbind(df, post$bgdp[,i]*post$sigma_gdp + post$mu_gdp)
}


# STEP 1 GET MEANS

set.seed(123)


mu  = rev(colMeans(df))

names(mu) <- rev(str_to_title(sectors))

# Calculate error bars for each column
errors <- apply(df, 2, PI, .9)[, ncol(df):1]

# Create dot chart
dotchart(mu,
         cex = .9,
         xlab = "Posterior estimate",
         pch = 16,
         xlim = c(-3., 3.))

# Add error bars
for (i in 1:ncol(df)) {
  #points(x = errors[,i ],
  #     y = rep(i, 2),
  #     pch = 19,
  #     col = i,
  #     xlim = range(df) + c(-1, 1) * max(errors),
  #     ylim = c(0.5, ncol(df) + 0.5))
  segments(x0 = errors[1,i],
           x1 = errors[2,i],
           y0 = rep(i, 1),
           y1 = rep(i, 1),
           lwd = 2,
           col = "black")
}
abline(v = 0, col = "black", lty = 2)