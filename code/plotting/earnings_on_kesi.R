#####################################################################################
################ Marginal Impact of Change in Sectoral Income on Illegal Harvests ###
#mm <- readRDS("data/sectors_all_moisture.RDS")
#post <- extract.samples(mm)


nsim  = 5000
npred = 10
K = ncol(data$y)
lambda = array(NA, dim=c(nsim, npred, data$K))
for (k in 1 : K){
  seq <- seq(min(data$y[,k]), max(data$y[,k]), length.out = 10)
  for (i in 1:length(seq)){
    y <- data$y
    y[,k]<-seq[i]
    mean_y=colMeans(y)
    m = matrix(NA, nrow = nrow(post$bgdp), ncol = data$K)
    for(j in 1:data$K){
      m[,j] = (post$bgdp[,j]*post$sigma_gdp + post$mu_gdp) * log1p(mean_y[j])
    }
    
    for (t in 1 : 1:26){
      m2 <- matrix(NA, nrow = length(post$a), ncol = data$DM)
      for(j in 1:data$DM){
        m2[,j] <- post$kDEnv[,j,data$denv_ind[j,t]]
      }
    }      
    lambda[,i,k]=with(post,
                      rnbinom(nsim, mu=log1p(exp( a +
                                                    rowSums(m) + 
                                                    rowMeans(kp) +
                                                    rowMeans(kX) +
                                                    btimber *mean(data$timber_prices) +
                                                    rowMeans(br_k) +
                                                    rowSums(m2))), omega)
    )
    
  }
} 


# pdf("figures/earnings_on_kesi.pdf", width = 7, height = 6)


i =1
#par(mfrow=c(3, 4), oma =c(3,3,3,0), mar = c(2,1,1,1))
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), ylim = c(0, 10), type = "l", lwd = lwidth, col = "#006c66", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[1]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq, col = col.alpha("#006c66", .2) )
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)



i=2
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10), lwd = lwidth, col = "#006c66", yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)



i=3
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10), lwd = lwidth, col = "#006c66", yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)



i=4
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10), lwd = lwidth, col = "#006c66", yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)



i=5
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10), lwd = lwidth, col = "#006c66", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)




i=6
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10), lwd = lwidth, col = "#006c66", yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)



i=7
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10), lwd = lwidth, col = "#006c66",  yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)


i=8
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10), lwd = lwidth, col = "#006c66",  yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)


i=9
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10),, lwd = lwidth, col = "#006c66", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)


i=10
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10),, lwd = lwidth, col = "#006c66", yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)


i=11
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10),, lwd = lwidth, col = "#006c66", yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)

i=12
seq<-seq(min(data$y[,i]), max(data$y[,i]), length.out = 10)
plot(seq, colMeans(lambda[,,i]), type = "l", ylim = c(0, 10),, lwd = lwidth, col = "#006c66", yaxt = "n", xlab= "", cex.axis = 1.2)
title(str_to_title(sectors[i]), line=-1.2, cex.main = 1.6)
shade(apply(lambda[,,i],2,PI), seq,, col = col.alpha("#006c66", .2))
#lines(c(1,1), c(0, 10), lty = 2)
lines(c(0,100), c(mean(lambda), mean(lambda)), lty = 2)

mtext("Earnings from sector (USD) per two week period", outer = T, side = 1, line = .5, cex = 1.3)
mtext("Predicted number of reported cases per two week period", outer = T, side = 2, line = .5, adj = .2, cex = 1.3)
# dev.off()