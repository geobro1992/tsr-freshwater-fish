# Bayesian MCMC to model slopes/intercepts of change-point LAA regression
# against temperature as test of the temperature-size rule in freshwater fish

# required packages
library(tidyverse)
library(R2jags)
library(runjags)
library(tidybayes)
library(modelr)
library(tidyr)
library(broom)
library(brms)
library(broom.mixed)
library(ggplot2)
library(ggdist)
library(ggpubr)
library(ggeffects)



# load in LAA data
load(here::here(sprintf(
  "laa_temp_data.RData")))

laa_temp_data <- 
  laa_temp_data |>
  mutate(
    nages = map_int(fish.data, \(x) length(unique(x$est.age))),
    nobs = map_int(fish.data, nrow)
  ) |>
  filter(nages >= 5, nobs >= 30) |>
  filter(species != "cisco")


# create a species index to use within JAGS model for spp mixed effects
spp_index <- 
  laa_temp_data |>
  distinct(species) |>
  mutate(spp.index = 1:n())

# create a lake index to eventually use within JAGS model for lake mixed effects
lake_index <- 
  laa_temp_data |>
  distinct(species, state, lake.id) |>
  arrange(species, lake.id) |>
  mutate(lake.index = 1:n())

# set up jags data
mod_data <- 
  laa_temp_data |>
  left_join(spp_index, by = "species") |>
  left_join(lake_index, by = c("species", "state", "lake.id")) |>
  arrange(spp.index, lake.index) |>
  mutate(
    lakeyear.index = 1:n(),
    nobs = map_int(fish.data, nrow),
    minage = map_int(fish.data, \(x) return(min(x$est.age, na.rm = TRUE))),
    age.95 = map_dbl(fish.data, \(x) return(quantile(x$est.age, probs = 0.95, na.rm = TRUE))),
    maxage = map_int(fish.data, \(x) return(max(x$est.age, na.rm = TRUE))),
    meanage = map_dbl(fish.data, \(x) return(mean(x$est.age, na.rm = TRUE)))
  )


mod_laa_data <- 
  mod_data |>
  select(spp.index, lakeyear.index, fish.data, alk.age.str) |>
  unnest(fish.data) |>
  select(spp.index, lakeyear.index, est.age, length, alk.age.str) |>
  filter(length > 0)


mod_temp_data <- 
  mod_data |>
  select(lakeyear.index, ten.year.mean.temp)


#################
# scale temp and size data
mod_temp_data$ten.year.mean.temp = scale(mod_temp_data$ten.year.mean.temp)
mod_laa_data$length = log(mod_laa_data$length)



# JAGS data
mod_data_list <- list(
  laa_data = as.matrix(mod_laa_data),
  temp_data = as.matrix(mod_temp_data),
  spp_lake_index = as.matrix(select(mod_data, spp.index, lake.index)),
  nspp = length(unique(mod_laa_data$spp.index)),
  nlakes = length(unique(mod_data$lake.index)),
  nobs = nrow(mod_laa_data),
  nlakeyears = max(mod_laa_data$lakeyear.index),
  minage = mod_data$minage,
  maxage = mod_data$maxage
)



jags1 <- "model {
  
  for (j in 1:nlakeyears){
  
  pred_alpha[j] <- a1[spp_lake_index[j,1]] + re_a[spp_lake_index[j,2]]
  pred_beta1[j] <- a2[spp_lake_index[j,1]] + (b1 * temp_data[j, 2]) + (rse_b1[spp_lake_index[j,1]] * temp_data[j, 2]) + re_b1[spp_lake_index[j,2]]
  pred_beta2[j] <- a3[spp_lake_index[j,1]] + (b2 * temp_data[j, 2]) + (rse_b2[spp_lake_index[j,1]] * temp_data[j, 2]) + re_b2[spp_lake_index[j,2]]
  pred_cp[j]    <- a4[spp_lake_index[j,1]] + (b3 * temp_data[j, 2]) + (rse_cp[spp_lake_index[j,1]] * temp_data[j, 2]) + re_cp[spp_lake_index[j,2]]
  
  alpha[j] ~ dnorm(pred_alpha[j], 1 / alpha_sig)T(0,10000)
  beta1[j] ~ dnorm(pred_beta1[j], 1 / beta1_sig)T(0, 10000)
  beta2[j] ~ dnorm(pred_beta2[j], 1 / beta2_sig)T(0, beta1[j])
  cp[j] ~ dnorm(pred_cp[j], 1 / cp_sig)T(minage[j], maxage[j])

  pred.sam[j] <- alpha[j] + (beta1[j] * cp[j])

}

# likelihoods
for (i in 1:nobs) {
  mu[i] <- log(
    alpha[laa_data[i,2]] + 
    (beta1[laa_data[i,2]] * min(laa_data[i,3], cp[laa_data[i,2]])) + 
    (beta2[laa_data[i,2]] * max(0,laa_data[i,3] - cp[laa_data[i,2]])) )
  laa_data[i,4] ~ dnorm(mu[i], 1/nu_laa)
}

# priors on regression coefficients
    for (s in 1:nspp){

      a1[s] ~ dunif(0, 30)
      a2[s] ~ dunif(0, 20)
      a3[s] ~ dunif(0, 10)
      a4[s] ~ dunif(0, 10)

  rse_b1[s] ~ dnorm(0, rse_b1_tau)
  rse_b2[s] ~ dnorm(0, rse_b2_tau)
  rse_cp[s] ~ dnorm(0, rse_cp_tau)

    }

# priors on temperature coefficients
b1 ~ dunif(-2, 2)
b2 ~ dunif(-2, 2)
b3 ~ dunif(-2, 2)

# priors on process error
alpha_sig ~ dgamma(1, 1)
beta1_sig ~ dgamma(1, 1)
beta2_sig ~ dgamma(1, 1)
cp_sig ~ dgamma(1, 1)
nu_laa ~ dgamma(1, 1)


# priors for random effects
re_a_sigma ~ dgamma(1, 1)
re_a_tau <- 1 / (re_a_sigma * re_a_sigma)
re_b1_sigma ~ dgamma(1, 1)
re_b1_tau <- 1 / (re_b1_sigma * re_b1_sigma)
re_b2_sigma ~ dgamma(1, 1)
re_b2_tau <- 1 / (re_b2_sigma * re_b2_sigma)
re_cp_sigma ~ dgamma(1, 1)
re_cp_tau <- 1 / (re_cp_sigma * re_cp_sigma)

rse_b1_sigma ~ dgamma(1, 1)
rse_b1_tau <- 1 / (rse_b1_sigma * rse_b1_sigma)
rse_b2_sigma ~ dgamma(1, 1)
rse_b2_tau <- 1 / (rse_b2_sigma * rse_b2_sigma)
rse_cp_sigma ~ dgamma(1, 1)
rse_cp_tau <- 1 / (rse_cp_sigma * rse_cp_sigma)

for (l in 1:nlakes) {
  re_a[l] ~ dnorm(0, re_a_tau)
  re_b1[l] ~ dnorm(0, re_b1_tau)
  re_b2[l] ~ dnorm(0, re_b2_tau)
  re_cp[l] ~ dnorm(0, re_cp_tau)
}


}"

writeLines(jags1, con="M1.txt")


jags_params <- c(
 paste0("a", 1:4), paste0("b", 1:3), "rse_b1", "rse_b2", "rse_cp", "cp", "pred.sam")


m <- jags.parallel(
  data = mod_data_list, 
  inits = NULL,
  parameters.to.save = jags_params,
  model.file = "M1.txt",
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 40000,
  n.thin = 100
)



save(m, file = "tsr_JAGS.RData")



#####################
# posterior summaries
#####################

# fixed effects

# posterior draws
mcmc = m$BUGSoutput$sims.list

# scaled temperature range
ts = seq(-2, 2, length.out = 100)


#
# temperature on juvenile growth
#

## Calculate the fitted values
newdata = data.frame(x = ts)
Xmat = model.matrix(~x, newdata)
coefs = cbind(rowMeans(mcmc[["a2"]]), mcmc[["b1"]])
fit = (coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), ts)

g1 = ggplot(newdata, aes(x = ts, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Temperature") + 
  ggtitle("Juvenile Growth Rate") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) 


#
# temperature on adult growth
#


## Calculate the fitted values
newdata = data.frame(x = ts)
Xmat = model.matrix(~x, newdata)
coefs = cbind(rowMeans(mcmc[["a3"]]), mcmc[["b2"]])
fit = (coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), ts)

g2 = ggplot(newdata, aes(x = ts, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Temperature") + 
  ggtitle("Adult Growth Rate") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) 


#
# temperature on age at maturity
#

## Calculate the fitted values
newdata = data.frame(x = ts)
Xmat = model.matrix(~x, newdata)
coefs = cbind(rowMeans(mcmc[["a4"]]), mcmc[["b3"]])
fit = (coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), ts)

g3 = ggplot(newdata, aes(x = ts, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Temperature") + 
  ggtitle("Age at Maturity") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) 


#
# temperature on size at maturity
#

## Calculate the fitted values
fit0 = rep(rowMeans(mcmc[["a1"]]), length(ts))

newdata = data.frame(x = ts)
Xmat = model.matrix(~x, newdata)
coefs = cbind(rowMeans(mcmc[["a4"]]), mcmc[["b3"]])
fit1 = (coefs %*% t(Xmat))

newdata = data.frame(x = ts)
Xmat = model.matrix(~x, newdata)
coefs = cbind(rowMeans(mcmc[["a2"]]), mcmc[["b1"]])
fit2 = (coefs %*% t(Xmat))

fit = fit0 + (fit1 * fit2) 

newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), ts)


g4 = ggplot(newdata, aes(x = ts, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Temperature") + 
  ggtitle("Size at Maturity") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

gg1 = ggarrange(g1, g2, g3, g4, nrow = 2, ncol = 2)



########################
# species random effects


p1 = as.data.frame(m$BUGSoutput$sims.list)


#
# temperature effect on juvenlie growth
#

pred = data.frame(b10 = p1$b1)
pred[,2] = (p1$b1) + (p1$rse_b1.1)
pred[,3] = (p1$b1) + (p1$rse_b1.2)
pred[,4] = (p1$b1) + (p1$rse_b1.3)
pred[,5] = (p1$b1) + (p1$rse_b1.4)
pred[,6] = (p1$b1) + (p1$rse_b1.5)
pred[,7] = (p1$b1) + (p1$rse_b1.6)
pred[,8] = (p1$b1) + (p1$rse_b1.7)

pred.b1 = gather(
  pred,
  key = "species",
  value = "value", -b10)


pred.b1$species = as.factor(pred.b1$species)
levels(pred.b1$species) = paste(spp_index$species)
pred.b1$species = factor(pred.b1$species, levels = rev(c("yellow_perch", "northern_pike", "walleye", "black_crappie", "smallmouth_bass", "largemouth_bass", "bluegill")))
levels(pred.b1$species) = rev(c("P. flavescens", "E. lucius", "S. vitreus", "P. nigromaculatus", "M. dolomieu", "M. salmoides", "L. macrochirus"))

aggregate(value ~ species, pred.b1, FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))

g1 = ggplot(pred.b1, aes(value, species)) +
  stat_halfeye() +
  xlab("Parameter estimate") + ylab("Temperature effect on juvenile growth") +
  theme_classic() + geom_vline(xintercept=0, linetype="dashed")+
  theme(axis.text.y = element_text(face = 'italic'))



#
# temperature effect on adult growth
#

pred = data.frame(b20 = p1$b2)
pred[,2] = (p1$b2) + (p1$rse_b2.1)
pred[,3] = (p1$b2) + (p1$rse_b2.2)
pred[,4] = (p1$b2) + (p1$rse_b2.3)
pred[,5] = (p1$b2) + (p1$rse_b2.4)
pred[,6] = (p1$b2) + (p1$rse_b2.5)
pred[,7] = (p1$b2) + (p1$rse_b2.6)
pred[,8] = (p1$b2) + (p1$rse_b2.7)

pred.b2 = gather(
  pred,
  key = "species",
  value = "value", -b20)


pred.b2$species = as.factor(pred.b2$species)
levels(pred.b2$species) = paste(spp_index$species)
pred.b2$species = factor(pred.b2$species, levels = rev(c("yellow_perch", "northern_pike", "walleye", "black_crappie", "smallmouth_bass", "largemouth_bass", "bluegill")))
levels(pred.b2$species) = rev(c("P. flavescens", "E. lucius", "S. vitreus", "P. nigromaculatus", "M. dolomieu", "M. salmoides", "L. macrochirus"))

aggregate(value ~ species, pred.b2, FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))


g2 = ggplot(pred.b2, aes(value, species)) +
  stat_halfeye() +
  xlab("Parameter estimate") + ylab("Temperature effect on adult growth") +
  theme_classic() + geom_vline(xintercept=0, linetype="dashed")+
  theme(axis.text.y = element_text(face = 'italic'))


#
# temperature effect on age at maturity
#

pred = data.frame(cp0 = p1$b3)
pred[,2] = (p1$b3) + (p1$rse_cp.1)
pred[,3] = (p1$b3) + (p1$rse_cp.2)
pred[,4] = (p1$b3) + (p1$rse_cp.3)
pred[,5] = (p1$b3) + (p1$rse_cp.4)
pred[,6] = (p1$b3) + (p1$rse_cp.5)
pred[,7] = (p1$b3) + (p1$rse_cp.6)
pred[,8] = (p1$b3) + (p1$rse_cp.7)

pred.cp = gather(
  pred,
  key = "species",
  value = "value", -cp0)


pred.cp$species = as.factor(pred.cp$species)
levels(pred.cp$species) = paste(spp_index$species)
pred.cp$species = factor(pred.cp$species, levels = rev(c("yellow_perch", "northern_pike", "walleye", "black_crappie", "smallmouth_bass", "largemouth_bass", "bluegill")))
levels(pred.cp$species) = rev(c("P. flavescens", "E. lucius", "S. vitreus", "P. nigromaculatus", "M. dolomieu", "M. salmoides", "L. macrochirus"))

aggregate(value ~ species, pred.cp, FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))


g3 = ggplot(pred.cp, aes(value, species)) +
  stat_halfeye() +
  xlab("Parameter estimate") + ylab("Temperature effect on age at maturity") +
  theme_classic() + geom_vline(xintercept=0, linetype="dashed") +
  theme(axis.text.y = element_text(face = 'italic'))




#
# temperature effect on size at maturity
#


ts = seq(-2, 2, length.out = 100)

pred = array(dim =c(length(p1$b3), length(ts), 7))

for (i in 1:length(ts)){
  pred[,i,1] = ((p1$b3) + (p1$rse_cp.1*ts[i])) * ((p1$b1) + (p1$rse_b1.1*ts[i]))
  pred[,i,2] = ((p1$b3) + (p1$rse_cp.2*ts[i])) * ((p1$b1) + (p1$rse_b1.2*ts[i]))
  pred[,i,3] = ((p1$b3) + (p1$rse_cp.3*ts[i])) * ((p1$b1) + (p1$rse_b1.3*ts[i]))
  pred[,i,4] = ((p1$b3) + (p1$rse_cp.4*ts[i])) * ((p1$b1) + (p1$rse_b1.4*ts[i]))
  pred[,i,5] = ((p1$b3) + (p1$rse_cp.5*ts[i])) * ((p1$b1) + (p1$rse_b1.5*ts[i]))
  pred[,i,6] = ((p1$b3) + (p1$rse_cp.6*ts[i])) * ((p1$b1) + (p1$rse_b1.6*ts[i]))
  pred[,i,7] = ((p1$b3) + (p1$rse_cp.7*ts[i])) * ((p1$b1) + (p1$rse_b1.7*ts[i]))
  
}

pred.cps = array(dim =c((length(ts)-1) * length(p1$b3), 7))

for(i in 1:7){
  pred.cps[,i] = as.vector(apply(pred[,,i], 1, FUN = diff))
}

pred.cps = as.data.frame(pred.cps)

pred.cps = gather(
  pred.cps,
  key = "species",
  value = "value")

pred.cps$value = pred.cps$value/(ts[2]-ts[1])

pred.cps$species = as.factor(pred.cps$species)
levels(pred.cps$species) = paste(spp_index$species)
pred.cps$species = factor(pred.cps$species, levels = rev(c("yellow_perch", "northern_pike", "walleye", "black_crappie", "smallmouth_bass", "largemouth_bass", "bluegill")))
levels(pred.cps$species) = rev(c("P. flavescens", "E. lucius", "S. vitreus", "P. nigromaculatus", "M. dolomieu", "M. salmoides", "L. macrochirus"))

aggregate(value ~ species, pred.cps, FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))


g4 = ggplot(pred.cps, aes(value, species)) +
  stat_halfeye() +
  xlab("Parameter estimate") + ylab("Temperature effect on size at maturity") +
  theme_classic() + geom_vline(xintercept=0, linetype="dashed") +
  theme(axis.text.y = element_text(face = 'italic')) + xlim(-1.5,1.5)


gg1 = ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2)



#########################
# quantile age regression
#########################

mod_data$species = factor(mod_data$species, levels = c("yellow_perch", "northern_pike", "walleye", "black_crappie", "smallmouth_bass", "largemouth_bass", "bluegill"))
levels(mod_data$species) = c("P. flavescens", "E. lucius", "S. vitreus", "P. nigromaculatus", "M. dolomieu", "M. salmoides", "L. macrochirus")


m1 = glm(age.95 ~ scale(ten.year.mean.temp)*species, data = mod_data, family = "poisson")
summary(m1)

# equivalent in Bayes
bayes.brms <- brm(age.95 ~ scale(ten.year.mean.temp)*species + (1 | lake.id), 
                  family = gaussian(),
                  data = mod_data,
                  chains = 3, # nb of chains
                  iter = 30000, # nb of iterations, including burnin
                  warmup = 20000, # burnin
                  thin = 10,  # thinning
                  set_prior("normal(0,10)", class = "b"))


# plot predictions
g1 = mod_data %>%
  group_by(species) %>%
  modelr::data_grid(ten.year.mean.temp = seq_range(ten.year.mean.temp, n = 51), lake.id = NA) %>%
  add_epred_draws(bayes.brms, re_formula = NA) %>%
  ggplot(aes(x = ten.year.mean.temp, y = age.95, color = species, fill = species)) +
  stat_lineribbon(aes(y = .epred), .width = c(.95), alpha = 0.9) +
  geom_point(data = mod_data, size = 0.5, alpha = 0.5) +
  theme_classic() + xlab("Temperature") + ylab("Maximum Age") + ylim(0, 20) +
  facet_wrap(~species) + theme(legend.position = "none") + 
  theme(strip.text = element_text(face = 'italic'))




################
# Growth curves
################


# max age predictions from GLM for plotting growth curves
max.age.preds = ggpredict(m1,c("ten.year.mean.temp [2098, 2636, 3174]", "species"))


p1 = m$BUGSoutput$mean

# scaled temperaure range
ts = c(-2, 0, 2)

# species order
ss = c(6,4,3,2,7,5,1)



par(cex.main = 1.5, mar = c(3.5, 4, 3.5, 0) + 0.1, mgp = c(2.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

par(mfrow = c(3,3))

count = 1

for(s in ss){

pred.alpha = p1$a1[s]
pred.beta1 = p1$a2[s]+(p1$b1*ts) + (p1$rse_b1[s]*ts)
pred.beta2 = p1$a3[s]+(p1$b2*ts) + (p1$rse_b2[s]*ts)
pred.cp = p1$a4[s]+(p1$b3*ts) + (p1$rse_cp[s]*ts)

pred.length.cold = data.frame(age = seq(from = 0, to = filter(max.age.preds, group == levels(mod_data$species)[count])$predicted[1], length.out = 100))
pred.length.mean = data.frame(age = seq(from = 0, to = filter(max.age.preds, group == levels(mod_data$species)[count])$predicted[2], length.out = 100))
pred.length.warm = data.frame(age = seq(from = 0, to = filter(max.age.preds, group == levels(mod_data$species)[count])$predicted[3], length.out = 100))

  for(i in 1:length(pred.length.mean[,1])){
    
    pred.length.cold[i,2] = log(
      pred.alpha + (pred.beta1[1] * pmin(pred.length.cold[i,1], pred.cp[1])) + 
        (pred.beta2[1] * pmax(0, (pred.length.cold[i,1] - pred.cp[1])))
    )    
    
    pred.length.mean[i,2] = log(
      pred.alpha + (pred.beta1[2] * pmin(pred.length.mean[i,1], pred.cp[2])) + 
      (pred.beta2[2] * pmax(0, (pred.length.mean[i,1] - pred.cp[2])))
    )  
    
    pred.length.warm[i,2] = log(
      pred.alpha + (pred.beta1[3] * pmin(pred.length.warm[i,1], pred.cp[3])) + 
        (pred.beta2[3] * pmax(0, (pred.length.warm[i,1] - pred.cp[3])))
    )    
    
  
  }
  

lab = levels(mod_data$species)[count]

plot(pred.length.mean$age, exp(pred.length.mean[,2]), type = "l", lwd = 2, 
     main = bquote(italic(.(paste(lab[1])))), xlab = "Age", ylab = "Length",
     xlim = c(0, max(filter(max.age.preds, group == levels(mod_data$species)[count])$predicted)+1),
     ylim = c(0, max(c(exp(pred.length.cold[,2]), exp(pred.length.mean[,2]), exp(pred.length.warm[,2]))))+2)
lines(pred.length.cold$age, exp(pred.length.cold[,2]), lwd = 2, col = "blue")
lines(pred.length.warm$age, exp(pred.length.warm[,2]), lwd = 2, col = "red")

count = count+1

}




