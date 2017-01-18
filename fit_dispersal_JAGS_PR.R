#===========================================================================================
# FIT HIERARCHICAL DISPERSAL MODELS TO DATA FROM MAPPED ARRAYS OF SEED TRAPS AND TREES
#===========================================================================================

# MODEL DESCRIPTION:
# Multiple islands -- parameters may differ by island or group of islands (FIXED effects)
# Multiple sites within each island -- parameters will always differ among sites (RANDOM effects)
# Multiple traps within each site -- likelihood of data (# seeds/trap) is negative binomial
# Fecundity per tree is proportional to basal area

options(device=windows)

library(R2jags)
library(matrixStats)
library(rgl)
library(lme4)
library(arm)
library(RColorBrewer)

source("k2d.predict.r")
source("k2d_radial_quantiles.r")
source("seeds.fit2sim.r")
source("minor.ticks.axis.r")
source("filled.contour.plt.r")
source("filled.contour.key.r")

#==================================================================
# READ IN DATA
# Premna (PR)
#==================================================================

species <- "PR"

traps <- read.table(paste(species, "trapseed.txt", sep="_"), sep="\t", header=T)
trees <- read.table(paste(species, "trees.txt", sep="_"), sep="\t", header=T)

if(species=="PS")
{
  trees$Total_dbh_summed <- trees$Total_dbh_summed/1e3  # convert from mm to m
  trees$Total_dbh_averaged <- trees$Total_dbh_averaged/1e3
  trees$ba.sum <- trees$ba.sum*1e4  # convert from m2 to cm2
  trees$ba.avg <- trees$ba.avg*1e4
} else
{
  trees$totdbh <- trees$totdbh/1e3
  trees$ba <- trees$ba*1e4
}

# Sort to ensure each island forms a contiguous block
traps <- traps[order(traps$site),]
trees <- trees[order(trees$site),]


#==================================================================
# ASSIGN DATA OBJECTS TO A LIST TO BE PASSED TO JAGS
#==================================================================

# Distance matrix (traps x trees)

trap.x <- matrix(rep(traps$x, nrow(trees)), nrow = nrow(traps))  # matrix of trap x-coords
trap.y <- matrix(rep(traps$y, nrow(trees)), nrow = nrow(traps))  # matrix of trap y-coords
tree.x <- matrix(rep(trees$x, each = nrow(traps)), nrow = nrow(traps))  # matrix of tree x-coords
tree.y <- matrix(rep(trees$y, each = nrow(traps)), nrow = nrow(traps))  # matrix of tree y-coords

dist <- sqrt((trap.x - tree.x)^2 + (trap.y - tree.y)^2)    # trap x tree distance matrix with sites in blocks

cc <- matrix(NA,nrow(traps),2)  # start and end cols of dist for each site
for(i in 1:nrow(traps))
{
  cc[i,1:2] <- range(which(trees$site==traps$site[i]))
}

# Assemble list of various data objects to be used in all models
# Note that JAGS doesn't know what to do with factors, so we convert all factors to numeric values

site.island <- as.vector(tapply(as.numeric(traps$island), as.numeric(traps$site), mean))  	# island membership of each site
trap.site <- as.numeric(traps$site)
tree.site <- as.numeric(trees$site)
trap.area <- 0.5					  # surface area of traps in m2
if(species=="PS")
{
  dbh <- trees$Total_dbh_averaged  # dbh in m
  basal.area <- trees$ba.avg       # basal area of trees in cm2; for multi-stem Psychotria use average
} else
{
  dbh <- trees$totdbh
  basal.area <- trees$ba
}
seeds <- traps$total			  # seed counts


#==================================================================
# REGRESSION TO PREDICT CANOPY DIAMETER FROM DBH
#==================================================================

crown.df <- read.csv("Crown Measurements_final_data_file.csv")
names(crown.df)[names(crown.df)=="species"] <- "spp"
names(crown.df)[names(crown.df)=="dbh"] <- "dbh.formula"
crown.df$total_dbh_multisum <- crown.df$total_dbh_multisum/1e3
crown.df$total_dbh_multiavg <- crown.df$total_dbh_multiavg/1e3
crown.df <- data.frame(crown.df, canopy_diam = crown.df$canopy_radius*2,
                       dbh = switch(species, PS=crown.df$total_dbh_multiavg, PR=crown.df$total_dbh_multisum))
summary(crown.df)

lm1.crown.diam <- lm(canopy_diam ~ dbh, data = crown.df, subset = spp==switch(species, PS="psychotria", PR="premna"))
lm2.crown.diam <- lm(canopy_diam ~ dbh*island, data = crown.df, subset = spp==switch(species, PS="psychotria", PR="premna"))
summary(lm1.crown.diam)
summary(lm2.crown.diam)
AICctab(lm1.crown.diam, lm2.crown.diam)

plot(lm1.crown.diam$model[,2], lm1.crown.diam$model[,1], xlab="dbh (m)", ylab="canopy diameter (m)")
abline(coef(lm1.crown.diam))

# predict crown diam; use one relationship for all islands
crown.diam <- predict(lm1.crown.diam, newdata = data.frame(dbh = dbh))


#==================================================================
# DEFINE INITIAL VALUES AND CALL JAGS TO FIT THE MODELS
# (Model files in native JAGS code are in separate R scripts)
#==================================================================

# VERY IMPORTANT!!!  The functions used below to generate initial values of kernel params will 
# dictate the number of island-specific fixed effects used in the model. For example,
# lshape.grp = rnorm(4,mean,sd) could be used to specify a distinct value for Guam vs. all other islands.  
# The actual correspondence between elements of parameter vectors and the values assigned
# to each island depends on the user-defined index vector kern.indx.

# Model structures
#
# all = all islands have same kernel params
# G = Guam kernel differs from others
# GRST = each island has its own kernel params

###############################
## POWER-EXPONENTIAL KERNELS ##
###############################

#------------------------------------------------------------------
# Power-exponential dispersal kernel
# Kernel identical for all islands
#------------------------------------------------------------------

# Index vector matching kernel parameters to islands
kern.indx <- c(1,1,1,1)

inits.pexp.all <- function()
{
  list(mu.lfec = runif(1,1,3), sigma.lfec = runif(1,0.5,1), lscale.grp=log(runif(max(kern.indx),0.05,0.1)), 
       r1 = runif(1,1,2), r2 = runif(1,0.5,1))
}


# Call JAGS to fit the model
jags.PR.pexp.all <- jags.parallel(data = c("site.island","trap.site","dist", "cc",
                                           "trap.area","basal.area","kern.indx","seeds"), 
                                  inits = inits.pexp.all, 
                                  parameters.to.save = c("ll.data","dev.data","dev.lfec","mu.lfec","sigma.lfec","lfec","lscale.grp","lshape.grp","ldmean.grp","ldmode.grp","r1","r2","fit"), 
                                  model.file = "HBdisp.pexp.r", 
                                  n.chains = 3, n.iter = 210000, n.burnin = 10000, n.thin = 200, 
                                  DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.PR.pexp.all,2);plot(jags.PR.pexp.all)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,!is.element(substring(dimnames(x)[[3]],1,3), c("fit","ll."))]), 1:3, list(x=jags.PR.pexp.all$BUGS$sims.array), SIMPLIFY=F)), ask=T)


#------------------------------------------------------------------
# Power-exponential dispersal kernel
# Kernel identical for all islands EXCEPT Guam
#------------------------------------------------------------------

# Index vectors matching kernel parameters to islands
kern.indx <- c(1,2,2,2)

inits.pexp.G <- function() 
{
  list(mu.lfec = runif(1,1,3), sigma.lfec = runif(1,0.5,1), lscale.grp=log(runif(max(kern.indx),0.05,0.1)), 
       r1 = runif(1,1,2), r2 = runif(1,0.5,1))
}


# Call JAGS to fit the model
jags.PR.pexp.G <- jags.parallel(data = c("site.island","trap.site","dist", "cc",
                                         "trap.area","basal.area","kern.indx","seeds"), 
                                inits = inits.pexp.G, 
                                parameters.to.save = c("ll.data","dev.data","dev.lfec","mu.lfec","sigma.lfec","lfec","lscale.grp","lshape.grp","ldmean.grp","ldmode.grp","r1","r2","fit"), 
                                model.file = "HBdisp.pexp.r", 
                                n.chains = 3, n.iter = 210000, n.burnin = 10000, n.thin = 200, 
                                DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.PR.pexp.G,2);plot(jags.PR.pexp.G)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,!is.element(substring(dimnames(x)[[3]],1,3), c("fit","ll."))]), 1:3, list(x=jags.PR.pexp.G$BUGS$sims.array), SIMPLIFY=F)), ask=T)


#------------------------------------------------------------------
# Power-exponential dispersal kernel
# Kernel distinct for each island
#------------------------------------------------------------------

# Index vectors matching fec and scale parameters to islands
kern.indx <- c(1,2,3,4)

inits.pexp.GRST <- function() 
{
  list(mu.lfec = runif(1,1,3), sigma.lfec = runif(1,0.5,1), lscale.grp=log(runif(max(kern.indx),0.05,0.1)), 
       r1 = runif(1,1,2), r2 = runif(1,0.5,1))
}


# Call JAGS to fit the model
jags.PR.pexp.GRST <- jags.parallel(data = c("site.island","trap.site","dist", "cc",
                                            "trap.area","basal.area","kern.indx","seeds"), 
                                   inits = inits.pexp.GRST, 
                                   parameters.to.save = c("ll.data","dev.data","dev.lfec","mu.lfec","sigma.lfec","lfec","lscale.grp","lshape.grp","ldmean.grp","ldmode.grp","r1","r2","fit"), 
                                   model.file = "HBdisp.pexp.r", 
                                   n.chains = 3, n.iter = 210000, n.burnin = 10000, n.thin = 200, 
                                   DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.PR.pexp.GRST,2);plot(jags.PR.pexp.GRST)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,!is.element(substring(dimnames(x)[[3]],1,3), c("fit","ll."))]), 1:3, list(x=jags.PR.pexp.GRST$BUGS$sims.array), SIMPLIFY=F)), ask=T)




#################
## 2Dt KERNELS ##
#################

#------------------------------------------------------------------
# 2Dt dispersal kernel
# Kernel identical for all islands
#------------------------------------------------------------------

# Index vector matching kernel parameters to islands
kern.indx <- c(1,1,1,1)

inits.2Dt.all <- function() 
{
  list(mu.lfec = runif(1,1,3), sigma.lfec = runif(1,0.5,1), lscale.grp=log(runif(max(kern.indx),0.05,0.1)), 
       r1 = runif(1,1,2), r2 = runif(1,0.5,1))
}

# Call JAGS to fit the model
jags.PR.2Dt.all <- jags.parallel(data = c("site.island","trap.site","dist", "cc",
                                          "trap.area","basal.area","kern.indx","seeds"), 
                                 inits = inits.2Dt.all, 
                                 parameters.to.save = c("ll.data","dev.data","dev.lfec","mu.lfec","sigma.lfec","lfec","lscale.grp","lshape.grp","ldmean.grp","ldmode.grp","r1","r2","fit"), 
                                 model.file = "HBdisp.2Dt.r", 
                                 n.chains = 3, n.iter = 110000, n.burnin = 10000, n.thin = 100, 
                                 DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.PR.2Dt.all,2);plot(jags.PR.2Dt.all)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,!is.element(substring(dimnames(x)[[3]],1,3), c("fit","ll."))]), 1:3, list(x=jags.PR.2Dt.all$BUGS$sims.array), SIMPLIFY=F)), ask=T)


#------------------------------------------------------------------
# 2Dt dispersal kernel
# Kernel identical for all islands EXCEPT Guam
#------------------------------------------------------------------

# Index vectors matching kernel parameters to islands
kern.indx <- c(1,2,2,2)

inits.2Dt.G <- function() 
{
  list(mu.lfec = runif(1,1,3), sigma.lfec = runif(1,0.5,1), lscale.grp=log(runif(max(kern.indx),0.05,0.1)), 
       r1 = runif(1,1,2), r2 = runif(1,0.5,1))
}


# Call JAGS to fit the model
jags.PR.2Dt.G <- jags.parallel(data = c("site.island","trap.site","dist", "cc",
                                        "trap.area","basal.area","kern.indx","seeds"), 
                               inits = inits.2Dt.G, 
                               parameters.to.save = c("ll.data","dev.data","dev.lfec","mu.lfec","sigma.lfec","lfec","lscale.grp","lshape.grp","ldmean.grp","ldmode.grp","r1","r2","fit"), 
                               model.file = "HBdisp.2Dt.r", 
                               n.chains = 3, n.iter = 110000, n.burnin = 10000, n.thin = 100, 
                               DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.PR.2Dt.G,2);plot(jags.PR.2Dt.G)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,!is.element(substring(dimnames(x)[[3]],1,3), c("fit","ll."))]), 1:3, list(x=jags.PR.2Dt.G$BUGS$sims.array), SIMPLIFY=F)), ask=T)


#------------------------------------------------------------------
# 2Dt dispersal kernel
# Kernel distinct for each island
#------------------------------------------------------------------

# Index vectors matching fec and scale parameters to islands
kern.indx <- c(1,2,3,4)

inits.2Dt.GRST <- function() 
{
  list(mu.lfec = runif(1,1,3), sigma.lfec = runif(1,0.5,1), lscale.grp=log(runif(max(kern.indx),0.05,0.1)), 
       r1 = runif(1,1,2), r2 = runif(1,0.5,1))
}

# Call JAGS to fit the model
jags.PR.2Dt.GRST <- jags.parallel(data = c("site.island","trap.site","dist", "cc",
                                           "trap.area","basal.area","kern.indx","seeds"), 
                                  inits = inits.2Dt.GRST, 
                                  parameters.to.save = c("ll.data","dev.data","dev.lfec","mu.lfec","sigma.lfec","lfec","lscale.grp","lshape.grp","ldmean.grp","ldmode.grp","r1","r2","fit"), 
                                  model.file = "HBdisp.2Dt.r", 
                                  n.chains = 3, n.iter = 110000, n.burnin = 10000, n.thin = 100, 
                                  DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.PR.2Dt.GRST,2);plot(jags.PR.2Dt.GRST)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,!is.element(substring(dimnames(x)[[3]],1,3), c("fit","ll."))]), 1:3, list(x=jags.PR.2Dt.GRST$BUGS$sims.array), SIMPLIFY=F)), ask=T)


#==================================================================
# MODEL COMPARISON
#==================================================================

# Parameter estimates and model selection by DIC and WAIC
mod.list <- list(jags.PR.pexp.all$BUGS, jags.PR.pexp.G$BUGS, jags.PR.pexp.GRST$BUGS,
                 jags.PR.2Dt.all$BUGS, jags.PR.2Dt.G$BUGS, jags.PR.2Dt.GRST$BUGS)
kern.list <- list("pexp","pexp","pexp","2Dt","2Dt","2Dt")
kern.indx.list <- list(c(1,1,1,1),c(1,2,2,2),c(1,2,3,4),c(1,1,1,1),c(1,2,2,2),c(1,2,3,4))
scale.names <- c("scale.G","scale.R","scale.S","scale.T")

mod.tab <- data.frame(kernel = unlist(kern.list), 
                      islands = c("all", "G", "GRST"),
                      mu.lfec = NA, sigma.lfec = NA,
                      scale.G = NA, scale.R = NA, scale.S = NA, scale.T = NA, r1 = NA, r2 = NA,
                      Dbar = NA, pV = NA, DIC = NA, dDIC = NA, wDIC = NA, 
                      lppd = NA, pWAIC = NA, WAIC = NA, dWAIC = NA, wWAIC = NA)

mod.tab$mu.lfec <- sapply(mod.list, function(x) {
  paste(round(x$mean$mu.lfec, 2), " (", round(quantile(x$sims.list$mu.lfec, 0.025), 2), ", ", 
        round(quantile(x$sims.list$mu.lfec, 0.975), 2), ")", sep="")})
mod.tab$sigma.lfec <- sapply(mod.list, function(x) {
  paste(round(x$mean$sigma.lfec, 2), " (", round(quantile(x$sims.list$sigma.lfec, 0.025), 2), ", ", 
        round(quantile(x$sims.list$sigma.lfec, 0.975), 2), ")", sep="")})
mod.tab$r1 <- sapply(mod.list, function(x) {
  paste(round(x$mean$r1, 2), " (", round(quantile(x$sims.list$r1, 0.025), 2), ", ", 
        round(quantile(x$sims.list$r1, 0.975), 2), ")", sep="")})
mod.tab$r2 <- sapply(mod.list, function(x) {
  paste(round(x$mean$r2, 2), " (", round(quantile(x$sims.list$r2, 0.025), 2), ", ", 
        round(quantile(x$sims.list$r2, 0.975), 2), ")", sep="")})
for(j in 1:4)
  mod.tab[,scale.names[j]] <- sapply(1:length(mod.list), function(i) {
    paste(round(mean(exp(mod.list[[i]]$sims.list$lscale.grp[,kern.indx.list[[i]][j]])), 2), " (", 
          round(quantile(exp(mod.list[[i]]$sims.list$lscale.grp[,kern.indx.list[[i]][j]]), 0.025), 2), ", ", 
          round(quantile(exp(mod.list[[i]]$sims.list$lscale.grp[,kern.indx.list[[i]][j]]), 0.975), 2), ")", sep="")})

mod.tab$Dbar <- sapply(mod.list, function(x) x$mean$dev.data)
mod.tab$pV <- sapply(mod.list, function(x) var(x$sims.list$dev.data)/2)
mod.tab$DIC <- mod.tab$Dbar + mod.tab$pV
mod.tab$dDIC <- mod.tab$DIC - min(mod.tab$DIC)
mod.tab$wDIC <- exp(-mod.tab$dDIC/2)/sum(exp(-mod.tab$dDIC/2))

# Calculate Watanabe-Akaike/widely available information criterion for each model
mod.tab$lppd <- sapply(mod.list, function(x) sum(colMeans(exp(x$sims.list$ll.data))))
mod.tab$pWAIC <- sapply(mod.list, function(x) sum(colVars(x$sims.list$ll.data)))
mod.tab$WAIC <- 2*(mod.tab$pWAIC - mod.tab$lppd)
mod.tab$dWAIC <- mod.tab$WAIC - min(mod.tab$WAIC)
mod.tab$wWAIC <- exp(-mod.tab$dWAIC/2)/sum(exp(-mod.tab$dWAIC/2))

mod.tab <- mod.tab[order(mod.tab$DIC),]
write.table(mod.tab, paste("mod.tab", species, "txt", sep="."), row.names=F, sep="\t")
rm(mod.list);rm(kern.list);rm(kern.indx.list);rm(scale.names)


# Compare mean and modal displacement and proportion of seeds falling < 2 m (PS) or 4 m (PR) from source
dist.tab <- data.frame(kernel = mod.tab$kern[which.min(mod.tab$DIC[mod.tab$islands=="G"])],
                       island = c("Guam","islands with birds"), dmean = NA, dmode = NA, prc = NA)
jags.fit <- get(paste("jags", species, dist.tab$kern[1], "G", sep="."))$BUGS
for(j in 1:2)
{
  dmean <- exp(jags.fit$sims.list$ldmean.grp[,j])
  dist.tab$dmean[j] <- paste(round(mean(dmean), 2), " (", 
                             round(quantile(dmean, 0.025), 2), ", ", 
                             round(quantile(dmean, 0.975), 2), ")", sep="")
  dmode <- exp(jags.fit$sims.list$ldmode.grp[,j])
  dist.tab$dmode[j] <- paste(round(mean(dmode), 2), " (", 
                             round(quantile(dmode, 0.025), 2), ", ", 
                             round(quantile(dmode, 0.975), 2), ")", sep="")
  pfunc <- switch(as.character(dist.tab$kernel[1]), pexp=ppexp, p2Dt)
  prc <- sapply(1:jags.fit$n.sims, function(i) pfunc(switch(species, PS=2, PR=4), 
                                                     shape=exp(jags.fit$sims.list$lshape.grp[i,j]), 
                                                     scale=exp(jags.fit$sims.list$lscale.grp[i,j]))$value)
  dist.tab$prc[j] <- paste(round(mean(prc), 2), " (", 
                           round(quantile(prc, 0.025), 2), ", ", 
                           round(quantile(prc, 0.975), 2), ")", sep="")
}
write.table(dist.tab, paste("dist.tab", species, "txt", sep="."), row.names=F, sep="\t")
rm(list=c("jags.fit","dmean","dmode","pfunc","prc"))


#-------------------------------------------------------------------------------------------
# Calculate the proportion of seed rain at each site that falls outside conspecific canopy
#-------------------------------------------------------------------------------------------

seed.canopy.tab <- data.frame(island=levels(trees$island)[site.island], site=levels(trees$site),
                              p.outside=NA, CI.025=NA, CI.975=NA)

kern <- mod.tab$kernel[1]
jags.fit <- get(paste("jags", species, kern, "G", sep="."))$BUGS$sims.list
npts <- 100
n.sims <- 1000

dev.new(width=14, height=7)
par(mfrow=c(3,5))

for(j in 1:nrow(seed.canopy.tab))
{ 
  site <- seed.canopy.tab$site[j]
  isl.grp <- ifelse(site.island[j]==1, 1, 2)
  map.traps <- traps[traps$site==seed.canopy.tab$site[j],]
  map.trees <- trees[trees$site==seed.canopy.tab$site[j],]
  if(species=="PS") names(map.trees)[names(map.trees)=="ba.avg"] <- "ba"
  
  map.grid <- expand.grid(x=seq(min(c(map.traps$x, map.trees$x)), max(c(map.traps$x, map.trees$x)), length=npts),
                          y=seq(min(c(map.traps$y, map.trees$y)), max(c(map.traps$y, map.trees$y)), length=npts))
  grid.x <- matrix(rep(map.grid$x, nrow(map.trees)), nrow = nrow(map.grid))  # matrix of grid x-coords
  grid.y <- matrix(rep(map.grid$y, nrow(map.trees)), nrow = nrow(map.grid))  # matrix of grid y-coords
  maptrees.x <- matrix(rep(map.trees$x, each = nrow(map.grid)), nrow = nrow(map.grid))  # matrix of tree x-coords
  maptrees.y <- matrix(rep(map.trees$y, each = nrow(map.grid)), nrow = nrow(map.grid))  # matrix of tree y-coords
  map.dist <- sqrt((grid.x - maptrees.x)^2 + (grid.y - maptrees.y)^2)    # grid cell x tree distance matrix
  outside <- apply(sweep(map.dist, 2, crown.diam[trees$site==site], ">"), 1, all)  # mask of cells outside any canopy
  
  map.params <- list(fec=matrix(1, n.sims, sum(tree.site==j)), 
                     shape=exp(as.matrix(jags.fit$lshape.grp)[1:n.sims,1]), 
                     scale=exp(as.matrix(jags.fit$lscale.grp)[1:n.sims,isl.grp]))
  seeds.surf <- k2d.predict(params=map.params, kernel=kern, traps=map.grid, trees=map.trees, trap.area=1)$total.seeds
  
  p.outside <- rowSums(sweep(seeds.surf, 2, outside, "*"))/rowSums(seeds.surf)
  
  seed.canopy.tab$p.outside[j] <- mean(p.outside)
  seed.canopy.tab$CI.025[j] <- quantile(p.outside, 0.025)
  seed.canopy.tab$CI.975[j] <- quantile(p.outside, 0.975)
  
  hist(p.outside, 50, prob=T, xlab="Proportion outside canopy", main=site) 
  
  rm(list=c("site","isl.grp","map.traps","map.trees","map.grid","grid.x","grid.y","outside",
            "maptrees.x","maptrees.y","map.dist","map.params","seeds.surf","p.outside"))
}

# write results
seed.canopy.tab <- seed.canopy.tab[order(seed.canopy.tab$island),]
write.table(seed.canopy.tab, paste("seed.canopy.tab", species, "txt", sep="."), sep="\t", row.names=F)

rm(list=c("kern","jags.fit","npts","n.sims"))


#----------------------------------------------------------------------------------------
# Calculate the proportion of mapped forest floor at each site that is predicted
# to receive at least 1 seed per year (accounts for NB distribution)
#----------------------------------------------------------------------------------------

seed.rain.tab <- data.frame(island=levels(trees$island)[site.island], site=levels(trees$site),
                            p.1plus=NA, CI.025=NA, CI.975=NA)

kern <- mod.tab$kernel[1]
jags.fit <- get(paste("jags", species, kern, "G", sep="."))$BUGS$sims.list
time.adj <- switch(species, PS = 12/4, PR = 12/10)
npts <- 100
n.sims <- 1000
prob.thresh <- 0.5

dev.new(width=14, height=7)
par(mfrow=c(3,5))

for(j in 1:nrow(seed.rain.tab))
{ 
  site <- seed.rain.tab$site[j]
  isl.grp <- ifelse(site.island[j]==1, 1, 2)
  map.traps <- traps[traps$site==seed.rain.tab$site[j],]
  map.trees <- trees[trees$site==seed.rain.tab$site[j],]
  if(species=="PS") names(map.trees)[names(map.trees)=="ba.avg"] <- "ba"
  
  map.grid <- expand.grid(x=seq(min(c(map.traps$x, map.trees$x)), max(c(map.traps$x, map.trees$x)), length=npts),
                          y=seq(min(c(map.traps$y, map.trees$y)), max(c(map.traps$y, map.trees$y)), length=npts))
  dx <- diff(sort(unique(map.grid$x)))[1] # grid cell size (for converting fitted density to count)
  dy <- diff(sort(unique(map.grid$y)))[1]
  grid.x <- matrix(rep(map.grid$x, nrow(map.trees)), nrow = nrow(map.grid))  # matrix of grid x-coords
  grid.y <- matrix(rep(map.grid$y, nrow(map.trees)), nrow = nrow(map.grid))  # matrix of grid y-coords
  maptrees.x <- matrix(rep(map.trees$x, each = nrow(map.grid)), nrow = nrow(map.grid))  # matrix of tree x-coords
  maptrees.y <- matrix(rep(map.trees$y, each = nrow(map.grid)), nrow = nrow(map.grid))  # matrix of tree y-coords
  map.dist <- sqrt((grid.x - maptrees.x)^2 + (grid.y - maptrees.y)^2)    # grid cell x tree distance matrix
  outs <- apply(map.dist, 1, min) > 20 # throw out grid cells that are >20 m from the nearest tree  
  map.grid <- data.frame(site=site, map.grid[!outs,])
  
  map.params <- list(fec=time.adj*exp(jags.fit$lfec[1:n.sims,tree.site[tree.site==j]]), 
                     shape=exp(as.matrix(jags.fit$lshape.grp)[1:n.sims,1]), 
                     scale=exp(as.matrix(jags.fit$lscale.grp)[1:n.sims,isl.grp]))
  seeds.surf <- k2d.predict(params=map.params, kernel=kern, traps=map.grid, trees=map.trees, trap.area=dx*dy)$total.seeds
  
  p.1plus <- vector("numeric", n.sims)
  
  for(i in 1:nrow(seeds.surf))
  {
    nb.p <- 1/(jags.fit$r1[i] + jags.fit$r2[i]*seeds.surf[i,])
    nb.r <- seeds.surf[i,]/(jags.fit$r1[i] - 1 + jags.fit$r2[i]*seeds.surf[i,])
    p.1plus[i] <- mean(1 - pnbinom(0, size = nb.r, prob = nb.p) >= prob.thresh)
  }
  
  seed.rain.tab$p.1plus[j] <- mean(p.1plus)
  seed.rain.tab$CI.025[j] <- quantile(p.1plus, 0.025)
  seed.rain.tab$CI.975[j] <- quantile(p.1plus, 0.975)
  
  hist(p.1plus, 50, prob=T, xlab=paste("Proportion with P(N > 0) >", prob.thresh), main=site) 
  
  rm(list=c("site","isl.grp","map.traps","map.trees","map.grid","dx","dy","grid.x","grid.y",
            "maptrees.x","maptrees.y","map.dist","outs","map.params","seeds.surf","p.1plus",
            "nb.p","nb.r"))
}

# write results
seed.rain.tab <- seed.rain.tab[order(seed.rain.tab$island),]
write.table(seed.rain.tab, paste("seed.rain.tab", species, "txt", sep="."), sep="\t", row.names=F)

rm(list=c("kern","jags.fit","time.adj","npts","n.sims","prob.thresh"))



#==================================================================
# FIGURES FOR PAPER
#==================================================================

#----------------------------------------------------------------------------------------
# Six-panel figure combining dispersal kernels and filled contour plots of seed density
#----------------------------------------------------------------------------------------

rm(list=ls())  ## Are you sure you want to do this?
spp <- c("PS","PR")

# png(filename=paste("kernels_and_seed_density", "png", sep="."), units="in",
#     width=5.6 + 5*2 + 0.5*2 + 5*0.25, height=5*2 + 0.5*2, res=300, pointsize=18, type="cairo-png")
dev.new(width=5.6 + 5*2 + 0.5*2 + 5*0.25, height=5*2 + 0.5*2, units="in")
layout(matrix(c(0,0,0,1,2,0,
                3,4,5,6,7,8,
                9,10,11,12,13,14,
                0,15,0,16,17,0), nrow=4, byrow=T), 
       widths = c(0.1,1.12,0.1,1,1,0.25), 
       heights = c(0.1,0.95,0.95,0.2))
par(mar=c(0,2.1,0,1))
plot(1, 1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
text(1, 1, "Islands with birds", cex=1.8/par("cex"), font=4)
par(mar=c(0,1,0,2.1))
plot(1, 1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
text(1, 1, "Guam", cex=1.8/par("cex"), font=4)

for(i in spp)
{
  load(paste("seedDispersal", i, "_prior_scale_LN.RData", sep=""))
  source("k2d_radial_quantiles.r")
  source("minor.ticks.axis.r")
  source("filled.contour.plt.r")
  source("filled.contour.key.r")

  par(mar=c(1.1,0,1,0))
  plot(1, 1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  text(0.95, 1, switch(species, PS="Psychotria", PR="Premna"), srt = 90, cex=1.8/par("cex"), font=4)
  abline(v=par("usr")[2], lwd=4)
  
  #----------------------------------------------------------------------------------------
  # Dispersal kernel (relative density per m2) for islands with vs. without birds
  #----------------------------------------------------------------------------------------
  
  kern <- mod.tab$kernel[1]
  kfunc <- switch(as.character(kern), pexp=kpexp, k2Dt)
  rmax <- 5
  
  jags.fit <- get(paste("jags", species, kern, "G", sep="."))$BUGS
  r <- seq(0,rmax,length=200)
  rr <- matrix(r, nrow=jags.fit$n.sims, ncol=length(r), byrow=T)
  scale.G <- matrix(exp(jags.fit$sims.list$lscale.grp[,1]),
                    nrow=jags.fit$n.sims, ncol=length(r))
  scale.RST <- matrix(exp(jags.fit$sims.list$lscale.grp[,2]),
                      nrow=jags.fit$n.sims, ncol=length(r))
  shape <- matrix(exp(jags.fit$sims.list$lshape.grp[,1]),
                  nrow=jags.fit$n.sims, ncol=length(r))
  dr.G <- log10(kfunc(rr, shape, scale.G))
  dr.RST <- log10(kfunc(rr, shape, scale.RST))
  cG <- col2rgb("#009E73")
  cG <- rgb(red=cG[1], green=cG[2], blue=cG[3],alpha=0.3*255, maxColorValue=255)
  cRST <- col2rgb("#E69F00")
  cRST <- rgb(red=cRST[1], green=cRST[2], blue=cRST[3],alpha=0.4*255, maxColorValue=255)
  
  par(mar = c(1.1,6,1,0))
  plot(r, colMeans(dr.G), type="l", lwd=3, col="#009E73", cex.axis=1.5, cex.lab=1.8, las=1, font.lab=2,
       xaxs="i", xaxt="n", yaxt="n", ylim=c(min(apply(dr.G,2,quantile,0.025)), max(apply(dr.G,2,quantile,0.975))), 
       xlab="",  ylab=expression(bold("Relative seed density ( " * m^-2 * ")")))
  polygon(c(r,rev(r)), c(apply(dr.G,2,quantile,0.025), rev(apply(dr.G,2,quantile,0.975))), col=cG, border=NA)
  lines(r, colMeans(dr.RST), lwd=3, col="#E69F00")
  polygon(c(r,rev(r)), c(apply(dr.RST,2,quantile,0.025), rev(apply(dr.RST,2,quantile,0.975))), col=cRST, border=NA)
  axis(1, labels=species=="PR", cex.axis=1.5)
  at <- pretty(colMeans(dr.G))
  axis(2, at, sapply(at,function(i) as.expression(bquote(10^ .(i)))), cex.axis=1.5, las=1)
  # minor.ticks.axis(2,9,0.5,las=1,cex.axis=1.5)
  if(species=="PS")
    legend("topright", c("Islands with birds","Guam"), lwd=3, pch="", col=c("#E69F00","#009E73"), text.font=2, cex=1.5)
  legend("topleft", "", title=switch(species, PS="a", PR="b"), cex=1.8, text.font=2, bty="n")

  #-----------------------------------------------------------------------
  # Filled contour plot of seed density for a selected site on Guam,
  # using "with bird" and "without bird" kernel parameters
  #-----------------------------------------------------------------------
  
  par(mar=c(1.1,0,1,0))
  plot(1, 1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  text(1, 1, "Location (m)", srt=90, cex=1.8, font=2)
  
  ### user inputs
  site <- which(levels(traps$site)=="RTN")
  #   kern <- mod.tab$kernel[1]
  npts <- 50
  ###
  
  #   jags.fit <- get(paste("jags", species, kern, "G", sep="."))$BUGS
  time.adj <- switch(species, PS = 12/4, PR = 12/10)
  
  map.trees <- trees[as.numeric(trees$site)==site,]
  map.trees$x <- map.trees$x - min(map.trees$x) + 5
  map.trees$y <- map.trees$y - min(map.trees$y) + 5
  if(species=="PS") 
  {
    indx <- map.trees$x <= 50 & map.trees$y <= 50
    names(map.trees)[names(map.trees)=="ba.avg"] <- "ba"
  }
  if(species=="PR") 
    indx <- map.trees$x > 3 & map.trees$x < 53 & map.trees$y > 25 & map.trees$y < 75
  map.trees <- map.trees[indx,]
  map.trees$x <- map.trees$x - min(map.trees$x) + 5
  map.trees$y <- map.trees$y - min(map.trees$y) + 5
  cr <- crown.diam[as.numeric(trees$site)==site]/2
  cr <- cr[indx]
  # if(species=="PS") names(map.trees)[names(map.trees)=="ba.avg"] <- "ba"
  # map.grid <- expand.grid(x=seq(min(map.trees$x) - 5, max(map.trees$x) + 5, length=npts),
  #                         y=seq(min(map.trees$y) - 5, max(map.trees$y) + 5, length=npts))
  map.grid <- expand.grid(x=seq(0, 55, length=npts),
                          y=seq(0, 55, length=npts))
  map.grid <- data.frame(site=map.trees$site[1], map.grid)
  
  map.paramsG <- list(fec=time.adj*exp(jags.fit$sims.list$lfec[1:500,tree.site[tree.site==site]]), 
                      shape=exp(as.matrix(jags.fit$sims.list$lshape.grp)[1:500,1]), 
                      scale=exp(as.matrix(jags.fit$sims.list$lscale.grp)[1:500,1]))
  seeds.surfG <- k2d.predict(params=map.paramsG, kernel=kern, traps=map.grid, trees=map.trees, trap.area=1)$total.seeds
  seeds.surfG <- colMeans(log10(seeds.surfG))
  
  map.paramsRST <- list(fec=time.adj*exp(jags.fit$sims.list$lfec[1:500,tree.site[tree.site==site]]), 
                        shape=exp(as.matrix(jags.fit$sims.list$lshape.grp)[1:500,2]), 
                        scale=exp(as.matrix(jags.fit$sims.list$lscale.grp)[1:500,2]))
  seeds.surfRST <- k2d.predict(params=map.paramsRST, kernel=kern, traps=map.grid, trees=map.trees, trap.area=1)$total.seeds
  seeds.surfRST <- colMeans(log10(seeds.surfRST))
  
  map.trees$x <- map.trees$x - min(map.grid$x)
  map.trees$y <- map.trees$y - min(map.grid$y)
  map.grid$x <- map.grid$x - min(map.grid$x)
  map.grid$y <- map.grid$y - min(map.grid$y)
  
  # Draw the first plot for the kernel with birds
  par(mar = c(1.1,2.1,1,1))
  filled.contour.plt(sort(unique(map.grid$x)), sort(unique(map.grid$y)), matrix(seeds.surfRST,npts,npts),
                     zlim=range(c(seeds.surfG, seeds.surfRST)), las=1,
                     #                      nlevels=50, color.palette=colorRampPalette(c("white","tan4","green")),
                     nlevels=50, color.palette=colorRampPalette(brewer.pal(9,"YlGn")),
                     #                      plot.title=title(main=expression("(A) Seed density with birds (" * m^-2 * ")"), 
                     xlab="", ylab="",
                     plot.axes={axis(1, labels=species=="PR", cex.axis=1.5); axis(2, cex.axis=1.5); 
                       symbols(map.trees$x, map.trees$y, circles=cr, inches=F, add=T)})
  legend("topleft", "", title=switch(species, PS="c", PR="d"), cex=1.8, text.font=2, bty="n")
  
  # Draw the second plot for the kernel without birds
  par(mar = c(1.1,1,1,2.1))
  filled.contour.plt(sort(unique(map.grid$x)), sort(unique(map.grid$y)), matrix(seeds.surfG,npts,npts),
                     zlim=range(c(seeds.surfG, seeds.surfRST)), las=1,
                     #                      nlevels=50, color.palette=colorRampPalette(c("white","tan4","green")),
                     nlevels=50, color.palette=colorRampPalette(brewer.pal(9,"YlGn")),
                     #                      plot.title=title(main=expression("(B) Seed density without birds (" * m^-2 * ")"), 
                     xlab="", ylab="",
                     plot.axes={axis(1, labels=species=="PR", cex.axis=1.5); axis(2, labels=F, cex.axis=1.5); 
                       symbols(map.trees$x, map.trees$y, circles=cr, inches=F, add=T)})
  legend("topleft", "", title=switch(species, PS="e", PR="f"), cex=1.8, text.font=2, bty="n")
  
  # Draw the key
  par(mar = c(1.1,0,1,4.1))
  key.lab <- pretty(range(c(seeds.surfG, seeds.surfRST)))
  filled.contour.key(sort(unique(map.grid$x)), sort(unique(map.grid$y)), matrix(seeds.surfRST,npts,npts),
                     zlim=range(c(seeds.surfG, seeds.surfRST)),
                     #                      nlevels=50, color.palette=colorRampPalette(c("white","tan4","green")),
                     nlevels=50, color.palette=colorRampPalette(brewer.pal(9,"YlGn")), border = NA,
                     key.axes=axis(4, key.lab, sapply(key.lab,function(i) as.expression(bquote(10^ .(i)))), cex.axis=1.5))
  
  rm(list=ls())
}
par(mar=c(0,6,0,0))
plot(1, 1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
text(1, 1, "Dispersal distance (m)", cex=1.8, font=2)
par(mar=c(0,2.1,0,1))
plot(1, 1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
text(1, 1, "Location (m)", cex=1.8, font=2)
par(mar=c(0,2.1,0,1))
plot(1, 1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
text(1, 1, "Location (m)", cex=1.8, font=2)
# dev.off()


#---------------------------------------------------------------------------------------------------
# Alternative version: Marginal PDF of radial displacement for islands with vs. without birds
#---------------------------------------------------------------------------------------------------

### user inputs
kern <- mod.tab$kernel[1]
rmax <- 20
###

jags.fit <- get(paste("jags", species, kern, "G", sep="."))$BUGS
r <- seq(1e-2,rmax,by=0.1)
rr <- matrix(r, nrow=jags.fit$n.sims, ncol=length(r), byrow=T)
scale.G <- matrix(exp(jags.fit$sims.list$lscale.grp[,1]),
                  nrow=jags.fit$n.sims, ncol=length(r))
scale.RST <- matrix(exp(jags.fit$sims.list$lscale.grp[,2]),
                    nrow=jags.fit$n.sims, ncol=length(r))
shape <- matrix(exp(jags.fit$sims.list$lshape.grp[,1]),
                nrow=jags.fit$n.sims, ncol=length(r))
dr.G <- log10(dpexp(rr, shape, scale.G))
dr.RST <- log10(dpexp(rr, shape, scale.RST))
cG <- col2rgb("green")
cG <- rgb(red=cG[1], green=cG[2], blue=cG[3],alpha=0.4*255, maxColorValue=255)
cRST <- col2rgb("blue")
cRST <- rgb(red=cRST[1], green=cRST[2], blue=cRST[3],alpha=0.4*255, maxColorValue=255)

# png(filename=paste("marg_pdf", species, "png", sep="."), width=7, height=7, units="in", res=200, type="cairo-png")
dev.new()
plot(r, colMeans(dr.G), type="l", lwd=2, col="green", cex.lab=1.5, las=1, yaxt="n",
     ylim=c(min(apply(dr.G,2,quantile,0.025)),max(apply(dr.G,2,quantile,0.975))), 
     xlab="Dispersal distance (m)",  ylab="Probability density", 
     main=switch(species, PR="Premna", PS="Psychotria"))
polygon(c(r,rev(r)), c(apply(dr.G,2,quantile,0.025), rev(apply(dr.G,2,quantile,0.975))), col=cG, border=NA)
lines(r, colMeans(dr.RST), lwd=2, col="blue")
polygon(c(r,rev(r)), c(apply(dr.RST,2,quantile,0.025), rev(apply(dr.RST,2,quantile,0.975))), col=cRST, border=NA)
minor.ticks.axis(2,9,0.5,las=1)
# dev.off()

rm(list=c("kern","jags.fit","r","rr","scale","shape","dr.G","dr.RST","cG","cRST"))

