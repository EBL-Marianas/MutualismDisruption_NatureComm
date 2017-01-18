#===========================================================================================
# POWER-EXPONENTIAL HIERARCHICAL DISPERSAL MODEL
#===========================================================================================

# Kernel is power-exponential parameterized as (fec, scale, shape), where shape is fixed at 
# a global value.
# 
# Seed density per m2 at a distance of dist in any direction is then
# 
#     dens = (fec*BA*shape/(2*pi*gamma(2/shape)*scale^2))*exp(-(dist/scale)^shape).
# 
# Fecundity per BA varies randomly among sites, so for all individual trees at site i
# 
#     log(fec[i]) ~ N(mu.lfec, sigma.lfec),
# 
# while scale is defined at the island level. Scale may vary among islands or groups of islands 
# according to the input vector kern.indx.  Log-mean fecundity does not vary across islands.
# 
# Additional summaries of the marginal distribution of radial displacement (mode and mean) 
# are calculated as
# 
#     log(dmode) = log(scale) + (-1/shape)*log(shape).
#     log(dmean) = log(scale + log(gamma(3/shape)) - gamma(2/shape)
# 
# (Note that the quantiles are not available in closed form, so they must be found numerically
#  using the MCMC output.)
# 
# The likelihood of observed seed counts is negative binomial, parameterized following 
# Linden and Mantyniemi (Ecology, 2011). Briefly, the variance scales with the mean as a polynomial
# of degree 2, with linear and quadratic coefficients r1 and r2, respectively.

data
{
  pi <- 3.1416
  shape.fix <- 0.5                       # fixed global kernel shape
  
  # vector of ones to ensure linear and quadratic NB variance coefs (sensu Linden and Mantyniemi 2011)
  for(i in length(seeds))
  {
    ones[i] <- 1
  }
}

model
{
  # Priors
  mu.lfec ~ dnorm(0,1e-6)                # hyper-mean(s) of log fecundity per basal area
  sigma.lfec ~ dunif(0,20)               # hyper-SD of log fecundity
  
  r1 ~ dunif(0,100)                      # linear NB variance coef
  r2 ~ dunif(0,100)                      # quadratic NB variance coef
  
  # Note that lscale may be a vector, so we loop over its elements
  # (i.e., groups of islands that differ in dispersal scale)
  for(i in 1:max(kern.indx))  
  {
    lscale.grp[i] ~ dnorm(0,1e-6)
    scale.grp[i] <- exp(lscale.grp[i])
    #     scale.grp[i] ~ dunif(0,1e6)
    #     lscale.grp[i] <- log(scale.grp[i])
    shape.grp[i] <- shape.fix
    lshape.grp[i] <- log(shape.fix)
    
    # log modal displacement
    ldmode.grp[i] <- lscale.grp[i] - lshape.grp[i]/shape.grp[i]
    
    # log mean displacement
    ldmean.grp[i] <- lscale.grp[i] + loggam(3/shape.grp[i]) - loggam(2/shape.grp[i])
  }
  
  # Site-specific fecundity
  for(j in 1:max(trap.site))
  {
    lfec[j] ~ dnorm(mu.lfec, pow(sigma.lfec,-2))
    ll.lfec[j] <- log(dnorm(lfec[j], mu.lfec, pow(sigma.lfec,-2)))
    fec[j] <- exp(lfec[j])
  }
  
  dev.lfec <- -2*sum(ll.lfec)            # deviance of log fecundity random effects
  
  # Calculate expected number of seeds contributed to each trap by each tree,
  # including only trap-tree combinations at the same site, then
  # add up seed contributions from all trees for each trap and evaluate likelihood
  for(i in 1:length(seeds))
  {
    # shape and scale for the site or island where trap i is located
    shape[i] <- shape.grp[kern.indx[site.island[trap.site[i]]]]
    scale[i] <- scale.grp[kern.indx[site.island[trap.site[i]]]]
    
    # necessary to loop over trees because exp function is not vectorized in JAGS
    for(j in cc[i,1]:cc[i,2])
    {
      fitmat[i,j] <- trap.area*fec[trap.site[i]]*basal.area[j]*shape[i]/(2*pi*exp(loggam(2/shape[i]))*scale[i]^2)*exp(-(dist[i,j]/scale[i])^shape[i])
#       fitmat[i,j] <- trap.area*fec[trap.site[i]]*shape[i]/(2*pi*exp(loggam(2/shape[i]))*scale[i]^2)*exp(-(dist[i,j]/scale[i])^shape[i])
    } 
    
    fit[i] <- sum(fitmat[i,cc[i,1]:cc[i,2]])
    p[i] <- 1/(r1 + r2*fit[i])                             # p parameter for negative binomial
    r[i] <- fit[i]/(r1 - 1 + r2*fit[i])                    # r parameter for negative binomial
    ones[i] ~ dbern(r[i] > 0)                              # ensure r > 0
    seeds[i] ~ dnegbin(p[i], max(r[i],1e-10))              # negative binomial likelihood of count
    ll.data[i] <- log(dnegbin(seeds[i], p[i], max(r[i],1e-10)))
  }
  
  dev.data <- -2*sum(ll.data)
}
