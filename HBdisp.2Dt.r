#===========================================================================================
# 2Dt HIERARCHICAL DISPERSAL MODEL
#===========================================================================================

# Kernel is 2Dt parameterized as (fec, scale, shape), where shape is fixed at a global value.
# 
# Seed density per m2 at a distance of dist in any direction is then
# 
#     dens = fec*BA*shape/(pi*scale^2*(1 + (dist/scale)^2)^(shape + 1)).
# 
# (Note that scale as defined here corresponds to sqrt(scale) in Clark et al. 1999.)
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
#     log(dmode) = log(scale) - 0.5*log(2*shape + 1).
#     log(dmean) = log(scale) + 0.5*log(pi) + log(gamma(shape - 1/2)) - log(2) - log(gamma(shape))
# 
# The likelihood of observed seed counts is negative binomial, parameterized following 
# Linden and Mantyniemi (Ecology, 2011). Briefly, the variance scales with the mean as a polynomial
# of degree 2, with linear and quadratic coefficients r1 and r2, respectively.

data
{
  pi <- 3.1416
  distsq <- dist^2                       # square of distance to source  
  shape.fix <- 1                         # fixed global kernel shape

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
  sigma.lfec ~ dunif(0,10)               # hyper-SD of log fecundity per basal area
  
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
    ldmode.grp[i] <- lscale.grp[i]/sqrt(2*shape.grp[i] + 1)
    
    # log mean displacement
    ldmean.grp[i] <- lscale.grp[i] + 0.5*log(pi) + loggam(shape.grp[i] - 0.5) - log(2) - loggam(shape.grp[i])
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
    # shape and scale for the island where trap i is located
    shape[i] <- shape.grp[kern.indx[site.island[trap.site[i]]]]
    scale[i] <- scale.grp[kern.indx[site.island[trap.site[i]]]]
    
    fit[i] <- sum(trap.area*fec[trap.site[i]]*basal.area[cc[i,1]:cc[i,2]]*shape[i]/(pi*scale[i]^2*pow(1 + distsq[i,cc[i,1]:cc[i,2]]/scale[i]^2, shape[i] + 1)))
    #     fit[i] <- sum(trap.area*fec[trap.site[i]]*shape[i]/(pi*scale[i]^2*pow(1 + distsq[i,cc[i,1]:cc[i,2]]/scale[i]^2, shape[i] + 1)))
    p[i] <- 1/(r1 + r2*fit[i])                             # p parameter for negative binomial
    r[i] <- fit[i]/(r1 - 1 + r2*fit[i])                    # r parameter for negative binomial
    ones[i] ~ dbern(r[i] > 0)                              # ensure r > 0
    seeds[i] ~ dnegbin(p[i], max(r[i],1e-10))              # negative binomial likelihood of count
    ll.data[i] <- log(dnegbin(seeds[i], p[i], max(r[i],1e-10)))
  }
  
  dev.data <- -2*sum(ll.data)
}
