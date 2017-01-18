seeds.fit2sim <- function(jags.fit, s)
{
  if(is.null(s))  s <- sample(jags.fit$BUGS$n.sims, 1)
  r1 <- jags.fit$BUGS$sims.list$r1[s]
  r2 <- jags.fit$BUGS$sims.list$r2[s]
  fit <- jags.fit$BUGS$sims.list$fit[s,]
  p <- 1/(r1 + r2*fit)
  r <- fit/(r1 - 1 + r2*fit)
  seeds.sim <- rnbinom(length(fit), prob=p, size=r)
    return(seeds.sim)
  #   sd <- sqrt(r*(1 - p))/p
  #   cv <- sd/fit
  #   return(list(seeds.sim=seeds.sim,fit=fit,r=r,cv=cv))
}
