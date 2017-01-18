#------------------------------------------------------------------------------------------------------
# Function to calculate fitted values (seeds per trap), given
# the data and one or more vectors of parameters

# Params is a list with elements named fec, shape, and scale, defined as appropriate 
# for the specified kernel, and (optionally) the negative binomial overdispersion parameter r.
# fec is an n.params x n.trees matrix with fecundity per BA for each tree.
# Each element describing the kernel params is either an n.params x n.traps matrix
# or a vector of length n.params.
# traps is a data frame with columns named x, y, and site, and (optionally) total.
# trees is a data frame with columns named x, y, site and ba (basal area).
# trap.area is a scalar.
# The fitted (posterior expectation) total seeds per trap (total.seeds) is returned 
# as a matrix (n.params x n.traps).
#------------------------------------------------------------------------------------------------------

k2d.predict <- function(params, kernel, traps, trees, trap.area)
{  
  # object dims
  n.traps <- nrow(traps)
  n.trees <- nrow(trees)
  if(is.matrix(params$shape)) n.params <- nrow(params$shape) else n.params <- length(params$shape)
  total.seeds <- matrix(0, nrow=n.params, ncol=n.traps)
  shape <- matrix(params$shape, nrow=n.params, ncol=n.traps)
  scale <- matrix(params$scale, nrow=n.params, ncol=n.traps)
  if(any(names(trees)=="ba")) ba <- trees$ba else ba <- trees$ba.avg
  
  # Distance matrix (traps x trees)  
  x <- expand.grid(traps.x=traps$x, trees.x=trees$x)  # matrices of trap x tree x-coords
  y <- expand.grid(traps.y=traps$y, trees.y=trees$y)  # matrices of trap x tree y-coords
  dist <- matrix(sqrt((x$traps.x - x$trees.x)^2 + (y$traps.y - y$trees.y)^2), 
                 nrow=n.traps, ncol=n.trees)          # trap x tree distance matrix
  #ba <- matrix(trees$ba, n.traps, n.trees, byrow=T)   # trap x tree matrix whose rows are ba by tree
  
  # Calculate fitted values under specified kernel, sum across trees
  for(i in 1:n.trees)
  {
    dist.i <- matrix(dist[,i], n.params, n.traps, byrow=T) 
    
    if(kernel=="2Dt")
    {
      fec <- matrix(params$fec[,i], nrow=n.params, ncol=n.traps)
      fit.seeds <- trap.area*fec*ba[i]*shape/(pi*scale^2*(1 + (dist.i/scale)^2)^(shape + 1))
      if(length(unique(traps$site)) > 1)
        fit.seeds[, traps$site != trees$site[i]] <- 0  # no contributions from trees at other sites
      total.seeds <- total.seeds + fit.seeds
    }
    
    else if(kernel=="pexp")
    {
      fec <- matrix(params$fec[,i], nrow=n.params, ncol=n.traps)
      fit.seeds <- trap.area*fec*ba[i]*(shape/(2*pi*gamma(2/shape)*scale^2))*exp(-(dist.i/scale)^shape)
      if(length(unique(traps$site)) > 1)
        fit.seeds[, traps$site != trees$site[i]] <- 0  # no contributions from trees at other sites
      total.seeds <- total.seeds + fit.seeds
    }
    
    else {cat("Error: unknown dispersal kernel specified"); break}
  }
    
  return(list(total.seeds=total.seeds))
}
