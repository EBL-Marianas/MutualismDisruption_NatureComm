#------------------------------------------------------------------
# Functions to find quantiles of marginal distribution of
# radial displacement under power-exponential or 2Dt kernel
#------------------------------------------------------------------

###############################
## POWER-EXPONENTIAL KERNELS ##
###############################

# 2D PDF (relative seed density)
kpexp <- function(r, shape, scale)
{
  (shape/(2*pi*gamma(2/shape)*scale^2))*exp(-(r/scale)^shape)
}

# marginal PDF of displacement
dpexp <- function(r, shape, scale)
{
  (r*shape/(gamma(2/shape)*scale^2))*exp(-(r/scale)^shape)
}

# marginal CDF of displacement
ppexp <- function(r, shape, scale)
{
  integrate(dpexp, lower = 0, upper = r, shape = shape, scale = scale)
}


#################
## 2Dt KERNELS ##
#################

# 2D PDF (relative seed density)
k2Dt <- function(r, shape, scale)
{
  shape/(pi*scale^2*(1 + (r/scale)^2)^(shape + 1))
}

# marginal PDF of displacement
d2Dt <- function(r, shape, scale)
{
  2*r*shape/(scale^2*(1 + (r/scale)^2)^(shape + 1))
}

# marginal CDF of displacement
p2Dt <- function(r, shape, scale)
{
  list(value = 1 - (scale^2/(scale^2 + r^2))^shape)
}


