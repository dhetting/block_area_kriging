library(sf)
library(LatticeKrig)

# load data
bc1060 <- st_read('~/Box Sync/2016-10-01_fpeam/data/bc1060_county_inmap.shp')

# subset to NE counties
ne = dplyr::filter(bc1060, grepl('^31', cnty_fips))
plot(st_geometry(ne))

# subset to counties not in NE
not_ne = dplyr::filter(bc1060, !grepl('^31', cnty_fips))

# get neighboring counties
ne_neighbors_mtrx = st_touches(not_ne, ne, sparse=FALSE)  # @TODO: slow, must not be using spatial index- docs say it creates spatial index on x on the fly
ne_neighbors_mtrx_idx = apply(ne_neighbors_mtrx, 1, any)
ne_neighbors = not_ne[ne_neighbors_mtrx_idx, ]
plot(st_geometry(ne_neighbors), col='blue', add=TRUE)

# merge neighbors with NE
ne_plus_neighbors = rbind(ne, ne_neighbors)
plot(ne_plus_neighbors['NH3'])

# make grid
system.time(
  {
    bc1060_grid_1km = st_make_grid(bc1060, cellsize=0.01) # 1km
  }
)
# user    system  elapsed 
# 627.457   4.408 632.743

saveRDS(bc1060_grid_1km, file='bc1060.grid_1km.RDS')

#plot(st_geometry(ne_plus_neighbors_grid_1km), add=TRUE)
# doesn't plot correctly, but data are all there

# LK on FPEAM data

system.time(
  {
    x = st_coordinates(st_centroid(bc1060))
  }
)

# user  system elapsed 
# 8.847   0.661  10.284 

y = bc1060[['NH3']]

system.time(
  {
    x_grid = st_coordinates(st_centroid(bc1060_grid_1km))
  }
)
# user  system elapsed 
# 404.139  71.761 521.457

saveRDS(x_grid, file='bc1060.x_grid.RDS')

nlevel=3

system.time(
  {
    LKinfo <- LKrigSetup(x_grid, a.wght=4.5, nlevel=nlevel, nu=1, NC=4, lambda=.01)
  }
)
# user  system elapsed 
# 12.874   2.546  17.962

# now observe a linear combination
system.time(
  {
    NNDist <- LKDist(x_grid, x_grid, delta=0.1, distance.type='GreatCircle') 
  }
)
# user  system elapsed 
# 561.307   0.983 562.128  < for NE + neighbors

saveRDS(NNDist, file='bc1060.NNDist.RDS')

A <- NNDist
A$ra <- exp(-NNDist$ra)
# A is a weight matrix based on neighbors close by and
# in spind sparse matrix format
# now convert to spam format
A <- spind2spam(A)

system.time(
  {
    TMatrix <- get(LKinfo$fixedFunction)(x = x_grid)
  }
)
# user  system elapsed 
# 0.128   0.003   0.130
# Tmatrix is a 3 column matrix of constant and the two spatial coordinates
#  i.e. a linear function of the spatial variables

system.time(
  {
    PHI <- LKrig.basis(x_grid, LKinfo)
  }
)

saveRDS(PHI, file='bc1060.PHI.RDS')

# user  system elapsed 
# 22.813   1.099  23.898
#  22.415   0.923  23.332
X <- A%*% PHI
U <- A%*%TMatrix

saveRDS(X, file='bc1060.X.RDS')
saveRDS(U, file='bc1060.U.RDS')

#yIndirect <- A%*%y
#
# X<-  A

# parallelized
#numCores <- detectCores()
#
#county_U_mp = matrix(nrow=nrow(bc1060), ncol=nlevel)
#
#fx <- function(i) {
#  print(i)
#  # get focus county
#  county = bc1060[i, ]
#  
#  # get 1st order neighbors
#  mtrx_n1 = st_touches(county, bc1060)
#  mtrx_n1_idx = apply(mtrx_n1, 1, any)
#  n1 = bc1060[mtrx_n1_idx, ]
#  
#  # get second order neighbors
#  mtrx_n2 = st_touches(county, bc1060)
#  mtrx_n2_idx = apply(mtrx_n2, 1, any)
#  n2 = bc1060[mtrx_n2_idx, ]
#  
#  n2_grid_1km = st_make_grid(n2, cellsize=0.01) # 1km
#  
#  x_grid = st_coordinates(st_centroid(n2_grid_1km))
#  
#  LKinfo <- LKrigSetup(x_grid, a.wght=4.5, nlevel=nlevel, nu=1, NC=4, lambda=.01)
#  
#  NNDist <- LKDist(x_grid, x_grid, delta=0.1, distance.type='GreatCircle') 
#  
#  TMatrix <- get(LKinfo$fixedFunction)(x = x_grid)
#  
#  PHI <- LKrig.basis(x_grid, LKinfo)
#  
#  A <- NNDist
#  A$ra <- exp(-NNDist$ra)
#  A <- spind2spam(A)
#  
#  X <- A %*% PHI
#  U <- A %*% TMatrix
#  
#  mtrx = st_intersects(n2_grid_1km, county, sparse=FALSE)
#  mtrx_idx = apply(mtrx, 1, any)
#
#  sums = colSums(U[mtrx_idx, ])
#
#  return(c(i, sums))
#}
#
#
#run_ids = 1:nrow(bc1060)
#



# parallelized
numCores <- detectCores()

county_U = matrix(nrow=nrow(bc1060), ncol=nlevel)

fx <- function(i) {
  print(i)
  county = bc1060[i, ]
  mtrx = st_intersects(bc1060_grid_1km, county, sparse=FALSE)  # @TODO: slow, must not be using spatial index- docs say it creates spatial index on x on the fly
  mtrx_idx = apply(mtrx, 1, any)
  sums = colSums(U[mtrx_idx, ])
  return(c(i, sums))
}


run_ids = 1:nrow(bc1060)

system.time(
  {
    results <- mclapply(run_ids, fx, mc.cores=numCores)
    
    for (result in results) {
      i = result[1]
      sums = result[2:4]
      county_U[i, ] = sums
    }
  }
)

saveRDS(county_U, file='bc1060.county_U.RDS')

# user   system  elapsed 
# 3828.126  321.593  348.205  <- ne_neighbors





#out0<- LatticeKrig(x,y, LKinfo=LKinfo)
#out1<- LatticeKrig(x,yIndirect, U=U, X=X, LKinfo=LKinfo)

# the predict function evaluates f in this case -- not the fitted values that
# are related to the
# observations
# partial agreement but not exact -- this is because the observational models
# assume different errors
#
#plot( predict(out0,x), predict( out1,x))

#out<- LKrig.sim.conditional( out1,M=5, nx=10, ny=10 )

