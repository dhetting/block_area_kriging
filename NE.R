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
    ne_plus_neighbors_grid_1km = st_make_grid(ne_plus_neighbors, cellsize=0.01) # 1km
  }
)
# user    system  elapsed 
# 19.030  0.040   19.068
# 23.245   0.128  23.368

#plot(st_geometry(ne_plus_neighbors_grid_1km), add=TRUE)
# doesn't plot correctly, but data are all there

# LK on FPEAM data

x = st_coordinates(st_centroid(ne_plus_neighbors))
y = ne_plus_neighbors[['NH3']]

system.time(
  {
    x_grid = st_coordinates(st_centroid(ne_plus_neighbors_grid_1km))
  }
)
# user  system elapsed 
# 9.384   0.206   9.592
# 11.101   0.249  11.351
# 11.250   0.179  11.437

nlevel=3

system.time(
  {
    LKinfo <- LKrigSetup(x_grid, a.wght=4.5, nlevel=nlevel, nu=1, NC=4, lambda=.01)
  }
)
# user  system elapsed 
# 0.508   0.005   0.514 
# 0.422   0.001   0.424
# 0.412   0.006   0.418

# now observe a linear combination
system.time(
  {
    NNDist <- LKDist(x_grid, x_grid, delta=0.1, distance.type='GreatCircle') 
  }
)
# user  system elapsed 
# 561.307   0.983 562.128

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
# user  system elapsed 
# 22.813   1.099  23.898
#  22.415   0.923  23.332
X <- A%*% PHI
U <- A%*%TMatrix

#yIndirect <- A%*%y
#
# X<-  A

# serial processing

system.time(
  {
    for (i in 1:nrow(ne_plus_neighbors)) {
      print(i)
      county = ne_plus_neighbors[i, ]
      mtrx = st_intersects(ne_plus_neighbors_grid_1km, county, sparse=FALSE)  # @TODO: slow, must not be using spatial index- docs say it creates spatial index on x on the fly
      mtrx_idx = apply(mtrx, 1, any)
#      cells = ne_plus_neighbors_grid_1km[mtrx_idx]
      sums = colSums(U[mtrx_idx, ])
      print(sums)
      county_U[i, ] = sums
    }
  }
)

# user   system  elapsed 
# 1997.546   21.736 2019.635
# 1916.270   18.714 1935.557

ne_plus_neighbors$L1= county_U[, 1]
ne_plus_neighbors$L2= county_U[, 2]
ne_plus_neighbors$L3= county_U[, 3]

par(mfrow=c(1, 3))
plot(ne_plus_neighbors['L1'])
plot(ne_plus_neighbors['L2'])
plot(ne_plus_neighbors['L3'])

# parallelized
numCores <- detectCores()

county_U_mp = matrix(nrow=nrow(ne_plus_neighbors), ncol=nlevel)

fx <- function(i) {
  print(i)
  county = ne_plus_neighbors[i, ]
  mtrx = st_intersects(ne_plus_neighbors_grid_1km, county, sparse=FALSE)  # @TODO: slow, must not be using spatial index- docs say it creates spatial index on x on the fly
  mtrx_idx = apply(mtrx, 1, any)
  sums = colSums(U[mtrx_idx, ])
  return(c(i, sums))
}


run_ids = 1:nrow(ne_plus_neighbors)

system.time(
  {
    results <- mclapply(run_ids, fx, mc.cores=numCores)
    
    for (result in results) {
      i = result[1]
      sums = result[2:4]
      county_U_mp[i, ] = sums
    }
  }
)

# user   system  elapsed 
# 3828.126  321.593  348.205



sum(county_U != county_U_mp)


out0<- LatticeKrig(x,y, LKinfo=LKinfo)
out1<- LatticeKrig(x,yIndirect, U=U, X=X, LKinfo=LKinfo)

# the predict function evaluates f in this case -- not the fitted values that
# are related to the
# observations
# partial agreement but not exact -- this is because the observational models
# assume different errors
#
plot( predict(out0,x), predict( out1,x))

out<- LKrig.sim.conditional( out1,M=5, nx=10, ny=10 )

