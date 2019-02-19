# Conditional Geostatistical Simulations

# Haddock, Division VIa, 2012 Quarter 1

# Script adapted from:
# Petitgas, P., Woillez, M., Rivoirard, J., Renard, D., and Bez, N. 2017. Handbook of 
# geostatistics in R for fisheries and marine ecology. ICES Cooperative Research Report
# No. 338. 177 pp.

##############################
### Load Packages and Data ###
##############################

library(RGeostats) # RGeostats_11.2.1
library(spatstat)  # spatstat_1.55-1

viapoly <- read.table("VIa_Polygon.txt", header=T)
dat <- read.table("haddock_example_2012_Q1.txt", header=T)

# EEZ polygons
out6a <- read.table("out6a.txt", header=T)
in6a <- read.table("in6a.txt", header=T)

#########################################
### Database management for RGeoStats ###
#########################################

# Create cod database
hadDb <- db.create(dat, flag.grid=F)

# Set coordinate locators
hadDb = db.locate(hadDb, "shootlong", "x",1)
hadDb = db.locate(hadDb, "shootlat", "x",2)

# Set swept-area density as simulation variable
hadDb = db.locate(hadDb, "density.wt","z")

# VIa area polygon
poly.data <- polygon.create(viapoly)
in6a <- polygon.create(in6a)
out6a <- polygon.create(out6a)

# Parameter for variogram transformation
n.H <- 50

# Approximate grid cell sizes at reference coordinates
projec.define(projection="mean", db=hadDb ,flag.update=T) 
mapCoords <- projec.operate(viapoly$lon, viapoly$lat)

# In nautical miles
nmix <- sum(abs(range(mapCoords$x)))
nmiy <- sum(abs(range(mapCoords$y)))

# Calculate number of km in x and y of area polygon; set as grid number
gnx <- round(nmix*1.852)
gny <- round(nmiy*1.852)

# Number of simulations to run
ns <- 500

###################
### Variography ###
###################

# Turn off any projections
projec.toggle(0)

# Plot dB
plot(hadDb, zmin=0.001, pch.low=3, cex.low=0.25, las=1,
     pch=21, col=1, inches=5, title="Haddock - 2012 Quarter 1", asp=1)

# Calculate influence surface for sample points
hadDb <- infl(hadDb, nodes=c(400,400), 
              origin=c(min(viapoly$lon),min(viapoly$lat)), 
              extend=c(sum(abs(range(viapoly$lon))),
                       max(viapoly$lat)-min(viapoly$lat)), 
              polygon=poly.data, plot=T, asp=1)

# Project database
projec.define(projection="mean", db=hadDb ,flag.update=T) 

# Plot projected data
plot(hadDb, zmin=0.001, pch.low=3, cex.low=0.25, las=1,
     pch=21, col=1, inches=5, 
     title="Haddock - 2012 Quarter 1 (Mean Projection)", asp=1)

# Compute lag & maximum distance
pvals <- projec.operate(hadDb@items$shootlong,
                        hadDb@items$shootlat)
llag <- mean(nndist(cbind(pvals$x,pvals$y)))
nlag <- round((max(pairdist(cbind(pvals$x,pvals$y)))/2)/llag)

# Plot raw data histogram
hist(hadDb@items$density.wt, main="Raw Swept-Area Density Data",
     xlab="", breaks=50, col="#7f7f7f", border="#7f7f7f")
mean(hadDb@items$density.wt)
var(hadDb@items$density.wt)

# Fit model anamorphosis
model.anam15 <- anam.fit(hadDb, type="emp",
                         draw=T, title="Anamorphosis")

# Create transformed database
db.data.trans <- anam.z2y(hadDb, anam=model.anam15)
db.data.trans <- db.rename(db.data.trans, name="Gaussian.density.wt", 
                           newname="Yp")
ycut3 <- qnorm(sum(db.extract(db.data.trans,"density.wt") == 0) / 
                 length(db.extract(db.data.trans, name="density.wt",
                                   flag.compress=T)))
Y3<- db.extract(db.data.trans,"Yp",flag.compress=T)

# Plot anamorphosis transformed data
hist(Y3, breaks = 50, main = "Anamorphosis Transformed", 
     col="#7f7f7f", border="#7f7f7f")

dens <- db.extract(db.data.trans,"density.wt")
Y3[dens == 0] <- ycut3
db.data.trans <- db.replace(db.data.trans,"Yp",Y3)

# Fitting variogram models
vario.Yp.15 <- vario.calc(db.data.trans, lag=llag, nlag=nlag)
vario.Y.15  <- vario.trans.cut(vario.Yp.15, ycut3, n.H)
plot(vario.Y.15)
model.vario.Y.15 <- model.auto(vario.Y.15, struc=melem.name(c(1,2)), draw=T)

# Define interval limits for the Gibbs sampler
Ymax <- db.extract(db.data.trans, name="Yp", flag.compress=F)
Ymin <- db.extract(db.data.trans, name="Yp", flag.compress=F)
Ymin[Ymin <= ycut3] <- -10
db.data.trans <- db.add(db.data.trans,Ymax)
db.data.trans <- db.locate(db.data.trans, db.data.trans$natt,"upper")
db.data.trans <- db.add(db.data.trans, Ymin)
db.data.trans <- db.locate(db.data.trans,db.data.trans$natt,"lower")
db.data.trans <- db.locate(db.data.trans, 9, "z")

# Gibbs sampler
db.data.trans <- gibbs(db = db.data.trans, model = model.vario.Y.15, seed = 1337, 
                       nboot = 1000, niter = 1000, flag.norm=FALSE, 
                       percent=0, toleps = 1, radix = "Gibbs",
                       modify.target = TRUE)
db.data.trans <- db.rename(db.data.trans, "Gibbs.G1", "Y2")

# Histograms to visualise data transformation
histYp <- hist(db.data.trans[,"Yp"],plot=F,breaks=seq(-4,4,.1))
hist(db.data.trans[,"Yp"], proba=T, breaks=seq(-4,4,.1), xlab="Y+", col=8, 
     main="", xlim=c(-4,4), ylim=c(0,ceiling(max(histYp$density)))) 

hist(db.extract(db.data.trans,"Y2",flag.compress=T), proba=T, breaks=seq(-4,4,.1),
     xlab="Y2", main="", col=8,
     xlim=c(-4,4), ylim=c(0,ceiling(max(histYp$density))))
lines(density(db.extract(db.data.trans,"Y2",flag.compress=T)))
lines(seq(-4,4,0.1), dnorm(seq(-4,4,0.1),0,1), col=2) 


###################
### Simulations ###
###################

# Grid for simulation
grid.simu <- db.grid.init(poly.data, nodes=c(gnx,gny))
neigh.simu <- neigh.create(ndim=2,type=2,nmini=1,nmaxi=100)
neigh.simu@dmax <- 100

# simulations
grid.simu <- simtub(dbin=db.data.trans, dbout=grid.simu, model=model.vario.Y.15, 
                    neigh=neigh.simu, uc="", mean=0, seed=1337, nbsimu=ns, 
                    nbtuba=1000, radix="Simu", modify.target=TRUE)

# db for mean realisation
grid.simu.mean <- db.compare(grid.simu, names="Simu*", fun="mean")
grid.simu.mean <- anam.y2z(grid.simu.mean, name="mean", anam=model.anam15)

# plotting mean of simulations
Raw.Simu.mean <- grid.simu.mean@items$Raw.mean
Raw.Simu.mean[round(Raw.Simu.mean,2)<=0.00] <- 0
grid.simu.mean <- db.add(grid.simu.mean, Raw.Simu.mean)
rm(Raw.Simu.mean)
grid.simu.121 <- db.polygon(grid.simu.mean, poly.data)
grid.simu.121 <- db.delete(grid.simu.121, "Simu.Y2*")
rm(grid.simu.mean)

plot(grid.simu.121, name="Raw.Simu.mean", pos.legend=0,
     flag.proj=F, xlim=c(-11,-4))   #(-12,-2)

# Backtransform all realisations
system.time(ysim2012q1 <- anam.y2z(grid.simu, names="Simu*", anam=model.anam15))
ysim2012q1 <- db.delete(ysim2012q1, "Simu.Y2*")
rm(grid.simu)

# Partition by EEZ and get distributions of total weights
# for all of VIa
ysim2012q1 <- db.polygon(ysim2012q1, poly.data)
vals <- ysim2012q1@items[ysim2012q1@items$Polygon==T,]
vals[vals <= 0] <- 0
via121 <- colSums(vals,na.rm = T)
via121 <- via121[grepl("Raw.Simu",names(via121))]
rm(vals)

# Subset the VIa simulation to that which falls outside UK EEZ
ysim2012q1 <- db.polygon(ysim2012q1, out6a, flag.replace = T)
vals <- ysim2012q1@items[ysim2012q1@items$Polygon==T,]
vals[vals <= 0] <- 0
eu121 <- colSums(vals,na.rm = T)
eu121 <- as.numeric(eu121[grepl("Raw.Simu",names(eu121))])
rm(vals)

# Subset the VIa simulation to that which falls within UK EEZ
ysim2012q1 <- db.polygon(ysim2012q1, in6a, flag.replace = T)
vals <- ysim2012q1@items[ysim2012q1@items$Polygon==T,]
vals[vals <= 0] <- 0
uk121 <- colSums(vals,na.rm = T)
uk121 <- as.numeric(uk121[grepl("Raw.Simu",names(uk121))])
rm(vals)

## Calculate areal proportions
boxplot(via121, eu121, uk121,
        border=c("#000000","#376ed3","#ec1336"))

# UK
uk5y <- uk121/via121
mean(uk5y); quantile(uk5y, probs = c(0.025,0.5,0.975))

# EU
eu5y <- eu121/via121
mean(eu5y);quantile(eu5y,probs = c(0.025,0.5,0.975))
