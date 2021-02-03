#
# Illumination Data
# Dec 28, 2020
# Authors: MGA & GAA
# R script to run Longitudinal Data Analytics

## Functions:
get.qi <- function(ti, mn=min(ti), mx=max(ti)) {
  # In order to use Legendre polynomials or other kinds of orthogonal polynomials, the time values 
  # (whole integer numbers) must be scaled to range from -1 to +1. The scaling formula is:
  
  -1 + 2*(ti - mn)/(mx - mn)
}

legendre.poly = function(no) {
  # Make the matrix LAMDA
  if(no > 9 ) no = 9
  nom = no - 1
  phi = matrix(0,nrow=9,ncol=9) 
  phi[1,1]=1
  phi[2,2]=1
  for(i in 2:nom){
    ia = i+1
    ib = ia - 1
    ic = ia - 2
    c = 2*(i-1) + 1 
    f=i-1
    c = c/i
    f = f/i
    for(j in 1:ia){
      if(j == 1) z = 0 
      else z = phi[ib,j-1] 
      phi[ia,j] = c*z - f*phi[ic,j]
    } 
  }
  for( m in 1:no){
    f = sqrt((2*(m-1)+1)/2) 
    phi[m, ] = phi[m, ]*f
  }
  return(phi[1:no,1:no])
}

time.matrix <- function(polyOrder, times.scaled) {
  # Build the matrix M
  n = length(times.scaled)
  M = matrix(0, nrow = n, ncol = polyOrder+1)
  M[,1] = 1
  if(polyOrder > 0)
    for(i in 1:polyOrder)
      M[,(i+1)] <- times.scaled^i
  M
}

#install.packages("mpoly")
library(mpoly)
library(polynom)
polyExp <- legendre(0:8, normalized=T)  # symbolic 
polyExp <- as.function(polyExp)  # make function
polyExp(q)  # = M %*% Lambda


# Read illumination data
data <- read.csv("/Users/marya/Documents/QS-II/Analysis/pixInt.csv")

# Get a starting covariance matrix from phenotypes
V.all <- cov(t(data[,-c(1,2)]))

round(V.all, 3)
# Goes up with increased illumination and then dwindles down

tper <- seq(2, 24, by = 3)
V <- V.all[tper, tper]  # 3-hourly V
invV <- solve(V)

ts <- get.qi(data$Time[tper], mn = 60, mx = 1475)  # Scaled times

M = time.matrix(polyOrder = 7, times.scaled = ts)  # Powered up scaled times
L = legendre.poly(8)  
PHI = M %*% t(L)
# OR, you can obtain PHI using the poly/polynom package
polyExp <- legendre(0:7, normalized=T)  # symbolic 
polyExp <- as.function(polyExp)  # make function
PHI = polyExp(ts)  # = M %*% Lambda

K = solve(PHI) %*% V %*% t(solve(PHI))

# Let's say that we want to get the cov's at times 90 and 150, 400
Mnew <- time.matrix(polyOrder = 7, times.scaled = get.qi(c(90, 150, 400), 60, 1475))
PHInew <- Mnew %*% t(L)
# V = PHI . K . PHI'
Vnew <- PHInew %*% K %*% t(PHInew)

# Fit legendre polynomials (of order p) within concentration
# -----------------------------------------------------------------
p = 7  # order of fit

# Restructure data
# ----------------
restData <- data.frame(time = rep(data$Time, 12), rep=rep(rep(1:2, 12), 12), 
                       strain = rep(1:12, each = 24), 
                       conc = c(rep('ctrl', 4*24), rep('mh15', 2*24), 
                                rep('mh30', 2*24), rep('mh60', 2*24),
                                rep('mh90', 2*24)), 
                       pxInt=unlist(c(data[,-(1:2)])))

restData$conc <- factor(restData$conc)
# Matrix X for the 'conc' factor has the same number of rows as in data and the 
#   same number of columns as in power+1 for each level. Let no = p + 1 
#   and N = nrow(data), dim(X) = N by nlevels(conc) * no
incX <- matrix(0, nrow=nrow(restData), ncol=(p+1)*nlevels(restData$conc))

# Construct polynomials for each time point in data
ts <- get.qi(restData$time, mn = 60, mx = 1475)  # Scaled times
M = time.matrix(polyOrder = p, times.scaled = ts)  # Powered-up scaled times
L = legendre.poly(no = (p+1))  
PHI = M %*% t(L)
# OR, you can obtain PHI using the poly/polynom package
polyExp <- as.function(legendre(0:p, normalized=T))  
PHI = polyExp(ts)  # = M %*% Lambda

# For each level in 'conc', add the corresponding row in PHI at the appropriate 
# starting column in the design matrix 
for(i in 1:nlevels(restData$conc)) {
  left = (i-1)*(p+1) + 1; right = i*(p+1)
  levelRows <- which(restData$conc == levels(restData$conc)[i])
  incX[levelRows, left:right] = PHI[levelRows,]
}
  
# Build and solve the normal equations 
a <- solve(t(incX) %*% incX) %*% t(incX) %*% restData$pxInt

# Equivalently, you can use lm()
d <- as.data.frame(cbind(restData$pxInt, incX))
b <- coef(lm(V1 ~ -1 + ., data = d))
cor(a, b)

summary(lm(V1 ~ -1 + ., data = d))
# p=2: R-squared:  0.9471,	Adjusted R-squared:  0.9442 
# p=5: R-squared:  0.9819,	Adjusted R-squared:  0.9798 
# p=6: R-squared:  0.9834,	Adjusted R-squared:  0.9812 
# p=7: R-squared:  0.9851,	Adjusted R-squared:  0.9827
# p=8: R-squared:  0.9852,	Adjusted R-squared:  0.9825 

# Get an estimate for px intensity w/t noise at each time point (12 pts)
phi <- unique(PHI)  # for 12 time pts

# obtain estimates for each time pt in the case of 'ctrl'
# 'ctrl' solution are the first p+1 estimates = a[1:(p+1)]
phi %*% a[1:(p+1)]  # 12 estimates for each time pt
mean(phi %*% a[1:(p+1)])  # an overall mean but much better than using a non longitudinal approach!

# Do the same for other concentrations

# Let's build a curve for each concentration (mn = 60, mx = 1475)
curve.times <- seq(60, 1475)
ts <- get.qi(curve.times, mn = 60, mx = 1475)  # Scaled times
# M = time.matrix(polyOrder = p, times.scaled = ts)  # Powered-up scaled times
# L = legendre.poly(no = (p+1))  
# PHI = M %*% t(L)
polyExp <- as.function(legendre(0:p, normalized=T))  
PHI = polyExp(ts) 

longt.results = data.frame(time=curve.times)
for(i in 1:nlevels(restData$conc)) {
  from = (i-1)*(p+1) + 1; to = i*(p+1)
  longt.results[,(i+1)] <- PHI %*% a[from:to]
  names(longt.results)[(i+1)] <- levels(restData$conc)[i]
}

# Illumination curves
par(mfrow=c(3,2))
graphLim = c(0, max(longt.results$ctrl)+1)
# ctrl
plot(longt.results$time, longt.results$ctrl, type='l', col='lightseagreen', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh15
plot(longt.results$time, longt.results$mh15, type='l', col='maroon4', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh30
plot(longt.results$time, longt.results$mh30, type='l', col='blue3', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh60
plot(longt.results$time, longt.results$mh60, type='l', col='purple', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh90
plot(longt.results$time, longt.results$mh90, type='l', col='orange', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)

# Plot all
# --------
plot(longt.results$time, longt.results$ctrl, type='l', col='lightseagreen', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh15
points(longt.results$time, longt.results$mh15, type='l', col='maroon4', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh30
points(longt.results$time, longt.results$mh30, type='l', col='blue3', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh60
points(longt.results$time, longt.results$mh60, type='l', col='purple', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)
# mh90
points(longt.results$time, longt.results$mh90, type='l', col='orange', 
     xlab='Time', ylab='Illumination Index', ylim=graphLim)

# Two barplots for maximum illumnation and peak time
par(mfrow=c(2,1))
peakIllumination <- c(max(longt.results$ctrl), max(longt.results$mh15), 
       max(longt.results$mh30), max(longt.results$mh60),
       max(longt.results$mh90))
barplot(peakIllumination, names.arg=c('Control','MH15', 'MH30', 'MH60', 'MH90'),
        xlab = 'Peak illumination intensity')
# Honey depresses the maximum possible illumination
peakTime <- c(longt.results[longt.results$ctrl == peakIllumination[1], 'time'],
              longt.results[longt.results$mh15 == peakIllumination[2], 'time'],
              longt.results[longt.results$mh30 == peakIllumination[3], 'time'],
              longt.results[longt.results$mh60 == peakIllumination[4], 'time'],
              longt.results[longt.results$mh90 == peakIllumination[5], 'time'])
barplot(peakTime, names.arg=c('Control','MH15', 'MH30', 'MH60', 'MH90'),
        xlab = 'Time to peak illumination (minutes)')
# Honey also extended the time for peak illumination (if ctrl is ignored)


