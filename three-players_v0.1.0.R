# three-player game
# PARAMETERS:
args <- as.numeric(commandArgs(trailingOnly = TRUE))
a <- args[1]
b <- args[2]
c <- args[3]

REPS <- 5
store.length <- 100 # number of steps that are stored at one time to check for convergence
step.sd <- 0.01

# PAYOFF FUNCTIONS:
# input: {t1,t2,t3} => output: {r1,r2,r3}
r <- function(T){
  if(T[1]+T[2]+T[3]==0) c(a/(a+b+c),b/(a+b+c),c/(a+b+c)) else {
    c(a*T[1]*(1-b*T[2])*(1-c*T[3])/
        (a*T[1]+b*T[2]+c*T[3]-2*(a*b*T[1]*T[2]+a*c*T[1]*T[3]+b*c*T[2]*T[3])+3*a*b*c*T[1]*T[2]*T[3]),
      b*T[2]*(1-a*T[1])*(1-c*T[3])/
        (a*T[1]+b*T[2]+c*T[3]-2*(a*b*T[1]*T[2]+a*c*T[1]*T[3]+b*c*T[2]*T[3])+3*a*b*c*T[1]*T[2]*T[3]),
      c*T[3]*(1-a*T[1])*(1-b*T[2])/
        (a*T[1]+b*T[2]+c*T[3]-2*(a*b*T[1]*T[2]+a*c*T[1]*T[3]+b*c*T[2]*T[3])+3*a*b*c*T[1]*T[2]*T[3]))
  }
}

# FUNCTION: calculate realized proportion of time spent displaying given strategies
# input: {u1,u2,v1,v3,w2,w3} => output: {t1,t2,t3} 
display.time <- function(X){
  c(max(0,1-max(X[1],X[2])-max(X[3],X[4])),
    max(0,1-max(X[1],X[2])-max(X[5],X[6])),
    max(0,1-max(X[3],X[4])-max(X[5],X[6])))
}

# STEP FUNCTION:
# input: all current strategies => output: new strategy
step <- function(X,UVW){
  new <- rnorm(1,X[UVW],step.sd)
  if(new<0) new <- 0
  if(new>1) new <- 1
  new
}

# create blank data frame for results
results <- data.frame("a" = rep(a,REPS),"b" = rep(b,REPS),"c" = rep(c,REPS),
                      "u1" = rep(NA,REPS),"u2" = rep(NA,REPS),"v1" = rep(NA,REPS),
                      "v3" = rep(NA,REPS),"w2" = rep(NA,REPS),"w3" = rep(NA,REPS),
                      "t1" = rep(NA,REPS),"t2" = rep(NA,REPS),"t3" = rep(NA,REPS))
# create blank matrix for path
path <- matrix(NA,store.length,3)

# START SIMULATION HERE:
for(i in 1:REPS){
  # set step counter
  STEPS <- 0
  
  # begin simulation with a randomly chosen set of strategies (seed)
  valid <- FALSE
  while(valid == FALSE){
    x <- runif(6)
    # verify that this is a valid seed
    valid <- x[1]+x[3]<=1 & x[2]+x[5]<=1 & x[4]+x[6]<=1
  }
  t <- display.time(x)
  
  # store current realized display times in path matrix
  path[1,] <- t
  
  # begin mcmc
  converge <- FALSE
  while(converge==FALSE){
    STEPS <- STEPS+1
    
    # select one strategy to change and make a new proposal
    uvw <- sample(1:6,1)
    valid <- FALSE
    while(valid==FALSE){
      prop <- x
      prop[uvw] <- step(x,uvw)
      # verify that this is a valid proposal
      valid <- prop[1]+prop[3]<=1 & prop[2]+prop[5]<=1 & prop[4]+prop[6]<=1
    }
    t.prop <- display.time(prop)
    
    # test new strategy for the relevant male and accept if r(t.prop) > r(t)
    male <- c(1,2,1,3,2,3)[uvw]
    if(r(t.prop)[male]>=r(t)[male]) x <- prop; t <- t.prop
    
    # store current realized display times in path matrix
    # but begin overwriting after completing "store.length" steps
    path[STEPS+1-floor(STEPS/store.length)*store.length,] <- t
    
    # see if path has converged
    converge <- length(unique(path))==3
    
    if(STEPS<10000) converge <- FALSE # run mcmc for at least 10,000 steps
    if(STEPS==5000000) converge <- TRUE # end mcmc after 5 million steps no matter what
  }
  # report the final positions of all six strategies
  results[i,4:9] <- x
  # round realized display times to 3 sigfigs and report path mode (equilibrium values)
  results[i,"t1"] <- as.numeric(names(sort(-table(round(path[,1],3)))))[1]
  results[i,"t2"] <- as.numeric(names(sort(-table(round(path[,2],3)))))[1]
  results[i,"t3"] <- as.numeric(names(sort(-table(round(path[,3],3)))))[1]
  
  results$meta_step.sd <- step.sd
  results$meta_store.length <- store.length
  results$meta_STEPS[i] <- STEPS
  results$meta_ess[i] <- length(unique(path))==3
}

# GENERATE OUTPUT:
out <- data.frame("a" = a,"b" = b,"c" = c)
out$t1 <- median(results[,"t1"])
out$t2 <- median(results[,"t2"])
out$t3 <- median(results[,"t3"])
out$u <- (1-out$t1-out$t2+out$t3)/2
out$v <- (1-out$t1-out$t3+out$t2)/2
out$w <- (1-out$t2-out$t3+out$t1)/2
out$r1 <- r(out[,4:6])[1]
out$r2 <- r(out[,4:6])[2]
out$r3 <- r(out[,4:6])[3]

out$meta_REPS <- REPS

# check that results converge within 0.01
out$converge <- sum(abs(results[,"t1"]-out$t1)<0.01,
                    abs(results[,"t2"]-out$t2)<0.01,
                    abs(results[,"t3"]-out$t3)<0.01)/(3*REPS)

write.csv(results,paste("./detail/",paste("detail",a,b,c,sep = "-"),".csv",sep = ""))
out <- as.matrix(out)
write.table(out,"data_out.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = FALSE)
