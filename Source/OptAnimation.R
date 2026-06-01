#optimisation

################################################################################
#                        setup + objective functions                           #

library(tidyverse,attach.required = TRUE)
library(gganimate)
library(gifski)
library(GA)
library(numDeriv)

scalefactor <- 1
minmax <- c(-5,5)*scalefactor
DIM <- 2
resolution <- 0.05*scalefactor

startbounds <- c(-5,-4)*scalefactor
iterlim <- 100


#                    objective function options
test_Rastrigin <- function(x){
  A<-10
  n<-length(x)
  f <- A*n + sum(x^2-A*cos(2*pi*x))
  return(f)
}
test_Ackley <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  f <- -20*exp(-0.2*sqrt(0.5*(x1^2+x2^2))) -
    exp(0.5*(cos(2*pi*x1)+cos(2*pi*x2))) +
    exp(1) + 20
  return(f)
} # 2 dim
test_sphere <- function(x){
  f <- t(x) %*% x
  return(f)
}
test_eggholder <- function(x){
  x1 <- x[1]*200 #rescaling
  x2 <- x[2]*200
  f  <- -(x2+47)*sin(sqrt(abs(x1/2+(x2+47))))-x1*sin(sqrt(abs(x1-(x2+47))))
  return(f)
} # 2 dimensions only
test_Rosenbrock <- function(x){
  n <- length(x)
  f <- sum(100*(x[-1] -(x[-n])^2)^2 + (1-x[-n])^2)
  return(f)
} # dim > 1

# drop ensures the 1x1 matrix becomes a scalar
#possibly useful?
obj.df <- function(X){ 
  x<-c(X) 
  t(x) %*% x
}

base_theme <- theme_minimal() +
  theme(aspect.ratio=1,legend.position = 'bottom'
        #panel.background = element_rect(fill='grey80')
  )
theme_set(base_theme)

xseq <- seq(minmax[1],minmax[2],resolution) #symmetric. 
# geom_tile is based on the centre of the tile.

# expand_grid for any amount of dimensions 
# (warning! grows by length(xseq) exponentially!)
# equivalent to: 
# expand.grid(rep(list(xseq),DIM)) 
# due to tidyr's different implementation of expand.grid
# or:
# expand_grid(xseq,xseq)   for DIM = 2

################################################################################

transform_val = "identity"
breaks_val = c(1,10,100,1000,20000)

#-#-# #-#-# #-#-# #-#-# #        TEST FUNCTION         # #-#-# #-#-# #-#-# #-#-#
obj <- function(x){ 
  return(test_Ackley(x))
}

################################################################################
#                        setup + optimisation functions                        #

ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                      .name_repair=~make.names(1:DIM))|>
  mutate(Y=pmap_vec(.l=pick(everything()),
                    .f=~ obj(c(...))))

#plot the objective surface
p_obj <- ggplot()+
  geom_raster(data=ggdata,aes(x=X1,y=X2,fill=Y)) + 
  scale_fill_viridis_c(option="E"
                       ,transform=transform_val,breaks=breaks_val
                       )
p_obj


#next we follow the path of the optimisation function
#set.seed(1)

start <- runif(n=DIM,startbounds[1],startbounds[2])
names(start)=names(ggdata)[1:DIM]
as_tibble(t(start))

opt_trace <- matrix(NA, nrow=iterlim, ncol=2*DIM+1,
                    dimnames= list(NULL,
                                   c('Iteration',
                                     names(ggdata)[1:DIM],
                                     paste0('Min',names(ggdata)[1:DIM]))
                                   )
                    )

nu_global<-0.05

# optimisation function options

step_gradient_descent  <- function(x,nu=0.2){
  next_x <- x-nu*grad(obj,x,method="simple")
  return(next_x)
}
step_random_search     <- function(x,nu=.2){
  u <- rnorm(n=DIM)
  r <- sqrt(t(u)%*%u)
  next_x <- x+nu*(u/r)
  return(next_x)
}
step_pure_rs           <- function(x,nu=.2,bounds=minmax){
  next_x <- runif(n=DIM,bounds[1],bounds[2]) 
  return(next_x)
}
step_newton_raphson    <- function(x,nu=0.5){
  next_x <- x - nu*(solve(hessian(obj,x))) %*% (grad(obj,x))
  return(next_x)
}
step_genetic_algorithm <- function(){
  #real valued genetic algorithm
  
}

################################################################################

#-#-# #-#-# #-#-# #-#-# #       UPDATE PROCEDURE       # #-#-# #-#-# #-#-# #-#-#
update_opt <- function(x, min_x, nu=nu_global){
  step_gradient_descent(x)
}

#                      run optimisation, plot, and animate

xi <- min_x <- start
min_y <- obj(min_x)

opt_trace[1,] <- c(1,xi,min_x)

for(i in 2:iterlim){
  
  xi <- update_opt(xi,min_x) 
  yi <- obj(xi) #objective being minimised
  
  if(yi<min_y){
    min_y <- yi
    min_x <- xi
  }
  
  opt_trace[i,] <- c(i,xi,min_x) 
  
}


p_static <- p_obj+
  geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
             colour='yellow',shape=4)+
  geom_point(data=opt_trace,aes(x=X1,y=X2,group=Iteration),
             colour="black",shape=4)+
  geom_point(data=opt_trace,aes(x=MinX1,y=MinX2,colour=Iteration))+
  scale_colour_viridis_c(option="H",begin=0.3,end=0.9) #+ xlim(minmax) + ylim(minmax)

p_static

p_anim <-
animate(
  p_obj+
  geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
             colour='yellow',shape=4)+
  geom_point(data=opt_trace,aes(x=X1,y=X2,group=Iteration),
             colour="white",shape=4)+
  geom_point(data=opt_trace,aes(x=MinX1,y=MinX2,colour=Iteration))+
  scale_colour_viridis_c(option="H",begin=0.3,end=0.9)+
  transition_states(states=Iteration,wrap=FALSE)+
  shadow_wake(wake_length=0.3)+
  labs(subtitle = "Iteration: {previous_state}.")  + xlim(minmax) + ylim(minmax)
)

p_anim

# anim_save(nframes= iterlim,
#           filename="rand_search_n0.2.gif",
#           path="../Render",
#           animation=p_anim)
