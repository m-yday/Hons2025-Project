#optimisation

################################################################################
#                        setup + objective functions                           #

library(tidyverse,attach.required = TRUE)
library(gganimate)
library(gifski)
library(GA)
library(numDeriv)
library(patchwork)

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
  min_point <<- rep(0,times=n)
  transform_val <- "identity"
  breaks_val    <- waiver()
  return(f)
}
test_Ackley <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  f <- -20*exp(-0.2*sqrt(0.5*(x1^2+x2^2))) -
    exp(0.5*(cos(2*pi*x1)+cos(2*pi*x2))) +
    exp(1) + 20
  min_point <<- c(0,0)
  transform_val <- "identity"
  breaks_val    <- waiver()
  return(f)
} # 2 dim
test_sphere <- function(x){
  n <- length(x)
  f <- t(x) %*% x
  min_point <<- rep(0,times=n)
  transform_val <- "identity"
  breaks_val    <- waiver()
  return(f)
}
test_eggholder <- function(x){
  x1 <- x[1]*200 #rescaling
  x2 <- x[2]*200
  f  <- -(x2+47)*sin(sqrt(abs(x1/2+(x2+47))))-x1*sin(sqrt(abs(x1-(x2+47))))
  min_point <<- c(512,404.2319)
  transform_val <- "identity"
  breaks_val    <- waiver()
  return(f)
} # 2 dimensions only
test_Rosenbrock <- function(x){
  n <- length(x)
  f <- sum(100*(x[-1] -(x[-n])^2)^2 + (1-x[-n])^2)
  min_point <<- rep(1,times=n)
  transform_val <- "log1p"
  breaks_val    <- c(1,10,100,1000,10000,50000)
  return(f)
} # dim > 1

# drop ensures the 1x1 matrix becomes a scalar
#possibly useful?
# obj.df <- function(X){ 
#   x<-c(X) 
#   t(x) %*% x
# }

base_theme <- theme_minimal() 
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

#-#-# #-#-# #-#-# #-#-# #        TEST FUNCTION         # #-#-# #-#-# #-#-# #-#-#
obj <- function(x){ 
  return(test_Rosenbrock(x))
}

# settings for colour scales:
# passed to scale_fill_viridis_c(), & eventually, continuous_scale()

# linear: 
# transform_val <- "identity"
# breaks_val    <- waiver()

# log:
transform_val <- "log1p"
breaks_val    <- c(1,10,100,1000,10000,50000)

################################################################################
#                            optimisation functions                            #

#                                setup surface

# gg is the surface information of the objective function. 
ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                      .name_repair=~make.names(1:DIM))|>
  mutate(Y=pmap_vec(.l=pick(everything()),
                    .f=~ obj(c(...))))

#plot the objective surface
p_obj <- ggplot()+
  geom_raster(data=ggdata,aes(x=X1,y=X2,fill=Y)) + 
  geom_point(data=as_tibble(t(min_point)), aes(x=V1,y=V2),
             colour='red',shape=4,alpha=0.5)+
  scale_fill_viridis_c(option="E",
                       transform=transform_val,breaks=breaks_val
  ) + coord_fixed()
p_obj

nu_global<-0.05

# optimisation function options
# these functions take in the current position (x), 
# returning the next x value based on their algorithm.

step_gradient_descent  <- function(x,nu=0.0002){
  next_x <- x-nu*grad(obj,x,method="simple")
  return(next_x)
}
# step_random_search must use MIN_X, not simply X when selecting update procedure.
# step_random_search (FSSRS) has a normalised (1) step size. nu acts as a multiplicative factor on this.
step_random_search     <- function(x,nu=.2){
  u <- rnorm(n=DIM)
  r <- sqrt(t(u)%*%u)
  next_x <- x+nu*(u/r)
  return(next_x)
}
# nu unused in step_pure_rs
step_pure_rs           <- function(x,nu=.2,bounds=minmax){
  next_x <- runif(n=DIM,bounds[1],bounds[2]) 
  return(next_x)
}
step_newton_raphson    <- function(x,nu=0.5){
  next_x <- x - nu*(solve(hessian(obj,x))) %*% (grad(obj,x))
  return(next_x)
}
step_genetic_algorithm <- function(x){ # might need to rethink algorithm selection
  #real valued genetic algorithm
  ga(type="real-valued",
     fitness= obj,
     
     )
  
}

################################################################################

#-#-# #-#-# #-#-# #-#-# #       UPDATE PROCEDURE       # #-#-# #-#-# #-#-# #-#-#
update_opt <- function(x, min_x, nu=nu_global){
  step_gradient_descent(x)
}

#-#-# #-#-# #-#-# #-#-# #       RUN OPTIMISATION       # #-#-# #-#-# #-#-# #-#-#

run_optim <- function(x_0, names){
  
}

#next we follow the path of the optimisation function
set.seed(111)

start <- runif(n=DIM,startbounds[1],startbounds[2])
names(start)=names(ggdata)[1:DIM]
#as_tibble(t(start)) # -4.41,-4.27

opt_trace <- matrix(NA, nrow=iterlim, ncol=2*(DIM+1)+1,
                    dimnames= list(NULL,
                                   c('Iteration',
                                     names(ggdata),
                                     paste0('Min',names(ggdata)))
                    )
)

#                              run optimisation


xi <- min_x <- start
yi <- min_y <- obj(min_x)

opt_trace[1,] <- c(1,xi,yi,min_x,min_y)

for(i in 2:iterlim){
  
  xi <- update_opt(xi,min_x) 
  yi <- obj(xi) #objective being minimised
  
  if(yi<min_y){
    min_y <- yi
    min_x <- xi
  }
  
  opt_trace[i,] <- c(i,xi,yi,min_x,min_y) 
  
}

#                              plot and animate


trace_segments <- 
  bind_rows(bind_cols(opt_trace[c(1,1),],facet=1), 
            #duplicated step 1 due to minimum count being 2.
            bind_cols(opt_trace[1:2   ,],facet=2),
            bind_cols(opt_trace[1:5   ,],facet=3),
            bind_cols(opt_trace[1:10  ,],facet=4),
            bind_cols(opt_trace[1:50  ,],facet=5),
            bind_cols(opt_trace[1:100 ,],facet=6))
facet_names <- c( `1` =   "1 Step",
                  `2` =   "2 Steps",
                  `3` =   "5 Steps",
                  `4` =  "10 Steps",
                  `5` =  "50 Steps",
                  `6` = "100 Steps")
p_static <- p_obj+
  geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
             colour='yellow',shape=4)+
  geom_point(data=trace_segments,aes(x=X1,y=X2,group=Iteration),
             colour="white",alpha=0.3,shape=4)+
  geom_point(data=trace_segments,aes(x=MinX1,y=MinX2,colour=Iteration))+
  scale_colour_viridis_c(option="H",begin=0.3,end=0.9) +
  facet_wrap(facets=~facet,nrow=3,
             labeller = as_labeller(facet_names))+ 
  xlim(minmax*1.01) + ylim(minmax*1.01)

p_static

p_performance <- ggplot(data=opt_trace, aes(x=Iteration)) +
  geom_line(aes(y=Y), colour="grey")+
  geom_line(aes(y=MinY, colour = Iteration),show.legend = FALSE)+
  scale_colour_viridis_c(option="H",begin=0.3,end=0.9) +
  labs(y="Y", x="Iteration",
       title="Behaviour of the gradient descent optimiser on a Rosenbrock surface.",
       subtitle="With nu = 0.0002.")+ 
  theme_bw() 

p_performance

p_full <- p_performance/p_static +
  plot_layout(heights=unit(c(20,1),c('mm','null')))

  
p_full
ggsave(plot=p_full,filename = "../Render/rosen-grad-desc-nu0.0002-%03d.svg",
       width=120, height=180, units="mm", dpi="retina")

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

anim_save(nframes= iterlim,
          filename="rosen-grad-desc-nu0.0002.gif",
          path="../Render",
          animation=p_anim)


# cannot train GA stepwise. The GA must be fitted first and then tested.

#minimisation problem must be changed to maximisation problem
ga_obj <- function(x){
  -obj(x)
}

ga_model <- ga(type="real-valued",
   fitness= ga_obj, 
   lower = rep(minmax[1],length.out=DIM), #length of lower/upper must = DIM
   upper = rep(minmax[2],length.out=DIM),
   keepBest = TRUE,
   monitor = plot
)

ga_model@population

# now, GA doesn't give us every single value tested.
# but it does give us the best value tested for each population iteration.

# we define a monitoring function so we can attempt to note all populations generated.
# let us do this with tidyverse
custom_ga_monitor <- function(model) {
  best <- which.max(model@fitness)
  pop <- model@population |> as_tibble() |> bind_cols(Iteration=model@iter)
  plot(p_obj+
         geom_point(data=pop,aes(x=V1,y=V2,group=Iteration),
                    colour="white",shape=4,alpha=0.2)+
         labs(subtitle=paste("iteration =", model@iter))+
         geom_point(data=pop[best,],aes(x=V1,y=V2,colour=Iteration))+
         scale_colour_viridis_c(option="H",begin=0.3,end=0.9)
    )
}
GA <- ga(type = "real-valued",
         fitness = ga_obj,
         lower = c(-5, -5), upper = c(5, 5), popSize = 10,
         maxiter = iterlim, 
         monitor = custom_ga_monitor)
summary(GA)