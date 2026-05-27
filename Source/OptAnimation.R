#optimisation

library(tidyverse,attach.required = TRUE)
library(gganimate)
library(gifski)
library(GA)
library(numDeriv)

scalefactor <- 1
minmax <- c(-5,5)*scalefactor
DIM <- 2
resolution <- 0.05*scalefactor

startbounds <- c(-1,-1)*scalefactor
iterlim <- 100

test_rastrigin <- function(x){
  A<-10
  n<-length(x)
  f <- A*n + sum(x^2-A*cos(0.5*pi*x))
  return(f)
}

# 2 dimensions only
test_eggholder <- function(x){
  x1<-x[1]
  x2<-x[2]
  f <- -(x2+47)*sin(sqrt(abs(x1/2+(x2+47))))-x1*sin(sqrt(abs(x1-(x2+47))))
}

obj <- function(x){ 
  return(test_rastrigin(x))
}
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

xseq <- seq(minmax[1],minmax[2],resolution) #symmetric. geom_tile is based on the centre of the tile.


#expand_grid for any amount of dimensions (warning! grows by length(xseq) exponentially!)
# equivalent to: 
# expand.grid(rep(list(xseq),DIM)) due to tidyr's different implementation of expand.grid
# or:
# expand_grid(xseq,xseq)   for DIM = 2

ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                      .name_repair=~make.names(1:DIM))|>
  mutate(Y=pmap_vec(.l=pick(everything()),
                    .f=~ obj(c(...))))

#plot the objective surface
pobj <- ggplot()+
  geom_raster(data=ggdata,aes(x=X1,y=X2,fill=Y)) + 
  scale_fill_viridis_c(option="E")
pobj


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

step_gradient_descent <- function(x,nu=0.2){
  next_x <- x-nu*grad(obj,x)
  return(next_x)
}

step_random_search <- function(x,nu=.2){
  u <- rnorm(n=DIM)
  r <- sqrt(t(u)%*%u)
  next_x <- x+nu*(u/r)
  return(next_x)
}

step_pure_rs <- function(x,nu=.2,bounds=minmax){
  next_x <- runif(n=DIM,bounds[1],bounds[2]) 
  return(next_x)
}

step_newton_raphson <- function(x,nu=0.5){
  next_x <- x - nu*(solve(hessian(obj,x))) %*% (grad(obj,x))
  return(next_x)
}

step_genetic_algorithm <- function(){
  #real valued genetic algorithm
  
}

# defines which updating equation is used
update_opt <- function(x, min_x, nu=nu_global){
  step_random_search(x)
}


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


p_static <- pobj+
  geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
             colour='yellow',shape=4)+
  geom_point(data=opt_trace,aes(x=X1,y=X2,group=Iteration),
             colour="white",shape=4)+
  geom_point(data=opt_trace,aes(x=MinX1,y=MinX2,colour=Iteration))+
  scale_colour_viridis_c(option="H",begin=0.3,end=0.9)

p_static

p_anim <-
animate(
  pobj+
  geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
             colour='yellow',shape=4)+
  geom_point(data=opt_trace,aes(x=X1,y=X2,group=Iteration),
             colour="white",shape=4)+
  geom_point(data=opt_trace,aes(x=MinX1,y=MinX2,colour=Iteration))+
  scale_colour_viridis_c(option="H",begin=0.3,end=0.9)+
  transition_states(states=Iteration,wrap=FALSE)+
  shadow_wake(wake_length=0.3)+
  labs(subtitle = "Iteration: {previous_state}.")
)

p_anim


# anim_save(nframes= iterlim,
#           filename="rand_search_n0.2.gif",
#           path="../Render",
#           animation=outrs)



# outgd
# outprs
# outrs
# outnr
#outga

# pnr
# pprs
# prs
# pgd
# pga

