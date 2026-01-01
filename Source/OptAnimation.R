#optimisation

library(tidyverse,attach.required = TRUE)
library(gganimate)
library(gifski)
library(GA)
library(numDeriv)

minmax <- c(-5,5)
DIM <- 2
resolution <- 0.05

startbounds <- c(-5,-4)
iterlim <- 100

obj <- function(x){ t(x)%*%x } #sum(x^2) 
# drop ensures the 1x1 matrix becomes a scalar

#possibly useful?
obj.df <- function(X){ 
  x<-c(X) 
  t(x) %*% x
}


base_theme <- theme_minimal() +
  theme(aspect.ratio=1,legend.position = 'bottom')
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

set.seed(1)
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

step_gradient_descent <- function(x,nu=nu_global){
  next_x <- x-nu*grad(obj,x)
  return(next_x)
}

step_random_search <- function(x,nu=nu_global,bounds=minmax){
  next_x <- runif(n=DIM,bounds[1],bounds[2])
  return(next_x)
}

step_newton_raphson <- function(x,nu=nu_global){

  next_x <- x - nu*(solve(hessian(obj,x))) %*% (grad(obj,x))
  return(next_x)
}

step_genetic_algorithm <- function(){
  
}

update_opt <- function(x, nu=nu_global){
  step_random_search(x,nu)
}


xi <- min_x <- start
min_y <- obj(min_x)

opt_trace[1,] <- c(1,xi,min_x)

for(i in 2:iterlim){
  
  xi <- update_opt(xi)
  yi <- obj(xi)
  
  if(yi<min_y){
    min_y <- yi
    min_x <- xi
  }
  
  opt_trace[i,] <- c(i,xi,min_x) 
  
}

#outrs <-
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
  labs(subtitle = "Iteration: {previous_state}."))

anim_save(nframes= iterlim,
          filename="newt_raph_1.gif",
          path="../Render",
          animation=outnr)

outgd
outrs
outnr
#outga



