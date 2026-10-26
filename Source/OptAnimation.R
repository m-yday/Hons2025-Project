################################################################################
#                        setup + objective functions                           #

library(tidyverse,attach.required = TRUE)
library(gganimate)
library(gifski)
library(GA)
library(numDeriv)
library(patchwork)
library(futurize)
library(progressify)

plan(multisession)
handlers(global=TRUE)


# Major effect on runtime. 
# 0.01: HIGH quality. 
# 0.05: Normal quality.
resolution <- 0.01 

# Plots work with 2 dimensions.
DIM <- 2


iterlim <- 100
popSize <- 10

trials <- 100

# NU: LEARNING RATE (fallback default, typically remains overwritten)
nu_global<-0.05

#                    objective function options (test surfaces)
test_Rastrigin <- function(x){
  A<-10
  n<-length(x)
  f <- A*n + sum(x^2-A*cos(2*pi*x))
  #min_point <- rep(0,times=n)
  return(f)
}
test_Ackley <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  f <- -20*exp(-0.2*sqrt(0.5*(x1^2+x2^2))) -
    exp(0.5*(cos(2*pi*x1)+cos(2*pi*x2))) +
    exp(1) + 20
  #min_point <<- c(0,0)
  return(f)
} # 2 dim
test_parabola <- function(x){
  #n <- length(x)
  f <- t(x) %*% x
  #min_point <<- rep(0,times=n)
  return(f)
}
test_eggholder <- function(x){
  x1 <- x[1]*200 #rescaling
  x2 <- x[2]*200
  f  <- -(x2+47)*sin(sqrt(abs(x1/2+(x2+47))))-x1*sin(sqrt(abs(x1-(x2+47))))
  #min_point <<- c(512,404.2319)
  return(f)
} # 2 dimensions only
test_Rosenbrock <- function(x){
  n <- length(x)
  f <- sum(100*(x[-1] -(x[-n])^2)^2 + (1-x[-n])^2)
  #min_point <<- rep(1,times=n)
  return(f)
} # dim > 1
test_Styblinski_Tang <- function (x){
  #n <- length(x)
  f <- 0.5 * sum( x^4 -16*x^2 +5*x )
  #min_point <<- rep(-2.903534,times=n)
  return(f)
}


# optimisation function options
# these functions take in the current position (x), 
# returning the next x value based on their algorithm.

step_gradient_descent  <- function(x,nu,OBJ){
  #OBJ: only func used
  next_x <- x-nu*grad(OBJ$func,x,method="simple")
  return(next_x)
}
# step_random_search must use MIN_X, not simply X when selecting update procedure.
# step_random_search (FSSRS) has a normalised (1) step size. nu acts as a multiplicative factor on this. 
step_random_search     <- function(x,nu,OBJ){
  # OBJ unused
  u <- rnorm(n=2)
  #hardcoded DIM
  r <- sqrt(t(u)%*%u)
  next_x <- x+nu*(u/r)
  return(next_x)
}

# nu unused in step_pure_rs
step_pure_rs           <- function(x,nu=NA,OBJ){
  #OBJ: only bounds used
  next_x <- runif(n=2,OBJ$bounds[1],OBJ$bounds[2]) 
  #hardcoded DIM
  return(next_x)
}

step_newton_raphson    <- function(x,nu,OBJ){
  #OBJ: only func used
  next_x <- x - nu*(solve(hessian(OBJ$func,x))) %*% (grad(OBJ$func,x))
  return(next_x)
}


# Objective functions
# these are specifically defined for 2 dimensions.

# S3 class constructor
new_objective_surface <- function(func, plot_name, file_prefix, min,
                                  startbounds, bounds, transform, breaks) {
  structure(
    list(func        = func,
         plot_name   = plot_name,
         file_prefix = file_prefix,
         min         = min,
         startbounds = startbounds,
         bounds      = bounds,
         transform   = transform,
         breaks      = breaks),
    class = "objective_surface"
  )
}
# func: name of test function.
# plot_name: names on ggplot graphs. "...on `plot_name` surface."
# file_prefix: rendered files saved with filename: `file_name`-`optim_name`.svg
# min: minimum point(s) of test function.
# startbounds: 1-dimensional bounds from which the initial position is randomly selected. This is recycled for both dimensions. 
# bounds: 1-dimensional bounds of the entire search space. This is recycled for both dimensions.
# transform, breaks: passed to ggplot scale functions to ensure correct scaling.

obj_Rastrigin <- new_objective_surface(
  func        = test_Rastrigin,
  plot_name   = "a Rastrigin",
  file_prefix = "rastrigin",
  min         = c(0,0),
  startbounds = c(-5,-4),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = waiver()
)

obj_Ackley <- new_objective_surface(
  func        = test_Ackley,
  plot_name   = "an Ackley",
  file_prefix = "ackley",
  min         = c(0,0),
  startbounds = c(-5,-4),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = waiver()
)

obj_parabola <- new_objective_surface(
  func        = test_parabola,
  plot_name   = "a parabolic",
  file_prefix = "para",
  min         = c(0,0),
  startbounds = c(-5,-4),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = waiver()
)

obj_eggholder <- new_objective_surface(
  func        = test_eggholder,
  plot_name   = "an eggholder",
  file_prefix = "egg",
  min         = c(512,404.2319),
  startbounds = c(-550,-450),
  bounds      = c(-550,550),
  transform   = "identity",
  breaks      = waiver()
)

obj_Rosenbrock <- new_objective_surface(
  func        = test_Rosenbrock,
  plot_name   = "a Rosenbrock",
  file_prefix = "rosen",
  min         = c(1,1),
  startbounds = c(-5,-4),
  bounds      = c(-5,5),
  transform   = "log1p",
  breaks      = c(1,10,100,1000,10000)
)

obj_Styblinski_Tang <- new_objective_surface(
  func        = test_Styblinski_Tang,
  plot_name   = "a Styblinski-Tang",
  file_prefix = "styb-tang",
  min         = c(-2.903534,-2.903534),
  startbounds = c(4,5),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = waiver()
)

# Update Procedures

# S3 class constructor
new_update_procedure <- function(func, use_min_x, optim_name, plot_name, nu) {
  structure(
    list(func       = func,
         use_min_x  = use_min_x,
         optim_name = optim_name,
         plot_name  = plot_name,
         nu         = nu),
    class = "update_procedure"
  )
}
#func: optimiser function
#use_min_x: boolean. selects `min_x` or `x` as parameter for `func`.
#optim_name: identifier used in filenames for this optimiser.
#plot_name: name used in ggplot graphics for this optimiser.
#overwrite_nu: 

up_gradient_descent <- new_update_procedure(
  func       = step_gradient_descent,
  use_min_x  = FALSE,
  optim_name = "grad-desc",
  plot_name  = "gradient descent",
  nu         = 0.2
)

up_random_search <- new_update_procedure(
  func       = step_random_search,
  use_min_x  = TRUE,
  optim_name = "rand-sear",
  plot_name  = "FSSRS",
  nu         = 0.2
)

up_newton_raphson <- new_update_procedure(
  func       = step_newton_raphson,
  use_min_x  = FALSE,
  optim_name = "newt-raph",
  plot_name  = "Newton-Raphson",
  nu         = 0.5
)

up_pure_rs <- new_update_procedure(
  func       = step_pure_rs,
  use_min_x  = FALSE,
  optim_name = "pure-rs",
  plot_name  = "PRS",
  nu         = NULL 
)

OBJECTIVE_SURFACES <- list(
  obj_parabola,
  obj_Ackley,
  obj_Rastrigin,
  obj_eggholder,
  obj_Rosenbrock,
  obj_Styblinski_Tang
  )

UPDATE_PROCEDURES <- list(
  up_gradient_descent,
  up_newton_raphson,
  up_pure_rs,
  up_random_search
  )

#-#-# #-#-# #-#-# #-#-# #          OPTIMISATION        # #-#-# #-#-# #-#-# #-#-#
 
# generates test surface
generate_surface <- function(OBJ,.resolution=resolution){
  base_theme <- theme_minimal() 
  theme_set(base_theme)
  
  xseq <- seq(OBJ$bounds[1],OBJ$bounds[2],.resolution) #symmetric. 
  # geom_tile is based on the centre of the tile.
  
  # ggdata is the surface information of the objective function. 
  ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                        .name_repair=~make.names(1:DIM))|>
    mutate(Y=pmap_vec(.l=pick(everything()),
                      .f=~ OBJ$func(c(...))))
  
  #plot the objective surface
  p_obj <- ggplot()+
    geom_raster(data=ggdata,aes(x=X1,y=X2,fill=Y)) + 
    geom_contour(data=ggdata,aes(x=X1,y=X2,z=Y),breaks=OBJ$breaks,
                 colour="white", linewidth = 0.2, alpha=0.8, linetype = 1)+
    geom_point(data=as_tibble(t(OBJ$min),
                              .name_repair=~make.names(1:DIM)),
               aes(x=X1,y=X2),
               colour='red',shape=4,alpha=0.8,show.legend = TRUE)+
    scale_fill_viridis_c(option="E",
                         transform=OBJ$transform,breaks=OBJ$breaks
    ) + coord_fixed() + 
    xlim(OBJ$bounds*1.01) +ylim(OBJ$bounds*1.01) # 1% expansion factor
  
  return(p_obj)
}

# runs optimisation. generates trace. generates plots. generates animation.
optim_main <- function(OBJ, UP, nu=NULL, trials=1,
                      .iterlim=iterlim, .resolution=resolution){
  
  DIM <- 2
  
  if(!is.null(nu)){
    UP$nu <- nu
  }
  
  base_theme <- theme_minimal() 
  theme_set(base_theme)
  
  xseq <- seq(OBJ$bounds[1],OBJ$bounds[2],.resolution) #symmetric. 
  # geom_tile is based on the centre of the tile.
  
  # ggdata is the surface information of the objective function. 
  ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                        .name_repair=~make.names(1:DIM))|>
    mutate(Y=pmap_vec(.l=pick(everything()),
                      .f=~ OBJ$func(c(...))))
  
  #plot the objective surface
  p_obj <- ggplot()+
    geom_raster(data=ggdata,aes(x=X1,y=X2,fill=Y)) + 
    geom_contour(data=ggdata,aes(x=X1,y=X2,z=Y),breaks=OBJ$breaks,
                 colour="white", linewidth = 0.2, alpha=0.8, linetype = 1)+
    geom_point(data=as_tibble(t(OBJ$min),
                              .name_repair=~make.names(1:DIM)),
               aes(x=X1,y=X2),
               colour='red',shape=4,alpha=0.8,show.legend = TRUE)+
    scale_fill_viridis_c(option="E",
                         transform=OBJ$transform,breaks=OBJ$breaks
    ) + coord_fixed() + 
    xlim(OBJ$bounds*1.01) +ylim(OBJ$bounds*1.01) # 1% expansion factor
  
  update_opt <- (
    if (UP$use_min_x) 
      function(x, min_x, nu=nu_global) UP$func(min_x, nu=UP$nu, OBJ) 
    else 
      function(x, min_x, nu=nu_global) UP$func(    x, nu=UP$nu, OBJ)
  )
  
  start <- runif(n=DIM,OBJ$startbounds[1],OBJ$startbounds[2])
  names(start)=names(ggdata)[1:DIM]
  
  opt_trace <- matrix(NA, nrow=.iterlim, ncol=7,
                      dimnames= list(NULL,
                                     c('Iteration',
                                       names(ggdata),
                                       paste0('Min',names(ggdata)))
                      )
  )
  # opt_trace <- array(NA, #trace matrix extended for multiple runs
  #                   dim = c(.iterlim,7,trials),
  #                   dimnames = list(iteration = 1:.iterlim,
  #                                   trace     = c("Iteration", "X1","X2","Y",
  #                                                 "MinX1","MinX2","MinY"),
  #                                   trial     = 1:trials))
  
  #                              run optimisation
  
  xi <- min_x <- start
  yi <- min_y <- OBJ$func(min_x)
  
  opt_trace[1,] <- c(1,xi,yi,min_x,min_y)
  
  #foreach(t = 1:trials) %dofuture%{
  for(i in 2:.iterlim){
    
    xi <- update_opt(xi,min_x) 
    yi <- OBJ$func(xi) #objective being minimised
    
    if(yi<min_y){
      min_y <- yi
      min_x <- xi
    }
    
    opt_trace[i,] <- c(i,xi,yi,min_x,min_y) 
    
  }
  #}
  #                              plot and animate
  
  # trace_segments: how much of the trace is used in each facet
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
               colour='aquamarine',shape=4)+
    geom_point(data=trace_segments,aes(x=X1,y=X2,group=Iteration),
               colour="white",alpha=0.3,shape=4)+
    geom_point(data=trace_segments,aes(x=MinX1,y=MinX2,colour=Iteration))+
    scale_colour_viridis_c(option="H",begin=0.3,end=0.9) +
    facet_wrap(facets=~facet,nrow=2,
               labeller = as_labeller(facet_names))
  
  p_performance <- ggplot(data=opt_trace, aes(x=Iteration)) +
    geom_line(aes(y=Y), colour="grey")+
    geom_line(aes(y=MinY, colour = Iteration),show.legend = FALSE)+
    scale_colour_viridis_c(option="H",begin=0.3,end=0.9) +
    labs(y="Y", x="Iteration",
         title= sprintf(
           "Behaviour of the %s optimiser on %s surface.",
           UP$plot_name, OBJ$plot_name))+ 
    theme_bw() 
  
  p_full <- p_performance/p_static +
    plot_layout(heights=unit(c(15,1),c('mm','null')))
  
  ggsave(plot=p_full,
         filename = sprintf("../Render/%s-%s.svg",
                            OBJ$file_prefix, UP$optim_name),
         width=180, height=150, units="mm", dpi=600)
  
  p_anim <-
    animate(
      p_obj+
        geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
                   colour='royalblue',shape=4)+
        geom_point(data=opt_trace,aes(x=X1,y=X2,group=Iteration),
                   colour="white",shape=4)+
        geom_point(data=opt_trace,aes(x=MinX1,y=MinX2,colour=Iteration))+
        scale_colour_viridis_c(option="H",begin=0.3,end=0.9)+
        transition_states(states=Iteration,wrap=FALSE)+
        shadow_wake(wake_length=0.3)+
        labs(subtitle = "Iteration: {previous_state}.")
      )
  
  anim_save(nframes= .iterlim,
            filename=sprintf("%s-%s.gif", OBJ$file_prefix, UP$optim_name),
            path="../Render",
            animation=p_anim)
  
  return(opt_trace)
}

set.seed(111)

p_Rosenbrock <- generate_surface(obj_Rosenbrock)
p_Rosenbrock + labs(title="Rosenbrock test surface.")
ggsave(filename = sprintf("../Render/%s.svg",obj_Rosenbrock$file_prefix),
      width=150, height=90, units="mm", dpi=600)

p_Styblinksi_Tang <- generate_surface(obj_Styblinski_Tang)
p_Styblinksi_Tang + labs(title="Styblinksi-Tang test surface.")
ggsave(filename = sprintf("../Render/%s.svg",obj_Styblinski_Tang$file_prefix),
       width=150, height=90, units="mm", dpi=600)

p_Rastrigin <- generate_surface(obj_Rastrigin)
p_Rastrigin + labs(title="Rastrigin test surface.")
ggsave(filename = sprintf("../Render/%s.svg",obj_Rastrigin$file_prefix),
       width=150, height=90, units="mm", dpi=600)

# trace_rosen_gd  <- optim_main(OBJ=obj_Rosenbrock, UP=up_gradient_descent, 
#                               nu=0.0002)
# trace_rosen_nr  <- optim_main(OBJ=obj_Rosenbrock, UP=up_newton_raphson) 
# trace_rosen_prs <- optim_main(OBJ=obj_Rosenbrock, UP=up_pure_rs) 
# trace_rosen_rs  <- optim_main(OBJ=obj_Rosenbrock, UP=up_random_search) 
# 
# trace_rast_gd  <- optim_main(OBJ=obj_Rastrigin, UP=up_gradient_descent)  
# trace_rast_nr  <- optim_main(OBJ=obj_Rastrigin, UP=up_newton_raphson) 
# trace_rast_prs <- optim_main(OBJ=obj_Rastrigin, UP=up_pure_rs) 
# trace_rast_rs  <- optim_main(OBJ=obj_Rastrigin, UP=up_random_search)
# 
# trace_st_gd  <- optim_main(OBJ=obj_Styblinski_Tang, UP=up_gradient_descent)
# trace_st_nr  <- optim_main(OBJ=obj_Styblinski_Tang, UP=up_newton_raphson)
# trace_st_prs <- optim_main(OBJ=obj_Styblinski_Tang, UP=up_pure_rs) 
# trace_st_rs  <- optim_main(OBJ=obj_Styblinski_Tang, UP=up_random_search) 
# 
# trace_rosen_gd
# trace_rosen_nr
# trace_rosen_prs
# trace_rosen_rs
# trace_rast_gd
# trace_rast_nr
# trace_rast_prs
# trace_rast_rs
# trace_st_gd
# trace_st_nr
# trace_st_prs
# trace_st_rs

# list(
#   list(OBJ=obj_Rosenbrock,      UP=up_gradient_descent, nu=0.0002),
#   list(OBJ=obj_Rosenbrock,      UP=up_newton_raphson,   nu=NULL),
#   list(OBJ=obj_Rosenbrock,      UP=up_pure_rs,          nu=NULL),
#   list(OBJ=obj_Rosenbrock,      UP=up_random_search,    nu=NULL),
#   list(OBJ=obj_Rastrigin,       UP=up_gradient_descent, nu=NULL),
#   list(OBJ=obj_Rastrigin,       UP=up_newton_raphson,   nu=NULL),
#   list(OBJ=obj_Rastrigin,       UP=up_pure_rs,          nu=NULL),
#   list(OBJ=obj_Rastrigin,       UP=up_random_search,    nu=NULL),
#   list(OBJ=obj_Styblinski_Tang, UP=up_gradient_descent, nu=NULL),
#   list(OBJ=obj_Styblinski_Tang, UP=up_newton_raphson,   nu=NULL),
#   list(OBJ=obj_Styblinski_Tang, UP=up_pure_rs,          nu=NULL),
#   list(OBJ=obj_Styblinski_Tang, UP=up_random_search,    nu=NULL)
# )


trace_list <- list(
  OBJ = list(obj_Rosenbrock, obj_Rosenbrock, obj_Rosenbrock, obj_Rosenbrock,
             obj_Rastrigin, obj_Rastrigin, obj_Rastrigin, obj_Rastrigin,
             obj_Styblinski_Tang, obj_Styblinski_Tang, obj_Styblinski_Tang,
             obj_Styblinski_Tang),
  UP  = list(up_gradient_descent, up_newton_raphson, up_pure_rs, 
             up_random_search, up_gradient_descent, up_newton_raphson, 
             up_pure_rs, up_random_search, up_gradient_descent, 
             up_newton_raphson, up_pure_rs, up_random_search),
  nu  = list(0.0002,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
  ) |>
  pmap(.f=optim_main) |> progressify() |> futurize(seed=TRUE)


# cannot train GA stepwise. The GA must be fitted first and then tested.

#minimisation problem must be changed to maximisation problem
OBJ <- obj_Rosenbrock

ga_obj <- function(x,OBJ=OBJ){
  -OBJ$func(x)
}

ga_model <- ga(type="real-valued",
   fitness= ga_obj, 
   lower = rep(OBJ$bounds[1],length.out=DIM), #length of lower/upper must = DIM
   upper = rep(OBJ$bounds[2],length.out=DIM),
   keepBest = TRUE,
   monitor = plot
)

ga_model@population

# now, GA doesn't give us every single value tested.
# but it does give us the best value tested for each population iteration.

# we define a monitoring function so we can attempt to note all populations generated.
# let us do this with tidyverse
plot_monitor <- function(model) {
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
         lower = rep(OBJ$bounds[1],times=2), 
         upper = rep(OBJ$bounds[2],times=2), 
         popSize = popSize,
         maxiter = iterlim/5, 
         monitor = plot_monitor)
summary(GA)



NOW <- now() #snapshot to label file generated

# custom monitor to print population to file
print_monitor <- function(ga_model){
  # print all populations to file, includes the iteration it was generated in
  bind_cols(Pop_Iteration=ga_model@iter,as_tibble(ga_model@population)) |> 
    write_csv(file=paste0("./population_generated/",NOW,".csv"),append=TRUE)
  print(paste("Population iteration:", ga_model@iter))
}
GA <- ga(type = "real-valued",
         fitness=run_batch,
         drones=d,
         obst=obst,
         lower = rep(-10, npars), upper = rep(10, npars),
         maxiter = GA_iterlim,
         monitor = print_monitor,
         keepBest = TRUE)
summary(GA)

fcalls <- read_csv(paste0("./population_generated/",NOW,".csv"),
                   col_names=FALSE)





# Shvet Maharaj 2026-10-26