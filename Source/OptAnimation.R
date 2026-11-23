################################################################################
#                        setup + objective functions                           #

library(tidyverse,attach.required = TRUE) # data tidying, visualisation + more tools
library(gganimate) # animates graphs, illustrating behaviour of update procedures
library(gifski)    # provides gif support for gganimate
library(GA)        # GA, Scrucca 2017
library(numDeriv)  # for derivative calculations for gradient descent & Newton-Raphson 
library(patchwork) # for merging multiple plots
library(futurize)  # simple parallel processing support
library(progressify) # simple progress reporting
library(tictoc)    # quick benchmarking
#library(ggnewscale)# allows for multiple colour scales on one ggplot

plan(multisession)
handlers(global=TRUE)

# Major effect on p_obj runtime and size of plot files. (1 is added to resolution to centralise).
# Exponential growth with DIM (resolution^DIM)
# resolution+1 = Number of values to be evaluated per dimension.
resolution <- 500
DIM <- 2               # Plots are defined for 2 dimensions.
iterlim <- 100 # amount of total iterations (function calls) allowed.
trials <- 1000  # amount of repeated trials for each objective and optimiser pairing

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
  #min_point <- c(0,0)
  return(f)
} # 2 dim
test_parabola <- function(x){
  #n <- length(x)
  f <- t(x) %*% x
  #min_point <- rep(0,times=n)
  return(f)
}
test_eggholder <- function(x){
  x1 <- x[1]*200 #rescaling
  x2 <- x[2]*200
  f  <- -(x2+47)*sin(sqrt(abs(x1/2+(x2+47))))-x1*sin(sqrt(abs(x1-(x2+47))))
  #min_point <- c(512,404.2319)
  return(f)
} # 2 dimensions only
test_Rosenbrock <- function(x){
  n <- length(x)
  f <- sum(100*(x[-1] -(x[-n])^2)^2 + (1-x[-n])^2)
  #min_point <- rep(1,times=n)
  return(f)
} # dim > 1
test_Styblinski_Tang <- function (x){
  #n <- length(x)
  f <- 0.5 * sum( x^4 -16*x^2 +5*x )
  #min_point <- rep(-2.903534,times=n)
  return(f +78.33234) #ensure never negative.
}


# optimisation function options
# these functions take in the current position (x), 
# returning the next x value based on their algorithm.

step_gradient_descent  <- function(x,nu,OBJ){
  #OBJ: only func used
  next_x <- x-nu*grad(OBJ$func,x,method="simple")
  return(next_x)
}
# step_random_search must use MinX, not simply X when selecting update procedure
# step_random_search (FSSRS) has a normalised (1) step size. 
# nu acts as a multiplicative factor on this. 
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
                                  bounds, transform, breaks) {
  structure(
    list(func        = func,
         plot_name   = plot_name,
         file_prefix = file_prefix,
         min         = min,
         bounds      = bounds,
         transform   = transform,
         breaks      = breaks),
    class = "objective_surface"
  )
}
# func: name of test function.
# plot_name: names on ggplot graphs. "...on `plot_name` surface."
# file_prefix: rendered files saved with filename: `file_name`-`optim_name`.pdf
# min: minimum point(s) of test function.
# bounds: 1-dimensional bounds of the entire search space. This is recycled for both dimensions.
# transform, breaks: passed to ggplot scale functions to ensure correct scaling.

# breaks were set by plotting the surfaces with "waiver()" at first, then manually filling the values given in the legend. This is how the contour lines match the ticks on the guide_colourbar().
is_objective_surface <- function(x){
  class(x) == "objective_surface"
}

obj_Rastrigin <- new_objective_surface(
  func        = test_Rastrigin,
  plot_name   = "a Rastrigin",
  file_prefix = "rastrigin",
  min         = c(0,0),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = c(0,20,40,60,80)
)

obj_Ackley <- new_objective_surface(
  func        = test_Ackley,
  plot_name   = "an Ackley",
  file_prefix = "ackley",
  min         = c(0,0),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = waiver()
)

obj_parabola <- new_objective_surface(
  func        = test_parabola,
  plot_name   = "a parabolic",
  file_prefix = "para",
  min         = c(0,0),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = c(0,10,20,30,40,50)
)

obj_eggholder <- new_objective_surface(
  func        = test_eggholder,
  plot_name   = "an eggholder",
  file_prefix = "egg",
  min         = c(512,404.2319),
  bounds      = c(-550,550),
  transform   = "identity",
  breaks      = waiver()
)

obj_Rosenbrock <- new_objective_surface(
  func        = test_Rosenbrock,
  plot_name   = "a Rosenbrock",
  file_prefix = "rosen",
  min         = c(1,1),
  bounds      = c(-5,5),
  transform   = "log1p",
  breaks      = c(1,10,100,1000,10000)
)

obj_Styblinski_Tang <- new_objective_surface(
  func        = test_Styblinski_Tang,
  plot_name   = "a Styblinski-Tang",
  file_prefix = "styb-tang",
  min         = c(-2.903534,-2.903534),
  bounds      = c(-5,5),
  transform   = "identity",
  breaks      = c(100,200,300)
)

# Update Procedures

# S3 class constructor
new_update_procedure <- function(func, use_min_x, optim_name, plot_name, nu, is_ga) {
  structure(
    list(func       = func,
         use_min_x  = use_min_x,
         optim_name = optim_name,
         plot_name  = plot_name,
         nu         = nu,
         is_ga      = is_ga),
    class = "update_procedure"
  )
}
#func:        optimiser function
#use_min_x:   boolean. selects `min_x` or `x` as parameter for `func`.
#optim_name:  identifier used in filenames for this optimiser.
#plot_name:   name used in ggplot graphics for this optimiser.
#nu:          default nu value for this optimiser. Can be overwritten.
#is_ga:       if TRUE, `nu` & `use_min_x` are unused, and `func` holds the crossover function.

is_update_procedure <- function(x){
  class(x) == "update_procedure"
}

up_gradient_descent <- new_update_procedure(
  func       = step_gradient_descent,
  use_min_x  = FALSE,
  optim_name = "grad-desc",
  plot_name  = "gradient descent",
  nu         = 0.2,
  is_ga      = FALSE
)

up_random_search <- new_update_procedure(
  func       = step_random_search,
  use_min_x  = TRUE,
  optim_name = "rand-sear",
  plot_name  = "FSSRS",
  nu         = 0.2,
  is_ga      = FALSE
)

up_newton_raphson <- new_update_procedure(
  func       = step_newton_raphson,
  use_min_x  = FALSE,
  optim_name = "newt-raph",
  plot_name  = "Newton-Raphson",
  nu         = 0.5,
  is_ga      = FALSE
)

up_pure_rs <- new_update_procedure(
  func       = step_pure_rs,
  use_min_x  = FALSE,
  optim_name = "pure-rs",
  plot_name  = "PRS",
  nu         = NULL ,
  is_ga      = FALSE
)

up_ga_none <- new_update_procedure(
  func       = NULL,
  use_min_x  = FALSE,
  optim_name = "ga-none",
  plot_name  = "GA with no crossover",
  nu         = NULL,
  is_ga      = TRUE
)

up_ga_wa <- new_update_procedure(
  func       = gareal_waCrossover,
  use_min_x  = FALSE,
  optim_name = "ga-wa",
  plot_name  = "GA with whole arithmetic crossover",
  nu         = NULL,
  is_ga      = TRUE
)

up_ga_la <- new_update_procedure(
  func       = gareal_laCrossover,
  use_min_x  = FALSE,
  optim_name = "ga-la",
  plot_name  = "GA with local arithmetic crossover",
  nu         = NULL,
  is_ga      = TRUE
)

up_ga_blx <- new_update_procedure(
  func       = gareal_blxCrossover,
  use_min_x  = FALSE,
  optim_name = "ga-blx",
  plot_name  = "GA with blend crossover",
  nu         = NULL,
  is_ga      = TRUE
)

up_ga_laplace <- new_update_procedure(
  func       = gareal_laplaceCrossover,
  use_min_x  = FALSE,
  optim_name = "ga-laplace",
  plot_name  = "GA with Laplace crossover",
  nu         = NULL,
  is_ga      = TRUE
)


# Ackley and eggholder not used in dissertation.
OBJECTIVE_SURFACES <- list(
  obj_parabola,
  #obj_Ackley,
  obj_Rastrigin,
  #obj_eggholder,
  obj_Rosenbrock,
  obj_Styblinski_Tang
  )

UPDATE_PROCEDURES <- list(
  up_gradient_descent,
  up_newton_raphson,
  up_pure_rs,
  up_random_search,
  up_ga_none,
  up_ga_wa,
  up_ga_la,
  up_ga_blx,
  up_ga_laplace
  )

#-#-# #-#-# #-#-# #-#-# #          OPTIMISATION        # #-#-# #-#-# #-#-# #-#-#
 
# generates test surface
generate_surface <- function(OBJ,.res=resolution){
  base_theme <- theme_minimal() 
  theme_set(base_theme)
  
  trace_names <- c("X1","X2","Y")
  
  xseq <- seq(OBJ$bounds[1],OBJ$bounds[2],length.out=.res+1) #symmetric. 
  # geom_tile is based on the centre of the tile.
  
  # ggdata is the surface information of the objective function. 
  ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                        .name_repair=~make.names(1:DIM))|>
    mutate(Y=pmap_vec(.l=pick(everything()),
                      .f=~ OBJ$func(c(...))))
  
  names(ggdata) <- trace_names
  
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
generate_surface_dark <- function(OBJ,.res=resolution){
  base_theme <- theme_minimal() 
  theme_set(base_theme)
  
  trace_names <- c("X1","X2","Y")
  
  xseq <- seq(OBJ$bounds[1],OBJ$bounds[2],length.out=.res+1) #symmetric. 
  # geom_tile is based on the centre of the tile.
  
  # ggdata is the surface information of the objective function. 
  ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                        .name_repair=~make.names(1:DIM))|>
    mutate(Y=pmap_vec(.l=pick(everything()),
                      .f=~ OBJ$func(c(...))))
  
  names(ggdata) <- trace_names
  
  #plot the objective surface
  p_obj <- ggplot()+
    geom_raster(data=ggdata,aes(x=X1,y=X2),fill='black',alpha=0.8)+
    geom_raster(data=ggdata,aes(x=X1,y=X2,fill=Y),alpha=0.5) + 
    geom_contour(data=ggdata,aes(x=X1,y=X2,z=Y),breaks=OBJ$breaks,
                 colour="white", linewidth = 0.2, alpha=0.4, linetype = 1)+
    geom_point(data=as_tibble(t(OBJ$min),
                              .name_repair=~make.names(1:DIM)),
               aes(x=X1,y=X2),
               colour='red',shape=4,alpha=0.5,show.legend = TRUE)+
    scale_fill_viridis_c(option="E",
                         transform=OBJ$transform,breaks=OBJ$breaks
    ) + coord_fixed() + 
    xlim(OBJ$bounds*1.01) +ylim(OBJ$bounds*1.01) # 1% expansion factor
  
  return(p_obj)
}

# runs optimisation. generates trace. generates plots. generates animation.
optim_main <- function(OBJ, UP, nu=NULL, trials=1,
                      .iterlim=iterlim, .res=resolution){
  
  DIM <- 2
  
  if(!is.null(nu)){
    UP$nu <- nu
  }
  
  base_theme <- theme_minimal() 
  theme_set(base_theme)
  
  xseq <- seq(OBJ$bounds[1],OBJ$bounds[2],length.out=.res+1) #symmetric. 
  # geom_tile is based on the centre of the tile.
  
  # ggdata is the surface information of the objective function. 
  ggdata <- expand_grid(!!! rep(list(xseq),DIM),
                        .name_repair=~make.names(1:DIM))|>
    mutate(Y=pmap_vec(.l=pick(everything()),
                      .f=~ OBJ$func(c(...))))
  
  trace_names <- c("X1","X2","Y")
  names(ggdata) <- trace_names
  
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
  
  # generate trace
  
  update_opt <- (
    if (UP$use_min_x) 
      function(x, min_x, nu=nu_global) UP$func(min_x, nu=UP$nu, OBJ) 
    else 
      function(x, min_x, nu=nu_global) UP$func(    x, nu=UP$nu, OBJ)
  )
  
  start <- runif(n=DIM,OBJ$bounds[1],OBJ$bounds[2])
  names(start)=names(ggdata)[1:DIM]
  
  opt_trace <- matrix(NA, nrow=.iterlim, ncol=7,
                      dimnames= list(NULL,
                                     c('Iteration',
                                       names(ggdata),
                                       paste0('Min',names(ggdata)))
                      )
  )
  
  #                              run optimisation
  
  xi <- min_x <- start
  yi <- min_y <- OBJ$func(min_x)
  
  opt_trace[1,] <- c(1,xi,yi,min_x,min_y)
  
  #message(sprintf("Beginning %s on %s surface.", UP$plot_name, OBJ$plot_name))
  
  #foreach(t = 1:trials) %do%{
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
    geom_point(data=trace_segments,aes(x=X1,y=X2,group=Iteration),
               colour="white",alpha=0.3,shape=4)+
    geom_point(data=trace_segments,aes(x=MinX1,y=MinX2,colour=Iteration))+
    scale_colour_viridis_c(option="H",begin=0.3,end=0.9) +
    geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
               colour='aquamarine',shape=4)+
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
    plot_layout(heights=unit(c(25,1),c('mm','null')))
  
  ggsave(plot=p_full,
         filename = sprintf("../Render/%s-%s.pdf",
                            OBJ$file_prefix, UP$optim_name),
         width=180, height=150, units="mm", dpi=400)
  
  #message(sprintf("Animating %s on %s surface.", UP$plot_name, OBJ$plot_name))
  
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
  
  #message(sprintf("Completed %s on %s surface.", UP$plot_name, OBJ$plot_name))
  
  return(opt_trace)
}

#uses `quantiles` to return data in the format desired by stat_summary
quartiles <- function(x){
  bind_cols(ymin=quantile(x,probs=0.25),ymax=quantile(x,probs=0.75))
}

# S3 class constructor
new_trace <- function(trace,up,obj,trials,iterations){
  trace <- structure(
    list("trace"=      trace,
         "up"=         up,
         "obj"=        obj,
         "trials"=     trials,
         "iterations"= iterations),
    class = "trace")
}
# trace:      tibble containing data from all trials.
# up:         `update_procedure` object of the update procedure used in the trials.
# obj:        `objective_surface` object of the objective surface tested in the trials.
# trials:     number of trials completed.
# iterations: number of iterations per trial.

is_trace <- function(x){
  class(x) == "trace"
}

# outputs `trace` object
# obj:          `objective_surface` object to run trials on.
# up:           `update_procedure` object to run trials on.
# nu:           (optional) value to overwrite `up$nu`.
# trials:       number of trials to run
# .iterlim:     iterations to complete.
generate_trace_base <- function(obj, up, nu=NULL, trials=1,
                           .iterlim=iterlim){
  DIM <- 2
  
  stopifnot("`up` must be an `update_procedure` object." =
              is_update_procedure(up),
            "`obj` must be an `objective_surface` object." =
              (is.null(obj)||is_objective_surface(obj))
  )
  
  if(!is.null(nu)){
    up$nu <- nu
  }
  
  update_opt <- (
    if (up$use_min_x) 
      function(x, min_x, nu=nu_global) up$func(min_x, nu=up$nu, obj) 
    else 
      function(x, min_x, nu=nu_global) up$func(    x, nu=up$nu, obj)
  )
  
  trace_names <- c("X1","X2","Y")
  
  trace <- 
    array(NA,
          dim = c(.iterlim,6,trials),
          dimnames = list(iteration = 1:.iterlim,
                          trace     = c(trace_names,paste0('Min',trace_names)),
                          trial     = 1:trials)
    )
  
  #                              run optimisation
  foreach(t = 1:trials) %do% {
    
    start <- runif(n=DIM,obj$bounds[1],obj$bounds[2])
    names(start) <- trace_names[1:2]
    
    xi <- min_x <- start
    yi <- min_y <- obj$func(min_x)
    
    trace[1,,t] <- c(xi,yi,min_x,min_y)
    
    for(i in 2:.iterlim){
      
      xi <- update_opt(xi,min_x) 
      yi <- obj$func(xi) #objective being minimised
      
      if(yi<min_y){
        min_y <- yi
        min_x <- xi
      }
      
      trace[i,,t] <- c(xi,yi,min_x,min_y) 
      
    }
  }
  
  traces <- trace |> data.table::as.data.table() |>
    as_tibble() |> 
    pivot_wider(names_from=trace,values_from=value) |> 
    mutate(trial =     as.integer(trial),
           iteration = as.integer(iteration)) |> 
    arrange(trial, iteration) |> 
    relocate(trial, iteration, X1, X2, Y)
  
  return(
    new_trace(trace=traces,
              up=up,
              obj=obj,
              trials=trials,
              iterations=.iterlim)
  )
}

# clears trace file.
init_print_monitor <- function(){
  unlink(x="../Populations/ga_trace.csv")
}

# custom monitor to print population to file
print_monitor <- function(ga_model){
  # print all populations to file, includes the iteration it was generated in
  bind_cols(Pop_Iteration=ga_model@iter,as_tibble(ga_model@population)) |> 
    write_csv(file="../Populations/ga_trace.csv",append=TRUE)
  #print(paste("Population iteration:", ga_model@iter))
}

# plot_monitor <- function(model,.surface=p_obj) {
#   best <- which.max(model@fitness)
#   pop <- model@population |> as_tibble() |> bind_cols(Iteration=model@iter)
#   plot(.surface+
#          geom_point(data=pop,aes(x=V1,y=V2,group=Iteration),
#                     colour="white",shape=4,alpha=0.2)+
#          labs(subtitle=paste("iteration =", model@iter))+
#          geom_point(data=pop[best,],aes(x=V1,y=V2,colour=Iteration))+
#          scale_colour_viridis_c(option="H",begin=0.3,end=0.9)
#     )
# }


# cannot train GA stepwise. The GA must be fitted first and then tested.

# outputs `trace` object
# obj:          `objective_surface` object to run trials on.
# up:           `update_procedure` object to run trials on.
# trials:       size of initial population.
# .iterlim:     maximum iterations (generations).
# pop:          population size. Ensure `trials` is divisible by `pop`.
generate_trace_ga <- function(obj, up, trials=100,
                              .iterlim=iterlim, pop=20){
  
  stopifnot("`up` must be an `update_procedure` object." =
              is_update_procedure(up),
            "`obj` must be an `objective_surface` object." =
              (is.null(obj)||is_objective_surface(obj))
  )

  trials <- trials/pop

  # GA attempts to maximise fitness. We must negate the objective function to minimise.
  ga_obj <- function(x){
    -obj$func(x)
  }
  
  crossover <- up$func
  
  
  # in order to set seed once at the top of the file, yet also remain reproducable, we shall generate a seed here
  #s <- sample.int(1000000,1)
  #message(paste("Seed:",s)) # for redundancy / to double check. 
  
  if(is.null(crossover)){
    pcross <- 0
    crossover <- gareal_waCrossover #won't do anything, but one needs to be set
  } else {
    pcross <- 0.8 #the default in ga
  }
  
  # gareal_spCrossover !?
  # gareal_waCrossover
  # gareal_laCrossover
  # gareal_blxCrossover
  # gareal_laplaceCrossover
  
  
  #skipping ga seed option. Unsure if it's required when seed is handled by future.
  
  trace <- foreach(t = 1:trials,.combine=bind_rows) %do% { 
    init_print_monitor()
    fitted <- ga(type = "real-valued",
                 fitness = ga_obj,
                 crossover = crossover,
                 pcrossover = pcross,
                 lower = rep(obj$bounds[1],times=2), 
                 upper = rep(obj$bounds[2],times=2), 
                 popSize = pop,
                 maxiter = .iterlim, 
                 monitor = print_monitor)
    
    trace_trial <- read_csv(file="../Populations/ga_trace.csv",col_names=c("iteration", "X1","X2")) |> suppressMessages() |> 
      bind_cols(trial = t) |> relocate(trial)
  }
  
  y <- select(trace, X1,X2) |> apply(MARGIN=1,FUN=obj$func)
  
  traces <- trace |> bind_cols(Y=y)
  
  trace_min <- traces |>  slice_min(order_by=Y,by=c(trial,iteration),with_ties = FALSE) |> rename(MinX1=X1,MinX2=X2,MinY=Y)
  
  full_traces <- traces |> left_join(trace_min,by=c("trial","iteration"))
  
  return(
    new_trace(trace=full_traces,
              up=up,
              obj=obj,
              trials=trials,
              iterations=.iterlim)
  )
}

# wrapper for `generate_trace_base` and `generate_trace_ga`.
generate_trace <- function(obj, up, nu=NULL, trials=100,
                           .iterlim=iterlim,pop=NULL){
  
  if(up$is_ga){
    if(is.null(pop)) {
      pop <- 20
    }
    trace <-generate_trace_ga(obj=obj,
                              up=up,
                              trials=trials,
                              .iterlim=.iterlim,
                              pop=pop)
  } else {
    trace <-generate_trace_base(obj, 
                                up, 
                                nu=nu, 
                                trials=trials,
                                .iterlim=.iterlim)
  }
  # ensure only non zero values 
  trace$trace <- trace$trace |> 
    mutate(Y= replace_when(   Y,    Y==0 ~ .Machine$double.eps),
        MinY= replace_when(MinY, MinY==0 ~ .Machine$double.eps))
  
  return(trace)
}

# creates tibble for use in the trace ggplot.
# Returns two components:
#   segments: the data split into segments for each facet.
#   names:    the names used for each facet.
form_faceting_segments <- function(x,segments,trial_plotted){
  
  n_segments <- length(segments)
  
  store       <- vector(mode="list",length=n_segments)
  end_store   <- vector(mode="list",length=n_segments)
  facet_names <- vector(mode="character",length=n_segments)
  
  for(i in 1:n_segments){
    
    store[[i]]<- x |>
      filter(trial==trial_plotted, iteration==1:segments[i]) |>
      bind_cols(facet=i)
    
    end_store[[i]] <- x |> 
      filter(iteration==segments[i]) |>
      bind_cols(facet=i)
    
    if(segments[i]==1){
      facet_names[i] <- "1 Step"
    } else {
      facet_names[i] <- sprintf("%d Steps",segments[i])
    }
    
  }
  
  trace_segments <- bind_rows(store)
  names(facet_names) <- 1:n_segments
  segment_ends <- bind_rows(end_store)
  
  return(list(segments=trace_segments,
              names=   facet_names,
              ends =   segment_ends))
}

# trace:          `trace` object from `generate_trace`
# surface:        (optional) a ggplot object containing the surface on which the trace was calculated.
# res:            the resolution desired for the plotted surface. Necessary if `surface` is NULL.
# segments:       the segments for the trace plot facets. 
# trial_plotted:  the number corresponding to the trial for which to trace. Meaningless for GA update procedures.
# animate:        whether to produce a .gif animation of the trace.
plot_trace <- function(trace,
                       res=resolution,
                       surface=NULL,
                       segments= c(1,5,10,25,50,100),
                       trial_plotted=NULL,
                       animate=FALSE){
  
  # is.integer does not check the value of a variable, simply its class.
  # x<-1 would not be considered an integer, but x<-1L would.
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol 
  # taken from `is.integer` examples.
  
  # input validation
  stopifnot("`res` must be a positive integer."=
              is.wholenumber(res)&&res>0,
            "`trace` must be a `trace` object." =
              is_trace(trace)
  )
  
  obj        <- trace$obj
  up         <- trace$up
  trace_data <- trace$trace
  
  maxiter <- trace_data |> summarise(maxiter = max(iteration)) |> pull(maxiter)
  if (max(segments) != maxiter) {
    warning("Segment max does not match iteration max.")
  }
  
  if(is.null(trial_plotted)){
    trace_range <- trace_data |> 
      summarise(max=max(MinY),min=min(MinY),.by=trial) |> 
      mutate(range=max-min) 
    trial_plotted <- trace_range |> 
      filter(range==max(range)) |> pull(trial)
    
    trial_plotted <- trace_data |> 
      slice_min(order_by=MinY,with_ties = FALSE) |> pull(trial)
  }
  
  n_trials <- max(trace_data$trial)
  
  # factors affected by GA / not.
  if(up$is_ga){
    n_pop <- trace_data |> summarise(n=n(),.by=c(iteration,trial)) |> 
      slice(1) |> pull(n) 
    perf_subtitle  <- sprintf("Population size: %d.",n_pop)
    perf_title     <- sprintf("Behaviour of %s on %s surface.",
                              up$plot_name, obj$plot_name)
  } else {
    perf_subtitle  <- sprintf("Total trials: %d.",n_trials)
    perf_title     <- sprintf("Behaviour of the %s optimiser on %s surface.",
                              up$plot_name, obj$plot_name)
  }
  
  trace_subtitle <- sprintf("Trial traced: %d.",trial_plotted)
  
  seg <- form_faceting_segments(x=trace_data,segments,trial_plotted) 
  # trace_segments: how much of the trace is used in each facet
  trace_segments <- seg$segments
  facet_names <-    seg$names 
  segment_ends <-   seg$ends
  
  if(is.null(surface)){
    surface <- generate_surface_dark(OBJ=obj,.res=res)
  }
  
  #                              plot and animate
  
  start <- filter(trace_data,iteration==1) |> bind_cols(facet=1) # show coverage of start points on first facet.
  
  p_static <- surface+
    # geom_point(data=trace_segments,aes(x=X1,y=X2),
    #            colour="white",alpha=0.2,shape=4)+
    scale_colour_viridis_c(option="H",begin=0.3,end=0.9)+
    # geom_point(data=start, aes(x=X1,y=X2),
    #            colour='aquamarine',shape=4,alpha=0.3)+
    geom_point(data=segment_ends,aes(x=MinX1,y=MinX2,colour=iteration),
               alpha=0.2,size=0.2)+
    geom_point(data=trace_segments,aes(x=MinX1,y=MinX2,colour=iteration),
               shape=21,size=2)+
    facet_wrap(facets=~facet,nrow=2,
               labeller = as_labeller(facet_names))+
    labs(colour="Iteration",
         subtitle=trace_subtitle)
  
  p_performance_simple <- ggplot(data=trace_data, aes(x=iteration, y=MinY)) +
    stat_summary(geom="ribbon",
                 fun.min=min,fun.max=max,
                 fill = 'black',alpha=0.5)+
    stat_summary(geom="ribbon",
                 fun.min=min,fun.max=max,
                 aes(fill = iteration),alpha=0.5,show.legend = FALSE)+
    stat_summary(geom="line",fun=min, colour='darkblue', alpha=0.5,aes(linetype="min"))+
    scale_y_log10()+
    scale_fill_viridis_c(option="H",begin=0.3,end=0.9) +
    scale_colour_viridis_c(option="H",begin=0.3,end=0.9) +
    labs(y="Y", x="Iteration",
         title= perf_title,
         subtitle= perf_subtitle)+ 
    #theme_bw() +
    guides(linetype=guide_legend(
      title="",
      theme=theme(legend.key=element_rect(fill="grey",colour="#00000000"))))+
    scale_linetype_manual(values=c("min"=1))
  
  p_performance <- ggplot(data=trace_data, aes(x=iteration, y=MinY)) +
    stat_summary(geom="ribbon",
                 fun.min=min,fun.max=max,
                 fill = 'black',alpha=0.5)+
    stat_summary(geom="ribbon",
                 fun.min=min,fun.max=max,
                 aes(fill = iteration),alpha=0.5,show.legend = FALSE)+
    stat_summary(geom="ribbon",
                 fun.data=quartiles,aes(fill = iteration),show.legend = FALSE)+
    stat_summary(geom="line",fun=min, colour="darkblue",alpha=0.5,aes(linetype='min'))+
    stat_summary(geom="line",fun=max, colour="darkblue",alpha=0.5,aes(linetype='max'))+
    stat_summary(geom="line",fun=quantile,
                 fun.args=list(probs=0.25),colour="white",aes(linetype='q1'))+
    stat_summary(geom="line",fun=quantile,
                 fun.args=list(probs=0.75),colour="white",aes(linetype='q3'))+
    stat_summary(geom="line",fun=median, colour="darkblue",aes(linetype='median'))+
    #geom_point(alpha=0.3,size=0.5,shape=1)+
    scale_y_log10()+
    scale_fill_viridis_c(option="H",begin=0.3,end=0.9) +
    labs(y="Y", x="Iteration",
         title   = perf_title,
         subtitle= perf_subtitle)+ 
    #theme_bw()+
    guides(linetype=guide_legend(
      title="",
      theme=theme(legend.key = element_rect(fill="grey",colour="#00000000"))))+
    scale_linetype_manual(breaks=c("max","q3","median","q1","min"),
                          values=c("max"   =1,
                                   "q3"    =3,
                                   "median"=2,
                                   "q1"    =3,
                                   "min"   =1))
  
  p_full_simple <- p_performance_simple/p_static +
    plot_layout(heights=unit(c(25,1),c('mm','null')))
  
  p_full <- p_performance/p_static +
    plot_layout(heights=unit(c(25,1),c('mm','null')))
  
  ggsave(plot=p_full_simple,
         filename = sprintf("../Render/%s-%s-simple.pdf",
                            obj$file_prefix, up$optim_name),
         width=180, height=150, units="mm", dpi=400)
  
  ggsave(plot=p_full,
         filename = sprintf("../Render/%s-%s.pdf",
                            obj$file_prefix, up$optim_name),
         width=180, height=150, units="mm", dpi=400)
  
  if(!animate){
    return(list(plot_simple= p_full_simple,
                plot=        p_full,
                perf_simple= p_performance_simple,
                perf=        p_performance))
  }
  
  p_anim <-
    animate(
      surface+
        geom_point(data=as_tibble(t(start)), aes(x=X1,y=X2),
                   colour='royalblue',shape=4)+
        geom_point(data=trace,aes(x=X1,y=X2,group=iteration),
                   colour="white",shape=4)+
        geom_point(data=trace,aes(x=MinX1,y=MinX2,colour=iteration))+
        scale_colour_viridis_d(option="H",begin=0.3,end=0.9)+
        transition_states(states=iteration,wrap=FALSE)+
        shadow_wake(wake_length=0.3)+
        labs(subtitle = "Iteration: {previous_state}.")
    )
  
  anim_save(nframes= nrow(trace),
            filename=sprintf("%s-%s.gif", obj$file_prefix, up$optim_name),
            path="../Render",
            animation=p_anim)
  
  return(list(plot_simple= p_full_simple,
              plot=        p_full,
              perf_simple= p_performance_simple,
              perf=        p_performance,
              anim=        p_anim))
}


# plot_trace <- function(trace,
#                        res=resolution,
#                        surface=NULL,
#                        segments= c(1,2,5,10,50,100),
#                        trial_plotted=NULL,
#                        animate=FALSE){
#   
#   if(trace$up$is_ga) {
#     return(plot_trace_ga(trace,res,surface,segments,trial_plotted,animate))
#   }
#   return(plot_trace_base(trace,res,surface,segments,trial_plotted,animate))
# }

#-#-# #-#-# #-#-# #-#-# #          RUN        # #-#-# #-#-# #-#-# #-#-#

set.seed(111)

#-#-# generate all surfaces into pdf images #-#-#

p_para <- generate_surface(obj_parabola)
p_para + labs(title="Parabolic test surface.")
ggsave(filename = sprintf("../Render/%s.pdf",obj_parabola$file_prefix),
       width=150, height=90, units="mm", dpi=400)

p_Rosenbrock <- generate_surface(obj_Rosenbrock)
p_Rosenbrock + labs(title="Rosenbrock test surface.")
ggsave(filename = sprintf("../Render/%s.pdf",obj_Rosenbrock$file_prefix),
       width=150, height=90, units="mm", dpi=400)

p_Styblinski_Tang <- generate_surface(obj_Styblinski_Tang)
p_Styblinski_Tang + labs(title="Styblinski-Tang test surface.")
ggsave(filename = sprintf("../Render/%s.pdf",obj_Styblinski_Tang$file_prefix),
       width=150, height=90, units="mm", dpi=400)

p_Rastrigin <- generate_surface(obj_Rastrigin)
p_Rastrigin + labs(title="Rastrigin test surface.")
ggsave(filename = sprintf("../Render/%s.pdf",obj_Rastrigin$file_prefix),
       width=150, height=90, units="mm", dpi=400)

N_OBJ <- length(OBJECTIVE_SURFACES)
N_UP  <- length(UPDATE_PROCEDURES)
N_TRACE <- N_UP*N_OBJ

N_GA <- 5
      #OBJ:  parabola, Rastrigin, Rosenbrock, Styblinski-Tang
.NU<- c(list(0.1,      0.04,      0.0001,     0.03, #UP: gradient descent
             1,        1,         1,          1,    #UP: Newton-Raphson
             NULL,     NULL,      NULL,       NULL, #UP: PRS (meaningless)
             0.4,      1,         0.4,        0.3), #UP: FSSRS
        rep(list(NULL),N_GA*N_OBJ)) # NULL repeated for each GA

traces <- foreach(i   = 1:N_TRACE,
                  OBJ = rep(OBJECTIVE_SURFACES,times =N_UP), 
                  UP  = rep(UPDATE_PROCEDURES, each  =N_OBJ),
                  NU  = .NU) %do% 
  { 
    tic(msg= sprintf("%d: %s on %s",i, UP$plot_name,OBJ$plot_name))
    trace <-generate_trace(obj = OBJ,up = UP, nu = NU, trials=trials)
    toc()
    return(trace)
    } |> progressify() # do not parallelise this. ga_trace.csv is not parallelise-friendly.

plots <- foreach(TR = traces) %do% {
  plot_trace(trace=TR)
  } |> progressify() |> futurize() # seed not important

all_traces <- foreach(i = 1:N_TRACE,.combine=bind_rows) %do% {
  bind_cols(trace= traces[[i]]$trace,
            obj  = traces[[i]]$obj$file_prefix,
            up   = traces[[i]]$up$optim_name)
}
minimums <- all_traces |> summarise(Y=quantile(MinY,probs=0.05),.by=c(obj,up)) |> 
  filter(up != "grad-desc",up != "newt-raph") 

all_traces |>
  filter(up != "grad-desc",up != "newt-raph") |>
  ggplot(aes(x=iteration, y=MinY)) +
  # stat_summary(geom="ribbon",
  #              fun.min=min,fun.max=max,
  #              fill = 'black',alpha=0.05)+
  # stat_summary(geom="ribbon",
  #              fun.min=min,fun.max=max,
  #              aes(fill = up),alpha=0)+
  stat_summary(geom="ribbon", aes(fill=up),
               fun.data=quartiles,alpha=0.15)+
  stat_summary(geom="line",fun=median,aes(y=MinY,colour=up))+
  geom_hline(data=minimums,aes(yintercept=Y,colour=up),linetype=2)+
  # stat_summary(geom="line",fun=max,alpha=0.5,aes(linetype='max',colour=up))+
  # stat_summary(geom="line",fun=quantile,
  #              fun.args=list(probs=0.25),aes(linetype='q1',colour=up))+
  # stat_summary(geom="line",fun=quantile,
  #              fun.args=list(probs=0.75), aes(linetype='q3',colour=up))+
  # stat_summary(geom="line",fun=median, aes(linetype='median',colour=up))+
  #geom_point(alpha=0.3,size=0.5,shape=1)+
  scale_y_log10()+
  # scale_fill_viridis_d(option="H",begin=0.3,end=0.9, aesthetics = c("fill","colour")) +
  labs(y="Y", x="Iteration",
       title   = "Optimiser performance on different objective surfaces.",
       subtitle= 
"Ribbon shows interquartile range. 
Dotted line shows 5th percentile of the best Y value.")+ 
  scale_colour_discrete(name="Optimiser",
                        breaks=waiver(),
                        labels = c("ga-blx"    ="GA Blend",
                                   "ga-la"     ="GA Local A",
                                   "ga-wa"     ="GA Whole A",
                                   "ga-laplace"="GA Laplace",
                                   "ga-none"   ="GA None",
                                   "pure-rs"   ="PRS",
                                   "rand-sear" ="FSSRS"),
                        aesthetics = c("fill","colour"))+
  #theme_bw()+
  # guides(linetype=guide_legend(
  #   title="",
  #   theme=theme(legend.key = element_rect(fill="grey",colour="#00000000"))))+
  scale_linetype_manual(breaks=c("max","q3","median","q1","min"),
                        values=c("max"   =1,
                                 "q3"    =3,
                                 "median"=2,
                                 "q1"    =3,
                                 "min"   =1)) +
  facet_wrap(vars(obj),scales="free_y",
             labeller = as_labeller(c(
               "styb-tang"="Styblinski-Tang",
               "rosen"    ="Rosenbrock",
               "rastrigin"="Rastrigin",
               "para"     ="Parabola"
             )))


