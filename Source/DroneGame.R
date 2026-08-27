# game params

library(tidyverse)
library(GA)

n_obst    <- 40 #J 
obst_size <- 0.02 #r (radius) 2*r = square side length
 
delta <- 0.05
iterlim <- 100 #k
GA_iterlim <- 100

bounds <- matrix(c(-1,1,  #X_1
                   -1,1), #X_2
                 nrow=2, byrow=FALSE)
colnames(bounds) <- c('X_1','X_2')
rownames(bounds) <- c('min','max')

start_bounds <- bounds #copy, so there's no need to rename again lol
start_bounds[1:2,1:2] <- matrix(c(-0.8,0.8,
                                  -1, -1
), ncol=2)
obst_bounds <- bounds
obst_bounds[1:2,1:2] <- matrix(c(  -1, 1,
                                   -0.5, 1), ncol=2)


#boundary points
corner_bl  <- bounds[1,] #minimums (bottom left)
corner_tr  <- bounds[2,] #maximums (top right)
corner_br  <- c(corner_tr[1],corner_bl[2]) #bottom right
corner_tl  <- c(corner_bl[1],corner_tr[2]) #top left

start_l    <- c(start_bounds[1,]) #left
start_r    <- c(start_bounds[2,]) #right

lines <- matrix(c(
  corner_bl, start_l,   #bottom_l
  start_l,   start_r,   #start
  start_r,   corner_br, #bottom_r
  corner_bl, corner_tl, #left
  corner_br, corner_tr, #right
  corner_tl, corner_tr  #top/end
), ncol=4, byrow=TRUE) |> as_tibble() 
colnames(lines) <- c('x_1','x_2','x_1.end','x_2.end')

lines <- lines |> add_column(
  "colour"=c("bound",
             "start",
             "bound",
             "bound",
             "bound",
             "finish")
)

generate_field <- function(colours,lines=lines){
  return(
    ggplot() + 
      geom_segment(data=lines,aes(x=x_1,
                                  y=x_2,
                                  xend=x_1.end,
                                  yend=x_2.end,
                                  colour=colour),
                   show.legend = FALSE) +
      scale_colour_manual(values=colours) +
      scale_x_continuous(sec.axis = sec_axis(transform = ~.* 1, 
                                             name = "North")) +
      scale_y_continuous(sec.axis = sec_axis(transform = ~.* 1, 
                                             name = "East")) +
      labs(x= "South", y="West")
  )
}

generate_obstacles <- function(obst, fillcolour='forestgreen'){ 
  return(
    geom_tile(data = obst,
              fill = fillcolour, 
              aes(x=X_1, y=X_2,
                  height=2*Size,
                  width =2*Size
                  #fill  = fillcolour
              )
    )
  )
  # chose to plot trees as squares instead of circles for plotting accuracy -
  # i.e., the plot matches exactly the reality of the game (except that the 
  # corners of the trees do not cause collision). 
  # The size of the trees are exactly what is shown.
}

detect_collision <- function(x,obst){
  if(abs(x[1]) >= bounds[2,1] || x[2] <= bounds[1,2]) #boundary detection
  {
    return(1)
  }
  
  # object detection
  J <- nrow(obst)
  dists <- abs((matrix(1,nrow=J,ncol=1) %*% x) - obst[,-3]) 
  
  col <- (dists > obst[,3]) #since comparisons are faster than float-pt mult
  # tests if outside Manhattan distance (outside the 'square tree'
  # i.e., as represented in the plots)
  
  # now, tests only the trees for which the drone is within the 'square tree'
  # testing the actual circle hitbox :o
  
  o <- !(col[,1]|col[,2])
  
  d <- rowSums(dists[o,,drop=FALSE]**2) #should be faster than matrix 
  # multiplication, then extracting diagonals?
  
  if (any(d<=obst[o,3]))
    return(1)
  
  return(0) # no collisions!
}

eval_game_state <- function(x,obst,step){
  if (step >= iterlim || detect_collision(x, obst))
    return(-1)  
  
  if (x[2] >= bounds[2,2])
    return(1)
  
  return(0)
}

#-----------------------------NEURAL NETWORK------------------------------------

create_IterativeNN <- function(theta, m=c(5,5),p=2,q=2,
                               npars=p*m[1] +m[1]*m[2] +m[2]*q +m[1]+m[2]+q){
  stopifnot(is.numeric(theta) || length(theta)==npars)
  
  # --- activation functions ---
  sig1 <- function(z) tanh(z)    # first m hidden layer
  sig2 <- function(z) tanh(z)    # second m hidden layer
  sig3 <- function(z) tanh(z)    # output layer (output constrained from -1 to 1)
  
  # allocating theta to weights and biases
  index <- 1:(p*m[1])
  W1 <- matrix(theta[index],p,m[1])    # p*m1
  index <- max(index)+ 1:(m[1]*m[2])
  W2 <- matrix(theta[index],m[1],m[2]) # m1*m2
  index <- max(index)+ 1:(m[2]*q)
  W3 <- matrix(theta[index],m,q)       # m2*q
  index <- max(index)+ 1:m[1]
  b1 <- matrix(theta[index],m[2],1)    # m1*1
  index <- max(index)+ 1:m[2]
  b2 <- matrix(theta[index],m[2],1)    # m2*1
  index <- max(index)+ 1:q
  b3 <- matrix(theta[index],q,1)       # q*1
  
  W <- list(W1,W2,W3)
  b <- list(b1,b2,b3)
  sig <- list(sig1,sig2,sig3)
  
  model_object <- structure(list(W=W,b=b,sig=sig),
                            class="IterativeNN",
                            hidden=length(m)) 
  
  model_func <- function(x){ #func takes position |> outputs control vector
    #note: since the coordinates are ONE observation, 
    # it must be a 1x2 matrix when it starts. 1 observation of 2 position vectors
    A0   <- t(x)                              # p*1
    A1   <- sig[[1]]( t(W[[1]])%*%A0+b[[1]] ) # m1*1
    A2   <- sig[[2]]( t(W[[2]])%*%A1+b[[2]] ) # m2*1
    A3   <- sig[[3]]( t(W[[3]])%*%A2+b[[3]] ) # q*1
    c    <- t(A3)  # 1*q
    
    return(c)
  }
  
  return(list(func=model_func, comps=model_object))
}

#--------------------------------RUN GAME--------------------------------------#

# drones: drones to run per batch
run_game <- function(iterative_func, drones = 10, obst){
  
  position <- array(NA, #position matrix extended for multiple runs
                    dim = c(iterlim,3,drones),
                    dimnames = list(step=paste0("k=",1:iterlim),
                                    position=c("X_1","X_2","State"),
                                    run=paste0("r=",1:drones)))
  
  #x & y position & state storage for all runs
  #dev.new(noRStudioGD = TRUE)
  
  for(r in 1:drones){ # 
    k <- 1
    x <- matrix(c(
      runif(n=1,min=start_bounds[1,1],max=start_bounds[2,1]),
      start_bounds[1,2]),
      ncol=2,
      dimnames=list(NULL, #row
                    c("X_1","X_2"))) #col
    sk <- -2 #denotes start point 
    
    #while(position[k,3,r]==0||k==1){ #only the 1st condition will be checked
    # for the most part (exceptional condition (thus do-while); shouldn't  
    # add overhead)
    while(sk==0||sk==-2){
      
      position[k,,r] <- c(x,sk) # store x in row k+1, for run r
      
      # print(field+geom_point(
      #   data=position[[r]],
      #   mapping=aes(x=X_1,y=X_2)))
      
      k <- k+1
      
      c <- iterative_func(x)
      new_x <- x + c*delta
      
      sk <- eval_game_state(x=new_x,obst=obst,step=k)
      
      if (sk == -1){
        position[k,,r] <- 
          c(new_x,sk) # someone needs to store the final position,...
      }
      if (sk == 1){
        position[k,,r] <- 
          c(new_x,sk) # ...it sure won't be the next loop
      }
      if (sk == 0){
        #all is well! GO ON YA DWEEB, NEXT LOOP
        x <- new_x
      }
    }
  }
  
  run_scores <- colSums(position[-1,3,],na.rm=TRUE)
  score <- sum(run_scores,na.rm=TRUE)
  
  data <- do.call(rbind,asplit(position,3)) |> 
    cbind("Run" = rep(1:drones,each=iterlim)) |>
    na.omit() |>
    as_tibble() |>
    mutate(State=as_factor(State),
           Run=as_factor(Run+1))
  
  return(list(data=data, scores=run_scores, obstacles=obst, raw=position))
}


# plotting output of run_game
plot_game <- function(game_obj){
  d <- dim(game_obj$raw)[3]
  run_colours <- viridisLite::viridis(n=d,
                                      alpha=0.7,begin=0.2,end=0.8,option="B")
  names(run_colours)<-c(1:d+1)
  
  colours <- c("bound" = 'firebrick',
               "finish"= 'limegreen',
               "start" = 'yellow3',
               "obst"  = 'forestgreen',
               "1"     = 'green',
               "0"     = '#00000000', #alpha of 0
               "-1"    = 'red',
               "-2"    = 'yellow4',
               run_colours)
  
  base_field <- generate_field(colours,lines=lines)
  
  
  blank_points <- rbind( #expand the frame slightly
    c(-1.1,-1.1),
    c( 1.1, 1.1)) |> as.data.frame()
  
  field_obst <- generate_obstacles(game_obj$obstacles, colours[4]) 
  field <- base_field+field_obst +geom_point(
    data=blank_points,
    aes(x=V1,y=V2),                    
    alpha=0,size=0.01)
  
  
  
  #dev.new(noRStudioGD = TRUE)
  
  # as to not create an excessive amount of layers,
  # the entire dataset is plotted each time
  
  plotted_field <- field + geom_point( #plot all runs as one layer
    data=game_obj$data,
    mapping=aes(x=X_1,y=X_2,fill=Run,colour=State),
    shape=21,
    size=1,
    stroke=1,
    show.legend = FALSE
  )
  plot(plotted_field)
}

# run and plot a game based on the theta values.
plot_theta <- function(theta,d,obst){
  
  run_colours <- viridisLite::viridis(n=d,
                                      alpha=0.7,begin=0.2,end=0.8,option="B")
  names(run_colours)<-c(1:d+1)
  
  colours <- c("bound" = 'firebrick',
               "finish"= 'limegreen',
               "start" = 'yellow3',
               "obst"  = 'forestgreen',
               "1"     = 'green',
               "0"     = '#00000000', #alpha of 0
               "-1"    = 'red',
               "-2"    = 'yellow4',
               run_colours)
  
  base_field <- generate_field(colours,lines=lines)
  
  
  blank_points <- rbind( #expand the frame slightly
    c(-1.1,-1.1),
    c( 1.1, 1.1)) |> as.data.frame()
  
  field_obst <- generate_obstacles(obst, colours[4]) 
  field <- base_field+field_obst +geom_point(
    data=blank_points,
    aes(x=V1,y=V2),                    
    alpha=0,size=0.01)
  
  # as to not create an excessive amount of layers,
  # the entire dataset is plotted each time
  
  game_obj<-run_game(iterative_func=create_IterativeNN(theta),drones=d,obst)
  
  plotted_field <- field + geom_point(
    data=game_obj$data,
    mapping=aes(x=X_1,y=X_2,fill=Run,colour=State),
    shape=21,
    size=1,
    stroke=1,
    show.legend = FALSE
  )
  plot(plotted_field)
}


# wrapper which simply outputs scores 
run_batch <- function(theta, drones = 10, obst,
                      plot=FALSE){   # plotting slows down the process
  
  #theta -> NN function
  model <- create_IterativeNN(theta)
  
  game_out <- run_game(model$func, drones=drones, obst=obst)
  if(plot){
    plot_game(game_out)
  }
  return(score=sum(game_out$scores))
}


# -----------------------RUNNING GAMES WITH GA---------------------------------

# --- --- --- --- how to run games --- --- --- ---
# score <- run_batch(theta.init, 
#                     drones=d, 
#                     obst=obst) #only returns total score
# 
# 
# game_out <- run_game(iterative_func = create_IterativeNN(theta.init)$func,
#                      drones=d,
#                      obst=obst) #returns model object and scores
# 
# plot_game(game_obj = game_out)


npars <- 57

# random generation here
theta.init <- runif(npars,-1,1)
obst <- cbind(X_1  = runif(n_obst, obst_bounds[1,1],obst_bounds[2,1]),
              X_2  = runif(n_obst, obst_bounds[1,2],obst_bounds[2,2]),
              Size = rep(obst_size, n_obst))
# end random generation

d <- 10


# we must vary theta after each run based on GA
obj<-run_batch(theta=theta.init,drones=d,obst=obst)

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

n_fcalls <- fcalls |> count()
distinct_fcalls <- fcalls |> distinct(pick(!X1),.keep_all = TRUE)
n_distinct <- distinct_fcalls |> count()


# ~ 81% distinct.
prop_distinct <- n_distinct/n_fcalls

 
duplicate_fcalls <- anti_join(fcalls,distinct_fcalls)
n_duplicated <- duplicate_fcalls |> count() #issue here. TODO, figure out. Missing entries.
# want to inspect for patterns in duplicated function calls.





n_fcalls <- fcalls |> count()

distinct_fcalls <- fcalls |> distinct(pick(!X1),.keep_all = TRUE)
n_distinct <- distinct_fcalls |> count()

duplicate_fcalls <- anti_join(fcalls,distinct_fcalls)
n_duplicated <- duplicate_fcalls |> count() 

n_duplicated + n_distinct


