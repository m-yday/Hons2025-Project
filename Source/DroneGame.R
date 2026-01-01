# game params

n_obst    <- 40 #J
obst_size <- 0.02 #r (radius) 2*r = square side length

delta <- 0.05
iterlim <- 100 #k

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


# theta needs to be defined such that all the game parameters are present, such that the model can 'see' the details of the game
# anything not present in the theta parameters is 'unknown' to the model or 
# "assumed" by the model (i.e., universal properties -> dimensionality = 2) 

# above does not seem to be true. that is based on a misunderstanding of theta.
# the params do not include object locations.

# theta structure:
# num obstacles / trees
# obstacle size (future: sizes?)
# obstacle locations (x-y pairs)
#

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

generate_field <- function(colours){
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
  if (step > iterlim || detect_collision(x, obst))
    return(-1)  
  
  if (x[2] >= bounds[2,2])
    return(1)
  
  return(0)
}

#--------------------------------RUN GAME--------------------------------------#

runs <- 10

run_colours <- viridisLite::viridis(n=runs,alpha=0.7,begin=0.2,end=0.8,option="B")
names(run_colours)<-c(1:runs+1)

colours <- c("bound" = 'firebrick',
             "finish"= 'limegreen',
             "start" = 'yellow3',
             "obst"  = 'forestgreen',
             "1"     = 'green',
             "0"     = '#00000000', #alpha of 0
             "-1"    = 'red',
             "-2"    = 'yellow4',
             run_colours)

base_field <- generate_field(colours)

set.seed(2)
obst <- cbind(X_1  = runif(n_obst, obst_bounds[1,1],obst_bounds[2,1]),
              X_2  = runif(n_obst, obst_bounds[1,2],obst_bounds[2,2]),
              Size = rep(obst_size, n_obst))

blank_points <- rbind( #expand the frame slightly
  c(-1.1,-1.1),
  c( 1.1, 1.1)) |> as.data.frame()

field_obst <- generate_obstacles(obst, colours[4]) 
field <- base_field+field_obst +geom_point(
  data=blank_points,
  aes(x=V1,y=V2),                    
  alpha=0,size=0.01)

run_game <- function(iterative_func){
  
  position <- array(NA, #position matrix extended for multiple runs
                    dim = c(iterlim,3,runs),
                    dimnames = list(step=paste0("k=",1:iterlim),
                                    position=c("X_1","X_2","State"),
                                    run=paste0("r=",1:runs)))
  
  #x & y position & state storage for all runs
  #dev.new(noRStudioGD = TRUE)
  
  for(r in 1:runs){
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
    cbind("Run" = rep(1:runs,each=iterlim)) |>
    na.omit() |>
    as_tibble() |>
    mutate(State=as_factor(State),
           Run=as_factor(Run+1))
  
  return(list(data=data, scores=run_scores, raw=position))
}


#initial theta
theta <- runif(npars,-1,1)

#initial model
model <- create_IterativeNN(theta)


for(seed in 1:10){
  set.seed(seed)
  game_out <- run_game(model$func)
  #dev.new(noRStudioGD = TRUE)
  this_field <- field
  
  this_field <- this_field + geom_point(
    data=game_out$data,
    mapping=aes(x=X_1,y=X_2,fill=Run,colour=State),
    shape=21,
    size=1,
    stroke=0.5,
    show.legend = FALSE
  ) +
    labs(title=sprintf("Seed = %d",seed))
  
  plot(this_field)
}