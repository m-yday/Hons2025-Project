m <- c(5,5) #hidden layer size
p <- 2  #predictors (x)
q <- 2  #responses  (c)
npars <- p*m[1] +m[1]*m[2] +m[2]*q +m[1]+m[2]+q # =57
# p*m[1] +(m[2]+1)*(m[1]+q) + m[2]

# creates IterativeNN from param list

# IterativeNN is a NN which only evaluates one observation at a time.
# The NN does not update between these observations,
# but the result does depend on the (position) row vector input to the NN.

# returns a NN as a function of x, and the model components in an object
# temporarily constrained to 2 hidden layer structure
# very hardcoded at the moment

create_IterativeNN <- function(theta, m=c(5,5),p=2,q=2,
                               npars=p*m[1] +m[1]*m[2] +m[2]*q +m[1]+m[2]+q){
  stopifnot(is.numeric(theta) || length(theta)==npars)
  
  # --- activation functions ------------------------------------------------
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
                            hidden=2) 
  
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


# # USAGE

# #prepare any parameter list:
# theta <- runif(npars,-1,1)
# 
# # [create the model]:p
# model <- create_IterativeNN(theta)
# 
# # [use the model]:
# # get position
# x <- c(0,-1) |> matrix(nrow=1)
# # ENSURE position is a 1*p matrix
# 
# # call function
# model$func(x)
# 
# # [extract model details]:
# model$comps

