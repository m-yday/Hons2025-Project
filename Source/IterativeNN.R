m <- c(5,5) #hidden layer size
p <- 2  #predictors (x)
q <- 2  #responses  (c)
npars <- p*m[1] +m[1]*m[2] +m[2]*q +m[1]+m[2]+q # =57

# creates IterativeNN from param list
# returns a NN as a function of x, and the model components in an object
# temporarily constrained to 2 hidden layer structure
# very hardcoded at the moment

create_IterativeNN <- function(theta){
  stopifnot(is.numeric(theta))
  
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
  
  model_func <-function(x){
    A0   <- t(x)                  # p*1
    A1   <- sig[[1]](t(W[[1]])%*%A0+b[[1]])   # m1*1
    A2   <- sig[[2]](t(W[[2]])%*%A1+b[[2]])   # m2*1
    A3   <- sig[[3]](t(W[[3]])%*%A2+b[[3]])   # q*1
    c    <- t(A3)                 # 1*q
  }
  
    return(list(func=model_func, comps=model_object))
}


#USAGE

#prepare any parameter list:
theta <- runif(npars,-1,1)

#create the model:
model <- create_IterativeNN(theta)

#use the model:
model$func(c(0,0))

#extract model details:
model$comps