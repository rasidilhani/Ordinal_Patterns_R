map <- function(k, x0, n){ 
  x <- vector(length=n)
  x[1] <- x0 
  for (i in 2:n){
    x[i]=k*x[i-1]*(1-x[i-1])
  }
return(x)  
}

(xfar1 <- map(3.9, 0.50001, 1000000)[1000000])
(xfar2 <- map(3.9, 0.50002, 1000000)[1000000])
(xfar3 <- map(3.9, 0.5000001, 1000000)[1000000])

my300observations <- map(3.9, 0.50002, 1000000)[500000:500300]

plot(my300observations, type="b")
