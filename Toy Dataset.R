t <- seq(0,1, length.out = 100)

w <- 25

a <- seq(0.25,3, length.out = 10)
b <- seq(0.25,1, length.out = 10)

y <- list()

for(i in 1:length(a)){
 y[[i]] <- a[i]*sin(w*t + b[i]) + a[i]*sin(w*t + b[i])^3 + a[i]*sin(w*t + b[i])^5 
}

# x11()
plot(t, y[[1]], type = "l", ylim = c(-10,10))

for(i in 1:length(a)){
  lines(t, y[[i]], col = i)
}

