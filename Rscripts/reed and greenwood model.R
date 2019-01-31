# Reed-Frost Model
gen <- 10
s1 = 99
i1 = 1

q = 0.99 # escape probability in one contact
S <- I <- rep(0, gen)# initialize S,I,R
S[1] <- s1
I[1] <- i1

for (j in 1:gen) {
      prob = 1 - (1 - q)^I[j]
      I[j+1] <- rbinom(1, S[j], prob)
      S[j+1] <- S[j] - I[j+1]
}

par(mfrow = c(2,1))
plot(c(0:gen), I, type="b", ylim=c(0,100), main = "Reed-Frost Model")


# Greenwood Model
gen <- 10
s1 = 99
i1 = 1

q=0.99 # escape probability in one contact
S <- I <- rep(0, gen) # initialize S,I,R
S[1] <- s1
I[1] <- i1

for (j in 1:gen) {
      prob = 1 - (1 - q)# Here comes the difference between Greenwood and Reed Model
      I[j+1] <- rbinom(1, S[j], prob)
      S[j+1] <- S[j] - I[j+1]
}

plot(c(0:gen), I, type="b", ylim=c(0,100), main="Greenwood Model")
