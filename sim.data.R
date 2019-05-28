#simulate the data!

#total number of patients
n <- 600

#number of patients that recieve the intervention
nI <- 300

#function for 'normal distribution' between 0 and 1
confounder.sim <- rbeta(n=n,shape1=2,shape2=2)

#intervention
intervention.sim <- c(rep(0,nI), rep(1,nI))

#set up data.frame
df <- data.frame(id = seq(1:n), intervention=intervention.sim, confounder=confounder.sim)

#we say the probability of the outcome with intervention=0 & confounder =0 is 0.3 (alpha)
#change in prob of outcome when intervention = 1 and confounder =0 is +0.1 (beta1)
#change in prob of outcome when intervention =0 and confounder =1 is +0.5 (beta2)
alpha <- 0.3
beta1 <- 0.1
beta2 <- 0.5

#vector of probabilities
probability  <- with(df, (alpha + beta1*intervention + beta2*confounder))

#simulate outcome
df$outcome <- rbinom(n=n, size=1, prob=probability)

