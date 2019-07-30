#regression adventure

#define the logistic (or inverse logit) function
logistic <- function(x){
  odds <- exp(x)
  prob <- odds/(odds+1)
  return(prob)
}

#step 1: get the data (from sim.data, data.frame=df)
data.list <- list(
  id=df$id,
  outcome=df$outcome,
  intervention= df$intervention,
  confounder=df$confounder
)

data.list$N <- nrow(df)

#step 2: get the package
require(rstan)

#step 3: specify the model
stan.1 <- '
data{
    int<lower=1> N;                           //number of observations
    int<lower=0, upper=1> outcome[N];         //binary outcome
    int<lower=0, upper=1> intervention[N];    //intervention
    real<lower=0, upper=1> confounder[N];     //confounder (proportion)
}

parameters{
    real alpha;
    real beta[2];
}

model{
    //vector of probabilities
    vector[N] p;
    //specify prior distributions
    alpha ~ normal(0, 10);
    beta ~ normal(0, 10);

    //set up regression
    for(i in 1:N){
        p[i] = inv_logit(alpha + beta[1]*intervention[i] + beta[2]*confounder[i]);
    }
    //likelihood function
    outcome ~ bernoulli(p);
}'

#step 4: run the model
mod.1.out <- stan(model_code=stan.1, data=data.list,
                  iter=5000, warmup=2000, chains=1, cores=1)

#step 5: check the model output
plot(mod.1.out) #plot parameters on log-odds scale
posterior <- extract(mod.1.out) #extract posterior distribution

#examine chains for convergence
plot(posterior$alpha)
lines(posterior$alpha)

#density plots of parameters (converted to probability scale)
plot(density(logistic(posterior$alpha)))
plot(density(logistic(posterior$alpha+posterior$beta[,1])))
plot(density(logistic(posterior$alpha+posterior$beta[,2])))

#retrieve initial values (when no confounding)
alpha.post <- logistic(posterior$alpha)
alphaBeta1 <- logistic(posterior$alpha + posterior$beta[,1])
alphaBeta2 <- logistic(posterior$alpha + posterior$beta[,2])
