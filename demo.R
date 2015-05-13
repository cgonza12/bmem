# install all the important stuff

install.packages('devtools')
install_github('rmcelreath/glmer2stan')
source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
install_rstan()
source_url("https://github.com/stan-dev/shinystan/raw/develop/install_shinystan.R")
install_shinystan()
install.packages('lmerTest')


library(devtools)
library(shinyStan)
library(glmer2stan)
library(rstan)
library(lmerTest)
# load my functions
source('mypostplot.R')
source('myglmer2stan.R')
source('panelfuns.R')
source('postPairs.R')
source('stanmer2.R')
source('propInRope.R')

# Outline
# 1) Intro to STAN
# 2) glmer2stan 
# 3) my functions
# 4) shinyStan




####### Intro to STAN ########
# STAN is... 
#   - MCMC sampler similar to JAGS and BUGS
#   - Runs in C++
#   - rstan is the interface between r and STAN
#   - write code -> compile model -> sample -> return output


height = rnorm(100)
earn = height + rnorm(100)

dat = data.frame(height,earn)

library(ggplot2)
qplot(height,earn)+geom_smooth(method=lm,se=F)

# Simple frequentist model
summary(lm(earn~height,data = dat))
confint(lm(earn~height,data = dat))

# Important annoying fact #1: STAN needs data as a list not a dataframe!!
# specify data as well as meta data (i.e. the number of groups)

earn_dat <- list(N = 100 , #specify number of observations as a scalar
                    earn = earn, # data vector
                    height = height # data vector (predictor) 
                    )

earn_code = 'data {
  // First we declare all of our variables in the data block
  int<lower=0> N;// Number of observations
  vector[N] earn; //Identify our predictor as a vector
  vector[N] height;  //Identify our outcome variable as a vector
}
parameters {
  vector[2] beta; //Our betas are a vector of length 2 (intercept and slope)
  real<lower=0> sigma; //error parameter
}
model {
  //Priors
  beta[1] ~ normal( 5 , .001); //intercept
  beta[2] ~ normal( 0 , 100 ); //slope
  sigma ~ uniform( 0 , 100 ); //error
  earn ~ normal(beta[1] + beta[2] * height, sigma);
}'


fit1 <- stan(model_code = earn_code, data = earn_dat,
             warmup = 100,
             iter = 1000, 
             chains = 4)

print(fit1)
# Output
# Mean: parameter estimate
# se_mean: standard error of mean estimate
# sd: standard deviation of posterior
# 2.5 - 97.5: quantiles of posterior 
# n_eff: effective sample size, a measure of autocorrelation among samples 
# Rhat: split-chain convergence diagnostic, >1.1 suggests poor convergance


# extract posterior samples for each parameter
fit1_samples = extract(fit1)
str(fit1_samples)
# subset just the betas
betas = fit1_samples[[1]]

qplot(betas[,1]) # intercept posterior samples
qplot(betas[,2]) # slope posterior samples

int = mean(betas[,1]) # posterior intercept estimate
slope = mean(betas[,2])# posterior slope estimate

# check fit 
qplot(height,earn)+geom_smooth(method=lm,se=F)+geom_abline(intercept=int,slope=slope,color='red')

# View trace plots of each parameter (fuzzier is better)
traceplot(fit1)

##### glmer2stan ####
# - allows user to write mixed effects formulas using lme4 syntax
# - converts formula and data to STAN friendly formats
# - returns stanmer object
# - computes WAIC (a bayesian model comparison statistic)


data(sleepstudy) # load data
# The average reaction time per day for subjects in a sleep deprivation study. 
# On day 0 the subjects had their normal amount of sleep.
# Starting that night they were restricted to 3 hours of sleep per night. 
# The observations represent the average reaction time on a series of tests 
# given each day to each subject

library(ggplot2)

# plot effects first!

# no random effects model
m1_lm <- lm(Reaction~Days,data=sleepstudy)
confint(m1_lm)
summary(m1_lm)

ggplot(sleepstudy,aes(x=Days,y=Reaction))+
  geom_point()+
  guides(color=F)+
  geom_smooth(method=lm,se = F)

# with nesting of intercepts and slopes within subjs
ggplot(sleepstudy,aes(x=Days,y=Reaction,color=Subject,group=Subject))+
  geom_point()+
  guides(color=F)+
  geom_smooth(method=lm,se = F)

# model in regular ol' lme4, random slopes and intercepts set
m1_lme4 <- lmer( Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE )
summary(m1_lme4)
confint(m1_lme4)

AIC(m1_lm,m1_lme4) # fixed effects the same, model is more parsimonious i.e. random effects = good

# Important annoying fact #2: STAN doesn't deal with non-numeric variables!!
# - factors must be converted into contrast codes (don't worry glmer2stan does that for you)
# - grouping variables,however, must be manually converted to integers

# convert factor of subject ids to sequential integers
sleepstudy$subject_index <- as.integer(as.factor(sleepstudy$Subject)) 

# basic MCMC parameter, should probably be a bit higher but we dont have all day!
nwarm = 100 # burn-in period, these samples are not included in estimation
niter = 500 # number of steps per chain, more is better (but takes longer)
chains = 4 # number of chains, usually at least 2

m1_g2s <- lmer2stan( Reaction ~ Days + (Days | subject_index), data=sleepstudy,
                     calcWAIC=T,
                     warmup=nwarm, 
                     iter = niter, 
                     chains=chains) 

print(m1_g2s) # standard stan output
stanmer(m1_g2s) # cleaned up stan output
plot(m1_g2s) # *looks like shit
traceplot(m1_g2s)

# Note: random effects increase processing time exponentially, see: nested for loops, no random slope

m1noslope_g2s <- lmer2stan( Reaction ~ Days + (1 | subject_index), data=sleepstudy,
                     calcWAIC=T,
                     warmup=nwarm, #burn-in period
                     iter = niter, # number of steps per chain
                     chains=chains) #number of chains

stanmer(m1noslope_g2s)


# output is an S4 class of object, access different chunks with the @ symbol
m1_g2s@stanmodel #stan model

# can also run the model with sampling turned off. This returns a list with the model code
# as well as the converted data object. Very helpful if you want to change priors

m1_g2s_noSamples <- lmer2stan( Reaction ~ Days + (Days | subject_index), data=sleepstudy,
                     calcWAIC=T,
                     warmup=nwarm, #burn-in period
                     iter = niter, # number of steps per chain
                     chains=chains,
                     sample=F #turns sampling off
                     )

str(m1_g2s_noSamples)

m1_g2s_noSamples$data # data object formatted to be STAN friendly

cat(m1_g2s_noSamples$model) # stan model as text

# use sink() to write output to a text file instead of the file

sink('example1.txt') #open connection
cat(m1_g2s_noSamples$model) # any output between these two functions gets written to that file
sink() # close connection

#### my functions ####
# myglmer2stan ####
# - fixes bug in formula code which crashes with interactions specified by ':' rather than '*'
# - allows user to use homebrewed model as a function argument 

# the problem with glmer2stan is that you can extract the model code no problem, but its not 
# straightforward at all to re run the model with edits.

# First we'll re run with no edits to the model
m2_g2s <- myglmer2stan( Reaction ~ Days + (Days | subject_index), data=sleepstudy,
                     calcWAIC=T,
                     warmup=nwarm, 
                     iter = niter,
                     chains=chains, 
                     mymodel=NULL # when this is NULL the model is built and run like normal
                     ) 


stanmer(m2_g2s) # output is the same as m1_g2s

cat(m1_g2s_noSamples$model) # print model

# glmer2stan defaults to diffuse normal priors for all fixed effects betas normal(0,100)
# by printing the model we can edit those to be more reasonable based on the scale of the data
# or informed from prior work etc.

# save model text to an object
m2_g2s_priors = "data{
  int N;
  real Reaction[N];
  real Days[N];
  int subject_index[N];
  int N_subject_index;
}

transformed data{
  vector[2] zeros_subject_index;
  for ( i in 1:2 ) zeros_subject_index[i] <- 0;
}

parameters{
  real Intercept;
  real beta_Days;
  real<lower=0> sigma;
  vector[2] vary_subject_index[N_subject_index];
  cov_matrix[2] Sigma_subject_index;
}

model{
  real vary[N];
  real glm[N];
  // Priors
  Intercept ~ normal( 0 , 100 );
  beta_Days ~ normal( 10 , 2 ); // Hey! Look! An informed prior!
  sigma ~ uniform( 0 , 100 );
  // Varying effects
  for ( j in 1:N_subject_index ) vary_subject_index[j] ~ multi_normal( zeros_subject_index , 
        Sigma_subject_index );
  // Fixed effects
  for ( i in 1:N ) {
    vary[i] <- vary_subject_index[subject_index[i],1]
    + vary_subject_index[subject_index[i],2] * Days[i];
    glm[i] <- vary[i] + Intercept
    + beta_Days * Days[i];
  }
  Reaction ~ normal( glm , sigma );
}

generated quantities{
  real dev;
  real vary[N];
  real glm[N];
  dev <- 0;
  for ( i in 1:N ) {
    vary[i] <- vary_subject_index[subject_index[i],1]
    + vary_subject_index[subject_index[i],2] * Days[i];
    glm[i] <- vary[i] + Intercept
    + beta_Days * Days[i];
    dev <- dev + (-2) * normal_log( Reaction[i] , glm[i] , sigma );
  }
}"


m2_g2s_withpriors <- myglmer2stan( Reaction ~ Days + (Days | subject_index), data=sleepstudy,
                        calcWAIC=T,
                        warmup=100, 
                        iter = 500,
                        chains=4, 
                        mymodel=m2_g2s_priors # include object with edited model text
) 

stanmer(m2_g2s_withpriors)

##### Posterior predictive check ####
# In addition to estimation it's usually a good idea to do some posterior predictive checks
# i.e. make sure the model output fits the data reasonably well


# overall model fit
ggplot(sleepstudy,aes(x=Days,y=Reaction))+
  geom_point()+
  guides(color=F)+
  geom_smooth(method=lm,se=F)+
  geom_abline(intercept = stanmer2(m1_g2s)[[1]][1],
              slope=stanmer2(m1_g2s)[[1]][2],color='red')+
  theme_bw()+
  theme(text=element_text(size=16))


# individual slopes and intercepts

m1_g2s_samples = extract(m1_g2s) # extract posterior samples for each parameter


# how many posterior samples for each participant?
nsim = 20
#create an empty array to fill with slopes and intercepts
check_subj_fit = array(0,c(max(sleepstudy$subject_index)*nsim,3)) 

for(i in 1:length(sleepstudy$subject_index)){
  # crap to fill the array right
  idx = seq(1,nsim*length(sleepstudy$subject_index),nsim)
  x = idx[i]:idx[i+1]-1
  start = idx[i]
  end = x[length(x)]
  
  check_subj_fit[start:end,1] = i # subject
  # N samples of intercepts for subject i (each is a deflection from the grand mean of ints)
  check_subj_fit[start:end,2] = stanmer2(m1_g2s)[[1]][1]+sample(
                                            m1_g2s_samples$vary_subject_index[,i,1],size = nsim)  
  # N samples of slopes for subject i (each is a deflection from the grand mean of slopes)
  check_subj_fit[start:end,3] = stanmer2(m1_g2s)[[1]][2]+sample(
                                            m1_g2s_samples$vary_subject_index[,i,2],size = nsim)
  
  
}
check_subj_fit = as.data.frame(check_subj_fit)
names(check_subj_fit) <- c('subject_index','int','slope')

ggplot(sleepstudy,aes(x=Days,y=Reaction))+
  geom_abline(data=check_subj_fit,aes(intercept=int,slope=slope),alpha=.5)+
  geom_point(shape=5,color='red',size=3)+
  geom_smooth(method=lm,se=F)+
  facet_wrap(~subject_index)+
  theme_bw()+
  theme(text=element_text(size=16))


# propInRope ####
# - computes proportion of posterior samples for fixed effects in specified ROPE
# - set ROPE to whatever makes sense

propInRope(m2_g2s,ropemin = -1, ropemax = 1)

# mypostplot ####
# - visualizes mixed effects model output
# - visualizes ROPE
# - plots X% HDIs for fixed effects
# - colored based on proportion in ROPE

# make a new variable to build complexity into the model
sleepstudy$Days2 = rnorm(length(sleepstudy$Days)) # just noise

m3_g2s <- myglmer2stan( Reaction ~ Days*Days2 + (Days | subject_index), data=sleepstudy,
                        calcWAIC=T,
                        warmup=nwarm, 
                        iter = niter, 
                        chains=chains 
) 

stanmer(m3_g2s)
mypostplot(m2_g2s,ropemin = -1,ropemax = 1)

# postPairs ####
# - look at posterior sample correlations in a scatterplot matrix
# - can also output text correlations
# - and save too!
# Helpful for diagnosing multicollinearity in your model

postPairs(m3_g2s)

##### Shiny Stan ####

launch_shinystan(m3_g2s)
