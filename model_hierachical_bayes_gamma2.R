#================================================================================== 
#                           Hierarchical model simulate
#==================================================================================

rm(list=ls(all=TRUE))
setwd("C:/Users/joaquin/Desktop/Desktop/Stan/ProjectAalto/hierarchical")

# Packages
library(TMB)
library(reshape2)
library(tmbstan)
library(shinystan)
library(rstan)
library(rstanarm)



options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')



#===================================================================================================================
#                                Firts, I must create the model in TMB
#===================================================================================================================
model_hierachical_bayes_gamma2 = '
#include <TMB.hpp>

//
//Custom likelihood functions, used be used in template
//below. These are not built into TMB like dnorm and dgamma are.

//log-normal likelihood
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
// Inverse Gaussian
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}


// dcaychy for hyperparameters
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{



//===================================================================================================================
// Data

  DATA_INTEGER(likelihood); 	                  // Likelihood flag (to compare differents models)
  DATA_VECTOR(y_is);                            // vector of the observations
  DATA_IVECTOR(year);                           // year as a factor
  DATA_IVECTOR(depth);                          // depth
  DATA_IVECTOR(trim);                           // trimester
  DATA_IVECTOR(destine);                        // destine
  DATA_IVECTOR(site);                           // random effects (site of fishing)
//===================================================================================================================

//===================================================================================================================
// Parameters fixed
   PARAMETER(intercept);
   PARAMETER_VECTOR(beta_year);                   //(beta for year)
   PARAMETER(beta_depth);                         //(beta for depth)
   PARAMETER_VECTOR(beta_trim);                   //(beta for trim)
   PARAMETER_VECTOR(beta_destine);                //(beta for destine)
   
// Sigmas in betas   
   PARAMETER(sigma_beta_year);
   PARAMETER(sigma_beta_depth);
   PARAMETER(sigma_beta_trim);
   PARAMETER(sigma_beta_destine);

   PARAMETER(logsigma);		                     // log of observation sd
   PARAMETER(logsigma_space);	                 // spatial variance



// Postive transformations, jacobians below
   Type yearSD2=exp(sigma_beta_year);
   Type depthSD2=exp(sigma_beta_depth);
   Type trimSD2=exp(sigma_beta_trim);
   Type destineSD2=exp(sigma_beta_destine);

   
   Type nlp=0.0; // negative log prior
   Type nll=0.0; // negative log likelihood

// Jacobian adjustment for variances
   nll -= sigma_beta_depth + sigma_beta_year + sigma_beta_trim + sigma_beta_destine;


// Parameter non-centered random effects
   PARAMETER_VECTOR(u);
   

//===================================================================================================================

//===================================================================================================================
// Priors

  // priors
  nlp-= dnorm(intercept,     Type(0.0), Type(5), true);
  nlp-= dnorm(beta_year,     Type(0.0), Type(5), true).sum();
  nlp-= dnorm(beta_depth,    Type(0.0), Type(5), true);
  nlp-= dnorm(beta_trim,     Type(0.0), Type(1.5), true).sum();
  nlp-= dnorm(beta_destine,  Type(0.0), Type(1.5), true).sum();
  
  nlp-= dcauchy(yearSD2,     Type(0), Type(0.5), true);
  nlp-= dcauchy(depthSD2,    Type(0), Type(0.2), true);
  nlp-= dcauchy(trimSD2,     Type(0), Type(0.2), true);
  nlp-= dcauchy(destineSD2,  Type(0), Type(0.2), true);
  
  nlp-= dcauchy(logsigma,          Type(0), Type(0.5), true);
  nlp-= dcauchy(logsigma_space,    Type(0), Type(0.2), true);

// random effects non-centered
// nlp-=dnorm(u, Type(0.0), Type(5.0)).sum();
//===================================================================================================================  



//===================================================================================================================

// Linear predictor
   Type sigma = exp(logsigma);
   Type sigma_space = exp(logsigma_space);
   int n = y_is.size();
   vector<Type> pred(n);
   for(int i=0; i<n; i++){
   pred(i) = exp(intercept + beta_year(year(i)) + beta_depth*depth(i) + beta_trim(trim(i)) + beta_destine(destine(i)) + u(site(i)));		// the spatial effect
}

//   


// Probability of data conditional on fixed effect values
  for(int i=0; i<n; i++){
      // Likelihood
      if(likelihood==1)                                                            // lognormal
	    nll -= dinvgauss(y_is(i), pred(i), sigma, true);
      else if(likelihood==2)                                                       // inverse gaussiana 
	    nll -= dlognorm(y_is(i),  log(pred(i)), sigma, true);
	    else if(likelihood==3)                                                       // gamma 
      nll -= dgamma(y_is(i), 1/pow(sigma,2), pred(i)*pow(sigma,2), true);
      else {
      	std::cout << "Invalid likelihood specified" << std::endl;
        return 0;
      }
}


// Probability of site means
  int n_site = u.size();
  for( int s = 0; s < n_site; s++) {
    nll-= dnorm(u(s), intercept, sigma_space, true);      // Comment: because sites (s) ~ N(intecerpet, sigma_space)
    
    
// Random effects; non-centered
   nll-=dnorm(u, Type(0.0), Type(1.0), true).sum();
  }

//===================================================================================================================


//===================================================================================================================
// Reporting
  REPORT(intercept);
  REPORT(beta_year);
  REPORT(beta_depth);
  REPORT(beta_trim);
  REPORT(beta_destine);
  REPORT(sigma_beta_year);
  REPORT(sigma_beta_depth);
  REPORT(sigma_beta_trim);
  REPORT(sigma_beta_destine);
  REPORT(sigma);
  REPORT(sigma_space);
  REPORT(u);
  REPORT(pred);
  
  ADREPORT(intercept);
  ADREPORT(beta_year);
  ADREPORT(beta_depth);
  ADREPORT(beta_trim);
  ADREPORT(beta_destine);
  ADREPORT(sigma_beta_year);
  ADREPORT(sigma_beta_depth);
  ADREPORT(sigma_beta_trim);
  ADREPORT(sigma_beta_destine);
  ADREPORT(u);
  
  Type nld=nll + nlp; // for obtain negative log density
  return(nld);
}
'
#====================================================================================================================

write(model_hierachical_bayes_gamma2, file = "model_hierachical_bayes_gamma2.cpp")
compile("model_hierachical_bayes_gamma2.cpp")



#====================================================================================================================
#                                        data (cpue)
#====================================================================================================================

# load data "data_cpue"
dat = read.csv("north2.csv", header=T)
dim(dat)
head(dat)
str(dat)

dat$year = as.factor(dat$year)
dat$site = as.factor(dat$site)
dat$trim = as.factor(dat$trim)
dat$destine = as.factor(dat$destine)

#dat <- subset(dat, year==c("2010", "2011", "2012", "2013",  "2014", "2015", "2016"))

# Create the data
data = list(likelihood = 3,                               # 1 = inverse gaussian, 2 = lognormal, 3 = gamma 
                  y_is = dat$cpue,
                  year = as.numeric(dat$year)-1,
                  depth = dat$depth,
                  trim =  as.numeric(dat$trim)-1, 
                  destine = as.numeric(dat$destine)-1,
                  site =  as.numeric(dat$site)-1)


parameters =  list(intercept = 0, 
                  beta_year = rep(0, length(levels(dat$year))),
                  beta_depth = 0, 
                  beta_trim = rep(0, length(levels(dat$trim))),
                  beta_destine = rep(0, length(levels(dat$destine))),
                  sigma_beta_year=0,
                  sigma_beta_depth=0.5,
                  sigma_beta_trim=0.5,
                  sigma_beta_destine=0.5,
                  logsigma = 1, 
                  logsigma_space = 0.5, 
                  u = rep(0,length(levels(dat$site))))




#======================================================================================================
#                                   load the model created in TMB 
                          dyn.load(dynlib("model_hierachical_bayes_gamma2"))
#======================================================================================================

#==============================
#      INVERSE GAUSSIAN
#==============================
obj = MakeADFun(data = data, parameters = parameters, random = "u", DLL="model_hierachical_bayes_gamma2")

# Optimize
obj$fn()
opt = with(obj, nlminb(par, fn, gr)) #restart 
opt$par

rep = sdreport(obj)
summary(rep, "random")                      ## Only random effects
summary(rep, "fixed", p.value = TRUE)       ## Only non-random effects
summary(rep, "report")                      ## Only report




#==============================
#           LOGNORMAL
#==============================
obj2 = MakeADFun(data = data, parameters = parameters, random = "u", DLL="model_hierachical_bayes_gamma2")

# Optimize
obj2$fn()
opt2 = with(obj2, nlminb(par, fn, gr)) #restart 
opt2$par

rep2 = sdreport(obj2)
summary(rep2, "random")                      ## Only random effects
summary(rep2, "fixed", p.value = TRUE)       ## Only non-random effects
summary(rep2, "report")                      ## Only report




#==============================
#           GAMMA
#==============================
obj3 = MakeADFun(data = data, parameters = parameters, random = "u", DLL="model_hierachical_bayes_gamma2")

# Optimize
obj3$fn()
opt3 = with(obj3, nlminb(par, fn, gr)) #restart
opt3$par

rep3 = sdreport(obj3)
summary(rep3, "random")                      ## Only random effects
summary(rep3, "fixed", p.value = TRUE)       ## Only non-random effects
summary(rep3, "report")                      ## Only report




## Calculate AIC
AIC1 <- 2*opt$objective +2*length(opt$par)
AIC2 <- 2*opt2$objective +2*length(opt2$par)
AIC3 <- 2*opt3$objective +2*length(opt3$par)
c(AIC1, AIC2, AIC3)




#==========================================================================
#                   obtain outputs from better model (GAMMA)
#==========================================================================
rep = sdreport(obj)
summary(rep, "random")                      ## Only random effects
summary(rep, "fixed", p.value = TRUE)       ## Only non-random effects
summary(rep, "report")                      ## Only report

obj$report()


#======================================================================================================
#                                      Types of optmization 
#======================================================================================================

# Using optimHess for finite-difference hessian using function only
Hess = optimHess(opt$par, fn=obj$fn )
SD = sdreport( obj, hessian.fixed=Hess )
SD
summary(SD, "random")

# Also fails using finite-difference hessian
sqrt(diag(solve(Hess)))


# Works using different finite-difference method
Hess = numDeriv::hessian( func=obj$fn, x=opt$par )
SD = sdreport(obj, hessian.fixed=Hess )
SD
summary(SD, "random")


#=====================================================================================================
#                                             tmbstan
#=====================================================================================================
library(tictoc)

tic("Time of estimation")
fit_mcmc = tmbstan(obj, chains=3, control = list(max_treedepth = 15, adapt_delta = 0.99), iter=3000, laplace = FALSE)
toc()

## Methods provided by 'rstan'
class(fit_mcmc)
methods(class="stanfit")

launch_shinystan(fit_mcmc)


# to obtain marginal posteriors of specyfy parameters
# 1) for each parameters
post <- as.matrix(fit_mcmc, pars = c("intercept", "beta_year", "u"))

# 2) all parameters
posterior = as.matrix(fit_mcmc)

x11()
par(mfrow=c(4,4))
plot(posterior[,'u[1]'], type = "l", col = "cadetblue")
plot(posterior[,'u[2]'], type = "l", col = "cadetblue")
plot(posterior[,'u[3]'], type = "l", col = "cadetblue")
plot(posterior[,'u[4]'], type = "l", col = "cadetblue")
plot(posterior[,'u[5]'], type = "l", col = "cadetblue")
plot(posterior[,'u[6]'], type = "l", col = "cadetblue")
plot(posterior[,'u[7]'], type = "l", col = "cadetblue")
plot(posterior[,'u[8]'], type = "l", col = "cadetblue")
plot(posterior[,'u[9]'], type = "l", col = "cadetblue")
plot(posterior[,'u[10]'], type = "l", col = "cadetblue")
plot(posterior[,'u[11]'], type = "l", col = "cadetblue")
plot(posterior[,'u[12]'], type = "l", col = "cadetblue")
plot(posterior[,'u[13]'], type = "l", col = "cadetblue")

x11()
par(mfrow=c(4,4))
hist(posterior[,'u[1]'], col = "darksalmon")
hist(posterior[,'u[2]'], col = "darksalmon")
hist(posterior[,'u[3]'], col = "darksalmon")
hist(posterior[,'u[4]'], col = "darksalmon")
hist(posterior[,'u[5]'], col = "darksalmon")
hist(posterior[,'u[6]'], col = "darksalmon")
hist(posterior[,'u[7]'], col = "darksalmon")
hist(posterior[,'u[8]'], col = "darksalmon")
hist(posterior[,'u[9]'], col = "darksalmon")
hist(posterior[,'u[10]'],col = "darksalmon")
hist(posterior[,'u[11]'],col = "darksalmon")
hist(posterior[,'u[12]'],col = "darksalmon")
hist(posterior[,'u[13]'],col = "darksalmon")






# bayesplot
library(bayesplot)
plot_title <- ggtitle("Posterior distributions", "with medians and 80% intervals")

x11()
mcmc_areas(posterior, pars = c("intercept", "beta_year[1]", "beta_trim[1]", "beta_destine[1]"), prob = 0.8) + plot_title


posterior2 <- extract(fit_mcmc, inc_warmup = TRUE, permuted = FALSE)

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(posterior2,  pars = c("beta_year[1]", "beta_trim[1]", "beta_destine[1]"), n_warmup = 300, facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)



y <- dat$cpue
yrep_mcmc = posterior_predict(fit_mcmc, draws = 500)
yrep_nb = posterior_predict(fit_nb, draws = 500)
dim(yrep_poisson)

get_posterior_mean(fit_mcmc)
plot(fit_mcmc)
log_posterior(fit_mcmc)




#=====================================================================================================
#             From now, we can do all the analysis like if we would have a fit model in stan
#=====================================================================================================
pairs(fit_mcmc, pars=names(obj$par))

## Trace plot
traceplot(fit_mcmc, pars=names(obj$par), inc_warmup=TRUE)

## 
# Extract posterior draws for later use
library(bayesplot)
available_mcmc(pattern = "_nuts_")

posterior_cp <- as.array(fit_mcmc)

lp_cp <- log_posterior(fit_mcmc)
head(lp_cp)

np_cp <- nuts_params(fit_mcmc)
head(np_cp)


color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp)

x11()
mcmc_pairs(posterior_cp, np = np_cp, pars = c("intercept", "logsigma", "logsigma_space"), off_diag_args = list(size = 0.75))


color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cp, pars = c("intercept", "logsigma", "logsigma_space") , np = np_cp) +
  xlab("Post-warmup iteration")


color_scheme_set("red")
mcmc_nuts_energy(np_cp)


ratios_cp <- neff_ratio(fit_mcmc)
print(ratios_cp)

mcmc_neff(ratios_cp, size = 2)


mcmc_acf(posterior_cp, pars = "intercept", lags = 10)






## Can extract marginal posteriors easily
post <- as.matrix(fit_mcmc)
hist(post[,'u[1]'])
hist(post[,'u[2]'])
hist(post[,'u[3]'])# random effect
hist(post[,'logsigma'])                   # fixed effect
hist(post[,'logsigma_space'])

## What if you want a posterior for derived quantities in the report? Just
## loop through each posterior sample (row) and call the report function
## which returns a list. The last column is the log-posterior density (lp__)
## and needs to be dropped
obj$report(post[1,-ncol(post)])         # sd0 is only element
sd0 <- rep(NA, len=nrow(post))
for(i in 1:nrow(post)){
  r <- obj$report(post[i,-ncol(post)])
  sd0[i] <- r$logsigma
}
hist(sd0)

names(as.data.frame(fit_mcmc))


