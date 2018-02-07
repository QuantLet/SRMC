#' @title Consider the mix noraml-laplace distribution to model CRIX log return
#' Call estimation.epsilon.mu.std(crix.log.return) o estimate the parameters 
#' Estimate the corresp. parameters by EM algorithm for contaminated data
#' @param crix.log.return, the log return of CRIX, read CRIX (data.frame) which is given by CRIX.r
#' @return mu, std and epsilon vector for normal and laplace subsequently
#' @usage the parameters used to give estimations of VaR and ES 
#' 
#' @title Estimate the VaR and ES with input parameters given by estimation.epsilon.mu.std 
#' @call estimation.VaR.ES(x, epsilon.vector, mu.vector, std.vector, alpha)
#' @param x, the data to estimate its VaR and ES. Here x is the log return of crix price, given by CRIX.r
#' @param epsilon.vector, vector containing the weight of all components
#' @param alpha, the probability level for VaR and ES
#' @param std.vector, vector containing the standard deviations of each component
#' @param mu.vector, vector containing the mean of each component
#' 
#' @return a vector of estimations of VaR and ES by three methods
#' @return q.HS, the historical simulation estimation of VaR at alpha level
#' @return q.laplace, estimation of VaR based on laplace estimation at alpha/epsilon
#' @return q.mix, estimation of VaR based on mix model at alpha level 
#' @return ES.HS, ES.lapace, ES.mix, the alternative ES wih the same understanding as for VaR
#'  
#' @references see Theorem 2.1

# clear variables
rm(list = ls(all = TRUE))

# install packages
libraries = c("ggplot2", "dplyr", "reshape2")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {install.packages(x)
  })
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

options(scipen = 999)
set.seed(1)

# quantile of Laplace with parameter mean and std with alpha < 0.5
qlaplace = function(alpha, mean, std) ifelse(alpha < 0.5, mean + std * 1/sqrt(2) * log(2 * 
    alpha), mean + std * (-1/sqrt(2) * log(2 * (1 - alpha))))

# double-sided exponential distribution with parameter mean and std
plaplace = function(x, mean, std) {
    x.star = (x - mean)/std
    Laplace = ifelse(x.star < 0, 1/2 * exp(sqrt(2) * x.star), 1 - sqrt(2) * exp(-sqrt(2) * 
        x.star))
    return(Laplace)
}

# The cumulative distribution function of mixed distribution of symmetry normal and
# Laplace distribution with weight, mean and standard deviation Use of uniroot to give VaR
pmix = function(x, epsilon.vector, mu.vector, std.vector) epsilon.vector[1] * pnorm((x - 
    mu.vector[1])/std.vector[1]) + epsilon.vector[2] * plaplace(x, mu.vector[2], std.vector[2])


#' density of Laplace distribution with mean and standard deviation 
#' @param mean the mean or the symmetry center of L
#' @param std the standard deviation or scale parameter of L
dlaplace = function(x, mean, std) 
    1/(std * sqrt(2)) * exp(-sqrt(2)/std * abs(x - mean))

# lower partial moment of normal-laplace mixture at threshold u
Lower.partial.moment.mix = function(u, epsilon.vector, mu.vector, std.vector) {
    u.star = (u - mu.vector)/std.vector
    Lpm = epsilon.vector[1] * (std.vector[1] * (-dnorm(u.star[1])) + mu.vector[1] * pnorm(u.star[1]))
    +epsilon.vector[2] * (std.vector[2]/sqrt(2) * 1/2 * (sqrt(2) * u.star[2] - 1) * exp(sqrt(2) * 
        u.star[2]) + mu.vector[2] * 1/2 * exp(sqrt(2) * u.star[2]))
    return(Lpm)
}


#' Expectation Step of the EM Algorithm
#'
#' Calculate the posterior probabilities (soft labels) that each component
#' has to each data point.
#'
#' @param std.vector Vector containing the standard deviations of each component
#' @param mu.vector Vector containing the mean of each component
#' @param epsilon.vector Vector containing the mixing weights  of each component
#' @return Named list containing the loglik and posterior.df
e_step = function(x, mu.vector, std.vector, epsilon.vector) {
    comp1.prod = dnorm(x, mu.vector[1], std.vector[1]) * epsilon.vector[1]
    comp2.prod = dlaplace(x, mu.vector[2], std.vector[2]) * epsilon.vector[2]
    sum.of.comps = comp1.prod + comp2.prod
    comp1.post = comp1.prod/sum.of.comps
    comp2.post = comp2.prod/sum.of.comps
    
    sum.of.comps.ln = log(sum.of.comps, base = exp(1))
    sum.of.comps.ln.sum = sum(sum.of.comps.ln)
    
    list(loglik = sum.of.comps.ln.sum, posterior.df = cbind(comp1.post, comp2.post))
}

#' Maximization Step of the EM Algorithm
#' Update the Component Parameters
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' @return Named list containing the mean (mu), standard deviation (std), and mixing
#'   weights (epsilon) for each component.
m_step = function(x, posterior.df) {
    comp1.n = sum(posterior.df[, 1])
    comp2.n = sum(posterior.df[, 2])
    
    comp1.mu = 1/comp1.n * sum(posterior.df[, 1] * x)
    comp2.mu = 1/comp2.n * sum(posterior.df[, 2] * x)
    
    comp1.sd = sqrt(sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n)
    comp2.sd = sqrt(sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n)
    
    comp1.epsilon = comp1.n/length(x)
    comp2.epsilon = comp2.n/length(x)
    
    list(mu = c(comp1.mu, comp2.mu), std = c(comp1.sd, comp2.sd), epsilon = c(comp1.epsilon, 
        comp2.epsilon))
}

#' estimate parameters involved in normal-Laplace mixture for a given data
#' input x Vector data
#' output a iist including the epsilon, mu and sigma vector for normal and Laplace component

estimation.epsilon.mu.std = function(x) {
    x.kmeans = kmeans(x, 2)
    x.kmeans.cluster = x.kmeans$cluster
    
    x.df = data_frame(x = x, cluster = x.kmeans.cluster)
    
    x.df %>% mutate(num = row_number()) %>% ggplot(aes(y = num, x = x, color = factor(cluster))) + 
        geom_point() + ylab("Values") + ylab("Data Point Number") + scale_color_discrete(name = "Cluster") + 
        ggtitle("K-means Clustering")
    
    x.summary.df = x.df %>% group_by(cluster) %>% summarize(mu = mean(x), 
        variance = var(x), std = sd(x), size = n())
    
    x.summary.df %>% select(cluster, mu, variance, std)
    
    x.summary.df = x.summary.df %>% mutate(epsilon = size/sum(size))
    
    x.summary.df %>% select(cluster, size, epsilon)
    
    comp1.prod = dnorm(x = x, mean = x.summary.df$mu[1], sd = x.summary.df$std[1]) * 
        x.summary.df$epsilon[1]
    
    comp2.prod = dlaplace(x = x, mean = x.summary.df$mu[2], 
        std = x.summary.df$std[2]) * x.summary.df$epsilon[2]
    
    normalizer = comp1.prod + comp2.prod
    
    comp1.post = comp1.prod/normalizer
    comp2.post = comp2.prod/normalizer
    
    comp1.n = sum(comp1.post)
    comp2.n = sum(comp2.post)
    
    comp1.mu = 1/comp1.n * sum(comp1.post * x)
    comp2.mu = 1/comp2.n * sum(comp2.post * x)
    
    comp1.std = sqrt(sum(comp1.post * (x - comp1.mu)^2) * 1/comp1.n)
    comp2.std = sqrt(sum(comp2.post * (x - comp2.mu)^2) * 1/comp2.n)
    
    comp1.epsilon = comp1.n/length(x)
    comp2.epsilon = comp2.n/length(x)
    
    comp.params.df = data.frame(comp = c("comp1", "comp2"), comp.mu = c(comp1.mu, comp2.mu), 
        comp.std = c(comp1.std, comp2.std), comp.epsilon = c(comp1.epsilon, comp2.epsilon), 
        comp.cal = c("self", "self"))
    sum.of.comps = comp1.prod + comp2.prod
    sum.of.comps.ln = log(sum.of.comps, base = exp(1))
    sum(sum.of.comps.ln)
    
    for (i in 1:50) {
        if (i == 1) {
            # Initialization
            e.step = e_step(x, x.summary.df[["mu"]], x.summary.df[["std"]], 
                x.summary.df[["epsilon"]])
            m.step = m_step(x, e.step[["posterior.df"]])
            cur.loglik = e.step[["loglik"]]
            loglik.vector = e.step[["loglik"]]
        } else {
            # Repeat E and M steps till convergence
            e.step = e_step(x, m.step[["mu"]], m.step[["std"]], m.step[["epsilon"]])
            m.step = m_step(x, e.step[["posterior.df"]])
            loglik.vector = c(loglik.vector, e.step[["loglik"]])
            loglik.diff = abs((cur.loglik - e.step[["loglik"]]))
            if (loglik.diff < 1e-06) {
                break
            } else {
                cur.loglik = e.step[["loglik"]]
            }
            m.step
        }
    }
    list(mu = m.step$mu, std = m.step$std, epsilon = m.step$epsilon)
}

#' estimate VaR and ES of data x at alpha level assuming that the pre-supposed model
#' is a normal-Laplace mixture with parameters epsilon.vector, mu.vector, std.vector
estimation.VaR.ES = function(alpha, x, epsilon.vector, mu.vector, std.vector) {
    q.HS = quantile(x, alpha)
    ES.HS = mean(x[which(x < q.HS | x == q.HS)])
    # Empirical estimation of VaR and ES
    
    u = quantile(x, alpha)
    u.star = (u - mu.vector)/std.vector
    u.star.laplace = u.star[2]
    std.laplace = std.vector[2]
    mu.laplace = mu.vector[2]
    scaled.alpha = alpha/epsilon.vector[2]
    # u.star, the standardized vector of threshold with u, historical simulation of VaR
    
    q.laplace = qlaplace(alpha = scaled.alpha, mean = mu.laplace, std = std.laplace)
    ES.laplace = (std.laplace/sqrt(2) * (u.star.laplace * sqrt(2) - 1)/2 * exp(u.star.laplace * 
        sqrt(2)) + mu.laplace * exp(u.star.laplace * sqrt(2))/2)/scaled.alpha
    # Estimation of VaR and ES based on Laplace at the scaled level alpha/epsilon For the
    # ES.laplace, we keep the VaR from the estimation by historical simulation
    
    ff = function(x, epsilon.vector, mu.vector, std.vector, alpha) pmix(x, epsilon.vector, 
        mu.vector, std.vector) - alpha
    
    q.mix = uniroot(ff, c(-10, 10), epsilon.vector, mu.vector, std.vector, alpha)$root
    ES.mix = 1/alpha * Lower.partial.moment.mix(u, epsilon.vector, mu.vector, std.vector)
    # Calculate the VaR and ES by mixed normal-laplace at alpha level
    
    estimations = c(q.HS, q.laplace, q.mix, ES.HS, ES.laplace, ES.mix)
    names(estimations) = c("q.HS", "q.laplace", "q.mix", "ES.HS", "ES.lapace", "ES.mix")
    return(round(estimations, 3))
}

# define the period of log returns of CRIX
date       = read.csv(file = 'SRMC_crix')$date 
sub.date.1 = 1: which(date == '2016-04-01')
sub.date.2 = which(date == '2016-04-01'): length(date)

# asign data samples for the wholeperiod, period.1 and period.2
crix.log.return   = read.csv(file = "SRMC_crix")$price
crix.log.return.1 = crix.log.return[sub.date.1] 
crix.log.return.2 = crix.log.return[sub.date.2]

# set seeds for estimations of parameters involved in normal-Laplace mixture
# consider the three periods
set.seed(10)
param   = estimation.epsilon.mu.std(crix.log.return)
param.1 = estimation.epsilon.mu.std(crix.log.return.1)
param.2 = estimation.epsilon.mu.std(crix.log.return.2)

# consider three probability levels for VaR and ES
alpha.vector   = c(0.005, 0.01, 0.05)

# asign estimated parameters for whole period 
epsilon.vector = param$epsilon
mu.vector      = param$mu
std.vector     = param$std

# asign estimated parameters for period.1
epsilon.vector.1 = param.1$epsilon
mu.vector.1      = param.1$mu
std.vector.1     = param.1$std

# asign estimated parameters for period.2
epsilon.vector.2 = param.2$epsilon
mu.vector.2     = param.2$mu
std.vector.2     = param.2$std

# give estimations of VaR, ES by historical simulations, based on Laplace, mixted distribution for three periods
estimation   = sapply(alpha.vector, estimation.VaR.ES, crix.log.return, epsilon.vector, mu.vector, std.vector)
estimation.1 = sapply(alpha.vector, estimation.VaR.ES, crix.log.return.1, epsilon.vector.1, mu.vector.1, std.vector.1)
estimation.2 = sapply(alpha.vector, estimation.VaR.ES, crix.log.return.2, epsilon.vector.2, mu.vector.2, std.vector.2)

# give estimations of epsilon, mu and sigma for three periods
es.epsilon  = round(c(epsilon.vector, epsilon.vector.1, epsilon.vector.2), 3) 
es.mu       = round(c(mu.vector, mu.vector.1, mu.vector.2), 3) 
es.std      = round(c(std.vector, std.vector.1, std.vector.2), 3) 

# print estimated parameter
es.parameter           = cbind(es.epsilon, es.mu, es.std)
rownames(es.parameter) = paste(c(rep("period", 2), rep("period.1", 2), rep("period.2", 2)), rep(c('.normal', '.Laplace'), 3))
print(es.parameter)

# print estimated VaR and ES
es.VaR.ES       = round(cbind(estimation, estimation.1, estimation.2), 3) 
alpha           = t(rep(alpha.vector, 3))
rownames(alpha) = 'alpha'
output          = t(rbind(alpha, es.VaR.ES))
print(output)

write.csv(es.parameter, file = "SRMC_crix_es.parameter")
write.csv(output, file = "SRMC_crix_es.VaR.ES", row.names = FALSE)



