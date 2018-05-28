[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SRMC_RelativeError_epsilon** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: SRMC_RelativeError_epsilon

Published in: Sensitivity_RiskMeasure_HuberContamination

Description: 'Shows small relative error of VaR, ES with the estimation given by the influence function and small epsilon'

Keywords: 'influence function, relative error, normal distribution, Laplace, heavy tail, Huber contamination'

Author: Chengxiu Ling

See also: 'SRMC_crix, SRMC_Estimation_crix'

Submitted: 07/02/2018

Input: 'SRMC_crix.csv'

Output: 'Table 1, Give VaR, ES and its relative error with alpha = rep(0.25, 3) and scale = c(2, 1.2, 1)'

```

### R Code
```r

#' @title Give VaR, ES and its relative error for given alpha and small varying epsilon
#' Here we consider some of the contamination model H being normal, laplace or Power-like mix distribution

#' @call IF.mix(epsilon, index, scale, alpha), the function gives data in Table 
#' @param epsilon, the contamination level should be small to show its infinestimal influence
#' @param index, the mixed type vector being any combination of 'Normal', 'Laplace', 'Power-like'
#' @param scale, the parameter vector of the contamination model, the same dim of index
#' @param alpha, the level vrctor of VaR and ES, could be any value in (0,1), the same dim of index
#' @return the input and VaR.RB, ES.RB for VaR and ES and its relative error as well
#' @references see Theorem 2.5
#' 
#' clear variables
rm(list = ls(all = TRUE))

# Standard double-sided exponential distribution with mean zero and variance 1
plaplace = function(x) {
    Laplace = ifelse(x < 0, 1/2 * exp(sqrt(2) * x), 1 - sqrt(2) * exp(-sqrt(2) * x))
    return(Laplace)
}

# Standard Power tail with \gamma = 1/2, with mean zero and infinitely variance, see
# Example 3.3
ppower_like = function(x) {
    Power = ifelse(x < 0, 1/2 * (1 + x/sqrt(4 + x^2)), 1 - 1/2 * (1 - x/sqrt(4 + x^2)))
    return(Power)
}

# cumulative distribution function of mixed model
mix_cdf = function(x, index, epsilon, scale) {
    if (index == "Normal") {
        mix_cdf = (1 - epsilon) * pnorm(x) + epsilon * pnorm(x/scale)
    } else if (index == "Laplace") {
        mix_cdf = (1 - epsilon) * pnorm(x) + epsilon * plaplace(x/scale)
    } else {
        mix_cdf = (1 - epsilon) * pnorm(x) + epsilon * ppower_like(x/scale)
    }
    return(mix_cdf)
}

# Lower partial moment function at threshold x
Lower.partial.moment = function(x, index, epsilon, scale) {
    if (index == "Normal") {
        Lower.partial.moment = (1 - epsilon) * (-dnorm(x)) + epsilon * (-scale * dnorm(x/scale))
    } else if (index == "Laplace") {
        Lower.partial.moment = (1 - epsilon) * (-dnorm(x)) + epsilon * (scale/sqrt(2) * 
            1/2 * (sqrt(2)/scale * x - 1) * exp(sqrt(2)/scale * x))
    } else {
        Lower.partial.moment = (1 - epsilon) * (-dnorm(x)) + epsilon * (-scale * 2/sqrt(4 + 
            (x/scale)^2))
    }
    return(Lower.partial.moment)
}

# Lower partial moment function at threshold x with possible different epsilon
Lower.partial.moment.epsilon = function(epsilon, x, index, scale) {
    if (index == "Normal") {
        Lower.partial.moment = (1 - epsilon) * (-dnorm(x)) + epsilon * (-scale * dnorm(x/scale))
    } else if (index == "Laplace") {
        Lower.partial.moment = (1 - epsilon) * (-dnorm(x)) + epsilon * (scale/sqrt(2) * 
            1/2 * (sqrt(2)/scale * x - 1) * exp(sqrt(2)/scale * x))
    } else {
        Lower.partial.moment = (1 - epsilon) * (-dnorm(x)) + epsilon * (-scale * 2/sqrt(4 + 
            (x/scale)^2))
    }
    return(Lower.partial.moment)
}

# VaR of mixed model
mix_VaR = function(alpha, index, epsilon, scale) {
    ff = function(x, alpha, index, epsilon, scale) mix_cdf(x, index, epsilon, scale) - 
        alpha
    VaR = uniroot(ff, c(-20, 10), alpha, index, epsilon, scale)$root  #the accurate alpha-th quantile
    return(VaR)
}

# expected shortfall of mixed model
mix_ES = function(alpha, index, epsilon, scale) {
    quantile.epsilon = mix_VaR(alpha, index, epsilon, scale)
    mix_ES = alpha^(-1) * Lower.partial.moment(quantile.epsilon, index, epsilon, scale)  #the accurate alpha-th Expected-Shortfall
    return(mix_ES)
}

# VaR of mixed model
mix_VaR_epsilon = function(epsilon, index, alpha, scale) {
    ff = function(x, alpha, index, epsilon, scale) mix_cdf(x, index, epsilon, scale) - 
        alpha
    VaR = uniroot(ff, c(-20, 10), alpha, index, epsilon, scale)$root  #the accurate alpha-th quantile
    return(VaR)
}

# expected shortfall for mixted model with epsilon etc 
mix_ES_epsilon = function(epsilon, index, alpha, scale) {
    quantile.epsilon = mix_VaR(alpha, index, epsilon, scale)
    mix_ES = alpha^(-1) * Lower.partial.moment(quantile.epsilon, index, epsilon, scale)  #the accurate alpha-th Expected-Shortfall
    return(mix_ES)
}

# Give the infulence function of VaR and ES for the contamination model (index) Save
# it as IF_VaR_ES, see Remark 2.6
IF.mix = function(epsilon, index, scale, alpha) {
    size.epsilon = length(epsilon)
    size.index = length(index)
    
    IF_VaR_ES = rep(1, size.index * 2)
    names(IF_VaR_ES) = paste(c(rep("IF of VaR", size.index), rep("IF of ES", size.index)), 
        rep(index, 2))
    VaR_Normal = qnorm(alpha)
    
    # Calculate influence function of VaR and ES for every contamination indicated by index 
    # at given alpha level
    for (i in 1:size.index) {
        IF_VaR_ES[i] = (alpha[i] - mix_cdf(VaR_Normal[i], index[i], 1, scale[i]))/dnorm(VaR_Normal[i])
        Lpm_H = Lower.partial.moment.epsilon(1, VaR_Normal[i], index[i], scale[i])
        Lpm_Normal = Lower.partial.moment.epsilon(0, VaR_Normal[i], index[i], 1)
        IF_VaR_ES[i + 3] = (VaR_Normal[i] * (alpha[i] - mix_cdf(VaR_Normal[i], index[i], 
            1, scale[i])) + Lpm_H - Lpm_Normal)/alpha[i]
    }
    
    # Give the true VaR and its relative bias (RB) for contamination model (indicated by
    # index) Save it in VaR.RB as a matrix of size.epsilon *(2*size.index) The same for
    # ES.RB
    VaR.RB = matrix(rep(1, size.epsilon * size.index * 2), nrow = size.epsilon)
    ES.RB = matrix(rep(1, size.epsilon * size.index * 2), nrow = size.epsilon)
    colnames(VaR.RB) = paste(c(rep("VaR of ", size.index), rep("RB(VaR) of", size.index)), 
        rep(index, 2))
    colnames(ES.RB) = paste(c(rep("ES of", size.index), rep("RB(ES) of", size.index)), 
        rep(index, 2))
    for (i in 1:size.index) {
        VaR.RB[, i] = sapply(epsilon, mix_VaR_epsilon, index[i], alpha[i], scale[i])
        VaR_appro = VaR_Normal[i] + IF_VaR_ES[i] * epsilon
        VaR.RB[, i + size.index] = (VaR_appro - VaR.RB[, i])/VaR.RB[, i]
        
        ES.RB[, i] = sapply(epsilon, mix_ES_epsilon, index[i], alpha[i], scale[i])
        Lpm_Normal = Lower.partial.moment.epsilon(0, VaR_Normal[i], index[i], 1)
        ES_appro = Lpm_Normal/alpha[i] + IF_VaR_ES[i + 3] * epsilon
        ES.RB[, i + size.index] = (ES_appro - ES.RB[, i])/ES.RB[, i]
    }
    list(index = index, alpha = alpha, scale = scale, epsilon = epsilon, IF_VaR_ES = IF_VaR_ES, 
        VaR.RB = VaR.RB, ES.RB = ES.RB)
}

#' asign specific input to the parameter and give the results of IF etc
epsilon = seq(from = 0.001, to = 0.1, by = 0.01)
index   = c("Normal", "Laplace", "Power-like")
alpha   = rep(0.25, 3)
scale   = c(2, 1.2, 1)

#' call IF.mix to give the VaR, ES and its relative error
output  = IF.mix(epsilon, index, scale, alpha)

#' give Table 1 including VaR, ES and its relative error for normal, Laplace
#' and Power-like contamination model with varying epsilon
table   = cbind(output$epsilon, output$VaR.RB, output$ES.RB)
write.csv(table, file = "SRMC_RelativeError_epsilon_table", row.names = FALSE)


```

automatically created on 2018-05-28