<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: SRMCappr_epsilon

Published in: unpublished

Description: Shows the approximations of VaR, ES and level ratio of expectile vs. quantile with varying epsilon and given alpha

Keywords: normal mixture, Laplace, heavy tail, approximation, Expected shortfall, VaR, expectile, 

Author: Chengxiu Ling

See also: SRMCappr_alpha,  SRMCappr_epsilon_alpha

Submitted: 07/02/2018

Input: 
- alpha: 'the probability level vector for VaR and ES
- index: 'any subset of (''Normal'', ''Laplace'', ''Power-like''), indicating the choice of contamination model
- scale: 'the scale vector of the index function
- epsilon: 'the contamination level vector

Output: Approximations and its theoretical values of VaR, ES and level ratio of expectile vs. quantile with varying alpha

Example: SRMCappr_epsilon.png, plot of approximations of VaR, ES and w_alpha / alpha with alpha = (0.5%, 0.5%, 0.5%),

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/SRMC/master/SRMCappr_epsilon/SRMCappr_epsilon.png" alt="Image" />
</div>

