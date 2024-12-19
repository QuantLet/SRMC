<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: SRMCappr_epsilon_alpha

Published in: unpublished

Description: 'Shows the approximations of VaR, ES with varying alpha and epsilon, a power function of alpha with exponent tau,

Keywords: 'normal mixture, Laplace, heavy tail, approximation, Expected shortfall, VaR, Huber contamination'

Author: Chengxiu Ling

See also: 'SRMCappr_alpha,  SRMCappr_epsilon'

Submitted: 07/02/2018

Input: 

- alpha: 'the probability level vector for VaR and ES'

- index: 'any subset of (''Normal'', ''Laplace'', ''Power-like''), indicating the choice of contamination model'

- scale: 'the scale vector of the contamination model indicated by index'

- tau: 'the power indice alpha to define epsilon, i.e., epsilon = alpha^\tau'

Output: 'Approximations and its theoretical values of VaR, ES of epsilon mixture model with epsilon = alpha^tau'

Example: 

- 1: 'SRMCappr_epsilon_alpha_1.png, approximations of VaR, ES  based on normal model with scale = c(1.1, 1.6, 1)

- 2: 'SRMCappr_epsilon_alpha_2.png, approximations of VaR, ES  based on contamination model with scale = (1.6, 0.95, 1), tau = (0.5, 0.1, 0.5)

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/SRMC/master/SRMCappr_epsilon_alpha/SRMCappr_epsilon_alpha_1.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/SRMC/master/SRMCappr_epsilon_alpha/SRMCappr_epsilon_alpha_2.png" alt="Image" />
</div>

