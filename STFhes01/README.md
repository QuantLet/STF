<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of QuantLet: STFhes01

Published in: Statistical Tools for Finance and Insurance

Description: 'Creates plots of the GBM, Heston spot and Heston volatility processes.'

Keywords: 'heston, simulation, geometric-brownian-motion, visualization, graphical representation, wiener-process, volatility'

See also: SFSbb, simGBM, simHeston

Author: Rafal Weron

Submitted: Tue, September 18 2012 by Dedy Dwi Prastyo

Input: 'SIGMA- volatility

Output: GBM and Heston dynamics

Example: 'User inputs the SFEWienerProcess parameters as sample input spot = 4 mu = .02 kappa = 2 theta = .04 sigma = .3 rho = -.05 days = 100 time = (0:days)/days STATE = [3621423255; 1292471671] randn("state",STATE) no = normrnd(0,1,length(time)-1,2); needed supplied functions are simGBM, simHeston.'

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STFhes01/plot.png" alt="Image" />
</div>

