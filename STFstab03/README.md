<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of QuantLet: STFstab03

Published in: Statistical Tools for Finance and Insurance

Description: STFstab03.R compares S0 and S parameterizations and plots both. STFstab03.m shows (top panels) semilog and loglog plots of symmetric 1.7-stable, symmetric tempered stable (TSD) with alpha = 1.7 abd lambda = 0.2, and Gaussian pdfs. Bottom panels show semilog and loglog plots of symmetric TSD pdfs with alpha = 1.7 and four truncation coefficients: lambda = 5, 0.5, 0.2, 0.01. Note, that for large lambdas the distribution approaches the Gaussian (though with a different scale) anf for small lambdas the stable law with the same shape parameter alpha. Requires "stabpdf_fft.m", "tstabpdf.m", and "tstabpdf3.m" functions.

Keywords: graphical representation, visualization, stable distribution, pdf, simulation

See also: STFstab01, STFstab02, stabpdf_fft, tstabpdf, tstabpdf3

Author: Zografia Anastasiadou, Rafal Weron

Submitted: Tue, September 18 2012 by Dedy Dwi Prastyo

Example: STFstab03.R produces the plots of the S and S0 parameterizations. stabpdf_fft.m, tstabpdf3.m, tstabpdf.m in STFstab03.m  produce the panel plots.

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STFstab03/plotR1.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STFstab03/plotR2.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STFstab03/plotm.png" alt="Image" />
</div>

