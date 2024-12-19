<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of QuantLet: STF2svm01

Published in: Statistical Tools for Finance and Insurance

Description: 'Plots the area of two different groups. It refers to the four plots in the SVM chapter.'

Keywords: SVM, risk management, kernel, default, classification

See also: MVAsvmSig05C200, MVAsvmSig100C1, MVAsvmSig2C1, SMSsvmbankrupt, SMSsvmspiral

Author: Linda Hoffmann, Awdesch Melzer

Submitted: Tue, October 15 2013 by Awdesch Melzer

Datafile: Training100by100noNA.txt

Input: 'G- Data Matrix, needs to include labels

Output: 'Four 2-DIm Plot of a svm classification. The RBF kernel and hinge loss are used. The triangles represent solvent (y = -1) and circles represent insolvent companies (y = +1). The shaded object represent the support vectors. The colored background corresponds to different score values $f$. The bluer the area, the higher the score (y= +1). Most successful companies (y=-1) are in the red area.'

Example: 'Plots for commands: stf2svm01(G,28,7,3,1,100), stf2svm01(G,28,7,3,1,2),  stf2svm01(G,28,7,3,1,1/2) and stf2svm01(G,28,7,3,200,2)'

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STF2svm01/plot1.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STF2svm01/plot2.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STF2svm01/plot3.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/STF/master/STF2svm01/plot4.png" alt="Image" />
</div>

