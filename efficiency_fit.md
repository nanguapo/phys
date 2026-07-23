Instructions for the efficiency_fit.py script

What it did.
This code was based of off Basile's code that he wrote for the efficiency calibration of the Zr FIPPS campaigns.

A few changes were made. 
Previously, he included the calculations for the efficiency of Eu and Al alongside this curve fitting process which only could be used for absolute efficiency. 
It also was set to fit a fourth degree polynomial. He used the uncertainties library for erro propagation which is very convenient.

What it does.
My changes

Now it can be used for both absolute and relative efficiency calculations being agnostic to the input data. Absolute and relative efficiency calculations and determination are handled in a google sheet for legibility/understanding. 
It fits the polynomial based on the data via an F-statistic (F for Fischer) test which determines whether a higher degree polynomial should be accepted based off of the chi^2. Note that here even though unlikely a 5th or 6th degree polynomial might cause some issues so be wary if you get to that point.

Required Inputs: 
  Energy, Efficiency or Counts per yield, Efficiency or Coutns per Yield absolute error
  Also change line 176 based off whether you're doing an absolute or relative efficiency calibration so that the graph's y axis reflects      this.
Outputs:
  A pretty lin lin plot so that you can see the energy values clearly
  The coefficients
  Degree of the polynomial
  Reduced chi^2
  Residuals
  Mu
  Sig (standard deviation)
  Uncertainties scaled by the square root of the reduced chi^2
  Equations
  

How to get the relative efficiency at certain lines:
  Navigate to this sheet: https://docs.google.com/spreadsheets/d/1lpTHxNY7b5FSSnZFpxfjWzzeSIVC-zpX/edit?usp=sharing&ouid=114456628968102514395&rtpof=true&sd=true and go to the NEW RELATIVE EFFICIENCY sheet
  Place the inputs into their respective cells. If you don't have an input like v10 etc that is because that input is for a higher order polynomial and thus the cell can be left blank and will not affect the result. The efficiency values should autopopulate based off of the formulae once the energy values are input.

Good luck.
  
  
  
