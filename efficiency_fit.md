Instructions for the efficiency_fit.py script

What it did.
This code was based of off Basile's code that he wrote for the efficiency calibration of the Zr FIPPS campaigns.

A few changes were made. 
Previously, he included the calculations for the efficiency of Eu and Al alongside this curve fitting process which only could be used for absolute efficiency. 
It also was set to fit a fourth degree polynomial.

What it does.
My changes

Now it can be used for both absolute and relative efficiency calculations being agnostic to the input data. Absolute and relative efficiency calculations and determination are handled in a google sheet for legibility/understanding. 
It fits the polynomial based on the data via an F-statistic (F for Fischer) test which determines whether a higher degree polynomial should be accepted based off of the chi^2. Note that here even though unlikely a 5th or 6th degree polynomial might cause some issues so be wary if you get to that point.
It outputs the degree of the polynomial as well as the chi^2, the residuals, mu, sig(standard deviation),various coefficients, and equations for calculations in the sheet.