import statsmodels.api as sm
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.iolib.table import (SimpleTable, default_txt_fmt)


#input data - paste from speadsheet
y = [14.8024786,
14.29593546,
14.10630147,
13.94759451,
13.86945773,
13.49573418,
13.36967609,
13.3422941,
13.2945571,
13.27706317,
13.25198225,
13.08119509]
x = [4.804266916,
5.500849964,
5.842152144,
6.019493053,
6.096387467,
6.658421745,
6.766018889,
6.871672039,
6.990615552,
6.994244853,
7.014481861,
7.250408375]

# Uncertainty - the default is relative uncertainty
sigma_y = np.array([0.001584706691,
0.003431322344,
0.006262742421,
0.008460420437,
0.007120964522,
0.003922159274,
0.007666515422,
0.003623833747,
0.004493234363,
0.01359314576,
0.003898803912,
0.003241638153])
w = 1 / (sigma_y**2)

#add a constant
X = sm.add_constant(x)


mod_wls = sm.WLS(y, X, weights=w)
res_wls = mod_wls.fit()

#override the scaling so that it reads the weights as absolute weights rather than relative weights and gives the absolute uncertainties of the slope and intercept
res_wls_absolute = mod_wls.fit(cov_type='fixed scale', cov_kwds={'scale': 1.0})

print("\n Absolute Weights (True Instrumental Propagation) ---")
print(f"Intercept Error: {res_wls_absolute.bse[0]}")
print(f"Slope Error: {res_wls_absolute.bse[1]}")

# Extract the full covariance matrix
cov_matrix = res_wls_absolute.cov_params()

# Print the covariance matrix 
print(cov_matrix)

#Uncertainty and covariance terms
var_A = cov_matrix[0, 0]  # This is (delta A)^2 at row 0 column 0
var_B = cov_matrix[1, 1]    # This is (delta B)^2 at row 1 column 1
cov_AB = cov_matrix[0, 1]  # This is Cov(A, B) at row 0 column 1

print(f"Intercept Variance: {var_A}")
print(f"Slope Variance:     {var_B}")
print(f"Covariance (A, B):  {cov_AB}")

print(res_wls_absolute.summary())

plt.figure(figsize=(8, 5))
plt.scatter(x, y, c=w, cmap='viridis', s=100, label='Data (color = weight)', edgecolor='k', zorder=3)
plt.plot(x, res_wls.fittedvalues, color='red', linewidth=2, label='WLS Fitted Line')

plt.title('Weighted Least Squares (WLS) Regression Fit for 7/03/26 Simulated Eu Gas Relative Efficiency')
plt.xlabel('Ln(Energy in keV)')
plt.ylabel('Ln(Counts Per Yield)')
plt.colorbar(label='Observation Weight')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()