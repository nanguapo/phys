import statsmodels.api as sm
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.iolib.table import (SimpleTable, default_txt_fmt)


#input data - paste from speadsheet
y = [-4.144195406,
-4.615789451,
-4.87771676,
-5.037329818,
-5.6008665,
-5.704750874,
-5.803858765,
-5.885170468,
-5.910485919,
-6.027966746,
-6.019148544,
-6.09361663]

x = [4.802230098,
5.500022754,
5.841450923,
6.018876143,
6.657887803,
6.765474868,
6.871173241,
6.990135867,
7.013982019,
7.100809039,
7.169457786,
7.249929798]

# Uncertainty - the default is relative uncertainty
sigma_y = np.array([2.57E-03,
6.33E-03,
3.86E-03,
1.44E-02,
7.93E-03,
1.46E-02,
8.26E-03,
1.03E-02,
9.01E-03,
2.96E-02,
2.76E-02,
7.96E-03])
w = 1 / (sigma_y**2)

#add a constant
X = sm.add_constant(x)


mod_wls = sm.WLS(y, X, weights=w)
res_wls = mod_wls.fit()

#override the scaling so that it reads the weights as absolute weights rather than relative weights and gives the absolute uncertainties of the slope and intercept
res_wls_absolute = mod_wls.fit(cov_type='fixed scale', cov_kwds={'scale': 1.0})

intercept = res_wls_absolute.params[0]
slope = res_wls_absolute.params[1]

print("\n Absolute Weights")
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

plt.title('Weighted Least Squares (WLS) Regression Fit for FIPPS Eu Relative Efficiency w/o 1805 Line')
plt.xlabel('Ln(Energy in keV)')
plt.ylabel('Ln(Counts Per Yield)')
plt.colorbar(label='Observation Weight')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()

#turning it back into regular graph for easy legibility of the energy
raw_energy = np.exp(x)
raw_counts = np.exp(y)

#smooth x array
raw_energy_smooth = np.linspace(min(raw_energy), max(raw_energy), 500)
#curve of best fit
raw_counts_fit = np.exp(intercept) * (raw_energy_smooth ** slope)

plt.figure(figsize=(8, 5))
# plot the raw data points
plt.scatter(raw_energy, raw_counts, c=w, cmap='viridis', s=100, label='Raw Data (color = weight)', edgecolor='k', zorder=3)
# plot the curved line of best fit
plt.plot(raw_energy_smooth, raw_counts_fit, color='red', linewidth=2, label='WLS Fitted Power-Law Curve')

plt.title('Linear Scale: FIPPS Eu Relative Efficiency w/o 1805 Line')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts Per Yield')
plt.colorbar(label='Observation Weight')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()
