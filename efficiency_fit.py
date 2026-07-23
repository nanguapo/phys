#author: Julian Herrera - jrherrera@udallas.edu . Previous author is Basile
#libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import uncertainties
from uncertainties import unumpy
import io
from scipy import stats

#paste your data here from the sheet as a triple-quoted string
pasted = """Energy, y, y_err
121.7817,2699142.203, 2842.84967
244.6975, 1623185.43, 10234.17028
344.2785, 1332957.879, 60203.61625
411.1163, 1139164.433, 10895.65954
778.904, 724631.4572, 4258.32842
867.378, 640591.5626, 6013.215981
964.079, 627291.7241, 3441.233682
1112.074, 578959.7315, 3428.096311
1212.948, 513947.7401, 8557.606572
1299.14, 522841.3962, 7046.04422
1408.006, 483285.8513, 2410.150622
"""
df = pd.read_csv(io.StringIO(pasted), skipinitialspace=True)
E, y, y_err = (df[c].values for c in ["Energy","y","y_err"])


#log by log prep
#center and scale the basis and save for later
t = np.log(E); mu, sig_t = t.mean(), t.std()
#fit in log space: sigma on ln(y) is the relative error on y
ly = np.log(y); ly_err = y_err / y

#fit plot range to the data, points for a smooth curve
xfit = np.linspace(E.min(),E.max(), 1000) 


#function for curve fitting to a polynomial; weighted fit of ln(eff) at fixed degree; log by log but outputs it in lin vs lin for legibility
def fitter(xfit, E, ly, eq, ly_err, degree):
    popt, pcov = curve_fit(eq, E, ly, p0=np.zeros(degree+1),
                           sigma=ly_err, absolute_sigma=True)
    
    fit = np.exp(eq(xfit, *popt))

    pars = uncertainties.correlated_values(popt, pcov)
    uncertain = unumpy.exp(eq(xfit, *pars))
    nom = unumpy.nominal_values(uncertain)
    std = unumpy.std_devs(uncertain)

    

    return fit, nom, std, popt, pcov, pars

#fit polynomial - new u as the pivot for 1st degree
def polyLog(E, *c):
    # u is a rescaled energy axis that is scale to be 0 at the middle of our data set (circa 635 keV) to minimize the covariance, one unit equals one standard deviation
    u = (np.log(E) - mu) / sig_t
    return sum(ck * u**k for k, ck in enumerate(c))

#choose the polynomial degree by nested F-test instead of guessing/hard-coding !!!
def select_degree(E, ly, ly_err, eq, max_degree=6, alpha=0.05):
    N = len(E)
    max_degree = min(max_degree, N - 3)          #keep at least 2 degrees of freedom
    def chi2_at(d):
        p, _ = curve_fit(eq, E, ly, p0=np.zeros(d+1),
                         sigma=ly_err, absolute_sigma=True)
        return float(np.sum(((ly - eq(E, *p)) / ly_err)**2))
    degree = 1
    chi2 = chi2_at(degree)
    while degree < max_degree:
        chi2_next = chi2_at(degree + 1)
        dof_next = N - (degree + 2)
        F = (chi2 - chi2_next) / (chi2_next / dof_next) # good old pal Fischer
        if stats.f.sf(F, 1, dof_next) < alpha:    #is the extra term justified?
            degree, chi2 = degree + 1, chi2_next
        else:
            break
    return degree


#outputs for the line of best fit - take these and input them into the google sheet (MAKE SURE YOU USE ALL DIGITS)
def report(E, mu, sig_t, popt, pcov, red_chi2, polyLog, inflate=True, #inflate =True scales covariance by the reduced chi^2 so that uncertainty represents observed scatter. It can be set to false if chi^2 is around 1
           checks=(200., 661.7, 1173.2)):

    C = pcov * (red_chi2 if inflate else 1.0)
    d = len(popt) - 1

    pars = uncertainties.correlated_values(popt, C)
    u_s = np.linspace(-2.5, 2.5, 2*d + 1)
    var_s = np.array([polyLog(np.exp(mu + u*sig_t), *pars).std_dev**2 for u in u_s])
    var_coeffs = np.polyfit(u_s, var_s, 2*d)[::-1]

    #making the outputs readable for transciption
    print("  ln(eff) = " + " + ".join(f"c{k}*u^{k}" for k in range(d+1)).replace("*u^0", ""))
    print(f"  u       = (LN(E) - mu)/sig        pivot = {np.exp(mu):.6g} keV")

    if d == 1:      # the special case of a power law in the instance of a first order polynomial
        print(f"\n  i.e.  eff = {np.exp(popt[0]):.6g} * (E/{np.exp(mu):.6g})^({popt[1]/sig_t:.5g})")
        print(f"        power-law exponent = {popt[1]/sig_t:+.5g} "
              f"+/- {np.sqrt(C[1,1])/sig_t:.4g}")

    print(f"\n  F1  mu  = {mu:.10g}")
    print(f"  F2  sig = {sig_t:.10g}")
    for k, c in enumerate(popt):
        print(f"  G{k+1}  c{k}  = {c:.10g}   +/- {np.sqrt(C[k,k]):.6g}")
    for m, vm in enumerate(var_coeffs):
        print(f"  H{m+1}  v{m}  = {vm:.10g}")

    poly_c = " + ".join(f"$G${k+1}*B1^{k}" for k in range(d+1)
                        ).replace("*B1^0", "").replace("B1^1", "B1")
    poly_v = " + ".join(f"$H${m+1}*B1^{m}" for m in range(2*d+1)
                        ).replace("*B1^0", "").replace("B1^1", "B1")
    print("\n  sheet formulas (energy in A1):")
    print("    B1  u     =(LN(A1)-$F$1)/$F$2")
    print(f"    C1  eff   =EXP({poly_c})")
    print(f"    D1  sigma =C1*SQRT({poly_v})")

    print(f"\n  valid {E.min():.2f} to {E.max():.2f} keV -- do not extrapolate")
    if inflate:
        print(f"  uncertainties scaled by sqrt(red.chi2) = {np.sqrt(red_chi2):.3f}")

    print("\n  CHECK -- sheet must reproduce:")
    print(f"  {'E (keV)':>9} {'eff':>15} {'1-sigma':>13} {'rel':>7}")
    for E0 in checks:
        if not (E.min() <= E0 <= E.max()):
            continue
        u = (np.log(E0) - mu) / sig_t
        eff = np.exp(sum(c*u**k for k, c in enumerate(popt)))
        rel = np.sqrt(sum(vm*u**m for m, vm in enumerate(var_coeffs)))
        print(f"  {E0:>9.2f} {eff:>15.6e} {eff*rel:>13.4e} {100*rel:>6.2f}%")


#call curve fitting function    
degree = select_degree(E, ly, ly_err, polyLog)                                     
eff_fit, eff_nom, eff_std, popteff, pcoveff, parseff = fitter(xfit, E, ly, polyLog, ly_err, degree)

#residuals and chi^2 in log space
resid = (ly - polyLog(E, *popteff)) / ly_err
red_chi2 = np.sum(resid**2) / (len(E) - len(popteff))
print(f"selected degree : {degree}")
print(f"reduced chi^2   : {red_chi2:.2f}")
if red_chi2 > 2 or red_chi2 < 0.3:
    print("  NOTE: far from 1. A TREND in the residuals means the shape is wrong;")
    print("  large random signs mean the error bars are incomplete; one or two")
    print("  big outliers mean point-specific physics (e.g. coincidence summing) so check those points.")

print("\nresiduals (sorted by chi^2 contribution):")
for i in np.argsort(-resid**2):
    print(f"  {E[i]:9.4f}  {resid[i]:+6.2f} sigma")

print()
report(E, mu, sig_t, popteff, pcoveff, red_chi2,polyLog)

#rebuild with the inflated covariance so plot, box, and report all agree
parseff = uncertainties.correlated_values(popteff, pcoveff * red_chi2)
band = unumpy.exp(polyLog(xfit, *parseff))
eff_nom, eff_std = unumpy.nominal_values(band), unumpy.std_devs(band)


#parameter box for graph
parameterseff = "\n".join(
    rf"$c_{{{k}}} = {p.nominal_value:.2e} \pm {p.std_dev:.2e}$"
    for k, p in enumerate(parseff)
)

#efficiency and 1-sigma at any energy
def efficiency_at(energy):
    ln_eff = polyLog(np.atleast_1d(np.asarray(energy, float)), *parseff)
    eff = unumpy.exp(ln_eff)
    return unumpy.nominal_values(eff), unumpy.std_devs(eff)



YLABEL = "Counts per yield"     # either absolute or relative efficiency change depending on whichever you're running

#plot stuff
fig, ax = plt.subplots(figsize=(20,10))
plt.rcParams.update({'font.size': 20})


ax.errorbar(E, y, yerr=y_err, fmt='o', c='r', capsize=4, label=r'Experimental data')
ax.plot(xfit, eff_fit, 'g', label=r" Fitted curve; reduced $\chi^2$ = {:.2f}".format(red_chi2))

ax.fill_between(xfit, eff_nom - eff_std, eff_nom + eff_std, alpha=0.3, color='g', label="±1σ")
ax.set_title('Efficiency fitted curve')
ax.set_xlim(E.min(), E.max())
ax.set_xlabel('Energy (keV)', fontsize=20)
ax.set_ylabel(rf'{YLABEL}', fontsize=20)
ax.grid()



ax.legend()
ax.text(0.5, 0.9, parameterseff, transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black'), va='top')
plt.tight_layout()
#plt.savefig('Efficiency.png')
plt.show()

