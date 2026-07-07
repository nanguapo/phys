import pandas as pd
import numpy as np
import os
import math


#This is the time interval that these data have been run for. If there are different steplen across then we need to define more than one time interval
time_interval = 50#seconds
cd_time_interval =2000#seconds
 #taken from spreadsheet


#Nickel Variables
og_Ni_abundance = 1.55E+21
Ni_decay_constant = 0.0000001132
ni_btw_EOI_Boc = 395382
Ni_target_isotope_abundance = 0.68077
Ni_peak_efficiency = 0.003404046236
Ni_branching_ratio = 0.9945
Ni_crosssectional_Area = 3.5E-26
Ni_abundance = og_Ni_abundance #Ni_target_isotope_abundance #- testing to see if sigphi does this calculation for us

# Gold-Cadmium-Covered variables
Au_branching_ratio = 0.9562
Au_Epithermal_crosssectional_area = 1.55E-21
Au_fast_crosssectional_area = 9.8E-26
Au_percent_abundance = 1
Au_MM = 196.967
Au_peak_efficiency = 0.00592809347158196
Au_decay_constant = 0.000002976
Cd_foil_mass = 0.02048
Cd_foil_purity = .0061
CdAu_abundance = (6.0221408e+23*Cd_foil_mass*Cd_foil_purity*Au_percent_abundance/Au_MM)
CdAu_decay_time = 1341457
Cd_livetime = 264047.81

#plain gold wire variables 

au_foil_mass = 0.02067
Au_decay_time = 386483
Au_livetime = 6076
Au_thermal_crosssectional_area = 9.865E-23
Au_abundance = (6.0221408e+23* au_foil_mass *Cd_foil_purity*Au_percent_abundance/Au_MM)


def calculate_uncertainty(input_file):
    uncertainty =np.load(input_file)
    activity_uncertainty = np.sqrt(np.sum(uncertainty**2))

    return float(activity_uncertainty)



def calculate_activity_and_ReactionRate(input_path, count_time, efficiency, decay_constant, branching_ratio, time_irradiated, abundance, decay_time):
    counts = np.load(input_path)
    
    Net_counts = math.fsum(counts)
    activity = Net_counts * decay_constant / ((1 - np.exp(-decay_constant * count_time)) * efficiency * branching_ratio)
    #np.save(output_path, activity)
    print("net counts", Net_counts)
    print("activity", activity)
    Reaction_Rate = (activity*time_irradiated) / (abundance *(np.exp(-decay_constant*decay_time))*(1-np.exp(-decay_constant*time_irradiated)))
    print("reaction rate is", Reaction_Rate)
    return float(Reaction_Rate) 




Ni_RR = calculate_activity_and_ReactionRate(
    input_path=r'/Volumes/SandySSD/Physics/Xenon/Activities_Counts_Sigmas/Rerun_of_Ni_measured_activity_50s.npy',
    count_time= 42884.54,
    efficiency=Ni_peak_efficiency,
    decay_constant= Ni_decay_constant,
    branching_ratio= Ni_branching_ratio,
    time_irradiated = 10,
    abundance = Ni_abundance,
    decay_time = ni_btw_EOI_Boc
    
)
ni_net_counts_uncertainty = calculate_uncertainty(r"/Volumes/SandySSD/Physics/Xenon/Activities_Counts_Sigmas/Ni_sigma_activity_rerun.npy")

print("NI counts Uncertainty", ni_net_counts_uncertainty)
Ni_flux = float(Ni_RR/(10*Ni_crosssectional_Area))
print("Fast flux",Ni_flux)

Cd_Au_RR = calculate_activity_and_ReactionRate(
    input_path=r'/Volumes/SandySSD/Physics/Xenon/Activities_Counts_Sigmas/CdAu_measured_activity_rerun.npy',
    count_time = Cd_livetime,
    efficiency = Au_peak_efficiency,
    decay_constant = Au_decay_constant,
    branching_ratio = Au_branching_ratio,
    time_irradiated = 10,
    abundance = CdAu_abundance,
    decay_time = CdAu_decay_time

)

CdAu_net_counts_uncertainty = calculate_uncertainty(r"/Volumes/SandySSD/Physics/Xenon/Activities_Counts_Sigmas/CdAu_sigma_activity_rerun.npy")
print("Cd-Au counts Uncertainty", CdAu_net_counts_uncertainty) 

CdAu_Flux = float(Cd_Au_RR/(10*Au_Epithermal_crosssectional_area))
print("Epithermal flux", CdAu_Flux)

Au_RR = calculate_activity_and_ReactionRate(
    input_path=r'/Volumes/SandySSD/Physics/Xenon/Activities_Counts_Sigmas/Rerun_of_Au_measured_activity_50s.npy',
    count_time = Au_livetime,
    efficiency = Au_peak_efficiency,
    decay_constant = Au_decay_constant,
    branching_ratio = Au_branching_ratio,
    time_irradiated = 10,
    abundance = Au_abundance,
    decay_time = Au_decay_time

)

Au_net_counts_uncertainty = calculate_uncertainty(r"/Volumes/SandySSD/Physics/Xenon/Activities_Counts_Sigmas/Au_sigma_activity_rerun.npy")
print("Au counts Uncertainty", Au_net_counts_uncertainty) 

thermal_flux = float((Au_RR - Cd_Au_RR)/(10*Au_thermal_crosssectional_area))
print("thermal flux",thermal_flux)






