import pandas as pd
import numpy as np

#lets find activity

#This is the time interval that these data have been run for. If there are different steplen across then we need to define more than one time interval
time_interval = 50#seconds
cd_time_interval =20000
Au_peak_efficiency = 0.00592809347158196 #taken from spreadsheet
Ni_peak_efficiency = 0.003404046236

def calculate_activity(input_path, output_path, interval, efficiency):
    counts = np.load(input_path)
    activity = counts / (interval * efficiency)
    np.save(output_path, activity)

# Au (plain)
calculate_activity(
    input_path=r'/Volumes/SandySSD/Physics/Xenon/Activities, Counts and Sigmas /Au_measured_activity.npy',
    output_path=f'/Volumes/SandySSD/Physics/Xenon/Plain_Gold_Activity_{time_interval}s',
    interval=time_interval,
    efficiency=Au_peak_efficiency
)

# Au (cadmium covered)
calculate_activity(
    input_path=r'/Volumes/SandySSD/Physics/Xenon/Activities, Counts and Sigmas /Ni_measured_activity.npy',
    output_path=f'/Volumes/SandySSD/Physics/Xenon/Cadmium_Gold_Activity_{time_interval}s',
    interval=cd_time_interval,
    efficiency=Au_peak_efficiency
)

# Ni
calculate_activity(
    input_path=r'/Volumes/SandySSD/Physics/Xenon/Activities, Counts and Sigmas /Ni_measured_activity.npy',
    output_path=f'/Volumes/SandySSD/Physics/Xenon/Nickel_Activity_{time_interval}s',
    interval=time_interval,
    efficiency=Ni_peak_efficiency
)





irrad_time = 15
time_btw_BOC_and_EOB = 345600
thermal_resonance_crosssection = 98.65E-24#Au198 cm^2
Resonance_region_crosssection = 1550E-24#Au198 cm^2
Au_f_crossection = 
og_sample_abundance =  3.8548810174E17#Au197Cd n
decay_constant= 2.977E-6 #Au198 n/s
#measured_activity = np.load(r'/Volumes/SandySSD/Physics/Xenon/Au_measured_activity.npy')


#Calculate plain gold activity based on counts and efficiency do we need to do error prop or what do we do with the sigma?
Au_peak_efficiency = 
measured_activity = measured_counts/Au_peak_efficiency
#print(measured_activity.shape)

#Cadmium-covered Gold knowns
cd_decay=2.977E-6
cd_irrad_time=15 #s
cd_time_btw_BOC_and_EOB= #s
#cd_measured_activity =  np.load(r'/Volumes/SandySSD/Physics/Xenon/Ni_measured_activity.npy')
cd_t_crosssection = 0
cd_r_crosssection = 1550E-24
cd_f_crossection =
cd_og_abundance = 3.81803676E17

#Calculate Cadmium-covered gold activity based on counts and efficiency
cd_measured_activity = cd_measured_counts / Au_peak_efficiency 

#Ni knowns
Ni_decay=1.1322E-7
Ni_irrad_time=15 #s
Ni_time_btw_BOC_and_EOB= 432000#s
#Ni_measured_activity =  np.load(r'/Volumes/SandySSD/Physics/Xenon/Ni_measured_activity.npy')
Ni_t_crosssection= # Nickel has essentially no t or r bc (n,p) is a threshold reaction???
Ni_r_crosssection = 4.620E-24 #1550E-24
Ni_f_crossection =
Ni_og_abundance= 1.55359946E21 #3.81803676E17

#calculate Ni activity based on counts and efficiency
Ni_peak_efficiency = 
Ni_measured_activity = Ni_measured_counts / Ni_peak_efficiency


#compacting variables for flux equation
dt=np.exp(-decay_constant*time_btw_BOC_and_EOB)
ti=1-np.exp(-decay_constant*irrad_time)


cd_dt=np.exp(-cd_decay*cd_time_btw_BOC_and_EOB)
cd_ti=np.exp(-cd_decay*cd_irrad_time)

Ni_dt = np.exp(-Ni_decay*Ni_time_btw_BOC_and_EOB)
Ni_ti = np.exp(-Ni_decay*cd_irrad_time)


#EpithermalCdFlux= activity/(dt*ti*Resonance_region_crosssection)

#for i in cd_r_crosssection and for j in cd_measured_activity:
   # A = np.array([cd_t_crosssection,cd_r_crosssection],[thermal_resonance_crosssection,Resonance_region_crosssection])
    
   #B = np.array([cd_measured_activity/(cd_og_abundance*cd_ti*cd_dt)],[measured_activity/(og_sample_abundance*ti*dt)])
    #X = np.dot(np.linalg.inv(A),B)

#We want to solve using X=B/A

#Get Number of Matrices/Solutions
N = len(measured_activity)

#the right hand side of each equation in 3 
B1 = cd_measured_activity / (cd_og_abundance*cd_ti*cd_dt) 
B2 = measured_activity / (og_sample_abundance*ti*dt)
B3 = Ni_measured_activity/(Ni_og_abundance*Ni_ti*Ni_dt)

#Making B stack to account for the fact that we have three solutions 
B_stack = np.array([B1,B2,B3])

#A array
A = np.array([[cd_t_crosssection,cd_r_crosssection,cd_f_crossection],[thermal_resonance_crosssection, Resonance_region_crosssection,Au_f_crossection],[Ni_t_crosssection,Ni_r_crosssection,Ni_f_crossection]])

#invert A then dot B to get B/A - using x0 to be the beginning flux and xf to be the final flux? or does staysl require activity?
X = np.dot(np.linalg.inv(A),B_stack)

X0=np.dot(np.linalg.inv(A[0:0]),B_stack[0:0])
Xf = np.dot(np.linalg.inv(A),B_stack)

#Now separate out thermal and epithermal flux from X matrix
thermal_flux = X[0,:]
epithermal_flux = X[1,:]
fast_flux = X[2,:]

print(thermal_flux + epithermal_flux + fast_flux)









