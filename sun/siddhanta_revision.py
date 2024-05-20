import numpy as np
from numpy import sin, cos, pi, deg2rad
import shelve
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

print('Starting Revision of parameters\n')

#####---------------- Sun's treatment according to Siddhanta ----------------#####
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from useful_functions import *
import siddhantaic_functions as sid
from functions import sun_theta_0
from functions import mandocca, r0byR

with shelve.open('./sun_observed_data/config_data') as config_data:
    steps       = config_data['steps']
    num_years   = config_data['num_years']
config_data.close()
with shelve.open('../modern_astronomical_constants/modern_astronomical_constants') as modern_astro_data:
    ANOM_YR     = modern_astro_data['ANOM_YR']
    SIDE_YR     = modern_astro_data['SIDE_YR']
modern_astro_data.close()
with shelve.open('./sun_observed_data/anomalistic_data') as anom_data:
    swelongs_an = anom_data['swelongs_an']
anom_data.close()

sun_revs = 4320000
civ_days_new = 1577917500

print('Imported modules and data')
#####----- Compute new parameters by comparing in one anomalistic cycle -----#####
show_all_plots = False
show_important_plots = True
#####-------------------- Effect of changing year length --------------------#####

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new, 0), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='civ_days = %d (KAPA)'%civ_days_new)

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, SIDE_YR*sun_revs, 0), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='civ_days = %d (SIDE)'%(SIDE_YR*sun_revs))

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, ANOM_YR*sun_revs), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='civ_days = %d (ANOM)'%(ANOM_YR*sun_revs))

plt.title('Effect of parameter: year length')
plt.ylabel('Difference between observed and siddhanta anomalistic longitudes (degrees)')
plt.xlabel('Step interval: ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/year_length.png')
if(show_all_plots): plt.show()
plt.clf()
print('Finished year length computations')

#####---------------------------- Modifying r0byR ---------------------------#####

##----- Make r0byR a  different constant -----##
civ_days_new=ANOM_YR*sun_revs
# Get ballpark range of new r0byR
r0byR=13.5/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))
r0byR=12/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, 12/360) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(12))

plt.title('Get ballpark range of new $r_0/R$')
plt.ylabel('Difference in observed and Siddhanta longitudes (degrees)')
plt.xlabel('Time since apogee: Step interval = ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0byR_original.png')
if(show_all_plots): plt.show()
plt.clf()
print('Finished ballpark $r_0/R$ computations')
##-------- Make r0byR vary with time --------##

# Demonstrate effect of varying r0byR
r0byR=12.1/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(sidlongs, swelongs_an), '-', label='$r_0/R$ = %f/360'%(r0byR*360))
r0byR=12/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, 12/360) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(sidlongs, swelongs_an), '-', label='$r_0/R$ = %f/360'%(r0byR*360))
r0byR=11.9/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc(sidlongs, swelongs_an), '-', label='$r_0/R$ = %f/360'%(r0byR*360))

plt.title('Effect of varying $r_0/R$')
plt.ylabel('$\\theta_{{ms}}-\\theta_{{jpl}}$, degrees')
plt.xlabel('Time since apogee: Step interval = ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0byR_variation.png')
if(show_important_plots): plt.show()
plt.clf()

# Get "actual" variation of r0byR(t) with `epsilon` error in anomalistic longitude

def get_r0byR(ag, long_an):
    if sin(deg2rad(sun_theta_0(ag, sun_revs, civ_days_new))) == 0:
        ag=ANOM_YR/4
        long_an=swelongs_an[int(steps/4)]
    _min=11.7/360
    _max = 12.3/360
    computed_r0 = sin(-deg2rad(long_an - sun_theta_0(ag, sun_revs, civ_days_new)))/sin(deg2rad(sun_theta_0(ag, sun_revs, civ_days_new)))
    if _min<computed_r0:
        if computed_r0<_max: return computed_r0
        else: return _max
    else: return _min

# Get "actual" r0s over half an anomalistic year
r0s = [360*get_r0byR(i*ANOM_YR/steps, swelongs_an[i]) for i in range(int(steps))]
plt.plot(range(int(steps)), r0s, '.', label='Actual $r_0[n]$')

# Define model for variation of r0 over one anomalistic cycle
def r0(ag, m, k):
    return m-k*cos(deg2rad(sun_theta_0(ag, sun_revs, civ_days_new)))
# Fit model to "actual" r0s
p0=[12, 0.18]
# plt.plot(range(int(steps)), [r0(i*ANOM_YR/steps, *p0) for i in range(int(steps))], '-', label='Initial Guess: %s'%p0)
sun_r0byR_popt, r0byR_pcov = curve_fit(r0, [i*ANOM_YR/steps for i in range(int(steps))], np.array(r0s), p0=p0)
plt.plot(range(int(steps)), [r0(i*ANOM_YR/steps, *sun_r0byR_popt) for i in range(int(steps))], '-', label='Fitted Curve: [m, k] = %s'%sun_r0byR_popt)

plt.ylim(11.7,12.3)
plt.title('Fitting model parameters for epicyclic radius $r_0$ over one anomalistic cycle')
plt.ylabel('$r_0(t)$')
plt.xlabel('Time (t) since apogee: Step interval: 1 unit = ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0_model_halfcycle.png')
if(show_important_plots): plt.show()
plt.clf()
print('Finished r0byR parameter fitting computations')
# Extend model r0 to full anomalistic cycle
from functions import r0byR_revised

plt.plot(range(int(steps)), [360*r0byR_revised(sun_theta_0(i*ANOM_YR/steps, sun_revs, civ_days_new), 0, *sun_r0byR_popt) for i in range(int(steps))], '-', label='Fitted Curve: %s'%sun_r0byR_popt)

plt.title('Variation of epicyclic radius $r_0$ over one anomalistic year')
plt.ylabel('$r_0(t)$')
plt.xlabel('Time (t) since apogee: Step interval: 1 unit = ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0_model_fullcycle.png')
if(show_important_plots): plt.show()
plt.clf()
print('Plotted r0 over full cycle')
##------- Evaluating model precision --------##
r0byR=13.5/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc_secs(sidlongs, swelongs_an), '-', label='Original $r_0/R$ = %f/360'%(r0byR*360))

r0byR=sun_r0byR_popt[0]/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc_secs(sidlongs, swelongs_an), '-', label='New constant $r_0/R$ = %f/360'%(r0byR*360))

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, civ_days_new), 0, r0byR_revised(sun_theta_0(ag, sun_revs, civ_days_new), 0, *sun_r0byR_popt)) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), subcirc_secs(sidlongs, swelongs_an), '-', label='New oscillating $r_0/R$ with [m, k] = %s'%str(sun_r0byR_popt))

plt.ylabel(' $\\theta_{{ms}}-\\theta_{{jpl}}$ (seconds)')
plt.xlabel('Time (t) since apogee: Step interval: 1 unit = ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0byR_final_comparison.png')
if(show_all_plots): plt.show()
plt.ylim(-0.001*3600,0.001*3600)
plt.savefig('./graphs/r0byR_final_comparison_zoomed.png')
if(show_important_plots): plt.show()
plt.clf()
print('Plotted r0byR comparisons')

#####--------------------------- Modifying mandocca -------------------------#####
print('\nStarting Mandocca computations:')

civ_days_new = SIDE_YR*sun_revs
mandocca_revs = sun_revs-(SIDE_YR*sun_revs)/ANOM_YR

with shelve.open('./sun_observed_data/sidereal_data') as sid_data:
    swelongs_si          = sid_data['swelongs_si']
    mandocca_true_mesadi = sid_data['mandocca_true_mesadi']
    mandocca_mean        = sid_data['mandocca_mean']
sid_data.close()

def fit_sun_dhruvas(_inputs, _sun_dhruva, _sun_mandocca_dhruva):
    _ag, _obs_long = _inputs
    mean_sun = sid.theta_0(_ag, sun_revs, _sun_dhruva, civ_days_new)
    sun_mandocca = sid.theta_0(_ag, mandocca_revs, _sun_mandocca_dhruva, civ_days_new)
    r0 = sid.r0_revised(mean_sun, sun_mandocca, *sun_r0byR_popt)
    return (sid.theta_ms(mean_sun, sun_mandocca, r0/360)-_obs_long+180)%360-180

def get_sun(_ag, _sun_dhruva, _sun_mandocca_dhruva):
    mean_sun = sid.theta_0(_ag, sun_revs, _sun_dhruva, civ_days_new)
    sun_mandocca = sid.theta_0(_ag, mandocca_revs, _sun_mandocca_dhruva, civ_days_new)
    r0 = sid.r0_revised(mean_sun, sun_mandocca, *sun_r0byR_popt)
    return sid.theta_ms(mean_sun, sun_mandocca, r0/360)

# Fit dhruvas to observed data
p0=[-1.9, 79]
# plt.plot(range(num_years*steps), subcirc([get_sun(i*SIDE_YR/steps, *p0) for i in range(int(num_years*steps))], swelongs_si), '-', label='Initial Guess: %s'%p0)
sun_dhruva_popt, sun_dhruva_pcov = curve_fit(fit_sun_dhruvas, ([i*SIDE_YR/steps for i in range(int(num_years*steps))], swelongs_si), np.zeros(len(swelongs_si)), p0=p0)
plt.plot(range(int(num_years*steps)), subcirc_secs([get_sun(i*SIDE_YR/steps, *sun_dhruva_popt) for i in range(int(num_years*steps))], swelongs_si), '-', label='Fitted Curve: %s'%sun_dhruva_popt)

plt.title('Finding best-fit dhruva for the Sun and its apogee')
plt.ylabel('$\\theta_{{ms}}-\\theta_{{jpl}}$ (seconds)')
plt.xlabel('Time (t) since apogee: Step interval: 1 unit = SIDE_YR/%d = %f days'%(steps, SIDE_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/sun_dhruva.png')
if show_all_plots: plt.show()
plt.clf()

sun_dhruva = sun_dhruva_popt[0]
mandocca_dhruva = sun_dhruva_popt[1]

def fit_sun_mean_mandocca(_inputs, _sun_mean_mandocca):
    _ag, _obs_long = _inputs
    mean_sun = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
    r0 = sid.r0_revised(mean_sun, _sun_mean_mandocca, *sun_r0byR_popt)
    return (sid.theta_ms(mean_sun, _sun_mean_mandocca, r0/360)-_obs_long+180)%360-180

def get_sun_mean_mandocca(_ag, _sun_mean_mandocca):
    mean_sun = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
    r0 = sid.r0_revised(mean_sun, _sun_mean_mandocca, *sun_r0byR_popt)
    return sid.theta_ms(mean_sun, _sun_mean_mandocca, r0/360)

# Fit mean_mandocca to observed data
p0=[79]
# plt.plot(range(num_years*steps), subcirc([get_sun(i*SIDE_YR/steps, *p0) for i in range(int(num_years*steps))], swelongs_si), '-', label='Initial Guess: %s'%p0)
sun_mean_mandocca_popt, sun_mean_mandocca_pcov = curve_fit(fit_sun_mean_mandocca, ([i*SIDE_YR/steps for i in range(int(num_years*steps))], swelongs_si), np.zeros(len(swelongs_si)), p0=p0)
plt.plot(range(int(num_years*steps)), subcirc_secs([get_sun_mean_mandocca(i*SIDE_YR/steps, *sun_mean_mandocca_popt) for i in range(int(num_years*steps))], swelongs_si), '-', label='Fitted Curve: %s'%sun_mean_mandocca_popt)

plt.title('Finding best-fit mean value for the Sun\'s apogee')
plt.ylabel('$\\theta_{{ms}}-\\theta_{{jpl}}$ (seconds)')
plt.xlabel('Time (t) since apogee: Step interval: 1 unit = SIDE_YR/%d = %f days'%(steps, SIDE_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/sun_mean_mandocca.png')
if True: plt.show()
plt.clf()

sidlongs_true_mandocca = []
sidlongs_const_mandocca_avg = []
sidlongs_const_mandocca_dhruva = []
for year in range(num_years):
    for i in range(steps):
        ag = year*SIDE_YR+i*SIDE_YR/steps
        mean_sun = sun_theta_0(ag, sun_revs, civ_days_new, sun_dhruva)
        true_mandocca = sun_theta_0(ag, mandocca_revs, civ_days_new, mandocca_dhruva)
        _r0 = sid.r0_revised(mean_sun, true_mandocca, *sun_r0byR_popt)
        sidlongs_true_mandocca.append(sid.theta_ms(mean_sun, true_mandocca, _r0/360))

        mandocca = mandocca_mean
        _r0 = sid.r0_revised(mean_sun, mandocca, *sun_r0byR_popt)
        sidlongs_const_mandocca_avg.append(sid.theta_ms(mean_sun, mandocca, _r0/360))

        mandocca = mandocca_true_mesadi
        _r0 = sid.r0_revised(mean_sun, mandocca, *sun_r0byR_popt)
        sidlongs_const_mandocca_dhruva.append(sid.theta_ms(mean_sun, mandocca, _r0/360))
plt.plot(range(num_years*steps), subcirc_secs(sidlongs_true_mandocca, swelongs_si), '-', label='Mandocca linearly increasing from %f'%mandocca_dhruva)
plt.plot(range(num_years*steps), subcirc_secs(sidlongs_const_mandocca_avg, swelongs_si), '-', label='Mandocca constant at %f'%mandocca_mean)
plt.plot(range(num_years*steps), subcirc_secs(sidlongs_const_mandocca_dhruva, swelongs_si), '-', label='Mandocca constant at %f'%mandocca_true_mesadi)

plt.ylabel('$\\theta_{{ms}}-\\theta_{{jpl}}$ (seconds)')
plt.xlabel('Time since Aries Transit: Step interval: SIDE_YR/%d = %f days'%(steps, SIDE_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/Mandocca_comparison.png')
if(True): plt.show()
plt.clf()
print('Finished mandocca computations')

#####--------------------------- Print parameters -------------------------#####

print("""\n-----New Parameters:-----\n\
      nC = %f\n\n\
      mandocca_revs = %f\n\
      mandocca_dhruva = %s\n\
      mandocca_mean = %s\n\n\
      sun_dhruva = See epoch.py\n\
      r0byR_popt = %s\
      """%(civ_days_new, mandocca_revs, todms(mandocca_dhruva), todms(mandocca_mean), sun_r0byR_popt))

params_file = shelve.open('./sun_computed_data/sun_modified_parameters')
params_file['sun_civ_days'] = civ_days_new
params_file['r0byR_popt'] = sun_r0byR_popt
params_file['mandocca_dhruva'] = mandocca_dhruva
params_file['mandocca_mean'] = mandocca_mean
params_file['mandocca_revs'] = mandocca_revs
params_file['mandocca_dhruva'] = mandocca_dhruva
params_file['sun_dhruva'] = sun_dhruva
params_file.close()