import numpy as np
from numpy import sin, cos, pi, deg2rad
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def subcirc(l1, l2):
    lres = []
    for i in range(len(l1)):
        _l = l1-l2
        if _l > 180:
            _l=_l-360
        if _l < -180:
            _l=_l+360
        lres.append(_l)
    return lres

print('Starting Revision of parameters\n')

#####---------------- Sun's treatment according to Siddhanta ----------------#####
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from modern_astro_params import ANOM_YR, SIDE_YR
from generate_data import steps, swelongs_an, swelongs_si, num_years
import siddhantaic_functions as sid
from functions import sun_theta_0
from functions import mandocca, r0byR

sun_revs = 4320000
sun_civ_days = 1577917500

print('Imported modules and data')
#####----- Compute new parameters by comparing in one anomalistic cycle -----#####

#####-------------------- Effect of changing year length --------------------#####

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days, 0), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='civ_days = %d (KAPA)'%sun_civ_days)

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, SIDE_YR*sun_revs, 0), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='civ_days = %d (SIDE)'%(SIDE_YR*sun_revs))

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, ANOM_YR*sun_revs), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='civ_days = %d (ANOM)'%(ANOM_YR*sun_revs))

plt.title('Effect of parameter: year length')
plt.ylabel('Difference between observed and siddhanta anomalistic longitudes (degrees)')
plt.xlabel('Step interval: ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/yearlength.png')
# plt.show()
plt.clf()
print('Finished year length computations')

#####---------------------------- Modifying r0byR ---------------------------#####

##----- Make r0byR a  different constant -----##
sun_civ_days=ANOM_YR*sun_revs
# Get ballpark range of new r0byR
r0byR=13.5/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))
r0byR=12/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, 12/360) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(12))

plt.title('Get ballpark range of new r0byR')
plt.ylabel('Difference in observed and Siddhanta longitudes (degrees)')
plt.xlabel('Time since apogee: Step interval = ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0byR_original.png')
# plt.show()
plt.clf()
print('Finished ballpark r0byR computations')
##-------- Make r0byR vary with time --------##

# Demonstrate effect of varying r0byR
r0byR=12.1/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))
r0byR=12/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, 12/360) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))
r0byR=11.9/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))

plt.title('Effect of varying r0byR')
plt.ylabel('Difference in observed and Siddhanta longitudes (degrees)')
plt.xlabel('Time since apogee: Step interval = ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0byR_variation.png')
# plt.show()
plt.clf()

# Get "actual" variation of r0byR(t) with `epsilon` error in anomalistic longitude

# def get_r0byR(ag, long_an):
#     if sin(deg2rad(theta_0(ag))) == 0: ag=(ANOM_YR/steps)/2
#     return sin(-deg2rad(long_an - theta_0(ag)))/sin(deg2rad(theta_0(ag)))

def get_r0byR(ag, long_an):
    timeout=0
    epsilon = 0.000001
    r0R=r0byR
    delta = long_an - sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0R)
    while abs(delta) >= epsilon:
        timeout=timeout+1
        if(timeout>=10**5): print('get_r0byR TIMEOUT')
        delta = long_an - sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0R)
        if abs(delta) >= epsilon:
            if long_an <= 0:
                long_an = long_an + 360
            r0R = r0R - sin(deg2rad(delta))/sin(deg2rad(sun_theta_0(ag, sun_revs, sun_civ_days)))
    return r0R

# Get "actual" r0s over half an anomalistic year
r0s = [360*get_r0byR(i*ANOM_YR/steps, swelongs_an[i]) for i in range(int(steps/2))]
plt.plot(range(int(steps/2)), r0s, '.', label='Actual r0(t)')

# Get minima for manual computation of params - IGNORE
ag0 = r0s.index(min(r0s))*ANOM_YR/steps
p1p2 = 1/sin(2*pi*(2*ag0)/(0.5*ANOM_YR))

# Define model for variation of r0 over anomalistic half-cycle
def r0(ag, m, k, p1, p2):
    return m-k*np.arctan(p1*sin(2*deg2rad(sun_theta_0(ag, sun_revs, sun_civ_days)))/(1-p2*cos(2*deg2rad(sun_theta_0(ag, sun_revs, sun_civ_days)))))

# Fit model to "actual" r0s
p0=[12, 0.18, 0.5, 0.5]
plt.plot(range(int(steps/2)), [r0(i*ANOM_YR/steps, *p0) for i in range(int(steps/2))], '-', label='Initial Guess: %s'%p0)
r0byR_popt, r0byR_pcov = curve_fit(r0, [i*ANOM_YR/steps for i in range(int(steps/2))], np.array(r0s), p0=p0)
plt.plot(range(int(steps/2)), [r0(i*ANOM_YR/steps, *r0byR_popt) for i in range(int(steps/2))], '-', label='Fitted Curve: %s'%r0byR_popt)

plt.title('Fitting model parameters for epicyclic radius r0 using data over half an anomalistic cycle')
plt.ylabel('r0(t)')
plt.xlabel('Time since apogee: Step interval: ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0_model_halfcycle.png')
# plt.show()
plt.clf()
print('Finished r0byR parameter fitting computations')
# Extend model r0 to full anomalistic cycle
from functions import r0byR_revised

plt.plot(range(int(steps)), [360*r0byR_revised(sun_theta_0(i*ANOM_YR/steps, sun_revs, sun_civ_days), 0, *r0byR_popt) for i in range(int(steps))], '-', label='Fitted Curve: %s'%r0byR_popt)

plt.title('Variation of epicyclic ratio r0 over one anomalistic year')
plt.ylabel('r0(t): epicyclic radius as a function of time')
plt.xlabel('Time since apogee: Step interval: ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0_model_fullcycle.png')
# plt.show()
plt.clf()
print('Plotted r0 over full cycle')
##------- Evaluating model precision --------##
r0byR=13.5/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='Original r0byR = %f/360'%(r0byR*360))

r0byR=r0byR_popt[0]/360
sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='New constant r0byR = %f/360'%(r0byR*360))

sidlongs = [sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), 0, r0byR_revised(sun_theta_0(ag, sun_revs, sun_civ_days), 0, *r0byR_popt)) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='Varying r0byR with parameters: %s'%str(r0byR_popt))

plt.ylabel('Difference between observed and siddhanta anomalistic longitudes (degrees)')
plt.xlabel('Time since apogee: Step interval: ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0byR_final_comparison.png')
# plt.show()
plt.ylim(-0.002,0.002)
plt.savefig('./graphs/r0byR_final_comparison_zoomed.png')
# plt.show()
plt.clf()
print('Plotted r0byR comparisons')
#####--------------------------- Modifying mandocca -------------------------#####
print('\nStarting Mandocca computations:')

sun_civ_days = SIDE_YR*sun_revs
mandocca_revs = sun_revs-(SIDE_YR*sun_revs)/ANOM_YR
from generate_data import mandocca_true_mesadi, mandocca_mean

def get_longitude_ag(target, ag):
    timeout=0
    epsilon = 1e-8
    delta = sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), sun_theta_0(ag, mandocca_revs, sun_civ_days, mandocca_true_mesadi), r0byR_revised(sun_theta_0(ag, sun_revs, sun_civ_days), sun_theta_0(ag, mandocca_revs, sun_civ_days, mandocca_true_mesadi),*r0byR_popt)) - target
    while abs(delta) >= epsilon:
        timeout=timeout+1
        if(timeout>=10): print('get_longitude_ag: TIMEOUT')
        lon_na = sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), sun_theta_0(ag, mandocca_revs, sun_civ_days, mandocca_true_mesadi), r0byR_revised(sun_theta_0(ag, sun_revs, sun_civ_days), sun_theta_0(ag, mandocca_revs, sun_civ_days, mandocca_true_mesadi),*r0byR_popt))
        delta = lon_na - target
        if abs(delta) >= epsilon:
            if lon_na >= 180:
                lon_na = lon_na - 360
            ag = ag - delta
    return ag

ag_start = get_longitude_ag(0, 0)#1873765

print('ag_start = %f, sid.theta_ms(ag_start) = %f'\
      %(ag_start, sid.theta_ms(sun_theta_0(ag_start, sun_revs, sun_civ_days), sun_theta_0(ag_start, mandocca_revs, sun_civ_days, mandocca_true_mesadi), r0byR_revised(sun_theta_0(ag_start, sun_revs, sun_civ_days), sun_theta_0(ag_start, mandocca_revs, sun_civ_days, mandocca_true_mesadi),*r0byR_popt))))

mandocca_dhruva = sun_theta_0(-ag_start, mandocca_revs, sun_civ_days, mandocca_true_mesadi)
sidlongs=[]
for year in range(num_years):
    for i in range(steps):
        ag = ag_start+year*SIDE_YR+i*SIDE_YR/steps
        sidlongs.append(sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), sun_theta_0(ag, mandocca_revs, sun_civ_days, mandocca_dhruva), r0byR_revised(sun_theta_0(ag, sun_revs, sun_civ_days), sun_theta_0(ag, mandocca_revs, sun_civ_days, mandocca_dhruva),*r0byR_popt)))
plt.plot(range(num_years*steps), np.subtract(swelongs_si,sidlongs), '-', label='Mandocca linearly increasing from %f'%mandocca_dhruva)

m = mandocca_mean
sidlongs=[]
for year in range(num_years):
    for i in range(steps):
        ag = ag_start+year*SIDE_YR+i*SIDE_YR/steps
        sidlongs.append(sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), m, r0byR_revised(sun_theta_0(ag, sun_revs, sun_civ_days), m, *r0byR_popt)))
plt.plot(range(num_years*steps), np.subtract(swelongs_si,sidlongs), '-', label='Mandocca constant at %f'%m)

m = mandocca_true_mesadi
sidlongs=[]
for year in range(num_years):
    for i in range(steps):
        ag = ag_start+year*SIDE_YR+i*SIDE_YR/steps
        sidlongs.append(sid.theta_ms(sun_theta_0(ag, sun_revs, sun_civ_days), m, r0byR_revised(sun_theta_0(ag, sun_revs, sun_civ_days), m, *r0byR_popt)))
plt.plot(range(num_years*steps), np.subtract(swelongs_si,sidlongs), '-', label='Mandocca constant at %f'%m)

plt.ylim(-0.015,0.015)
plt.ylabel('Difference between observed and siddhanta sidereal longitudes (degrees)')
plt.xlabel('Time since Aries Transit: Step interval: SIDE_YR/%d = %f days'%(steps, SIDE_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/Mandocca_comparison.png')
plt.show()
plt.clf()
print('Finished mandocca computations')

#####--------------------------- Modifying dhruva -------------------------#####
#####--------------------------- Print parameters -------------------------#####

import shelve

print('\n-----New Parameters:-----\n\
      Sun\'s mean mandocca = %f\n\
      nC = %f\n\
      mandocca_revs = %f'%(mandocca_mean, sun_civ_days, mandocca_revs))

params_file = shelve.open('./tables/sun_modified_parameters')
params_file['sun_civ_days'] = sun_civ_days
params_file['r0byR_popt'] = r0byR_popt
params_file['mandocca_dhruva'] = mandocca_dhruva
params_file['mandocca_mean'] = mandocca_mean
params_file['mandocca_revs'] = mandocca_revs
params_file.close()