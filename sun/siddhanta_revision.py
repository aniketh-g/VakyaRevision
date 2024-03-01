import numpy as np
from numpy import sqrt, sin, cos, arcsin, pi, rad2deg, deg2rad
from generate_data import steps, swelongs_an, swelongs_si, mandocca_true, ANOM_YR, SIDE_YR
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

# Madhyagraha
sun_revs = 4320000
civ_days = 1577917500

theta_0 = lambda ahargana: ((ahargana*sun_revs/civ_days)%1)*360

# Mandasamskara
mandocca = 78
r0byR = 13.5/360

theta_ms = lambda _theta_0, _mandocca, _r0byR: (_theta_0 - rad2deg(arcsin(_r0byR*sin(deg2rad(_theta_0-_mandocca)))))%360

#####----- Compute new parameters by comparing in one anomalistic cycle -----#####

#####-------------------- Effect of changing year length --------------------#####
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='civ_days = %d (KAPA)'%civ_days)

civ_days=SIDE_YR*sun_revs
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='civ_days = %d (SIDE)'%civ_days)

civ_days=ANOM_YR*sun_revs
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='civ_days = %d (ANOM)'%civ_days)

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

##----- Make r0byR a different constant -----##

# Get ballpark range of new r0byR
r0byR=13.5/360
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))
r0byR=12/360
sidlongs = [theta_ms(theta_0(ag), 0, 12/360) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
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
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))
r0byR=12/360
sidlongs = [theta_ms(theta_0(ag), 0, 12/360) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='r0byR = %f/360'%(r0byR*360))
r0byR=11.9/360
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
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
    delta = long_an - theta_ms(theta_0(ag), 0, r0R)
    while abs(delta) >= epsilon:
        timeout=timeout+1
        if(timeout>=10**5): print('get_r0byR TIMEOUT')
        delta = long_an - theta_ms(theta_0(ag), 0, r0R)
        if abs(delta) >= epsilon:
            if long_an <= 0:
                long_an = long_an + 360
            r0R = r0R - sin(deg2rad(delta))/sin(deg2rad(theta_0(ag)))
    return r0R

# Get "actual" r0s over half an anomalistic year
r0s = [360*get_r0byR(i*ANOM_YR/steps, swelongs_an[i]) for i in range(int(steps/2))]
plt.plot(range(int(steps/2)), r0s, '.', label='Actual r0(t)')

# Get minima for manual computation of params - IGNORE
ag0 = r0s.index(min(r0s))*ANOM_YR/steps
p1p2 = 1/sin(2*pi*(2*ag0)/(0.5*ANOM_YR))

# Define model for variation of r0 over anomalistic half-cycle
def r0(ag, m, k, p1, p2):
    return m-k*np.arctan(p1*sin(2*deg2rad(theta_0(ag)))/(1-p2*cos(2*deg2rad(theta_0(ag)))))

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
def r0byR_revised(t0, tm, *popt):
    [m, k, p1, p2] = popt
    delt = t0-tm
    if delt<0: delt=delt+360
    if delt<180: return (m-k*np.arctan(p1*sin(2*deg2rad(delt))/(1-p2*cos(2*deg2rad(delt)))))/360
    else: return (m+k*np.arctan(p1*sin(2*deg2rad(delt))/(1-p2*cos(2*deg2rad(delt)))))/360
plt.plot(range(int(steps)), [360*r0byR_revised(theta_0(i*ANOM_YR/steps), 0, *r0byR_popt) for i in range(int(steps))], '-', label='Fitted Curve: %s'%r0byR_popt)

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
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='Original r0byR = %f/360'%(r0byR*360))

r0byR=12.00969326/360
sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_an,sidlongs), '-', label='New constant r0byR = %f/360'%(r0byR*360))

sidlongs = [theta_ms(theta_0(ag), 0, r0byR_revised(theta_0(ag), 0, *r0byR_popt)) for ag in np.arange(0, ANOM_YR, ANOM_YR/steps)]
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
print('Starting Mandocca computations:')

civ_days = SIDE_YR*sun_revs
mandocca_revs = sun_revs-civ_days/ANOM_YR
theta_m = lambda ahargana, mean: mean+ ((ahargana*mandocca_revs/civ_days)%1)*360

theta_ms_new = lambda _theta_0, _theta_m, _r0byR_revised: (_theta_0 - rad2deg(arcsin(_r0byR_revised*sin(deg2rad(_theta_0-_theta_m)))))%360

def get_longitude_ag(target, ag):
    timeout=0
    epsilon = 0.000000001
    delta = theta_ms_new(theta_0(ag), theta_m(ag, mandocca_true), r0byR_revised(theta_0(ag),theta_m(ag, mandocca_true),*r0byR_popt)) - target
    while abs(delta) >= epsilon:
        timeout=timeout+1
        if(timeout>=10): print('get_longitude_ag: TIMEOUT')
        lon_na = theta_ms_new(theta_0(ag), theta_m(ag, mandocca_true), r0byR_revised(theta_0(ag),theta_m(ag, mandocca_true),*r0byR_popt))
        delta = lon_na - target
        if abs(delta) >= epsilon:
            if lon_na >= 180:
                lon_na = lon_na - 360
            ag = ag - delta
    return ag

ag_start = get_longitude_ag(0, 0)#1873765

print('Starting at ag = %f, theta_ms = %f'\
      %(ag_start, theta_ms_new(theta_0(ag_start), theta_m(ag_start, mandocca_true), r0byR_revised(theta_0(ag_start),theta_m(ag_start, mandocca_true),*r0byR_popt))))

m = mandocca_true
sidlongs = [theta_ms_new(theta_0(ag), theta_m(ag, m), r0byR_revised(theta_0(ag),theta_m(ag, m),*r0byR_popt))\
             for ag in np.arange(ag_start, ag_start+SIDE_YR, SIDE_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_si,sidlongs), '-', label='Mandocca linearly increasing from %f'%m)

m = mandocca_true
sidlongs = [theta_ms(theta_0(ag), m, r0byR_revised(theta_0(ag), m, *r0byR_popt))\
             for ag in np.arange(ag_start, ag_start+SIDE_YR, SIDE_YR/steps)]
plt.plot(range(steps), np.subtract(swelongs_si,sidlongs), '-', label='Mandocca constant at %f'%m)

plt.ylabel('Difference between observed and siddhanta sidereal longitudes (degrees)')
plt.xlabel('Time since Aries Transit: Step interval: SIDE_YR/%d = %f days'%(steps, SIDE_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/Mandocca_comparison.png')
# plt.show()
plt.clf()
print('Finished mandocca computations')

#####--------------------------- Modifying dhruva -------------------------#####
#####--------------------------- Print parameters -------------------------#####

print('\n-----New Parameters:-----\n\
      theta_m = %f\n\
      nC = %f\n\
      mandocca_revs = %f'%(mandocca_true, civ_days, mandocca_revs))
