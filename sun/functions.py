from numpy import sqrt, sin, cos, arcsin, pi, rad2deg, deg2rad, arctan
import shelve
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')

#####--------------- Siddhanta Model ---------------#####

# Madhyagraha

sun_theta_0 = lambda ahargana, _revs, _sun_civ_days, _dhruva=0: (_dhruva + ((ahargana*_revs/_sun_civ_days)%1)*360)%360

# Mandasamskara
mandocca = 78
r0byR = 13.5/360

from siddhantaic_functions import theta_ms

#####--------------- Modified Model ---------------#####

def r0byR_revised(t0, tm, *popt):
    [m, k] = popt
    delt = t0-tm
    return (m-k*cos(deg2rad(delt)))/360