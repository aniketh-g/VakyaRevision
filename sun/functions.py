from numpy import sqrt, sin, cos, arcsin, pi, rad2deg, deg2rad, arctan
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from modern_astro_params import SIDE_YR, ANOM_YR

#####--------------- Siddhanta Model ---------------#####

# Madhyagraha

sun_theta_0 = lambda ahargana, _revs, _sun_civ_days, _dhruva=0: (_dhruva + ((ahargana*_revs/_sun_civ_days)%1)*360)%360

# Mandasamskara
mandocca = 78
r0byR = 13.5/360

from siddhantaic_functions import theta_ms

#####--------------- Modified Model ---------------#####

# Modified r0_by_R
def r0byR_revised(t0, tm, *popt):
    [m, k, p1, p2] = popt
    delt = t0-tm
    if delt<0: delt=delt+360
    if delt<180: return (m-k*arctan(p1*sin(2*deg2rad(delt))/(1-p2*cos(2*deg2rad(delt)))))/360
    else: return (m+k*arctan(p1*sin(2*deg2rad(delt))/(1-p2*cos(2*deg2rad(delt)))))/360

# Modified Mandocca
    # See Main code