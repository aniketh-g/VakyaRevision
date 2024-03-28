# This file generates a list of `steps` elements, containing the
# longitudes of the Sun with respect to its apogee (mandocca), 
# over one anomalistic cycle, at  intervals of `ANOM_YR`/`steps`,
# averaged over `numyears` years

from numpy import abs
import sys
from scipy.stats import circmean
import swisseph as swe

swe.set_ephe_path('~/swisseph/ephe') # Set Ephemeris Path
swe.set_topo(75.7683, 0.0, 0.0) # Positions as observed at Lanka

#####------------------------ Astronomical Parameters -----------------------#####
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from modern_astro_params import ANOM_YR, SIDE_YR

#####-------- Set Initial year, averaging window, and number of steps -------#####
startyear = 2029 # AD
num_years = 60    # Averaging window
steps = 1000     # Number of readings per year
target_year = 2029 # Sidereal single year

# Configure Parameters
COMPUTE_PRESENT_AN = 'T'#sys.argv[1]
COMPUTE_PRESENT_SI = 'T'#sys.argv[2]
COMPUTE_HISTORICAL_AN = 'F'#sys.argv[3]

#####--------- Get anomalistic longitudes over one anomalistic cycle --------#####
def get_mandocca_date(jd):
    epsilon = 0.000000001
    delta = epsilon + 1
    while abs(delta) >= epsilon:
        lon_sa, _ = swe.calc(jd, swe.SUN)
        mandocca_sa = swe.nod_aps(jd, swe.SUN)[3][0]
        delta = lon_sa[0] - (mandocca_sa-0.0069586069359039325) #Set such that average is close to 0
        if abs(delta) >= epsilon:
            if lon_sa[0] >= 180:
                lon_sa[0] = lon_sa[0] - 360
            if lon_sa[0] <= 0:
                lon_sa[0] = lon_sa[0] + 360
            jd = jd - delta
    return jd

def get_mean_anom_longs(startyear):
    jd_ut = swe.julday(startyear, 5, 15, 0)
    jd = get_mandocca_date(jd_ut + swe.deltat(jd_ut))

    datas = [[swe.calc(jd+year*ANOM_YR+i*ANOM_YR/steps, swe.SUN, swe.FLG_SIDEREAL)[0][0] \
            - (swe.nod_aps(jd+year*ANOM_YR+i*ANOM_YR/steps, swe.SUN)[3][0] \
            - swe.get_ayanamsa_ex(tjdet=jd+year*ANOM_YR+i*ANOM_YR/steps, flags=0)[1]) \
                for i in range(steps)] for year in range(num_years)]

    return [circmean([data[i] for data in datas], high=360, low=0) for i in range(steps)]

if COMPUTE_PRESENT_AN == 'T':
    print('Computing average anomalistic longitudes from AD %d to %d'%(startyear, startyear+num_years))
    swelongs_an = get_mean_anom_longs(startyear)
    print('Computed anomalistic longitudes! Longitude on first day = %f\n'%(swelongs_an[0]))

#####----------------------- Get variation of mandocca ----------------------#####

#####----------- Get sidereal longitudes over one sidereal cycle ------------#####
swe.set_sid_mode(mode = swe.SIDM_LAHIRI) # Ayanamsa

def get_longitude_date(target, jd):
    timeout = 0
    epsilon = 1e-9
    delta = swe.calc(jd, swe.SUN, swe.FLG_SIDEREAL)[0][0] - target
    while abs(delta) >= epsilon:
        timeout=timeout+1
        if(timeout>=10**4): print('get_longitude_date: TIMEOUT')
        (lon_na, lat, dist, vlon, vlat, vdist), _ = swe.calc(jd, swe.SUN, swe.FLG_SIDEREAL)
        delta = lon_na - target
        if abs(delta) >= epsilon:
            if lon_na >= 180:
                lon_na = lon_na - 360
            jd = jd - delta
    return jd

if COMPUTE_PRESENT_SI == 'T':
    print('Computing sidereal longitudes on AD %d'%(startyear))

    jd_ut = swe.julday(startyear, 5, 15, 0)
    jd = get_longitude_date(0, jd_ut + swe.deltat(jd_ut)) # 359.997719 Set so that average is 0

    mandocca_true_mesadi = swe.nod_aps(jd, swe.SUN)[3][0] - swe.get_ayanamsa_ex(tjdet=jd, flags=0)[1]
    print('Mandocca on jd %f = %f degrees\n'%(jd, mandocca_true_mesadi))
    swelongs_si = []# [swe.calc(jd+i*SIDE_YR/steps, swe.SUN, swe.FLG_SIDEREAL)[0][0] for i in range(steps)]
    for year in range(num_years):
        for i in range(steps):
            swelongs_si.append(swe.calc(jd+year*SIDE_YR+i*SIDE_YR/steps, swe.SUN, swe.FLG_SIDEREAL)[0][0])

    mandocca_mean = circmean([swe.nod_aps(jd+i*SIDE_YR/steps, swe.SUN)[3][0] - swe.get_ayanamsa_ex(tjdet=jd+i*SIDE_YR/steps, flags=0)[1] for i in range(steps*num_years)], high=360, low=0)
    print('Mandocca averaged =%f\n'%mandocca_mean)
    # datas = [[swe.calc(jd+year*SIDE_YR+i*ANOM_YR/steps, swe.SUN, swe.FLG_SIDEREAL)[0][0] for i in range(steps)] for year in range(numyears)]
    # swelongs_si = [circmean([data[i] for data in datas], high=360, low=0) for i in range(steps)]
    print('Computed swelongs_si starts on jd = %f (%s) at sidereal longitude = %f'%(jd, swe.revjul(jd), swelongs_si[0]))

#####--- Get historical anomalistic longitudes over one anomalistic cycle ---#####
if COMPUTE_HISTORICAL_AN == 'T':
    steps = 360
    swelongs_an_hi=[]
    for year in range(-3000, 2001, 500):
        print('Computing average anomalistic longitudes from AD %d to %d, step interval = ANOM_YR/%d'%(year, year+num_years, steps))
        jd = swe.julday(year, 5, 15, 0) + swe.deltat(swe.julday(year, 5, 15, 0))
        mandocca = (swe.nod_aps(jd, swe.SUN)[3][0] - swe.get_ayanamsa_ex(tjdet=jd, flags=0)[1])%360
        swelongs_an_hi.append([year, get_mean_anom_longs(year), mandocca])