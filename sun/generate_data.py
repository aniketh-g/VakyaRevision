# This file generates a list of `steps` elements, containing the
# longitudes of the Sun with respect to its apogee (mandocca), 
# over one anomalistic cycle, at  intervals of `ANOM_YR`/`steps`,
# averaged over `numyears` years

from numpy import abs
import shelve
import sys
from scipy.stats import circmean
import swisseph as swe

swe.set_ephe_path('~/swisseph/ephe') # Set Ephemeris Path
swe.set_topo(75.7683, 0.0, 0.0) # Positions as observed at Lanka

#####------------------------ Astronomical Parameters -----------------------#####
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
with shelve.open('../modern_astronomical_constants/modern_astronomical_constants') as modern_astro_data:
    ANOM_YR     = modern_astro_data['ANOM_YR']
    SIDE_YR     = modern_astro_data['SIDE_YR']
modern_astro_data.close()

#####-------- Set Initial year, averaging window, and number of steps -------#####
startyear = 2024 # AD
num_years = 100    # Averaging window
steps = 1000     # Number of readings per year
target_year = 2029 # Sidereal single year

#####--------------------------- Print config parameters -------------------------#####

print('\n-----Configuration Data:-----\n\
      steps = %f\n\
      num_years = %f'%(steps, num_years))

data_file = shelve.open('./sun_observed_data/config_data')
data_file['steps'] = steps
data_file['num_years'] = num_years
data_file.close()

# Configure Parameters
COMPUTE_PRESENT_AN = 'T'#sys.argv[1]
COMPUTE_PRESENT_SI = 'T'#sys.argv[2]
COMPUTE_HISTORICAL_AN = 'F'#sys.argv[3]
COMPUTE_DHRUVA = 'T'#sys.argv[3]

#####--------- Get anomalistic longitudes over one anomalistic cycle --------#####
def get_mandocca_date(jd):
    epsilon = 0.000000001
    delta = epsilon + 1
    while abs(delta) >= epsilon:
        lon_sa, _ = swe.calc(jd, swe.SUN)
        mandocca_sa = swe.nod_aps(jd, swe.SUN)[3][0]
        delta = lon_sa[0] - (mandocca_sa-0) #Set such that average is close to 0 # 0.0069586069359039325
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
    print('\n-----Computing Anomalistic Data:-----\n')
    print('Computing average anomalistic longitudes from AD %d to %d'%(startyear, startyear+num_years))
    swelongs_an = get_mean_anom_longs(startyear)
    print('Computed anomalistic longitudes! Longitude on first day %s = %f'%(swe.revjul(get_mandocca_date(swe.julday(startyear, 5, 15, 0))), swelongs_an[0]))

    #####---------------------- Print anomalistic parameters --------------------#####
    print('swelongs_an[0] = %f\n\
\n-----Generated Anomalistic Data:-----\n'%(swelongs_an[0]))

    data_file = shelve.open('./sun_observed_data/anomalistic_data')
    data_file['swelongs_an'] = swelongs_an
    data_file.close()

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
    print('\n-----Computing sidereal longitudes on AD %d-----\n'%(startyear))

    jd_ut = swe.julday(startyear, 4, 10, 0)
    jdtt_mesadi = swe.solcross(0, jd_ut + swe.deltat(jd_ut), swe.FLG_SIDEREAL) # 359.997719 Set so that average is 0

    mandocca_true_mesadi = swe.nod_aps(jdtt_mesadi, swe.SUN)[3][0] - swe.get_ayanamsa_ex(tjdet=jdtt_mesadi, flags=0)[1]
    swelongs_si = []
    for year in range(num_years):
        for i in range(steps):
            swelongs_si.append(swe.calc(jdtt_mesadi+year*SIDE_YR+i*SIDE_YR/steps, swe.SUN, swe.FLG_SIDEREAL)[0][0])

    mandocca_mean = circmean([swe.nod_aps(jdtt_mesadi+i*SIDE_YR/steps, swe.SUN)[3][0] - swe.get_ayanamsa_ex(tjdet=jdtt_mesadi+i*SIDE_YR/steps, flags=0)[1] for i in range(steps*num_years)], high=360, low=0)
    print('Computed swelongs_si starts on jd = %f %s at sidereal longitude = %f'%(jdtt_mesadi, swe.revjul(jdtt_mesadi), swelongs_si[0]))
    print('Mandocca on jd %f = %f degrees'%(jdtt_mesadi, mandocca_true_mesadi))
    print('Mandocca averaged =%f'%mandocca_mean)

    #####----------------------- Print sidereal parameters ----------------------#####
    print('\
mandocca_mean = %f\n\
mandocca_true_mesadi = %f\n\
swelongs_si[0] = %f\
\n\n-----Generated Sidereal Data-----\n\n'%(mandocca_mean, mandocca_true_mesadi, swelongs_si[0]))

    data_file = shelve.open('./sun_observed_data/sidereal_data')
    data_file['mandocca_mean'] = mandocca_mean
    data_file['mandocca_true_mesadi'] = mandocca_true_mesadi
    data_file['swelongs_si'] = swelongs_si
    data_file['swelongs_si_start_jdtt'] = jdtt_mesadi
    data_file.close()

#####--- Get historical anomalistic longitudes over one anomalistic cycle ---#####
if COMPUTE_HISTORICAL_AN == 'T':
    steps = 360
    swelongs_an_hi=[]
    for year in range(-3000, 2001, 500):
        print('Computing average anomalistic longitudes from AD %d to %d, step interval = ANOM_YR/%d'%(year, year+num_years, steps))
        jd = swe.julday(year, 5, 15, 0) + swe.deltat(swe.julday(year, 5, 15, 0))
        mandocca = (swe.nod_aps(jd, swe.SUN)[3][0] - swe.get_ayanamsa_ex(tjdet=jd, flags=0)[1])%360
        swelongs_an_hi.append([year, get_mean_anom_longs(year), mandocca])

    #####---------------------- Print historical parameters ---------------------#####
    print('\n-----Generated Historical Data:-----\n\
            swelongs_an_hi[0] = %f'%(swelongs_an_hi))

    data_file = shelve.open('./sun_observed_data/historical_data')
    data_file['swelongs_an_hi'] = swelongs_an_hi
    data_file.close()