from numpy import abs
import sys
from scipy.stats import circmean
import swisseph as swe

import matplotlib.pyplot as plt

swe.set_ephe_path('~/swisseph/ephe') # Set Ephemeris Path
swe.set_topo(75.7683, 0.0, 0.0) # Positions as observed at Lanka

#####------------------------ Astronomical Parameters -----------------------#####
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')

#####------------------------ Import modern constants ------------------------#####
print('#####------------------------ Import modern constants ------------------------#####')
import shelve
with shelve.open('../modern_astronomical_constants/modern_astronomical_constants') as astro_const:
    ANOM_MT = astro_const['ANOM_MT']
    SIDE_MT = astro_const['SIDE_MT']

import shelve
params_file = shelve.open('./moon_computed_data/generated_data_short')
#####-------- Set Initial year, averaging window, and number of steps -------#####
startyear = 2024 # AD
num_months = 24
num_months_avg = 200    # Averaging window
steps = 1000     # Number of readings per month
target_year = 2029 # Sidereal single year

COMPUTE_SUN = False
COMPUTE_AVG_AN = False
COMPUTE_AN = False

print('----------Computing Swiss Ephemeris Lunar data----------')
#####-------- Get Anomalistic Longitudes --------#####
def get_mandocca_date(jd):
    timeout = 0
    epsilon = 1e-8
    delta = epsilon + 1

    lon_sa = swe.calc(jd, swe.MOON)[0][0]
    mandocca_sa = swe.nod_aps(jd, swe.MOON)[3][0]
    delta = lon_sa - mandocca_sa

    if delta>0:
        jd = jd+ANOM_MT*(360-delta)/360
    if delta<0:
        jd = jd+ANOM_MT*(-delta)/360

    while abs(delta) >= epsilon:
        lon_sa = swe.calc(jd, swe.MOON)[0][0]
        mandocca_sa = swe.nod_aps(jd, swe.MOON)[3][0]
        delta = lon_sa - mandocca_sa
        delta = (delta+180)%360-180 # Bring delta between -180 and 180
        if abs(delta) >= epsilon:
            jd = jd - delta*ANOM_MT/360

        timeout=timeout+1
        if(timeout>=10**4): print('get_mandocca_date: TIMEOUT')
    return jd

jd_start_ut = swe.julday(startyear, 5, 15, 0)
jd_start_tt = jd_start_ut + swe.deltat(jd_start_ut)
jd_lunar_apogee = get_mandocca_date(jd_start_tt)

def get_mean_anom_longs(startyear):
    datas = [[(swe.calc(jd_lunar_apogee+mon*ANOM_MT+i*ANOM_MT/steps, swe.MOON)[0][0] \
            - swe.nod_aps(jd_lunar_apogee+mon*ANOM_MT+i*ANOM_MT/steps, swe.MOON)[3][0])%360 \
                for i in range(steps)] for mon in range(num_months_avg)]
    return [circmean([data[i] for data in datas], high=360, low=0) for i in range(steps)]

if COMPUTE_AVG_AN:
    print('Computing average anomalistic longitudes over %f anomalistic months, %f steps per month...'%(num_months_avg, steps))
    swelongs_mean_an = get_mean_anom_longs(startyear)

#####-------- --------#####
swe.set_sid_mode(mode = swe.SIDM_LAHIRI) # Ayanamsa

jd_start_ut = swe.julday(startyear, 5, 15, 0)
jd_start_tt = jd_start_ut + swe.deltat(jd_start_ut)
jd_lunar_mesadi = swe.mooncross(0, jd_start_tt, swe.FLG_SIDEREAL)

print('Computing sidereal and anomalistic longitudes over %f anomalistic months, %f steps per month...'%(num_months, steps))
print("Starting on jd %f, Moon longitude %f"%(jd_lunar_mesadi, swe.calc(jd_lunar_mesadi, swe.MOON, swe.FLG_SIDEREAL)[0][0]))
swelongs_si = []
swemands_si=[]
swelongs_an=[]
for mon in range(num_months):
    for i in range(steps):
        jd_from_lunar_mesadi = jd_lunar_mesadi+mon*SIDE_MT+i*SIDE_MT/steps
        jd_from_lunar_apogee = jd_lunar_apogee+mon*SIDE_MT+i*SIDE_MT/steps
        swelongs_si.append(swe.calc(jd_from_lunar_mesadi, swe.MOON, swe.FLG_SIDEREAL)[0][0])
        swemands_si.append((swe.nod_aps(jd_from_lunar_mesadi, swe.MOON)[3][0] - swe.get_ayanamsa_ex(tjdet=jd_from_lunar_mesadi, flags=0)[1])%360)
        if COMPUTE_AN: swelongs_an.append((swe.calc(jd_from_lunar_apogee, swe.MOON)[0][0] - swe.nod_aps(jd_from_lunar_apogee, swe.MOON)[3][0])%360)

if COMPUTE_SUN:
    print('Computing Solar longitudes...')
    swesunls_si = []
    for mon in range(num_months):
        for i in range(steps):
            swesunls_si.append(swe.calc(jd_lunar_mesadi+mon*SIDE_MT+i*SIDE_MT/steps, swe.SUN, swe.FLG_SIDEREAL)[0][0])

moon_apogee_at_lunar_mesadi = (swe.nod_aps(jd_lunar_mesadi, swe.MOON)[3][0] - swe.get_ayanamsa_ex(tjdet=jd_lunar_mesadi, flags=0)[1])%360
moon_apogee_at_lunar_apogee = (swe.nod_aps(jd_lunar_apogee, swe.MOON)[3][0] - swe.get_ayanamsa_ex(tjdet=jd_lunar_apogee, flags=0)[1])%360
sun_lunar_mesadi = swe.calc(jd_lunar_mesadi, swe.SUN, swe.FLG_SIDEREAL)[0][0]
sun_lunar_apogee = swe.calc(jd_lunar_apogee, swe.SUN, swe.FLG_SIDEREAL)[0][0]
moon_node_at_lunar_mesadi = (swe.nod_aps(jd_lunar_mesadi, swe.MOON)[0][0] - swe.get_ayanamsa_ex(tjdet=jd_lunar_mesadi, flags=0)[1])%360
moon_node_at_lunar_apogee = (swe.nod_aps(jd_lunar_apogee, swe.MOON)[0][0] - swe.get_ayanamsa_ex(tjdet=jd_lunar_apogee, flags=0)[1])%360

jd_start_ut = swe.julday(startyear, 5, 1, 0)
jd_start_tt = jd_start_ut + swe.deltat(jd_start_ut)
jd_solar_mesadi = swe.solcross(0, jd_start_tt, swe.FLG_SIDEREAL)
moon_solar_mesadi = swe.calc(jd_solar_mesadi, swe.MOON, swe.FLG_SIDEREAL)[0][0]
moon_apogee_solar_mesadi = (swe.nod_aps(jd_solar_mesadi, swe.MOON)[3][0] - swe.get_ayanamsa_ex(tjdet=jd_solar_mesadi, flags=0)[1])%360

print('Computed All Data!')

params_file['steps'] = steps
params_file['num_months'] = num_months

# Longitude Lists
params_file['swelongs_si'] = swelongs_si
params_file['swemands_si'] = swemands_si
if COMPUTE_AVG_AN: params_file['swelongs_mean_an'] = swelongs_mean_an
if COMPUTE_AVG_AN: params_file['swelongs_an'] = swelongs_an
if COMPUTE_SUN: params_file['swesunls_si'] = swesunls_si

# Lunar Mesadi Dhruvas
params_file['swelongs_si_start_jdtt'] = jd_lunar_mesadi
params_file['moon_apogee_at_lunar_mesadi'] = moon_apogee_at_lunar_mesadi
params_file['moon_node_at_lunar_mesadi'] = moon_node_at_lunar_mesadi
params_file['sun_lunar_mesadi'] = sun_lunar_mesadi

# Lunar Apogee Dhruvas
params_file['moon_apogee_at_lunar_apogee'] = moon_apogee_at_lunar_apogee
params_file['moon_node_at_lunar_apogee'] = moon_node_at_lunar_apogee
params_file['sun_lunar_apogee'] = sun_lunar_apogee

# Solar Mesadi Dhruvas
params_file['moon_solar_mesadi'] = moon_solar_mesadi
params_file['moon_apogee_solar_mesadi'] = moon_apogee_solar_mesadi

params_file.close()