from numpy import abs, pi, cos, deg2rad
import shelve
import sys
from scipy.stats import circmean
import matplotlib.pyplot as plt
import swisseph as swe

swe.set_ephe_path('~/swisseph/ephe') # Set Ephemeris Path
swe.set_topo(75.7683, 0.0, 0.0) # Positions as observed at Lanka
swe.set_sid_mode(mode = swe.SIDM_LAHIRI) # Ayanamsa

import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from useful_functions import *
import siddhantaic_functions as sid
import epoch as ep

with shelve.open('./modern_astronomical_constants/modern_astronomical_constants') as modern_astro_data:
    
    ANOM_YR      = modern_astro_data['ANOM_YR']
    SIDE_YR      = modern_astro_data['SIDE_YR']

    sun_revs     = modern_astro_data['sun_revs']
    civ_days_new = modern_astro_data['civ_days_new']
    sun_mandocca_revs_new = modern_astro_data['sun_mandocca_revs_new']

    SIDE_MT      = modern_astro_data['SIDE_MT']
    ANOM_MT      = modern_astro_data['ANOM_MT']
    SYNO_MT      = modern_astro_data['SYNO_MT']

    moon_revs_new = modern_astro_data['moon_revs_new']
    moon_mandocca_revs_new = modern_astro_data['moon_mandocca_revs_new']

modern_astro_data.close()

with shelve.open('./common_data/epochal_values') as dhruvas:
    sun_dhruva              = dhruvas['sun_dhruva']
    sun_mandocca_dhruva     = dhruvas['sun_mandocca_dhruva']
    moon_dhruva             = dhruvas['moon_dhruva']
    moon_mandocca_dhruva    = dhruvas['moon_mandocca_dhruva']
dhruvas.close()

with shelve.open('./sun/sun_computed_data/sun_modified_parameters') as sunp:
    sun_manda_r_popt             = sunp['r0byR_popt']
    sun_mandocca_mean            = sunp['mandocca_mean']
sunp.close()
with shelve.open('./moon/moon_computed_data/moon_revised_parameters') as moonp:
    moon_manda_r_popt           = [7, moonp['manda_radius_oscillation']]
    moon_tungantara_r_by_mean_R = moonp['dvitiya_radius']
    moon_neg_paksika_r_by_tungantara_K = moonp['trtiya_radius']
    moon_digamsa_factor         = moonp['annual_equation_factor']
moonp.close()

print(moon_manda_r_popt, moon_tungantara_r_by_mean_R, moon_neg_paksika_r_by_tungantara_K, moon_digamsa_factor)

#####-------- Set Initial year, averaging window, and number of steps -------#####
COMPUTE_DATA = True
if COMPUTE_DATA:
    print("-----Computing Data-----")
    startyear = 1946  # Saka
    num_years = 400   # Number of years to compare
    steps = 100       # Number of readings per year

    ag_start = ep.get_ahargana(3179+startyear, 0, 0, "tue")
    jdtt_start = ep.ag_audayika2jdtt(ag_start)

    swe_suns_si_long_term = []
    swe_moons_si_long_term = []
    for year in range(num_years):
        for i in range(steps):
            swe_suns_si_long_term.append(swe.calc(jdtt_start+year*SIDE_YR+i*SIDE_YR/steps, swe.SUN, swe.FLG_SIDEREAL)[0][0])
            swe_suns_si_long_term.append(swe.calc(jdtt_start+year*SIDE_YR+i*SIDE_YR/steps, swe.MOON, swe.FLG_SIDEREAL)[0][0])

    long_term_datas = shelve.open('./common_data/long_term_data')
    
    long_term_datas['ag_start'] = ag_start
    long_term_datas['jdtt_start'] = jdtt_start
    long_term_datas['num_years'] = num_years
    long_term_datas['steps'] = steps

    long_term_datas['swesuns_si_600y_100s'] = swe_suns_si_long_term
    long_term_datas['swemoons_si_600y_100s'] = swe_moons_si_long_term
    long_term_datas.close()

else:
    with shelve.open('./common_data/long_term_data') as long_term_datas:
        ag_start = long_term_datas['ag_start']
        jdtt_start = long_term_datas['jdtt_start']
        num_years = long_term_datas['num_years']
        steps = long_term_datas['steps']

        swesuns_si_600y_100s  = long_term_datas['swesuns_si_600y_100s']
        swemoons_si_600y_100s = long_term_datas['swemoons_si_600y_100s']
    long_term_datas.close()

def sun(_ag):
    mean_sun = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
    sun_mandocca = sid.theta_0(_ag, sun_mandocca_revs_new, sun_mandocca_dhruva, civ_days_new)
    manda_r = sid.r0_revised(mean_sun, sun_mandocca, *sun_manda_r_popt)
    return sid.theta_ms(mean_sun, sun_mandocca, manda_r/360)

def approx_sun(_ag):
    mean_sun = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
    sun_mandocca = sun_mandocca_mean
    manda_r = sid.r0_revised(mean_sun, sun_mandocca, *sun_manda_r_popt)
    return sid.theta_ms(mean_sun, sun_mandocca, manda_r/360)

def moon(_ag):
    mean_moon = sid.theta_0(_ag, moon_revs_new, moon_dhruva, civ_days_new)
    mean_R = 216000/(2*pi)

    moon_mandocca = sid.theta_0(_ag, moon_mandocca_revs_new, moon_mandocca_dhruva, civ_days_new)
    manda_r = sid.r0_revised(mean_moon, moon_mandocca, *moon_manda_r_popt)
    manda_moon = sid.theta_ms(mean_moon, moon_mandocca, manda_r/80)
    manda_K = 10*mean_R*sid.KbyR(mean_moon, moon_mandocca, manda_r/80)

    moon_dvitiyocca = approx_sun(_ag)
    tungantara_r = (mean_R*moon_tungantara_r_by_mean_R)*cos(deg2rad(moon_dvitiyocca-moon_mandocca))
    tungantara_moon = sid.theta_ms(manda_moon, moon_dvitiyocca, tungantara_r/manda_K)
    tungantara_K = sid.KbyR(manda_moon, moon_dvitiyocca, tungantara_r/manda_K)

    moon_trtiyocca = 2*moon_dvitiyocca-tungantara_moon
    paksika_r = -moon_neg_paksika_r_by_tungantara_K*tungantara_K
    paksika_moon = sid.theta_ms(tungantara_moon, moon_trtiyocca, paksika_r/tungantara_K)

    # return paksika_moon
    mean_sun = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
    digamsa = moon_digamsa_factor*((mean_sun-sun(_ag)+180)%360-180)

    return paksika_moon+digamsa

print('----- Computing ags and jds -----')
aharganas = [ag_start+i*SIDE_YR/steps for i in range(num_years*steps)]
juldays   = [jdtt_start+i*SIDE_YR/steps for i in range(num_years*steps)]

if True:
    print('\n----- Plot Solar Longitudes -----\n')
    print(sun(ag_start)-swe.calc(jdtt_start, swe.SUN, swe.FLG_SIDEREAL)[0][0])
    plt.plot(aharganas, subcirc_secs([approx_sun(ag) for ag in aharganas], [swe.calc(jd, swe.SUN, swe.FLG_SIDEREAL)[0][0] for jd in juldays]), '.-', label='Constant Mean Mandocca')
    plt.plot(aharganas, subcirc_secs([sun(ag) for ag in aharganas], [swe.calc(jd, swe.SUN, swe.FLG_SIDEREAL)[0][0] for jd in juldays]), '.-', label='Varying Mandocca')
    plt.title('SUN')
    plt.ylabel('Difference between Siddhantaic and observed longitudes (seconds)')
    plt.xlabel('Ahargana')
    plt.legend()
    plt.grid()
    plt.gcf().set_size_inches(18, 10)
    plt.savefig('./sun/graphs/sun_long_term.png')
    plt.show()
    plt.clf()
    print("Moving mandocca:")
    print("Error less than 36'' till",swe.revjul(ep.ag_audayika2jdut(1925000)))
    print("Error less than 1' till",swe.revjul(ep.ag_audayika2jdut(1967000)))
    print("Stationary mandocca:")
    print("Error less than 36'' till",swe.revjul(ep.ag_audayika2jdut(1906000)))
    print("Error less than 1' till",swe.revjul(ep.ag_audayika2jdut(1927000)))


print('\n----- Plot Lunar Longitudes -----\n')

plt.plot(aharganas, subcirc_mins([moon(ag) for ag in aharganas], [swe.calc(jd, swe.MOON, swe.FLG_SIDEREAL)[0][0] for jd in juldays]), '.-')
plt.title('MOON')
plt.ylabel('Difference between Siddhantaic and observed longitudes (minutes)')
plt.xlabel('Ahargana')
plt.legend()
plt.grid()
plt.gcf().set_size_inches(18, 10)
# plt.savefig('./moon/graphs/moon_long_term.png')
plt.show()
plt.clf()