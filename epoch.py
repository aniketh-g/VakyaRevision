import numpy as np
import swisseph as swe
swe.set_ephe_path('~/swisseph/ephe') # Set Ephemeris Path
swe.set_sid_mode(mode = swe.SIDM_LAHIRI) # Ayanamsa

ujjain_long = 75+46/60+6/3600 # 75+45/60+0/3600 # in K Chandra Hari
gmt_ujj_delt = 24*ujjain_long/360

jd_kali_epoch_audayika=swe.julday(-3101, 2, 18, 6-gmt_ujj_delt, 0) # Proleptic Julian calendar
jd_kali_epoch_ardharatra=swe.julday(-3101, 2, 17, 12-gmt_ujj_delt, 0) # Proleptic Julian calendar
if False:
    print("JD on Audayika Ahargana = 0: UT %.10f | TT %.10f"%(jd_kali_epoch_audayika, jd_kali_epoch_audayika+swe.deltat(jd_kali_epoch_audayika)))
    print("DeltaT on this day, hours", swe.deltat(jd_kali_epoch_audayika)*24)
    print("JD on Ardharatra Ahargana = 0: UT %.10f | TT %.10f"%(jd_kali_epoch_ardharatra, jd_kali_epoch_ardharatra+swe.deltat(jd_kali_epoch_ardharatra)))
    print("DeltaT on this day, hours", swe.deltat(jd_kali_epoch_ardharatra)*24)
    print('\n')

import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from useful_functions import subcirc, todms, limsmall
import siddhantaic_functions as sid

#####------------------------ Import modern constants ------------------------#####
import shelve
with shelve.open('./modern_astronomical_constants/modern_astronomical_constants') as modern_astro_params_file:
    SIDE_YR = modern_astro_params_file['SIDE_YR']
    SIDE_MT = modern_astro_params_file['SIDE_MT']
    
    sun_revs = modern_astro_params_file['sun_revs']
    civ_days_new = modern_astro_params_file['civ_days_new']
    sun_mandocca_revs_new = modern_astro_params_file['sun_mandocca_revs_new']

    SIDE_MT = modern_astro_params_file['SIDE_MT']
    ANOM_MT = modern_astro_params_file['ANOM_MT']
    SYNO_MT = modern_astro_params_file['SYNO_MT']

    moon_revs_new = modern_astro_params_file['moon_revs_new']
    moon_mandocca_revs_new = modern_astro_params_file['moon_mandocca_revs_new']
modern_astro_params_file.close()

#####------------------------- Import observed data --------------------------#####
with shelve.open('./sun/sun_observed_data/sidereal_data') as sun_side_data:
    jdtt_sun_curvefit_start = sun_side_data['swelongs_si_start_jdtt']
sun_side_data.close()
with shelve.open('./sun/sun_computed_data/sun_modified_parameters') as sunp:
    mean_sun_sun_curvefit_start = sunp['sun_dhruva']
    sun_mandocca_sun_curvefit_start = sunp['mandocca_dhruva']
sunp.close()

with shelve.open('./moon/moon_computed_data/generated_data') as moon_gen_dat:
    jdtt_moon_curvefit_start = moon_gen_dat['swelongs_si_start_jdtt']
moon_gen_dat.close()
with shelve.open('./moon/moon_computed_data/moon_revised_parameters') as moonp:
    mean_moon_curvefit_start = moonp['moon_dhruva']
    moon_mandocca_curvefit_start = moonp['moon_mandocca_dhruva']
moonp.close()
#####------------------------ Save Computed Dhruvas --------------------------#####

#####-------------------------- Computing Ahargana ---------------------------#####
# total_civ_days = 1577917500
# moon_revs = 57753320
total_civ_days = civ_days_new
moon_revs = moon_revs_new

total_solar_months = 12*sun_revs
total_lunar_months = moon_revs - sun_revs
total_adhimasas = total_lunar_months - total_solar_months
total_tithis = 30*total_lunar_months
total_ksayatithis = total_tithis - total_civ_days

def get_ahargana(_elapsed_kali_years, _elapsed_lunar_months_present_kali, _elapsed_tithis_present_month, _weekday="null", verbose = False):

    day2num = {"fri":0, "sat":1, "sun":2, "mon":3, "tue":4, "wed":5, "thu":6}
    num2day = {0:"fri", 1:"sat", 2:"sun", 3:"mon", 4:"tue", 5:"wed", 6:"thu"}

    _elapsed_solar_months_from_kali = 12*_elapsed_kali_years+_elapsed_lunar_months_present_kali # m
    _elapsed_adhimasas_from_kali = _elapsed_solar_months_from_kali*total_adhimasas/total_solar_months # a
    _elapsed_lunar_months_from_kali = _elapsed_solar_months_from_kali + int(_elapsed_adhimasas_from_kali) # l
    _elapsed_tithis_from_kali = 30*_elapsed_lunar_months_from_kali + _elapsed_tithis_present_month # A'
    _elapsed_ksayatithis_from_kali = _elapsed_tithis_from_kali*total_ksayatithis/total_tithis # k
    
    _ahargana = _elapsed_tithis_from_kali - int(_elapsed_ksayatithis_from_kali)
    if verbose:
        print("\nComputing for Kali Year %d, Month %d, elapsed tithis %d"%(_elapsed_kali_years, _elapsed_lunar_months_present_kali, _elapsed_tithis_present_month))
        print("_elapsed_adhimasas_from_kali,\ta ", _elapsed_adhimasas_from_kali)
        print("_elapsed_lunar_months_from_kali,l ", _elapsed_lunar_months_from_kali)
        print("_elapsed_lunar_months_from_kali,A'", _elapsed_tithis_from_kali)
        print("_elapsed_ksayatithis_from_kali,\tk ", _elapsed_ksayatithis_from_kali)
        print("_ahargana without corrections,\tA ", _ahargana)
        print("computed day of week ", num2day[int(_ahargana)%7])

    if _weekday != "null":
        _weekday_num = day2num[_weekday]
        count = 0
        while np.abs(int(_ahargana)%7 - _weekday_num) != 0:
            if count == 0:
                if verbose: print("Rounding up k")
                _ahargana = _elapsed_tithis_from_kali - (int(_elapsed_ksayatithis_from_kali) + 1)
                if verbose: print("Computed weekday %s"%num2day[int(_ahargana)%7])
            if count == 1:
                if verbose: print("Rounding up a")
                _elapsed_lunar_months_from_kali = _elapsed_solar_months_from_kali + (int(_elapsed_adhimasas_from_kali) + 1) # l
                _elapsed_tithis_from_kali = 30*_elapsed_lunar_months_from_kali + _elapsed_tithis_present_month # A'
                _elapsed_ksayatithis_from_kali = _elapsed_tithis_from_kali*total_ksayatithis/total_tithis # k
                _ahargana = _elapsed_tithis_from_kali - int(_elapsed_ksayatithis_from_kali)
                if verbose: print("Computed weekday %s"%num2day[int(_ahargana)%7])
            if count == 2:
                if verbose: print("Rounding up both k and a")
                _elapsed_lunar_months_from_kali = _elapsed_solar_months_from_kali + int(_elapsed_adhimasas_from_kali+1) # l
                _elapsed_tithis_from_kali = 30*_elapsed_lunar_months_from_kali + _elapsed_tithis_present_month # A'
                _elapsed_ksayatithis_from_kali = _elapsed_tithis_from_kali*total_ksayatithis/total_tithis # k
                _ahargana = _elapsed_tithis_from_kali - (int(_elapsed_ksayatithis_from_kali) + 1)
                if verbose: print("Computed weekday %s"%num2day[int(_ahargana)%7])
            if count >= 5:
                print("WARNING: WRONG WEEKDAY PROVIDED FOR LUNAR DATE")
                break
            count = count+1
        if verbose: print("ahargana after correction: %d, %s"%(_ahargana, num2day[int(_ahargana)%7]))
    return int(_ahargana)

# ag_recent_epoch = get_ahargana(3179+1945, 11, 29, "mon")
ag_recent_epoch = get_ahargana(3179+1946, 0, 0, "tue")

jdut_recent_epoch_audayika=swe.julday(2024, 4, 9, 6-gmt_ujj_delt) # Greogorian calendar
jdut_recent_epoch_ardharatra=swe.julday(2024, 4, 9, 12-gmt_ujj_delt) # Gregorian calendar
if False:
    print("JD on Audayika Recent Epoch: UT %.10f | TT %.10f"%(jdut_recent_epoch_audayika, jdut_recent_epoch_audayika+swe.deltat(jdut_recent_epoch_audayika)))
    print("DeltaT on this day, hours", swe.deltat(jdut_recent_epoch_audayika)*24)
    print("JD on Ardharatra Recent Epoch: UT %.10f | TT %.10f"%(jdut_recent_epoch_ardharatra, jdut_recent_epoch_ardharatra+swe.deltat(jdut_recent_epoch_ardharatra)))
    print("DeltaT on this day, hours", swe.deltat(jdut_recent_epoch_ardharatra)*24)
    print('\n')

jdtt_recent_epoch_audayika = jdut_recent_epoch_audayika + swe.deltat(jdut_recent_epoch_audayika)
jdtt_recent_epoch_ardharatra = jdut_recent_epoch_ardharatra + swe.deltat(jdut_recent_epoch_ardharatra)

jdut2ag_audayika = lambda jdut: jdut - (jdut_recent_epoch_audayika-ag_recent_epoch)
jdtt2ag_audayika = lambda jdtt: jdtt - (jdtt_recent_epoch_audayika-ag_recent_epoch)
jdut2ag_ardhratra = lambda jdut: jdut2ag_audayika(jdut) - 0.25
jdtt2ag_ardharatra = lambda jdtt: jdtt2ag_audayika(jdtt) - 0.25

ag_audayika2jdtt = lambda ag_audayika: ag_audayika + (jdtt_recent_epoch_audayika-ag_recent_epoch)
def ag_audayika2jdut(ag_audayika):
    jdtt = ag_audayika2jdtt(ag_audayika)
    jdut = jdtt
    delta = jdtt-(jdut+swe.deltat(jdut))
    while np.abs(delta) >= 1e-5:
        jdut = jdut + delta
        delta = jdtt-(jdut+swe.deltat(jdut))
    return jdut

if __name__ == '__main__':
    epoch_file = shelve.open('./common_data/epochal_values')
    print("Epoch: %s = %.10f JDTT = %.10f JDUT = %f AG, Delta T = %f"%\
        (swe.revjul(jdut_recent_epoch_audayika), jdtt_recent_epoch_audayika, jdut_recent_epoch_audayika, ag_recent_epoch, swe.deltat(jdut_recent_epoch_audayika)))
    print("Conversion: AG = JDTT - %.12f"%(jdtt_recent_epoch_audayika-ag_recent_epoch))
    epoch_file['jdtt_recent_epoch_audayika_minus_ag_recent_epoch'] = jdtt_recent_epoch_audayika-ag_recent_epoch

if __name__ == '__main__':
    #####------------------------- Computing Sun Dhruva --------------------------#####
    print('\n#####------------------------- Computing Sun Dhruva --------------------------#####\n')

    ag_sun_curvefit_start = jdtt2ag_audayika(jdtt_sun_curvefit_start)

    print("Sun CurveFit Starts: %s = %.10f JDTT = %.10f AG"%(swe.revjul(ag_audayika2jdut(ag_sun_curvefit_start)), jdtt_sun_curvefit_start, ag_sun_curvefit_start))
    print(swe.calc(jdtt_sun_curvefit_start, swe.SUN, swe.FLG_SIDEREAL)[0][0])

    sun_dhruva = (mean_sun_sun_curvefit_start - sid.theta_0(ag_sun_curvefit_start, sun_revs, 0, civ_days_new))%360
    sun_mandocca_dhruva = sun_mandocca_sun_curvefit_start - sid.theta_0(ag_sun_curvefit_start, sun_mandocca_revs_new, 0, civ_days_new)

    epoch_file['sun_dhruva'] = sun_dhruva
    epoch_file['sun_mandocca_dhruva'] = sun_mandocca_dhruva

    print('sun_dhruva = -', todms(360-sun_dhruva))
    print('sun_mandocca_dhruva =', todms(sun_mandocca_dhruva))

    #####------------------------ Computing Moon Dhruva --------------------------#####
    print('\n#####------------------------ Computing Moon Dhruva --------------------------#####\n')

    ag_moon_curvefit_start = jdtt2ag_audayika(jdtt_moon_curvefit_start)

    print("Moon CurveFit Starts: %s = %.10f JDTT = %.10f AG"%(swe.revjul(jdtt_moon_curvefit_start), jdtt_moon_curvefit_start, ag_moon_curvefit_start))
    print(swe.calc(jdtt_moon_curvefit_start, swe.MOON, swe.FLG_SIDEREAL)[0][0])

    moon_dhruva = mean_moon_curvefit_start
    moon_mandocca_dhruva = moon_mandocca_curvefit_start
    print(moon_dhruva)
    epoch_file['moon_dhruva'] = moon_dhruva
    epoch_file['moon_mandocca_dhruva'] = moon_mandocca_dhruva

    print('moon_dhruva = -', todms(-moon_dhruva))
    print('moon_mandocca_dhruva = ', todms(moon_mandocca_dhruva))

    #####------------------------ Computing Moon Sodhya --------------------------#####
    print('\n#####------------------------ Computing Moon Sodhya --------------------------#####\n')

    #####-------- Get Anomalistic Longitudes --------#####
    def get_mean_mandocca_date(ag):
        timeout = 0
        epsilon = 1e-8
        delta = epsilon + 1

        mean_moon = sid.theta_0(ag, moon_revs_new, moon_dhruva, civ_days_new)
        moon_mandocca = sid.theta_0(ag, moon_mandocca_revs_new, moon_mandocca_dhruva, civ_days_new)
        delta = mean_moon - moon_mandocca

        if delta>0:
            ag = ag+ANOM_MT*(360-delta)/360
        if delta<0:
            ag = ag+ANOM_MT*(-delta)/360

        while abs(delta) >= epsilon:
            mean_moon = sid.theta_0(ag, moon_revs_new, moon_dhruva, civ_days_new)
            moon_mandocca = sid.theta_0(ag, moon_mandocca_revs_new, moon_mandocca_dhruva, civ_days_new)
            delta = mean_moon - moon_mandocca
            delta = (delta+180)%360-180 # Bring delta between -180 and 180
            if abs(delta) >= epsilon:
                ag = ag - delta*ANOM_MT/360

            timeout=timeout+1
            if(timeout>=10**4): print('get_mandocca_date: TIMEOUT')
        return ag

    anom_period = 1/(moon_revs_new/civ_days_new - moon_mandocca_revs_new/civ_days_new)
    ag_recent_moon_apogee = get_mean_mandocca_date(ag_recent_epoch-500*anom_period)
    for i in range(600):
        ag = ag_recent_moon_apogee + i*anom_period
        mean_moon = sid.theta_0(ag, moon_revs_new, moon_dhruva, civ_days_new)
        moon_apogee = sid.theta_0(ag, moon_mandocca_revs_new, moon_mandocca_dhruva, civ_days_new)
        mean_sun = sid.theta_0(ag, sun_revs, sun_dhruva, civ_days_new)
        if(abs(mean_sun-mean_moon)>=1 and abs(mean_sun-mean_moon)<=3):
            print(ag, swe.revjul(ag_audayika2jdut(ag)), mean_moon, moon_apogee, mean_sun, mean_sun-mean_moon)

    epoch_file['moon_sodhyadina_ag'] = 1865535.4569303067
    epoch_file.close()