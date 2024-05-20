import numpy as np
from numpy import sqrt, sin, cos, arcsin, pi, rad2deg, deg2rad
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import statistics as stat
import shelve

import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from useful_functions import *
import siddhantaic_functions as sid
import sun.functions as sun

#####------------------------ Import modern constants ------------------------#####
print('#####------------------------ Import modern constants ------------------------#####')
with shelve.open('../modern_astronomical_constants/modern_astronomical_constants') as astro_const:
    SIDE_MT = astro_const['SIDE_MT']
    ANOM_MT = astro_const['ANOM_MT']
    SIDE_YR = astro_const['SIDE_YR']

    civ_days_new = astro_const['civ_days_new']
    sun_revs = astro_const['sun_revs']
    moon_revs_new = astro_const['moon_revs_new']
    moon_apogee_revs_new = astro_const['moon_mandocca_revs_new']
astro_const.close()
#####----------------------- Import Revised Parameters -----------------------#####
print('#####----------------------- Import Revised Parameters -----------------------#####')
with shelve.open('../sun/sun_computed_data/sun_modified_parameters') as sunp:
    sun_civ_days = sunp['sun_civ_days']
    sun_mandocca_mean = sunp['mandocca_mean']
    sun_mandocca_revs = sunp['mandocca_revs']
    sun_r0byR_popt = sunp['r0byR_popt']
sunp.close()
#####------------------------- Import observed data --------------------------#####
print('------------------------- Import observed data --------------------------')
with shelve.open('../moon/moon_computed_data/generated_data') as moon_obs_data:
    steps = moon_obs_data['steps']
    num_months = moon_obs_data['num_months']
    
    jdtt_moon_curvefit_start = moon_obs_data['swelongs_si_start_jdtt']
    # swelongs_an = gen_dat['swelongs_an']
    moons_swe_si = moon_obs_data['swelongs_si']
    moon_apogee_at_lunar_mesadi = moon_obs_data['moon_apogee_at_lunar_mesadi']
    sun_lunar_mesadi = moon_obs_data['sun_lunar_mesadi']
    swemands_si = moon_obs_data['swemands_si']
    swesunls_si = moon_obs_data['swesunls_si']
    
    moon_apogee_at_lunar_apogee = moon_obs_data['moon_apogee_at_lunar_apogee']

    moon_solar_mesadi = moon_obs_data['moon_solar_mesadi']
    moon_apogee_solar_mesadi = moon_obs_data['moon_apogee_solar_mesadi']
    
    moon_node_at_lunar_mesadi = moon_obs_data['moon_node_at_lunar_mesadi']
    moon_node_at_lunar_apogee = moon_obs_data['moon_node_at_lunar_apogee']
moon_obs_data.close()
#####------------------------------ Plot Settings ----------------------------#####
plt.rcParams['axes.grid'] = True

#####------------------- Start Revision of Moon parameters -------------------#####
print('\n------------------- Start Revision of Moon parameters -------------------\n')
print('\nWARNING: USING DATA FROM FILE. Recompile generate_data.py to get latest changes\n')
#####---------------- Moon's treatment according to Siddhanta ----------------#####

# Old values
moon_revs        = 57753320
moon_apogee_revs = 488122
moon_node_revs   = 232300
civ_days = 1577917500
moonR = 10*216000/(2*pi)
# Mandasamskara
#####----- Compute new parameters by comparing in one anomalistic cycle -----#####
anom_rev_diff = civ_days/ANOM_MT
#####-------------------- Effect of changing month length -------------------#####
print('#####-------------------- Effect of changing month length -------------------##### -- !!!SKIPPED!!!')
if False:
    sidlongs = [sid.theta_ms(sid.theta_0(ag, moon_revs - moon_apogee_revs), 0, 7/80) for ag in np.arange(0, ANOM_MT, ANOM_MT/steps)]
    plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='moon_revs = %f (KAPA)'%(moon_revs - moon_apogee_revs))
    print('moon_revs - moon_apogee_revs: old = %f'%(moon_revs - moon_apogee_revs))

    sidlongs = [sid.theta_ms(sid.theta_0(ag,anom_rev_diff), 0, 7/80) for ag in np.arange(0, ANOM_MT, ANOM_MT/steps)]
    plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='moon_revs = %f (1961AD)'%anom_rev_diff)
    print('moon_revs - moon_apogee_revs: new = %f'%(anom_rev_diff))

    plt.title('Effect of changing length of anomalistic month')
    plt.ylabel('Difference between observed and siddhanta anomalistic longitudes (degrees)')
    plt.xlabel('Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.legend()
    plt.gcf().set_size_inches(8, 4)
    plt.savefig('./graphs/month_length.png')
    # plt.show()
    plt.clf()

#####---------------------------- Comparing Mandocca ---------------------------#####
print('#####---------------------------- Comparing Mandocca ---------------------------#####')

# moon_mandasphuta = lambda _ag: sid.theta_ms(sid.theta_0(_ag, moon_revs),  sid.theta_0(_ag, moon_apogee_revs, mandocca_true_apogee), sid.r0_revised(sid.theta_0(_ag, moon_revs), sid.theta_0(_ag, moon_apogee_revs, mandocca_true_apogee), *r0byR_popt)/80)
moon_dhruva = 0
mandocca_dhruva = moon_apogee_at_lunar_mesadi
moon_mandasphuta_old = lambda _ag: sid.theta_ms(sid.theta_0(_ag, moon_revs, moon_dhruva),  sid.theta_0(_ag, moon_apogee_revs, mandocca_dhruva), 7/80)

ag_start = sid.get_ag(moon_mandasphuta_old, 0, 0, (SIDE_MT/360)/5)
mandocca_dhruva = sid.theta_0(-ag_start, moon_apogee_revs, moon_apogee_at_lunar_mesadi)

if False:
    moon_apogees_madhya = []
    for mon in range(num_months):
        for i in range(steps):
            ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
            lu = sid.theta_0(ag, moon_apogee_revs, mandocca_dhruva)
            moon_apogees_madhya.append(lu)

    plt.plot(range(steps*num_months), subcirc(moon_apogees_madhya, swemands_si), '-', label='Mandocca Longitude Difference')
    plt.ylabel('Difference between observed and siddhanta sidereal longitudes of apogee (degrees)')
    plt.xlabel('Time since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.legend()

    plt.gcf().set_size_inches(8, 4)
    plt.savefig('./graphs/Mandocca_comparison.png')
    # plt.show()
    plt.clf()

#####---------------------- Ucca Analysis with madhyagraha ---------------------#####
print('---------------------- Ucca Analysis with madhyagraha ---------------------Skipped!')
moons_madhya = []
moons_manda = []
moon_mandaKbyRs = []
if False:
    for mon in range(num_months):
        for i in range(steps):
            ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
            moons_madhya.append(sid.theta_0(ag, moon_revs, 0))
            moons_manda.append(moon_mandasphuta_old(ag))
            # lu = moon_apogees_madhya[steps*mon+i]
            lu = sid.theta_0(ag, moon_apogee_revs, mandocca_dhruva)
            moon_mandaKbyRs.append(sid.KbyR(sid.theta_0(ag, moon_revs), lu, 7/80))
if False:
    ## Get zeroes of an array
    def get_zeroes(errors_vs_time):
        zeroes = []
        for mon in range(num_months):
            for i in range(steps):
                ag_present = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
                present_index = steps*mon+i
                ag_next = ag_start+(civ_days/moon_revs)*mon+(i+1)*(civ_days/moon_revs)/steps
                
                if(present_index < len(errors_vs_time)-1):
                    if(errors_vs_time[present_index]*errors_vs_time[present_index+1] <= 0):
                        zeroes.append((ag_present*errors_vs_time[present_index+1] - ag_next*errors_vs_time[present_index])/(errors_vs_time[present_index+1]-errors_vs_time[present_index]))
        return zeroes

    moon_madhya_errors = subcirc(moons_madhya, moons_swe_si)
    print('avg of madhya errors', stat.mean(moon_madhya_errors))
    moon_dhruva = -stat.mean(moon_madhya_errors)

    madhya_minus_swe_zero_ags = get_zeroes(subcirc(moons_madhya, moons_swe_si))
    total_ags = [ag_start+i*(civ_days/moon_revs)/steps for i in range(steps*num_months)]

    moon_alt_180 = []
    for i in range(len(madhya_minus_swe_zero_ags)):
        ag = madhya_minus_swe_zero_ags[i]
        if i%2 == 0:
            moon_alt_180.append((sid.theta_0(ag, moon_revs, 0)-180)%360)
        else:
            moon_alt_180.append(sid.theta_0(ag, moon_revs, 0)%360)

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,9), gridspec_kw={'height_ratios': [3, 1]})

    ax2.plot(total_ags, subcirc(moons_madhya, moons_swe_si), '-', label='theta_0 - theta_swe')
    ax2.plot(madhya_minus_swe_zero_ags, [0 for ag in madhya_minus_swe_zero_ags], '.', label = 'zeroes computed')
    ax2.plot(total_ags, [(sid.theta_0(ag, moon_revs, 0)-sid.theta_ms(sid.theta_0(ag, moon_revs, 0),  lu, 7/80)+180)%360-180 for ag in total_ags], '-', label = 'mandaphala')
    ax2.legend()

    ax1.plot(madhya_minus_swe_zero_ags, [sid.theta_0(ag, moon_revs, 0)%360 for ag in madhya_minus_swe_zero_ags], '.-', label = 'Mean Moon at zeroes computed')
    ax1.plot(madhya_minus_swe_zero_ags, moon_alt_180, '.-', label = 'Mean Moon at zeroes computed; 180 added at alernate places')
    ax1.plot(total_ags, [sid.theta_0(ag, moon_apogee_revs, mandocca_dhruva) for ag in total_ags], '-', label = 'Actual Mandocca')
    ax1.legend()

    plt.ylabel('Difference between observed and madhya sidereal longitudes of Moon (degrees)')
    plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.savefig('./graphs/Mean_comparison.png')
    # plt.show()
    plt.clf()

    plt.plot(total_ags, subcirc(moons_manda, moons_swe_si), '-', label='mandasphuta Longitude Difference')
    plt.ylabel('Difference between observed and mandasphuta sidereal longitudes of Moon (degrees)')
    plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.legend()
    plt.gcf().set_size_inches(8, 4)
    plt.savefig('./graphs/Mandasphuta_comparison.png')
    # plt.show()
    plt.clf()

    print('\n')
    print('Siddhanta Mandocca starts at %f deg'%moon_apogees_madhya[0])
    print('SwissEphe Mandocca starts at %f deg'%swemands_si[0])
    print('Moon at ag %f = %f deg'%(ag_start, sid.theta_ms(sid.theta_0(ag, moon_revs),  sid.theta_0(ag, moon_apogee_revs, mandocca_dhruva), 7/80)))
    print('Mandocca at ag %f = %f deg'%(ag_start, sid.theta_0(ag_start, moon_apogee_revs, mandocca_dhruva)))
    print('Mandocca Dhruva = %f deg'%mandocca_dhruva)

#####---------------------------- Comparing Evection ---------------------------#####
print('\n---------------------------- Comparing Evection ---------------------------\n')

Rem = 10*moonR
dvitiyasphutas_mandamoon_mandasun = []
evections_mandamoon_mandasun = []

sun_mandasphuta_func = lambda _ag: sid.theta_ms(sid.theta_0(_ag, sun_revs, 0, civ_days_new),  sun_mandocca_mean, sid.r0_revised(sid.theta_0(_ag, sun_revs, 0, civ_days_new), sun_mandocca_mean, *sun_r0byR_popt)/360)
moon_mandocca = lambda _ag: sid.theta_0(_ag, moon_apogee_revs_new, 0, civ_days_new)

sun_dhruva = sid.theta_0(sid.get_ag(sun_mandasphuta_func,sun_lunar_mesadi, 0, 1)-ag_start, sun_revs, 0, civ_days_new)
if False:
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    ax1.plot(range(steps*num_months), subcirc(moons_manda, moons_swe_si), '-', label='$\\theta_{{MS}}-\\theta_{{jpl}}$')
    ax1.legend()

def sun_mandasphuta(_ag):
    t0 = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
    return sid.theta_ms(t0,  sun_mandocca_mean, sid.r0_revised(t0, sun_mandocca_mean, *sun_r0byR_popt)/360)

if False:
    dvitiyasphutas_mandamoon_mandasun = []
    dvitiyasphutas_mean = []
    dvitiyasphutas_meanmoon_meansun = []
    dvitiyasphutas_mean_mandocca = []
    dvitiyasphutas_meanmoon = []

    evections_mandamoon_mandasun = []
    evections_mandamoon_meansun = []
    evections_meanmoon_meansun = []
    evections_mean_mandocca = []
    evections_meanmoon = []

    mandasuns = []
    madhyasuns = []

    f=10
    for mon in range(num_months):
        for i in range(steps):
            ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps

            lm_mean = moons_madhya[steps*mon+i]
            lu = swemands_si[steps*mon+i]

            lm_manda = moons_manda[steps*mon+i]
            mandaKbyR = moon_mandaKbyRs[steps*mon+i]

            ls_manda = sun_mandasphuta(ag)
            ls_mean = sid.theta_0(ag, civ_days/SIDE_YR, sun_dhruva)

            r = (moonR/2)*cos(deg2rad(ls_manda-lu))
            theta_ds = sid.theta_ms(lm_manda, ls_manda, r/(f*moonR*mandaKbyR))
            dvitiyasphutas_mandamoon_mandasun.append(theta_ds)
            evections_mandamoon_mandasun.append((lm_manda-theta_ds+180)%360-180)

            # theta_ds = sid.theta_ms(lm, ls, r/(f*moonR*mandaKbyR))
            # dvitiyasphutas_mean.append(theta_ds)
            # evections_mandamoon_meansun.append((lm_manda-theta_ds+180)%360-180)

            r = (moonR/2)*cos(deg2rad(ls_mean-lu))
            theta_ds = -lm_mean + lm_manda + sid.theta_ms(lm_mean, ls_mean, r/(f*moonR*mandaKbyR))
            dvitiyasphutas_meanmoon_meansun.append(theta_ds)
            evections_meanmoon_meansun.append((lm_manda-theta_ds+180)%360-180)

            # lu = moon_apogees_madhya[steps*mon+i]
            # theta_ds = sid.theta_ms(lm, ls, r/(f*moonR*mandaKbyR))
            # dvitiyasphutas_mean_mandocca.append(theta_ds)
            # evections_mean_mandocca.append((lm-theta_ds+180)%360-180)

            # lm = moons_madhya[steps*mon+i]
            # theta_ds = sid.theta_ms(lm, ls, r/(f*moonR))
            # dvitiyasphutas_meanmoon.append(moons_manda[steps*mon+i]+theta_ds-lm)
            # evections_meanmoon.append((lm-theta_ds+180)%360-180)

    if False:
        # ax2.plot(range(steps*num_months), evections, '-', label='Dvitiyaphala with True Moon, True Sun Rem = %f*R'%f)
        ax2.plot(range(steps*num_months), evections_mandamoon_mandasun, '-', label='$\\theta_{{MS}}-\\theta_{{DS}}$ with Manda Moon, Manda Sun')
        ax2.plot(range(steps*num_months), evections_meanmoon_meansun, '-', label='$\\theta_{{MS}}-\\theta_{{DS}}$ with Mean Moon, Mean Sun')
        # ax2.plot(range(steps*num_months), evections_manda, '-', label='Dvitiyaphala with Manda Moon, Manda Sun')
        # ax2.plot(range(steps*num_months), evections_mean_mandocca, '-', label='Dvitiyaphala with Mean mandocca, Rem = %f*R'%f)
        # ax2.plot(range(steps*num_months), evections_meanmoon, '-', label='Dvitiyaphala with Mean Moon, Manda Sun Rem = %f*R'%f)
        ax2.legend()

        # ax3.plot(range(steps*num_months), subcirc(moons_dvitiya, moons_swe_si), '-', label='Dvitiyasphuta, True Sun')
        ax3.plot(range(steps*num_months), subcirc(dvitiyasphutas_mandamoon_mandasun, moons_swe_si), '-', label='$\\theta_{{DS}}-\\theta_{{jpl}}$ with Manda Moon, Manda Sun')
        ax3.plot(range(steps*num_months), subcirc(dvitiyasphutas_meanmoon_meansun, moons_swe_si), '-', label='$\\theta_{{DS}}-\\theta_{{jpl}}$ with Mean Moon, Mean Sun')
        # ax3.plot(range(steps*num_months), subcirc(dvitiyasphutas_mean_mandocca, moons_swe_si), '-', label='Dvitiyasphuta, Mean mandocca')
        # ax3.plot(range(steps*num_months), subcirc(dvitiyasphutas_meanmoon, moons_swe_si), '-', label='Dvitiyasphuta, Mean moon')
        ax3.legend()

        plt.ylabel('Degrees')
        plt.xlabel('Time since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
        plt.gcf().set_size_inches(12, 6)
        plt.savefig('./graphs/Dvitiya_sphuta_comparison.png')
        plt.show()
        plt.clf()

#####---------------------------- Modifying Evection ---------------------------#####
if False:
    for f in [8, 10, 12]:
        dvitiyasphutas_mandamoon_mandasun = []
        dvitiyasphutas_mean = []
        dvitiyasphutas_meanmoon_meansun = []

        evections_mandamoon_mandasun = []
        evections_mandamoon_meansun = []
        evections_meanmoon_meansun = []

        mandasuns = []
        madhyasuns = []
        for mon in range(num_months):
            for i in range(steps):
                ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps

                lm_manda = moons_manda[steps*mon+i]
                lu = swemands_si[steps*mon+i]
                ls_manda = sun_mandasphuta(ag)
                
                mandaKbyR = moon_mandaKbyRs[steps*mon+i]
                r = (moonR/2)*cos(deg2rad(ls_manda-lu))
                theta_ds = sid.theta_ms(lm_manda, ls_manda, r/(f*moonR*mandaKbyR))
                dvitiyasphutas_mandamoon_mandasun.append(theta_ds)
                evections_mandamoon_mandasun.append((lm_manda-theta_ds+180)%360-180)

        plt.plot(range(steps*num_months), subcirc(dvitiyasphutas_mandamoon_mandasun, moons_swe_si), '-', label='Dvitiyasphuta, True Sun, Rem = %f*R'%f)

    plt.ylabel('Difference between observed and Dvitiyasphuta sidereal longitudes of apogee (degrees)')
    plt.xlabel('Time since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.legend()
    plt.gcf().set_size_inches(8, 4)
    plt.savefig('./graphs/Dvitiya_sphuta_comparison_varying_R.png')
    # plt.show()
    plt.clf()

#####-------------------------- Verifying Dvitiyocca --------------------------#####
print('\n-------------------------- Verifying Dvitiyocca --------------------------Skipped!\n')
if False:
    def moon_mandasphuta_old(_ag):
        lu = sid.theta_0(_ag, moon_apogee_revs, mandocca_dhruva)
        return sid.theta_ms(sid.theta_0(_ag, moon_revs, moon_dhruva),  lu, 7/80)

    def moon_dvitiyasphuta(_ag):
        _f=10
        lu = sid.theta_0(_ag, moon_apogee_revs, mandocca_dhruva)

        lm = sid.theta_ms(sid.theta_0(_ag, moon_revs, moon_dhruva),  lu, 7/80)
        mandaKbyR = sid.KbyR(sid.theta_0(_ag, moon_revs, moon_dhruva), lu, 7/80)
        
        ls = sun_mandasphuta(_ag)

        r = (moonR/2)*cos(deg2rad(ls-lu))
        return sid.theta_ms(lm, ls, r/(_f*moonR*mandaKbyR))

    moons_manda = []
    for mon in range(num_months):
        for i in range(steps):
            ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
            moons_manda.append(moon_mandasphuta_old(ag))

    manda_minus_swe_zero_ags = get_zeroes(subcirc(moons_manda, moons_swe_si))
    total_ags = [ag_start+i*(civ_days/moon_revs)/steps for i in range(steps*num_months)]

    moon_shift_180_closest = []
    for i in range(len(manda_minus_swe_zero_ags)):
        ag = manda_minus_swe_zero_ags[i]
        if i==0:
            moon_shift_180_closest.append(moon_mandasphuta_old(ag))
        else:
            r1 = abs((moon_mandasphuta_old(ag) - moon_shift_180_closest[i-1])%360)
            r2 = abs((moon_mandasphuta_old(ag) +180)%360 - moon_shift_180_closest[i-1])
            if(r1<r2): moon_shift_180_closest.append(moon_mandasphuta_old(ag))
            else: moon_shift_180_closest.append((moon_mandasphuta_old(ag) +180)%360)

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,9), gridspec_kw={'height_ratios': [3, 1]})

    ax2.plot(total_ags, subcirc(moons_manda, moons_swe_si), '-', label='theta_ms - theta_swe')
    ax2.plot(manda_minus_swe_zero_ags, [0 for ag in manda_minus_swe_zero_ags], '.', label = 'zeroes computed')
    ax2.plot(total_ags, [(moon_mandasphuta_old(ag)-moon_dvitiyasphuta(ag)+180)%360-180 for ag in total_ags], '-', label = 'Dvitiyaphala')
    ax2.legend()

    ax1.plot(manda_minus_swe_zero_ags, [moon_mandasphuta_old(ag) for ag in manda_minus_swe_zero_ags], '.', label = 'Mandasphuta Moon at zeroes computed')
    # plt.plot(manda_minus_swe_zero_ags, moon_shift_180_closest, '.', label = 'Mean Moon at zeroes computed; 180 added arbitrarily')
    # plt.plot(manda_minus_swe_zero_ags, [(moon_mandasphuta(ag)+180)%360 for ag in manda_minus_swe_zero_ags], '.', label = 'Mandasphuta Moon + 180 at zeroes computed')
    # plt.plot(total_ags, [moon_mandasphuta(ag) for ag in total_ags], '-', label = 'Mandasphuta Moon at every point')
    ax1.plot(total_ags, [sun_mandasphuta(ag) for ag in total_ags], '-', label = 'Actual Dvitiyocca')
    ax1.legend()

    plt.ylabel('Difference between observed and madhya sidereal longitudes of Moon (degrees)')
    plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.savefig('./graphs/Manda_comparison.png')
    # plt.show()
    plt.clf()

#####--------------------------- Gravitational Model -------------------------#####
print('\n--------------------------- Gravitational Model -------------------------\n')
if False:
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(12,6), gridspec_kw={'height_ratios': [1, 1, 1, 1]})
    total_ags = [ag_start+i*(civ_days/moon_revs)/steps for i in range(steps*num_months)]
    for x in [1]:
        moon_dhruva = 1
        moons_eqn_of_center = []
        moons_evection = []
        moons_evection_prime = []
        moons_variation = []
        moons_annual_equation = []
        moons_parallax = []
        moons_reduced_to_ecliptic = []
        for mon in range(num_months):
            for i in range(steps):
                ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps

                lm_mean = sid.theta_0(ag, moon_revs, moon_dhruva)
                lu = sid.theta_0(ag, moon_apogee_revs, mandocca_dhruva)
                ln = sid.theta_0(ag, -232100, moon_node_at_lunar_mesadi)

                lm_eqn_of_center = + (22639/3600) * sin(deg2rad(lm_mean-lu+180)) + x * (769/3600) * sin(2*deg2rad(lm_mean-lu+180)) #+ (36/3600) * sin(3*deg2rad(lm_mean-lu+180))
                
                lu_dvitiya = sun_mandasphuta(ag)
                
                # evection_prime = + (4586/3600) * sin(deg2rad(2*(lm_mean-lu_dvitiya) - (lm_mean-lu+180))) # Original Keplerian form

                evection = + (4586/3600) * ( 2* sin(deg2rad(lm_mean-lu_dvitiya)) * cos(deg2rad(lu_dvitiya-lu+180)) - sin(deg2rad(lm_mean-lu+180)))

                variation = + (2370/3600) * sin(deg2rad(2*(lm_mean-lu_dvitiya)))

                parallax = 0* - (125/3600) * sin(deg2rad(lm_mean-lu_dvitiya))

                annual_equation = -(668/3600)*sin(deg2rad(lu_dvitiya-sun_mandocca_mean+180))
                
                reduction_to_ecliptic = 0 *- (412/3600) * sin(2*deg2rad(lm_mean-ln))

                moons_eqn_of_center.append(lm_mean + lm_eqn_of_center)
                moons_evection.append(lm_mean + lm_eqn_of_center + evection)
                # moons_evection_prime.append(lm_mean + lm_eqn_of_center + evection_prime)
                moons_variation.append(lm_mean + lm_eqn_of_center + evection + variation)
                moons_annual_equation.append(lm_mean + lm_eqn_of_center + evection + variation + annual_equation)
                # moons_parallax.append(lm_mean + lm_eqn_of_center + evection + variation + annual_equation + parallax)
                # moons_reduced_to_ecliptic.append(lm_mean + lm_eqn_of_center + evection + variation + annual_equation + parallax + reduction_to_ecliptic)

        ax1.plot(total_ags, subcirc(moons_eqn_of_center, moons_swe_si), '-', label='Equation of Center')
        ax1.set_ylim([-3,3])
        ax1.set_ylabel('Degrees')
        ax1.legend()
        ax2.plot(total_ags, subcirc(moons_evection, moons_swe_si), '-', label='Evection')
        ax2.set_ylim([-1,1])
        ax2.set_ylabel('Degrees')
        ax2.legend()
        # ax2.plot(total_ags, subcirc(moons_evection_prime, moons_swe_si), '-', label='Evection_prime')
        # ax2.set_ylim([-1,1])
        # ax2.legend()
        ax3.plot(total_ags, subcirc_mins(moons_variation, moons_swe_si), '-', label='Variation')
        ax3.set_ylim([-25,25])
        ax3.set_ylabel('Minutes')
        ax3.legend()
        ax4.plot(total_ags, subcirc_mins(moons_annual_equation, moons_swe_si), '-', label='Annual Equation')
        ax4.set_ylim([-25,25])
        ax4.set_ylabel('Minutes')
        ax4.legend()
        # ax5.plot(range(steps*num_months), subcirc(moons_parallax, moons_swe_si), '-', label='Parallax x = %d'%x)
        # ax5.set_ylim([-0.3,0.3])
        # ax5.legend()
        # ax6.plot(range(steps*num_months), subcirc(moons_reduced_to_ecliptic, moons_swe_si), '-', label='Reduction x = %d'%x)
        # ax6.set_ylim([-0.25,0.25])
        # ax6.legend()

    # plt.ylabel('Difference between observed and Dvitiyasphuta sidereal longitudes of apogee (degrees)')
    ax1.set_title('Newtonian Gravitational Model')
    plt.xlabel('Time since Lunar Mesadi (days)')
    plt.gcf().set_size_inches(6, 7)
    plt.savefig('./graphs/gravitational_model.png')
    plt.show()
    plt.clf()

#####--------------------------- Modifying manda r0 -------------------------#####
print('\n--------------------------- Modifying manda r0 -------------------------\n')

def r0_manda(_anomaly, _p, _q):
    return _p-_q*cos(deg2rad(_anomaly))

if False:
    fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(16,9), gridspec_kw={'height_ratios': [1, 1, 1]})
    for q in [0.40, 0.45, 0.50]: #0.5965147532 # 0.01149008424 # 0.4446711083
        moons_manda = []
        dvitiyasphutas_mandamoon_mandasun = []
        moons_trtiya = []
        for mon in range(num_months):
            for i in range(steps):
                ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
                # print(ag, moons_swe_si[i])

                lm_mean = sid.theta_0(ag, moon_revs, moon_dhruva)
                K_mean = 80

                lu = sid.theta_0(ag, moon_apogee_revs, mandocca_dhruva)
                r_manda=r0_manda(lm_mean-lu, 7, 0)
                lm_manda = sid.theta_ms(lm_mean, lu, r_manda/K_mean)
                K_manda = 10*moonR*sid.KbyR(lm_mean, lu, r_manda/K_mean)

                lu_dvitiya = sun_mandasphuta(ag)
                r_dvitiya = (moonR*q)*cos(deg2rad(lu_dvitiya-lu))
                lm_dvitiya = sid.theta_ms(lm_manda, lu_dvitiya, r_dvitiya/K_manda) # 
                K_dvitiya = sid.KbyR(lm_manda, lu_dvitiya, r_dvitiya/K_manda)

                lu_trtiya = 2*lu_dvitiya-lm_dvitiya
                r_trtiya = -0.01149008424*K_dvitiya
                lm_trtiya = sid.theta_ms(lm_dvitiya, lu_trtiya, r_trtiya/K_dvitiya)
                
                moons_manda.append(lm_manda)
                dvitiyasphutas_mandamoon_mandasun.append(lm_dvitiya)
                moons_trtiya.append(lm_trtiya)

        ax1.plot(range(steps*num_months), subcirc(moons_manda, moons_swe_si), '-', label='Mandasphuta, q = %f'%q)
        ax2.plot(range(steps*num_months), subcirc(dvitiyasphutas_mandamoon_mandasun, moons_swe_si), '-', label='Dvitiyasphuta, 0.4446')
        ax3.plot(range(steps*num_months), subcirc(moons_trtiya, moons_swe_si), '-', label='Trtiyasphuta, 0.01149')

    plt.ylabel('Difference between observed and Dvitiyasphuta sidereal longitudes of apogee (degrees)')
    plt.xlabel('Time since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    ax1.legend()
    ax2.legend()
    ax3.legend()
    plt.gcf().set_size_inches(8, 4)
    plt.savefig('./graphs/Manda_r0_variation.png')
    plt.show()
    plt.clf()

#####-------------------------- Fitting Parameters --------------------------#####
print('\n#####-------------------------- Fitting Parameters --------------------------#####\n')
with shelve.open('../common_data/epochal_values') as epochfile:
    sun_dhruva = epochfile['sun_dhruva']
    jdtt_recent_epoch_audayika_minus_ag_recent_epoch = epochfile['jdtt_recent_epoch_audayika_minus_ag_recent_epoch']
ag_start = jdtt_moon_curvefit_start-jdtt_recent_epoch_audayika_minus_ag_recent_epoch#1871998.8010278898
print(ag_start)
def sid_swe_diff(_inputs, _p, _q, _r, _moon_dhruva, _k, _mandocca_dhruva):
    i, _swelong = _inputs
    _ag = ag_start+i*(SIDE_MT)/steps
    # print(_ag, moons_swe_si[i])

    lm_mean = sid.theta_0(_ag, moon_revs_new, _moon_dhruva, civ_days_new)
    K_mean = 80

    lu = sid.theta_0(_ag, moon_apogee_revs, _mandocca_dhruva, civ_days_new)
    r_manda=r0_manda(lm_mean-lu, 7, _p)
    lm_manda = sid.theta_ms(lm_mean, lu, r_manda/K_mean)
    K_manda = 10*moonR*sid.KbyR(lm_mean, lu, r_manda/K_mean)

    lu_dvitiya = sun_mandasphuta(_ag)
    r_dvitiya = (moonR*_q)*cos(deg2rad(lu_dvitiya-lu))
    lm_dvitiya = sid.theta_ms(lm_manda, lu_dvitiya, r_dvitiya/K_manda) # 
    K_dvitiya = sid.KbyR(lm_manda, lu_dvitiya, r_dvitiya/K_manda)

    lu_trtiya = 2*lu_dvitiya-lm_dvitiya
    r_trtiya = -_r*K_dvitiya
    lm_trtiya = sid.theta_ms(lm_dvitiya, lu_trtiya, r_trtiya/K_dvitiya)

    mean_sun = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
    lm_annual_equation = _k* ((mean_sun-sun_mandasphuta(_ag)+180)%360-180)

    return 60*(((lm_trtiya + lm_annual_equation-_swelong+180)%360)-180)

total_ags = [ag_start+i*(civ_days/moon_revs)/steps for i in range(steps*num_months)]
p0 = [0.5965147532/2, 0.5, 0.01149008424, -7.39, 0.1, 117]
# p0 = [0.273835802, 0.457860009, 0.00956421378, 1, 0.0926511722,]
plt.plot(total_ags, [sid_swe_diff((i, moons_swe_si[i]), *p0) for i in range(len(total_ags))], '-', label='Initial Guess: %s'%p0)
moon_popt, moon_pcov = curve_fit(sid_swe_diff, (range(len(total_ags)), moons_swe_si), np.zeros(len(total_ags)), p0=p0)
plt.plot(total_ags, [sid_swe_diff((i, moons_swe_si[i]), *moon_popt) for i in range(len(total_ags))], '-', label='Fitted Curve: %s'%moon_popt)
plt.ylabel('$\\theta_{{sid}}-\\theta_{{jpl}}$ (minutes)')
plt.xlabel('Time since Lunar Aries Transit (days)')
plt.legend()
plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/New_Siddhanta_Parameters_2024.png')
plt.show()
plt.clf()

#####-------------------------- Plot Revised Model --------------------------#####
print('\n#####-------------------------- Plot Revised Model --------------------------#####\n')

if False:
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(16,9), gridspec_kw={'height_ratios': [1, 1, 1, 1]})

    moons_manda = []
    moons_tungantara = []
    moons_paksika = []
    moons_digamsa = []
    total_ags = [ag_start+i*(civ_days/moon_revs)/steps for i in range(steps*num_months)]
    for mon in range(num_months):
        for i in range(steps):
            _ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
            [_p, _q, _r, _moon_dhruva, _k, _mandocca_dhruva] = moon_popt
            lm_mean = sid.theta_0(_ag, moon_revs_new, _moon_dhruva, civ_days_new)
            K_mean = 80

            lu = sid.theta_0(_ag, moon_apogee_revs, _mandocca_dhruva, civ_days_new)
            r_manda=r0_manda(lm_mean-lu, 7, _p)
            lm_manda = sid.theta_ms(lm_mean, lu, r_manda/K_mean)
            K_manda = 10*moonR*sid.KbyR(lm_mean, lu, r_manda/K_mean)
            moons_manda.append(lm_manda)

            lu_dvitiya = sun_mandasphuta(_ag)
            r_dvitiya = (moonR*_q)*cos(deg2rad(lu_dvitiya-lu))
            lm_dvitiya = sid.theta_ms(lm_manda, lu_dvitiya, r_dvitiya/K_manda) # 
            K_dvitiya = sid.KbyR(lm_manda, lu_dvitiya, r_dvitiya/K_manda)
            moons_tungantara.append(lm_dvitiya)

            lu_trtiya = 2*lu_dvitiya-lm_dvitiya
            r_trtiya = -_r*K_dvitiya
            lm_trtiya = sid.theta_ms(lm_dvitiya, lu_trtiya, r_trtiya/K_dvitiya)
            moons_paksika.append(lm_trtiya)

            mean_sun = sid.theta_0(_ag, sun_revs, sun_dhruva, civ_days_new)
            lm_annual_equation = _k* ((mean_sun-sun_mandasphuta(_ag)+180)%360-180)
            moons_digamsa.append(lm_trtiya+lm_annual_equation)

    ax1.plot(total_ags, subcirc(moons_manda, moons_swe_si), '-', label='Manda')
    ax1.set_ylim([-3,3])
    ax1.set_ylabel('Degrees')
    ax1.legend()
    ax2.plot(total_ags, subcirc(moons_tungantara, moons_swe_si), '-', label='Tungantara')
    ax2.set_ylim([-1,1])
    ax2.set_ylabel('Degrees')
    ax2.legend()
    ax3.plot(total_ags, subcirc_mins(moons_paksika, moons_swe_si), '-', label='Paksika')
    ax3.set_ylim([-25,25])
    ax3.set_ylabel('Minutes')
    ax3.legend()
    ax4.plot(total_ags, subcirc_mins(moons_digamsa, moons_swe_si), '-', label='Digamsa')
    ax4.set_ylim([-25,25])
    ax4.set_ylabel('Minutes')
    ax4.legend()

    ax1.set_title('Revised Siddhantaic Model')
    plt.xlabel('Time since Lunar Mesadi (days)')
    plt.gcf().set_size_inches(6, 7)
    plt.savefig('./graphs/revised_siddhantaic_model.png')
    plt.show()
    plt.clf()

ag = 1871998.8010278898
# print(ag, sid_moon(ag, *moon_popt))

#####-------------------------- Writing Parameters --------------------------#####
print('\n#####-------------------------- Writing Parameters --------------------------#####\n')

import shelve
moon_revised_params = shelve.open('./moon_computed_data/moon_revised_parameters')

moon_revised_params['moon_apogee_revs'] = moon_apogee_revs

moon_revised_params['moon_revs'] = moon_revs
moon_revised_params['moon_dhruva'] = moon_popt[3]
moon_revised_params['moon_mandocca_dhruva'] = moon_popt[5]

moon_revised_params['manda_radius_oscillation'] = moon_popt[0]
moon_revised_params['dvitiya_radius'] = moon_popt[1]
moon_revised_params['trtiya_radius'] = moon_popt[2]
moon_revised_params['annual_equation_factor'] = moon_popt[4]

moon_revised_params.close()
#####-------------------------- Introducing Trtiya --------------------------#####
# print('\n-------------------------- Introducing Trtiya --------------------------Skipped!\n')
if False:
    moon_ds_errors = subcirc(dvitiyasphutas_mandamoon_mandasun, moons_swe_si)
    print('avg of dvitiyasphutas', stat.mean(moon_ds_errors))
    moon_dhruva = -stat.mean(moon_ds_errors)

    def trtiyaphala(_dvitiyasphuta, _ls, _rbyK):
        return rad2deg(arcsin(_rbyK*sin(2*deg2rad(_ls-_dvitiyasphuta))))

    dvitiyasphutas_mandamoon_mandasun = []
    for mon in range(num_months):
        for i in range(steps):
            ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
            dvitiyasphutas_mandamoon_mandasun.append(moon_dvitiyasphuta(ag))
    moon_ds_errors = subcirc(dvitiyasphutas_mandamoon_mandasun, moons_swe_si)

    dvitiyasphuta_minus_swe_zero_ags = get_zeroes(moon_ds_errors)
    total_ags = [ag_start+i*(civ_days/moon_revs)/steps for i in range(steps*num_months)]

    moon_alt_180 = []
    for i in range(len(dvitiyasphuta_minus_swe_zero_ags)):
        ag = dvitiyasphuta_minus_swe_zero_ags[i]
        if i%2 == 0:
            moon_alt_180.append((moon_dvitiyasphuta(ag)+180)%360)
        else:
            moon_alt_180.append((moon_dvitiyasphuta(ag)))

    moon_shift_180_closest = []
    for i in range(len(dvitiyasphuta_minus_swe_zero_ags)):
        ag = dvitiyasphuta_minus_swe_zero_ags[i]
        if i==0:
            moon_shift_180_closest.append(moon_dvitiyasphuta(ag))
        else:
            r1 = abs((moon_dvitiyasphuta(ag) - moon_shift_180_closest[i-1])%360)
            r2 = abs((moon_dvitiyasphuta(ag) + 180)%360 - moon_shift_180_closest[i-1])
            if(r1<r2): moon_shift_180_closest.append(moon_dvitiyasphuta(ag))
            else: moon_shift_180_closest.append((moon_dvitiyasphuta(ag) + 180)%360)

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,9), gridspec_kw={'height_ratios': [3, 1]})

    ax2.plot(total_ags, moon_ds_errors, '-', label='theta_ds - theta_swe')
    ax2.plot(dvitiyasphuta_minus_swe_zero_ags, [0 for ag in dvitiyasphuta_minus_swe_zero_ags], '.', label = 'zeroes computed')

    # for fr in [0.001,0.005,0.01]:
    #     trtiyaphalas = []
    #     for mon in range(num_months):
    #         for i in range(steps):
    #             ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
    #             ls = sun_mandasphuta(ag)
    #             lm = moon_dvitiyasphuta(ag)
    #             trtiyaphalas.append(trtiyaphala(lm, ls, fr))
    #     ax2.plot(total_ags, trtiyaphalas, '-', label='trtiyaphala %f'%fr)
    ax2.legend()

    ax1.plot(dvitiyasphuta_minus_swe_zero_ags, [moon_dvitiyasphuta(ag) for ag in dvitiyasphuta_minus_swe_zero_ags], '.--', label = 'Dvitiyasphuta Moon at zeroes computed')
    # plt.plot(total_ags, [moon_dvitiyasphuta(ag) for ag in total_ags], '-', label = 'Dvitiyasphuta Moon at every point')
    # plt.plot(dvitiyasphuta_minus_swe_zero_ags, moon_alt_180, '.', label = 'Dvitiyasphuta Moon at zeroes computed; 180 added at alernate places')
    # plt.plot(dvitiyasphuta_minus_swe_zero_ags, moon_shift_180_closest, '.', label = 'Dvitiyasphuta Moon at zeroes computed; 180 added arbitrarily')

    # ax1.plot(dvitiyasphuta_minus_swe_zero_ags, [(moon_dvitiyasphuta(ag)+180)%360 for ag in dvitiyasphuta_minus_swe_zero_ags], '.', label = 'Dvitiyasphuta Moon + 180 at zeroes computed')
    ax1.plot(total_ags, [sun_mandasphuta(ag) for ag in total_ags], '-', label = 'Sun')
    ax1.plot(total_ags, [(sun_mandasphuta(ag)+90)%360 for ag in total_ags], '-', label = 'Sun+90')
    ax1.plot(total_ags, [(sun_mandasphuta(ag)+180)%360 for ag in total_ags], '-', label = 'Sun+180')
    ax1.plot(total_ags, [(sun_mandasphuta(ag)+270)%360 for ag in total_ags], '-', label = 'Sun+270')

    # plt.plot(dvitiyasphuta_minus_swe_zero_ags, [(2*sun_mandasphuta(ag))%360 for ag in dvitiyasphuta_minus_swe_zero_ags], 'o', label = '2*Sun')
    # plt.plot(total_ags, [(moons_dvitiyasphuta(ag)-moon_trtiyasphtua(ag)+180)%360-180 for ag in total_ags], '-', label = 'Trtiyaphala')
    ax1.legend()

    plt.ylabel('Difference between sid ans swe sidereal longitudes of Moon (degrees)')
    plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.savefig('./graphs/Identifying_trtiyocca.png')
    # plt.show()
    plt.clf()

    def get_r0byR(_theta_prev, _theta_true, _theta_ucca):
        if sin(deg2rad(_theta_prev - _theta_ucca)) == 0: return 0
        return sin(-deg2rad(_theta_true - _theta_prev))/sin(deg2rad(_theta_prev - _theta_ucca))

    def phala_model(_inputs, _p, _q, _a, _b):
        _theta_prev, _ag = _inputs
        _theta_ucca = 2*sun_mandasphuta(_ag)-_theta_prev
        _kendra = _theta_prev-_theta_ucca
        _rbyK = _a
        return rad2deg(arcsin(_rbyK*sin(deg2rad(_kendra))))

    fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(16,9), gridspec_kw={'height_ratios': [3, 1, 1]})
    p0 = [2, # _p
        0.025, # _q
        -0.01, # _a
        -0.001 #_b
        ]
    ax1.plot(total_ags, [phala_model((dvitiyasphutas_mandamoon_mandasun[i], total_ags[i]), *p0) for i in range(len(total_ags))], '--', label='Initial Guess: %s'%(p0))
    ax2.plot(total_ags, [p0[2]+p0[3]*cos(p0[0]*deg2rad((sun_mandasphuta(_ag)))) for _ag in total_ags], '--', label='rbyK: %s'%(p0))
    ax3.plot(total_ags, [cos(p0[0]*deg2rad((sun_mandasphuta(_ag)))) for _ag in total_ags], '--', label='rbyK argument: %s'%(p0))
    # trtiyaphala_popt, trtiyaphala_pcov = curve_fit(phala_model, (moons_dvitiya, total_ags), np.array(moon_ds_errors), p0=p0)
    # plt.plot(total_ags, [phala_model((moons_dvitiya[i], total_ags[i]), *trtiyaphala_popt) for i in range(len(total_ags))], '-', label='Fitted Curve: %s'%trtiyaphala_popt)
    ax1.plot(total_ags, moon_ds_errors, '-', label='theta_ds - theta_swe')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    plt.ylabel('Difference between sid and swe sidereal longitudes of Moon (degrees)')
    plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.savefig('./graphs/Trtiya_modelling.png')
    plt.show()
    plt.clf()

    plt.plot(total_ags, subcirc([dvitiyasphutas_mandamoon_mandasun[i] - phala_model((dvitiyasphutas_mandamoon_mandasun[i], total_ags[i]), *p0) for i in range(len(total_ags))], moons_swe_si), '-', label='Trtiyaphuta - swe with trtiyaparams %s'%p0)

    plt.legend()
    plt.ylim(-0.5,0.5)
    plt.ylabel('Difference between sid and swe sidereal longitudes of Moon (degrees)')
    plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.savefig('./graphs/Trtiya_modelling_result.png')
    plt.show()
    plt.clf()

    print('\nWARNING: DATA TAKEN FROM FILE. `make generate_data` to recompile generate_data.py')