import numpy as np
from numpy import sqrt, sin, cos, arcsin, pi, rad2deg, deg2rad
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import statistics as stat

import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
from modern_astro_params import ANOM_MT, SIDE_MT, subcirc
# from generate_data import steps, num_months, swelongs_an, swelongs_si,\
# mandocca_true_apogee, sun_lunar_mesadi, swemands_si, swesunls_si
import siddhantaic_functions as sid
import sun.functions as sun
import moon.functions as moon

#####------------------------- Import computed data --------------------------#####
print('------------------------- Import computed data --------------------------')
import shelve
with shelve.open('../sun/sun_computed_data/sun_modified_parameters') as sunp:
    sun_civ_days = sunp['sun_civ_days']
    sun_mandocca_dhruva = sunp['mandocca_dhruva']
    sun_mandocca_revs = sunp['mandocca_revs']
    sun_r0byR_popt = sunp['r0byR_popt']
with shelve.open('../moon/moon_computed_data/generated_data') as gen_dat:
    steps = gen_dat['steps']
    num_months = gen_dat['num_months']
    swelongs_an = gen_dat['swelongs_an']
    moons_swe_si = gen_dat['swelongs_si']
    mandocca_true_apogee = gen_dat['mandocca_true_apogee']
    sun_lunar_mesadi = gen_dat['sun_lunar_mesadi']
    swemands_si = gen_dat['swemands_si']
    swesunls_si = gen_dat['swesunls_si']
    
    moon_solar_mesadi = gen_dat['moon_solar_mesadi']
    moon_apogee_solar_mesadi = gen_dat['moon_apogee_solar_mesadi']
    SIDE_YR = gen_dat['SIDE_YR']

#####------------------------------ Plot Settings ----------------------------#####
plt.rcParams['axes.grid'] = True

#####------------------- Start Revision of Moon parameters -------------------#####
print('\n------------------- Start Revision of Moon parameters -------------------\n')
print('\nWARNING: USING DATA FROM FILE. Recompile generate_data.py to get latest changes\n')
#####---------------- Moon's treatment according to Siddhanta ----------------#####

# Mandasamskara
from functions import moon_revs, civ_days,\
    moon_apogee_revs, r0byR, R


#####----- Compute new parameters by comparing in one anomalistic cycle -----#####
anom_rev_diff = civ_days/ANOM_MT
#####-------------------- Effect of changing month length -------------------#####
sidlongs = [sid.theta_ms(sid.theta_0(ag, moon_revs - moon_apogee_revs), 0, r0byR) for ag in np.arange(0, ANOM_MT, ANOM_MT/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='moon_revs = %f (KAPA)'%(moon_revs - moon_apogee_revs))
print('moon_revs - moon_apogee_revs: old = %f'%(moon_revs - moon_apogee_revs))

sidlongs = [sid.theta_ms(sid.theta_0(ag,anom_rev_diff), 0, r0byR) for ag in np.arange(0, ANOM_MT, ANOM_MT/steps)]
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

# moon_mandasphuta = lambda _ag: sid.theta_ms(sid.theta_0(_ag, moon_revs),  sid.theta_0(_ag, moon_apogee_revs, mandocca_true_apogee), sid.r0_revised(sid.theta_0(_ag, moon_revs), sid.theta_0(_ag, moon_apogee_revs, mandocca_true_apogee), *r0byR_popt)/80)
moon_dhruva = 0
mandocca_dhruva = mandocca_true_apogee
moon_mandasphuta = lambda _ag: sid.theta_ms(sid.theta_0(_ag, moon_revs, moon_dhruva),  sid.theta_0(_ag, moon_apogee_revs, mandocca_dhruva), 7/80)

ag_start = sid.get_longitude_ag(moon_mandasphuta, 0, 0, (SIDE_MT/360)/5)
mandocca_dhruva = sid.theta_0(-ag_start, moon_apogee_revs, mandocca_true_apogee)

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
print('---------------------- Ucca Analysis with madhyagraha ---------------------')

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

moons_madhya = []
moons_manda = []
moon_mandaKbyRs = []
for mon in range(num_months):
    for i in range(steps):
        moons_madhya.append(sid.theta_0(ag, moon_revs, 0))
        moons_manda.append(moon_mandasphuta(ag))
        lu = moon_apogees_madhya[steps*mon+i]
        moon_mandaKbyRs.append(sid.iterated_KbyR(sid.theta_0(ag, moon_revs), lu, 7/80))

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

Rem = 10*moon.R
moons_dvitiya = []
evections = []

sun_mandasphuta_func = lambda _ag: sid.theta_ms(sid.theta_0(_ag, civ_days/SIDE_YR),  sun_mandocca_dhruva, sid.r0_revised(sid.theta_0(_ag, civ_days/SIDE_YR), sun_mandocca_dhruva, *sun_r0byR_popt)/360)
moon_mandocca = lambda _ag: sid.theta_0(_ag, moon_apogee_revs, 0)

sun_dhruva = sid.theta_0(sid.get_longitude_ag(sun_mandasphuta_func,sun_lunar_mesadi, 0, 1)-ag_start,civ_days/SIDE_YR,0)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
ax1.plot(range(steps*num_months), subcirc(moons_manda, moons_swe_si), '-', label='mandasphuta Longitude Difference')
ax1.legend()

def sun_mandasphuta(_ag):
    t0 = sid.theta_0(_ag, civ_days/SIDE_YR, sun_dhruva)
    return sid.theta_ms(t0,  sun_mandocca_dhruva, sid.r0_revised(t0, sun_mandocca_dhruva, *sun_r0byR_popt)/360)


moons_dvitiya = []
dvitiyasphutas_mean = []
dvitiyasphutas_manda = []
dvitiyasphutas_mean_mandocca = []

evections = []
evections_mean = []
evections_manda = []
evections_mean_mandocca = []

mandasuns = []
madhyasuns = []

f=10
for mon in range(num_months):
    for i in range(steps):
        ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps

        lm = moons_manda[steps*mon+i]
        lu = swemands_si[steps*mon+i]
        ls = swesunls_si[steps*mon+i]
        mandaKbyR = moon_mandaKbyRs[steps*mon+i]

        r = (moon.R/2)*cos(deg2rad(ls-lu))
        theta_ds = sid.theta_ms(lm, ls, r/(f*moon.R*mandaKbyR))
        moons_dvitiya.append(theta_ds)
        evections.append((lm-theta_ds+180)%360-180)

        ls = sid.theta_0(ag, civ_days/SIDE_YR, sun_dhruva)
        theta_ds = sid.theta_ms(lm, ls, r/(f*moon.R*mandaKbyR))
        dvitiyasphutas_mean.append(theta_ds)
        evections_mean.append((lm-theta_ds+180)%360-180)

        ls = sun_mandasphuta(ag)
        theta_ds = sid.theta_ms(lm, ls, r/(f*moon.R*mandaKbyR))
        dvitiyasphutas_manda.append(theta_ds)
        evections_manda.append((lm-theta_ds+180)%360-180)

        lu = moon_apogees_madhya[steps*mon+i]
        theta_ds = sid.theta_ms(lm, ls, r/(f*moon.R*mandaKbyR))
        dvitiyasphutas_mean_mandocca.append(theta_ds)
        evections_mean_mandocca.append((lm-theta_ds+180)%360-180)


ax2.plot(range(steps*num_months), evections, '-', label='Evections with True Sun Rem = %f*R'%f)
ax2.plot(range(steps*num_months), evections_mean, '-', label='Evections with Mean Sun, Rem = %f*R'%f)
ax2.plot(range(steps*num_months), evections_manda, '-', label='Evections with Manda Sun, Rem = %f*R'%f)
ax2.plot(range(steps*num_months), evections_mean_mandocca, '-', label='Evections with Mean mandocca, Rem = %f*R'%f)
ax2.legend()

ax3.plot(range(steps*num_months), subcirc(moons_dvitiya, moons_swe_si), '-', label='Dvitiyasphuta, True Sun')
ax3.plot(range(steps*num_months), subcirc(dvitiyasphutas_mean, moons_swe_si), '-', label='Dvitiyasphuta, Mean Sun')
ax3.plot(range(steps*num_months), subcirc(dvitiyasphutas_manda, moons_swe_si), '-', label='Dvitiyasphuta, Manda Sun')
ax3.plot(range(steps*num_months), subcirc(dvitiyasphutas_mean_mandocca, moons_swe_si), '-', label='Dvitiyasphuta, Mean mandocca')
ax3.legend()

plt.ylabel('Difference between observed and swe sidereal longitudes (degrees)')
plt.xlabel('Time since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
plt.gcf().set_size_inches(18, 10)
plt.savefig('./graphs/Dvitiya_sphuta_comparison.png')
# plt.show()
plt.clf()

#####---------------------------- Modifying Evection ---------------------------#####
if False:
    for f in [8, 10, 12]:
        moons_dvitiya = []
        dvitiyasphutas_mean = []
        dvitiyasphutas_manda = []

        evections = []
        evections_mean = []
        evections_manda = []

        mandasuns = []
        madhyasuns = []
        for mon in range(num_months):
            for i in range(steps):
                ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps

                lm = moons_manda[steps*mon+i]
                lu = swemands_si[steps*mon+i]
                ls = sun_mandasphuta(ag)
                
                mandaKbyR = moon_mandaKbyRs[steps*mon+i]
                r = (moon.R/2)*cos(deg2rad(ls-lu))
                theta_ds = sid.theta_ms(lm, ls, r/(f*moon.R*mandaKbyR))
                moons_dvitiya.append(theta_ds)
                evections.append((lm-theta_ds+180)%360-180)

        plt.plot(range(steps*num_months), subcirc(moons_dvitiya, moons_swe_si), '-', label='Dvitiyasphuta, True Sun, Rem = %f*R'%f)

    plt.ylabel('Difference between observed and Dvitiyasphuta sidereal longitudes of apogee (degrees)')
    plt.xlabel('Time since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
    plt.legend()
    plt.gcf().set_size_inches(18, 10)
    plt.savefig('./graphs/Dvitiya_sphuta_comparison_varying_R.png')
    # plt.show()
    plt.clf()

#####-------------------------- Wrapping up Evection --------------------------#####
print('\n-------------------------- Wrapping up Evection --------------------------\n')

def moon_mandasphuta(_ag):
    lu = sid.theta_0(_ag, moon_apogee_revs, mandocca_dhruva)
    return sid.theta_ms(sid.theta_0(_ag, moon_revs, moon_dhruva),  lu, 7/80)

def moon_dvitiyasphuta(_ag):
    _f=10
    lu = sid.theta_0(_ag, moon_apogee_revs, mandocca_dhruva)

    lm = sid.theta_ms(sid.theta_0(_ag, moon_revs, moon_dhruva),  lu, 7/80)
    mandaKbyR = sid.iterated_KbyR(sid.theta_0(_ag, moon_revs, moon_dhruva), lu, 7/80)
    
    ls = sun_mandasphuta(_ag)

    r = (moon.R/2)*cos(deg2rad(ls-lu))
    return sid.theta_ms(lm, ls, r/(_f*moon.R*mandaKbyR))

moons_manda = []
for mon in range(num_months):
    for i in range(steps):
        ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
        moons_manda.append(moon_mandasphuta(ag))

manda_minus_swe_zero_ags = get_zeroes(subcirc(moons_manda, moons_swe_si))
total_ags = [ag_start+i*(civ_days/moon_revs)/steps for i in range(steps*num_months)]

moon_alt_180 = []
for i in range(len(manda_minus_swe_zero_ags)):
    ag = manda_minus_swe_zero_ags[i]
    if i%2 == 0:
        moon_alt_180.append((moon_mandasphuta(ag)+180)%360)
    else:
        moon_alt_180.append((moon_mandasphuta(ag)))

moon_shift_180_closest = []
for i in range(len(manda_minus_swe_zero_ags)):
    ag = manda_minus_swe_zero_ags[i]
    if i==0:
        moon_shift_180_closest.append(moon_mandasphuta(ag))
    else:
        r1 = abs((moon_mandasphuta(ag) - moon_shift_180_closest[i-1])%360)
        r2 = abs((moon_mandasphuta(ag) +180)%360 - moon_shift_180_closest[i-1])
        if(r1<r2): moon_shift_180_closest.append(moon_mandasphuta(ag))
        else: moon_shift_180_closest.append((moon_mandasphuta(ag) +180)%360)

fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,9), gridspec_kw={'height_ratios': [3, 1]})

ax2.plot(total_ags, subcirc(moons_manda, moons_swe_si), '-', label='theta_ms - theta_swe')
ax2.plot(manda_minus_swe_zero_ags, [0 for ag in manda_minus_swe_zero_ags], '.', label = 'zeroes computed')
ax2.plot(total_ags, [(moon_mandasphuta(ag)-moon_dvitiyasphuta(ag)+180)%360-180 for ag in total_ags], '-', label = 'Dvitiyaphala')
ax2.legend()

ax1.plot(manda_minus_swe_zero_ags, [moon_mandasphuta(ag) for ag in manda_minus_swe_zero_ags], '.', label = 'Mandasphuta Moon at zeroes computed')
# plt.plot(manda_minus_swe_zero_ags, moon_alt_180, '.', label = 'Mean Moon at zeroes computed; 180 added at alernate places')
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

#####-------------------------- Introducing Variation --------------------------#####
print('\n-------------------------- Introducing Variation --------------------------\n')
moon_ds_errors = subcirc(moons_dvitiya, moons_swe_si)
print('avg of dvitiyasphutas', stat.mean(moon_ds_errors))
moon_dhruva = -stat.mean(moon_ds_errors)

def trtiyaphala(_dvitiyasphuta, _ls, _rbyK):
    return rad2deg(arcsin(_rbyK*sin(2*deg2rad(_ls-_dvitiyasphuta))))

moons_dvitiya = []
for mon in range(num_months):
    for i in range(steps):
        ag = ag_start+(civ_days/moon_revs)*mon+i*(civ_days/moon_revs)/steps
        moons_dvitiya.append(moon_dvitiyasphuta(ag))
moon_ds_errors = subcirc(moons_dvitiya, moons_swe_si)

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
plt.savefig('./graphs/Manda_comparison.png')
# plt.show()
plt.clf()

def get_r0byR(_theta_prev, _theta_true, _theta_ucca):
    if sin(deg2rad(_theta_prev - _theta_ucca)) == 0: return 0
    return sin(-deg2rad(_theta_true - _theta_prev))/sin(deg2rad(_theta_prev - _theta_ucca))

def phala_model(_inputs, _p, _q, _a, _b):
    _theta_prev, _ag = _inputs
    _theta_ucca = 2*sun_mandasphuta(_ag)-_theta_prev
    _kendra = _theta_prev-_theta_ucca
    _rbyK_argument = _p*_theta_ucca
    _rbyK = _a+_b*cos(deg2rad(_rbyK_argument))
    return rad2deg(arcsin(_rbyK*sin(deg2rad(_kendra)))) + rad2deg(arcsin(_q*sin(deg2rad(sun_mandasphuta(_ag)+90))))

# def phala_model(_inputs, _k, _a, _b):
#     _theta_prev, _ag = _inputs
#     _theta_ucca = _k * _theta_prev
#     _kendra = _theta_prev-_theta_ucca
#     _rbyK = _a+_b*cos(deg2rad(_kendra))
#     return rad2deg(arcsin(_rbyK*sin(deg2rad(_kendra))))

# p0 = [-4.9e7, -76, 0.009, 0.01]
# plt.plot(total_ags, [phala_model((moons_dvitiya[i], total_ags[i]), *p0) for i in range(len(total_ags))], '-', label='Initial Guess: %s'%p0)
p0 = [1, 0.025, -0.01, -0.001]
plt.plot(total_ags, [phala_model((moons_dvitiya[i], total_ags[i]), *p0) for i in range(len(total_ags))], '--', label='Initial Guess: %s'%(p0))
trtiyaphala_popt, trtiyaphala_pcov = curve_fit(phala_model, (moons_dvitiya, total_ags), np.array(moon_ds_errors), p0=p0)
plt.plot(total_ags, [phala_model((moons_dvitiya[i], total_ags[i]), *trtiyaphala_popt) for i in range(len(total_ags))], '-', label='Fitted Curve: %s'%trtiyaphala_popt)
plt.plot(total_ags, moon_ds_errors, '-', label='theta_ds - theta_swe')

plt.legend()
plt.ylabel('Difference between sid and swe sidereal longitudes of Moon (degrees)')
plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
plt.savefig('./graphs/Trtiya_modelling.png')
plt.show()
plt.clf()

plt.plot(total_ags, subcirc([moons_dvitiya[i] - phala_model((moons_dvitiya[i], total_ags[i]), *trtiyaphala_popt) for i in range(len(total_ags))], moons_swe_si), '-', label='Trtiyaphuta - swe with trtiyaparams %s'%trtiyaphala_popt)

plt.legend()
plt.ylabel('Difference between sid and swe sidereal longitudes of Moon (degrees)')
plt.xlabel('Days since Lunar Mesadi: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
plt.savefig('./graphs/Trtiya_modelling_result.png')
plt.show()
plt.clf()

#####--------------------------- Modifying manda r0 -------------------------#####
print('\n--------------------------- Modifying manda r0 -------------------------\n')

# Get "actual" variation of r0byR(t)

def get_manda_r0byR_moon(ag, long_an):
    if sin(deg2rad(sid.theta_0(ag, anom_rev_diff))) == 0:
        ag=0.25*ANOM_MT
        long_an=swelongs_an[int(steps/4)]
    return sin(-deg2rad(long_an - sid.theta_0(ag, anom_rev_diff)))/sin(deg2rad(sid.theta_0(ag, anom_rev_diff)))

# Get "actual" r0s over half an anomalistic year
r0s = [80*get_manda_r0byR_moon(i*ANOM_MT/steps, swelongs_an[i]) for i in range(int(steps))]
plt.plot(range(int(steps)), r0s, '.', label='Actual r0(t)')

def r0(ag, m, k):
    return m-k*cos(deg2rad(sid.theta_0(ag, anom_rev_diff)))

# Fit model to "actual" r0s
p0=[9, 0.5]
# plt.plot(range(int(steps)), [r0(i*ANOM_MT/steps, *p0) for i in range(int(steps))], '-', label='Initial Guess: %s'%p0)
r0byR_popt, r0byR_pcov = curve_fit(r0, [i*ANOM_MT/steps for i in range(int(steps))], np.array(r0s[0:int(steps)]), p0=p0)
plt.plot(range(int(steps)), [r0(i*ANOM_MT/steps, *r0byR_popt) for i in range(int(steps))], '-', label='Fitted Curve: %s'%r0byR_popt)

plt.ylim(8,10)
plt.title('Fitting model parameters for epicyclic radius r0 using data over half an anomalistic cycle')
plt.ylabel('r0(t)')
plt.xlabel('Time since apogee: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
plt.legend()

plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0_model_fitting_curve.png')
# plt.show()
plt.clf()

print('Finished r0byR parameter fitting computations')

plt.plot(range(int(steps)), [sid.r0_revised(sid.theta_0(i*ANOM_MT/steps, anom_rev_diff), 0, *r0byR_popt) for i in range(int(steps))], '-', label='Fitted Curve: %s'%r0byR_popt)

plt.title('Variation of epicyclic ratio r0 over one anomalistic year')
plt.ylabel('r0(t): epicyclic radius as a function of time')
plt.xlabel('Time since apogee: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))

plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0_model.png')
# plt.show()
plt.clf()
print('Plotted r0 over full cycle')

##------- Evaluating model precision --------##
r0byR=7/80
sidlongs = [sid.theta_ms(sid.theta_0(ag, anom_rev_diff), 0, r0byR) for ag in np.arange(0, ANOM_MT, ANOM_MT/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='Original r0byR = %f/80'%(r0byR*80))

r0byR=r0byR_popt[0]/80
sidlongs = [sid.theta_ms(sid.theta_0(ag, anom_rev_diff), 0, r0byR) for ag in np.arange(0, ANOM_MT, ANOM_MT/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='New constant r0byR = %f/80'%(r0byR*80))

sidlongs = [sid.theta_ms(sid.theta_0(ag, anom_rev_diff), 0, sid.r0_revised(sid.theta_0(ag, anom_rev_diff), 0, *r0byR_popt)/80) for ag in np.arange(0, ANOM_MT, ANOM_MT/steps)]
plt.plot(range(steps), subcirc(swelongs_an,sidlongs), '-', label='Varying r0byR with parameters: %s'%str(r0byR_popt))

plt.ylabel('Difference between observed and siddhanta anomalistic longitudes (degrees)')
plt.xlabel('Time since apogee: Step interval: ANOM_MT/%d = %f days'%(steps, ANOM_MT/steps))
plt.legend()

plt.gcf().set_size_inches(8, 4)
plt.savefig('./graphs/r0byR_final_comparison.png')
# plt.show()

plt.ylim(-0.1,0.1)
plt.savefig('./graphs/r0byR_final_comparison_zoomed.png')
# plt.show()
plt.clf()
print('Finished r0byR computations')


print('\nWARNING: DATA TAKEN FROM FILE. `make generate_data` to recompile generate_data.py')