from numpy import sin, cos, deg2rad, rad2deg, abs
import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')
import siddhantaic_functions as sid
import moon.functions as moon
from useful_functions import todms, intToRoman

#####------------------------ Import modern constants ------------------------#####
print('#####------------------------ Import modern constants ------------------------#####')
import shelve
with shelve.open('../modern_astronomical_constants/modern_astronomical_constants') as astro_const:
    ANOM_MT = astro_const['ANOM_MT']
    SIDE_MT = astro_const['SIDE_MT']
    SYNO_MT = astro_const['SYNO_MT']

    sun_revs = astro_const['sun_revs']
    civ_days_new = astro_const['civ_days_new']
    moon_revs_new = astro_const['moon_revs_new']
    moon_mandocca_revs_new = astro_const['moon_mandocca_revs_new']

#####----------------- Import revised Sun and Moon Parameters ----------------#####
print('#####----------------- Import revised Sun and Moon Parameters ----------------#####')

with shelve.open('./moon_computed_data/moon_revised_parameters') as moon_revised_params:
    manda_radius_oscillation = moon_revised_params['manda_radius_oscillation']
    dvitiya_radius = moon_revised_params['dvitiya_radius']
    trtiya_radius = moon_revised_params['trtiya_radius']
    annual_equation_factor = moon_revised_params['annual_equation_factor']
moon_manda_r0_popt = [7, manda_radius_oscillation]

with shelve.open('../sun/sun_computed_data/sun_modified_parameters') as sun_revised_params:
    sun_mandocca_mean = sun_revised_params['mandocca_mean']
    sun_mandocca_revs = sun_revised_params['mandocca_revs']
    sun_manda_r0_popt = sun_revised_params['r0byR_popt']

from useful_functions import subcirc
import matplotlib.pyplot as plt

print("Starting Candravakya")
# Compute Mandala-s

COMPUTE_MANDALAS = True
SYNODIC_CYCLES = False
if COMPUTE_MANDALAS:
    print("#####---------------- Computing Mandalas ----------------#####")
    if SYNODIC_CYCLES:
        syn_period = 1/(moon_revs_new/civ_days_new - sun_revs/civ_days_new)
        print("Synodic Period used %f vs Actual %f\n"%(syn_period, SYNO_MT))
        print("SynodicCycle\tMoonApogee\tMeanMoon\tMeanSun\t\tMandoccaDhruva")
        for i in range(20):
            ag = syn_period*i
            mandocca_dhruva = sid.theta_0(ag, moon_mandocca_revs_new, 0, civ_days_new)-sid.theta_0(ag, moon_revs_new, 0, civ_days_new)
            print("%d\t\t%f\t%f\t%f\t%f"%(i, sid.theta_0(ag, moon_mandocca_revs_new, 0, civ_days_new), sid.theta_0(ag, moon_revs_new, 0, civ_days_new), sid.theta_0(ag, sun_revs, 0, civ_days_new), mandocca_dhruva))

    anom_period = 1/(moon_revs_new/civ_days_new - moon_mandocca_revs_new/civ_days_new)
    print("Anomalistic Period used %f vs Actual %f\n"%(anom_period, ANOM_MT))
    print("AnomalyCycle\tDays\t\tMoonApogee\tMeanMoon\tMeanSun\t\tSunDhruva")
    for i in range(20):
        ag = anom_period*i
        sun_dhruva = sid.theta_0(ag, sun_revs, 0, civ_days_new)-sid.theta_0(ag, moon_revs_new, 0, civ_days_new)
        if(abs(sun_dhruva) <= 360): #0.011586
            print("%d\t\t%f\t%f\t%f\t%f\t%f"%(i, i*anom_period, sid.theta_0(ag, moon_mandocca_revs_new, 0, civ_days_new), sid.theta_0(ag, moon_revs_new, 0, civ_days_new), sid.theta_0(ag, sun_revs, 0, civ_days_new), sun_dhruva))

def trtiya_moon(_mean_moon, _moon_mandocca, _mean_sun):
    lm_mean = _mean_moon
    K_mean = 80

    lu = _moon_mandocca
    r_manda = sid.r0_revised(lm_mean, lu, *moon_manda_r0_popt)
    lm_manda = sid.theta_ms(lm_mean, lu, r_manda/K_mean)
    K_manda = 10*moon.R*sid.KbyR(lm_mean, lu, r_manda/K_mean)

    # lu_dvitiya = sid.theta_ms(_mean_sun, sun_mandocca_mean, sid.r0_revised(_mean_sun, sun_mandocca_mean, *sun_manda_r0_popt)/360)
    lu_dvitiya = _mean_sun
    r_dvitiya = (moon.R*dvitiya_radius)*cos(deg2rad(lu_dvitiya-lu))
    lm_dvitiya = sid.theta_ms(lm_manda, lu_dvitiya, r_dvitiya/K_manda)
    K_dvitiya = sid.KbyR(lm_manda, lu_dvitiya, r_dvitiya/K_manda)

    lu_trtiya = 2*lu_dvitiya-lm_dvitiya
    r_trtiya = -trtiya_radius*K_dvitiya
    lm_trtiya = sid.theta_ms(lm_dvitiya, lu_trtiya, r_trtiya/K_dvitiya)

    return lm_trtiya

def manda_sun(_mean_sun):
    r0 = sid.r0_revised(_mean_sun, sun_mandocca_mean, *sun_manda_r0_popt)
    manda_sun = sid.theta_ms(_mean_sun, sun_mandocca_mean, r0/360)
    return manda_sun

def format_candravakya(v, c, vakya_num):
    lras = v[0]
    ldeg = v[1]
    lmin = v[2]+round(v[3]/60)

    cmin = c[2]
    csec = c[2]*60+c[3]+round(c[4]/60)
    return "\
\def\cvlong{vnum}{{${lras:02d}^r{ldeg:02d}\degree{lmin:02d}'$}}\n\
\def\cvlong{vnum}vakya{{\\textit{{\skta }}}}\n\
\def\cvcorr{vnum}{{${csec:03d}''$}}\n\
\def\cvcorr{vnum}vakya{{\\textit{{\skta }}}}\n\n"\
.format(vnum=intToRoman(vakya_num), lras=lras, ldeg=ldeg, lmin=lmin, cmin=cmin, csec=csec)

def format_old_candravakya(v, vakya_num):
    lras = v[0]
    ldeg = v[1]
    lmin = v[2]+round(v[3]/60)

    if vakya_num <= 248:
        return "\
\def\oldcvlong{vnum}{{${lras:02d}^r{ldeg:02d}\degree{lmin:02d}'$}}\n\
\def\oldcvlong{vnum}vakya{{\\textit{{\skta }}}}\n\n\
".format(vnum=intToRoman(vakya_num), lras=lras, ldeg=ldeg, lmin=lmin)
    else:
        return "\
\def\oldcvlong{vnum}{{$*{lras:02d}^r{ldeg:02d}\degree{lmin:02d}'$}}\n\
\def\oldcvlong{vnum}vakya{{\\textit{{\skta }}}}\n\n\
".format(vnum=intToRoman(vakya_num), lras=lras, ldeg=ldeg, lmin=lmin)


cv_table_format = './tables/cv_table_format.tex'
with open(cv_table_format, 'w') as table_format: table_format.write("% Candravakyas Comparison Table Format\n\n")

new_candravakya_file = './tables/candra_vakyas_mean_sun.tex'
with open(new_candravakya_file, 'w') as new_vakya_file: new_vakya_file.write("% Candravakyas\n\n")
old_candravakya_file = './tables/old_candra_vakyas.tex'
with open(old_candravakya_file, 'w') as old_vakya_file: old_vakya_file.write("% Old Candravakyas\n\n")

old_vakyas = []
vakyas = []
correction_terms = []
for ag in range(round(ANOM_MT*15+2)):
    mean_moon = sid.theta_0(ag, moon_revs_new, 0, civ_days_new)
    mean_sun = sid.theta_0(ag, sun_revs, 0, civ_days_new)
    moon_mandocca = sid.theta_0(ag, moon_mandocca_revs_new, 0, civ_days_new)
    v = trtiya_moon(mean_moon, moon_mandocca, mean_sun)
    vakyas.append(v)

    dx = 1e-4
    c = (trtiya_moon(mean_moon, moon_mandocca, mean_sun+dx/2)-trtiya_moon(mean_moon, moon_mandocca, mean_sun-dx/2))/dx
    correction_terms.append(c)
    # print(todms(v), todms(c))
    # print(v, trtiya_moon(mean_moon, moon_mandocca, mean_sun+180), c, (trtiya_moon(mean_moon, moon_mandocca, mean_sun+180+dx/2)-trtiya_moon(mean_moon, moon_mandocca, mean_sun+180-dx/2))/dx)

    with open(cv_table_format, 'a') as table_format:
        table_format.write(
"""${vnum}$ & \oldcvlong{romannum} & \cvlong{romannum} & \cvcorr{romannum}\\\\\n""".format(vnum = ag, romannum = intToRoman(ag))
        )

    with open(new_candravakya_file, 'a') as new_vakya_file:
        new_vakya_file.write(format_candravakya(todms(v), todms(c), ag))
    with open(old_candravakya_file, 'a') as old_vakya_file:
        if ag <=248:
            old_moon = sid.theta_ms(mean_moon, moon_mandocca, 7/80)
            old_vakyas.append(old_moon)
            old_vakya_file.write(format_old_candravakya(todms(old_moon), ag))
        else:
            old_vakya_file.write(format_old_candravakya(todms(old_vakyas[248]+old_vakyas[ag-248]), ag))
print("Computed Candravakyas:\t\t%s"%new_candravakya_file)

# Test

for sun_dhruva_test in [11.52, -11.52]:
    actual_longs = []
    vakya_corrected_longs = []
    for ag in range(round(ANOM_MT*15+1)):
        mean_moon = sid.theta_0(ag, moon_revs_new, 0, civ_days_new)
        mean_sun = sid.theta_0(ag, sun_revs, sun_dhruva_test, civ_days_new)
        moon_mandocca = sid.theta_0(ag, moon_mandocca_revs_new, 0, civ_days_new)
        actual_longs.append(trtiya_moon(mean_moon, moon_mandocca, mean_sun))
        vakya_corrected_longs.append(vakyas[ag] + correction_terms[ag]*sun_dhruva_test)

    # plt.plot(range(405), [60*((actual_longs[ag]-vakyas[ag]+180)%360-180) for ag in range(405)], label='uncorrected vakya, sun_dhruva_test = %f'%sun_dhruva_test)
    plt.plot(range(405), [60*((actual_longs[ag]-vakya_corrected_longs[ag]+180)%360-180) for ag in range(405)], label='corrected vakya, $Sun\'s\ Dhruva = %f$'%sun_dhruva_test)
plt.ylabel('$Corrected\ Vakya - Revised\ Epicyclic\ Longitude$ (min)')
plt.xlabel('$Vakya\ number$')
plt.legend()
plt.grid()
plt.legend()
plt.gcf().set_size_inches(12, 6)
plt.savefig('./graphs/vakya_max_error.png')
plt.show()
plt.clf()