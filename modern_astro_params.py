import shelve
modern_astro_params_file = shelve.open('./modern_astronomical_constants/modern_astronomical_constants')

# All values taken from Astronomical Almanac for the year 2014

#####-------------- Solar Constants --------------#####

SIDE_YR =  365.256363 # Sidereal year
ANOM_YR =  365.259636 # Anomalistic year

modern_astro_params_file['SIDE_YR'] = SIDE_YR
modern_astro_params_file['ANOM_YR'] = ANOM_YR

sun_revs = 4320000
civ_days_new = round(sun_revs*SIDE_YR)
sun_mandocca_revs_new = round(sun_revs-civ_days_new/ANOM_YR)

modern_astro_params_file['sun_revs'] = sun_revs
modern_astro_params_file['civ_days_new'] = civ_days_new
modern_astro_params_file['sun_mandocca_revs_new'] = sun_mandocca_revs_new

print('civ_days_new', civ_days_new)
print('sun_mandocca_revs_new', sun_mandocca_revs_new)

#####-------------- Lunar Constants --------------#####

SYNO_MT = 29.530589 # Synodic Month
SIDE_MT = 27.321662 # Sidereal Month
ANOM_MT = 27.554550 # Anomalistic Month

modern_astro_params_file['SIDE_MT'] = SIDE_MT
modern_astro_params_file['ANOM_MT'] = ANOM_MT
modern_astro_params_file['SYNO_MT'] = SYNO_MT

moon_revs_new = round(civ_days_new/SIDE_MT)
moon_mandocca_revs_new = round(moon_revs_new-civ_days_new/ANOM_MT)

modern_astro_params_file['moon_revs_new'] = moon_revs_new
modern_astro_params_file['moon_mandocca_revs_new'] = moon_mandocca_revs_new

print('moon_revs_new', moon_revs_new)
print('moon_mandocca_revs_new', moon_mandocca_revs_new)