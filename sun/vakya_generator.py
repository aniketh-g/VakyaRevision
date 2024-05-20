from numpy import rad2deg, arcsin, sin, deg2rad

import sys
sys.path.insert(0, '/home/aniketh/VakyaRevision')

import shelve
with shelve.open('../sun/sun_computed_data/sun_modified_parameters') as sunp:
    sun_civ_days = sunp['sun_civ_days']
    sun_mandocca_mean = sunp['mandocca_mean']
    sun_r0byR_popt = sunp['r0byR_popt']
sunp.close()

sun_revs = 4320000

from useful_functions import *
import siddhantaic_functions as sid

def sun_mandasphuta(_ag):
    _t0 = sid.theta_0(_ag, sun_revs, 0, sun_civ_days)
    _tm = sun_mandocca_mean
    _r0byR = sid.r0_revised(_t0, _tm, *sun_r0byR_popt)/360
    return sid.theta_ms(_t0, _tm, _r0byR)

def sun_mandasphuta_old(_ag):
    _t0 = sid.theta_0(_ag, sun_revs)
    _tm = 78
    _r0byR = 13.5/360
    return sid.theta_ms(_t0, _tm, _r0byR)

def intToRoman(num):
    # Storing roman values of digits from 0-9
    # when placed at different places
    m = ["", "M", "MM", "MMM"]
    c = ["", "C", "CC", "CCC", "CD", "D", "DC", "DCC", "DCCC", "CM"]
    x = ["", "X", "XX", "XXX", "XL", "L", "LX", "LXX", "LXXX", "XC"]
    i = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"]
    # Converting to roman
    thousands = m[num // 1000]
    hundreds = c[(num % 1000) // 100]
    tens = x[(num % 100) // 10]
    ones = i[num % 10]
 
    ans = (thousands + hundreds + tens + ones) 
    return ans

months={
    1:'Mesa',
    2:'Vrsabha',
    3:'Mithuna',
    4:'Karkata',
    5:'Simha',
    6:'Kanya',
    7:'Tula',
    8:'Vrscika',
    9:'Dhanus',
    10:'Makara',
    11:'Kumbha',
    12:'Mina',
}

naksatras={
    1:'Asvini',
    2:'Bharani',
    3:'Krttika',
    4:'Rohini',
    5:'Mrgasira',
    6:'Ardra',
    7:'Punarvasu',
    8:'Pusya',
    9:'Aslesa',
    10:'Magha',
    11:'Purvaphalguni',
    12:'Uttaraphalguni',
    13:'Hasta',
    14:'Citra',
    15:'Svati',
    16:'Vaisakha',
    17:'Anuradha',
    18:'Jyestha',
    19:'Mula',
    20:'Purvasadha',
    21:'Uttarasadha',
    22:'Sravana',
    23:'Dhanistha',
    24:'Satabhisa',
    25:'Purvabhadrapada',
    26:'Uttarabhadrapada',
    27:'Revati',
}

def format_yogyadi_vakya(x, month_name, vakya_num):
    if x>0:
        return "\def\yv{month}{vnum}sgn{{$+$}}\n\
\def\yv{month}{vnum}{{${min:02d}'{sec:02d}''$}}\n\
\def\yv{month}{vnum}vakya{{\\textit{{}}}}\n\n"\
.format(month=month_name, vnum=vakya_num, min=int((60*(x-x//1))),sec=round(60*((60*(x-x//1))-(60*(x-x//1))//1)))
    else:
        x=-x
        return "\def\yv{month}{vnum}sgn{{$-$}}\n\
\def\yv{month}{vnum}{{${min:02d}'{sec:02d}''$}}\n\
\def\yv{month}{vnum}vakya{{\\textit{{}}}}\n\n"\
.format(month=month_name, vnum=vakya_num, min=int((60*(x-x//1))),sec=round(60*((60*(x-x//1))-(60*(x-x//1))//1)))

def format_rashi_vakyas(x, month_name):
    days = x%7
    nadis = 60*(days-int(days))
    vinadis = 60*(nadis-int(nadis))
    return "\def\mv{month}Value{{{masa_v}}}\n\
\def\mv{month}Vakya{{\\textit{{\skta }}}}\n\n\
\def\sv{month}Value{{${sank_v_d}^d{sank_v_n:02d}^n{sank_v_v:02d}^v$}}\n\
\def\sv{month}Vakya{{\\textit{{\skta }}}}\n\n"\
.format(month=month_name, masa_v=round(x), sank_v_d=int(days), sank_v_n=int(nadis), sank_v_v=round(vinadis))

def format_naksatra_vakyas(x, naks_name):
    return "\def\\nv{naks_name}Value{{${naks_v_d}^d{naks_v_n:02d}^n$}}\n\
\def\\nv{naks_name}Vakya{{\\textit{{\skta }}}}\n\n"\
.format(naks_name=naks_name, masa_v=round(x), naks_v_d=int(x%7), naks_v_n=round(60*(x%7-int(x%7))))

def format_bhupajnadi_vakya(x, vakya_num):
    return "\def\\bv{vakya_num}Value{{${min}'{sec}''$}}\n\
\def\\bv{vakya_num}Vakya{{\\textit{{\skta }}}}\n\n"\
.format(vakya_num=intToRoman(vakya_num), min=int(60*x), sec=round(60*((60*(x-x//1))-(60*(x-x//1))//1)))

print("-----Beginning Vakya Computations-----")

if False:
    with shelve.open('../modern_astronomical_constants/modern_astronomical_constants') as modern_astronomical_constants:
        SIDE_YR = modern_astronomical_constants['SIDE_YR']
    modern_astronomical_constants.close()
    ag0=sid.get_ag(sun_mandasphuta, 0, 0, 1)
    import matplotlib.pyplot as plt
    sidlongs=[]
    sidlongs_old=[]
    for year in range(num_years):
        for i in range(steps):
            ag = ag0+year*SIDE_YR+i*SIDE_YR/steps
            sidlongs.append(sun_mandasphuta(ag))
            sidlongs_old.append(sun_mandasphuta_old(ag))
    plt.plot(range(num_years*steps), subcirc(sidlongs, swelongs_si), '-', label='Modified')
    plt.plot(range(num_years*steps), subcirc(sidlongs_old, swelongs_si), '-', label='Old')

    plt.ylim(-0.015,0.015)
    plt.ylabel('Difference between observed and siddhanta sidereal longitudes (degrees)')
    plt.xlabel('Time since Aries Transit: Step interval: SIDE_YR/%d = %f days'%(steps, SIDE_YR/steps))
    plt.legend()
    plt.grid()
    plt.gcf().set_size_inches(8, 4)
    plt.show()
    plt.clf()

def print_yogyadi_vakyas(ag, filepath, month_name, sun_algorithm):
    longs = [0] * 48
    longs[0] = sun_algorithm(ag) # Assuming initial lambda_na
    vakyas = [0] * 48
    vakyas[0] = 0
    for i in range(1, 5):
        ag += 8
        longs[i] = sun_algorithm(ag)
        vakyas[i] = (longs[i] - longs[i - 1])%360 - 8
        with open(filepath, 'a') as vakya_file:
            vakya_file.write(format_yogyadi_vakya(vakyas[i], month_name, intToRoman(i)))

yogyadi_file = './tables/yogyadi_vakyas.txt'
rashi_file = './tables/rashi_vakyas.txt'
with open(yogyadi_file, 'w') as vakya_file: vakya_file.write("% Yogyadivakyas\n\n")
with open(rashi_file, 'w') as vakya_file: vakya_file.write("% Rashivakyas\n\n")

ag=0
ag0=sid.get_ag(sun_mandasphuta, 0, ag, 1)
for i in range(12):
    ag = sid.get_ag(sun_mandasphuta, i*30, ag, 1)
    print_yogyadi_vakyas(ag, yogyadi_file, months[i+1], sun_mandasphuta)

    if i == 0:
        delta_masa = sun_civ_days/sun_revs
    else:
        delta_masa = ag-ag0
    with open(rashi_file, 'a') as vakya_file:
        vakya_file.write(format_rashi_vakyas(delta_masa, months[i+1]))
print("Computed Yogyadi Vakyas:\t\t%s"%yogyadi_file)
print("Computed Masa and Sankranti Vakyas:\t%s"%rashi_file)

naksatra_file = './tables/naksatra_vakyas.txt'
with open(naksatra_file, 'w') as vakya_file: vakya_file.write("% Naksatravakyas\n\n")

ag=0
ag0=sid.get_ag(sun_mandasphuta, 0, ag, 1)
for i in range(27):
    ag = sid.get_ag(sun_mandasphuta, i*360/27, ag, 1)
    if i == 0:
        delta_naks = sun_civ_days/sun_revs
    else:
        delta_naks = ag-ag0
    with open(naksatra_file, 'a') as vakya_file:
        vakya_file.write(format_naksatra_vakyas(delta_naks, naksatras[i+1]))
print("Computed Naksatra Vakyas:\t\t%s"%naksatra_file)

bhupajnadi_file = './tables/bhupajnadi_vakyas.txt'
with open(bhupajnadi_file, 'w') as vakya_file: vakya_file.write("% Bhupajnadivakyas\n\n")
ag=0
ag0=sid.get_ag(sun_mandasphuta, 0, ag, 1)
for i in range(1, 38):
    ag = ag0+10*i
    v = 10*i-sun_mandasphuta(ag)
    with open(bhupajnadi_file, 'a') as vakya_file:
        vakya_file.write(format_bhupajnadi_vakya(v, i))
print("Computed Bhupajnadi Vakyas:\t\t%s"%bhupajnadi_file)

print("-----Finished Vakya Computations-----")