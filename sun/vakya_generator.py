from generate_data import SIDE_YR, ANOM_YR, mandocca_true
from numpy import rad2deg, arcsin, sin, deg2rad
from siddhanta_revision import r0byR_revised,\
    get_longitude_ag, theta_ms_new, r0byR_popt,\
    theta_m, theta_0, civ_days, sun_revs
import csv

romannum={1:'I', 2:'II', 3:'III', 4:'IV'}
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
\def\yv{month}{vnum}{{${min}'{sec}''$}}\n\
\def\yv{month}{vnum}vakya{{\\textit{{}}}}\n\n"\
.format(month=month_name, vnum=vakya_num, min=int((60*(x-x//1))),sec=round(60*((60*(x-x//1))-(60*(x-x//1))//1)))
    else:
        x=-x
        return "\def\yv{month}{vnum}sgn{{$-$}}\n\
\def\yv{month}{vnum}{{${min}'{sec}''$}}\n\
\def\yv{month}{vnum}vakya{{\\textit{{}}}}\n\n"\
.format(month=month_name, vnum=vakya_num, min=int((60*(x-x//1))),sec=round(60*((60*(x-x//1))-(60*(x-x//1))//1)))

def format_rashi_vakyas(x, month_name):
    return "\def\mv{month}Value{{{masa_v}}}\n\
\def\mv{month}Vakya{{\\textit{{\skta }}}}\n\n\
\def\sv{month}Value{{${sank_v_d}^d{sank_v_n}^n$}}\n\
\def\sv{month}Vakya{{\\textit{{\skta }}}}\n\n"\
.format(month=month_name, masa_v=round(x), sank_v_d=int(x%7), sank_v_n=round(60*(x%7-int(x%7))))

def format_naksatra_vakyas(x, naks_name):
    return "\def\\nv{month}Value{{${naks_v_d}^d{naks_v_n}^n$}}\n\
\def\\nv{month}Vakya{{\\textit{{\skta }}}}\n\n"\
.format(month=naks_name, masa_v=round(x), naks_v_d=int(x%7), naks_v_n=round(60*(x%7-int(x%7))))


def print_vakyas(ag, filepath, month_name):
    longs = [0] * 48
    longs[0] = theta_ms_new(theta_0(ag), mandocca_true, r0byR_revised(theta_0(ag), mandocca_true, *r0byR_popt)) # Assuming initial lambda_na
    vakyas = [0] * 48
    vakyas[0] = 0
    for i in range(1, 5):
        ag += 8
        longs[i] = theta_ms_new(theta_0(ag), mandocca_true, r0byR_revised(theta_0(ag), mandocca_true, *r0byR_popt))
        vakyas[i] = (longs[i] - longs[i - 1])%360 - 8
        with open(filepath, 'a') as vakya_file:
            vakya_file.write(format_yogyadi_vakya(vakyas[i], month_name, romannum[i]))

yogyadi_file = './tables/yogyadi_vakyas.txt'
masa_file = './tables/masa_vakyas.txt'
with open(yogyadi_file, 'w') as vakya_file: vakya_file.write("% Yogyadivakyas\n\n")
with open(masa_file, 'w') as vakya_file: vakya_file.write("% Masavakyas\n\n")

ag=0
ag0=get_longitude_ag(0, ag)
for i in range(12):
    ag = get_longitude_ag(i*30, ag)
    print_vakyas(ag, yogyadi_file, months[i+1])

    if i == 0:
        delta_masa = civ_days/sun_revs
    else:
        delta_masa = ((ag-ag0)/360)*civ_days/sun_revs
    with open(masa_file, 'a') as vakya_file:
        vakya_file.write(format_rashi_vakyas(delta_masa, months[i+1]))

naksatra_file = './tables/naksatra_vakyas.txt'
with open(naksatra_file, 'w') as vakya_file: vakya_file.write("% Naksatravakyas\n\n")
# naksatra_table = './tables/naksatra_table.txt'
# with open(naksatra_table, 'w') as vakya_file: vakya_file.write("")

ag=0
ag0=get_longitude_ag(0, ag)
for i in range(27):
    ag = get_longitude_ag(i*360/27, ag)
    if i == 0:
        delta_naks = civ_days/sun_revs
    else:
        delta_naks = ((ag-ag0)/360)*civ_days/sun_revs
    with open(naksatra_file, 'a') as vakya_file:
        vakya_file.write(format_naksatra_vakyas(delta_naks, naksatras[i+1]))
    # with open(naksatra_table, 'a') as vakya_file:
    #     vakya_file.write("\\textit{{{naks_name}}} &$ &\\nv{naks_name}Value\\\\\n\
# &{{\\skta }}&\\nv{naks_name}Vakya\\\\\n".format(naks_name=naksatras[i+1]))