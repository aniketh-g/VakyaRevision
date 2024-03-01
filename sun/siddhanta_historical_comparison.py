from numpy import sqrt, sin, cos, arcsin, pi, rad2deg, deg2rad, subtract, arange
from generate_data import swelongs_an_hi, ANOM_YR, steps
import matplotlib.pyplot as plt
#####---------------- Sun's treatment according to Siddhanta ----------------#####

# Madhyagraha
sun_revs = 4320000
civ_days = 1577917500

theta_0 = lambda ahargana: ((ahargana*sun_revs/civ_days)%1)*360

# Mandasamskara
mandocca = 78
r0byR = 13.5/360

theta_ms = lambda _theta_0, _mandocca, _r0byR: (_theta_0 - rad2deg(arcsin(_r0byR*sin(deg2rad(_theta_0-_mandocca)))))%360

#####-------------------- Effect of changing year length --------------------#####
for swelongs_an in swelongs_an_hi:
    sidlongs = [theta_ms(theta_0(ag), 0, r0byR) for ag in arange(0, ANOM_YR, ANOM_YR/steps)]
    plt.plot(range(steps), subtract(swelongs_an[1],sidlongs), '-', label='Year = %d AD; Mandocca = %d'%(swelongs_an[0], swelongs_an[2]))

plt.ylim(-0.3,0.3)
plt.title('Longitudes computed from Siddhanta: Plotted for multiple years')
plt.ylabel('Difference between observed and siddhanta anomalistic longitudes (degrees)')
plt.xlabel('Step interval: ANOM_YR/%d = %f days'%(steps, ANOM_YR/steps))
plt.legend()
plt.grid()
plt.gcf().set_size_inches(12, 7)
plt.savefig('./graphs/siddhanta_historical_comparison.png')
plt.show()
plt.clf()
print('Finished year length computations')