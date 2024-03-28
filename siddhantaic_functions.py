from numpy import sin, cos, arcsin, rad2deg, deg2rad, sqrt, arctan

#####-------------- Traditional Model --------------#####

# Madhyagraha
civ_days = 1577917500
theta_0 = lambda _ag, _revs, _dhruva=0: (_dhruva + ((_ag*_revs/civ_days)%1)*360)%360

# Mandasamskara
theta_ms = lambda _theta_0, _mandocca, _r0byR: (_theta_0 - rad2deg(arcsin(_r0byR*sin(deg2rad(_theta_0-_mandocca)))))%360

RvbyR = lambda _t0, _tm, _r0byR: sqrt(1-_r0byR**2 * (sin(_t0-_tm))**2) - _r0byR*cos(_t0-_tm)
viparita_K = lambda _R, _t0, _tm, _r0byR: _R/RvbyR(_t0, _tm, _r0byR)

def iterated_K(_R, _t0, _tm, _r0byR):
    epsilon = 0.0000001
    rbyK = _r0byR
    KbyR = [0,1]
    while(abs(KbyR[1]-KbyR[0])>=epsilon):
        KbyR = [KbyR[1],sqrt((1+rbyK*cos(deg2rad(_t0-_tm)))**2+(rbyK*sin(deg2rad(_t0-_tm)))**2)]
        rbyK=KbyR[1]*rbyK
    return KbyR*_R

#####--------------- Modified Model ---------------#####

# Modified r0_by_R
def r0_revised(t0, tm, *popt):
    [m, k, p1, p2] = popt
    delt = t0-tm
    if delt<0: delt=delt+360
    if delt<180: return (m-k*arctan(p1*sin(2*deg2rad(delt))/(1-p2*cos(2*deg2rad(delt)))))
    else: return (m+k*arctan(p1*sin(2*deg2rad(delt))/(1-p2*cos(2*deg2rad(delt)))))