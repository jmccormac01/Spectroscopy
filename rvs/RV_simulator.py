import math
import numpy as np
import matplotlib.pyplot as plt

rv_data = np.array([ -15.3670,
                -19.1260,
                -28.5270,
                -27.0880,
                -15.4430,
                -15.5170,
                -19.2930,
                -26.8680,
                -26.9280,
                -28.5540,
                -28.3790,
                -27.0540,
                -27.0160,
                -23.0050,
                -22.9080,
                -20.8120,
                -20.8250,
                -18.5430,
                -12.7070,
                -12.2190])

bjd = np.array([2456609.151870,
                2456610.146860,
                2456613.123780,
                2456616.165370,
                2456643.087190,
                2456643.109860,
                2456644.078960,
                2456646.085930,
                2456646.108650,
                2456647.090270,
                2456649.080180,
                2456650.078290,
                2456650.101180,
                2456669.091720,
                2456669.114510,
                2456670.085100,
                2456670.107830,
                2456671.088170,
                2456674.097610,
                2456675.087320])

def meanAnomoly(ta, e):
    """
    ta = true anomoly
    e = eccentricity
    the = tan half eccentricity
    ea = eccentric anomoly
    """
    the = math.tan(ta/2.)*math.sqrt((1.-e)/(1.+e))
    ea = 2.*math.atan(the)
    ma = ea - e*math.sin(ea)
    return ma

def trueAnomoly(ma, e):
    """
    ma = mean anomoly
    e = eccentricity
    ta = true anomoly
    guess = initial guess at ea
    ea = eccentric anomoly
    thnu = tan half nu
    """
    guess = ma
    diff = 1e9
    c = 0
    while diff > 1./206265:
        ea = ma + e*math.sin(guess)
        diff = abs(ea-guess)
        guess = ea
        if c > 100:
            print('true anomoly not goverged, breaking...')
            return None
        c += 1
    thnu = math.sqrt((1.+e)/(1.-e))*math.tan(ea/2.)
    ta = math.atan(thnu)*2.
    return ta

def kepler(bjd, tperi, e, K, w, gamma):
    """
    Calculate the radial velocites as per ACC
    mcmc code
    """
    vrad = np.empty(len(bjd))
    for j in range(0,len(bjd)):
        ma = 2.*math.pi*((bjd[j]-tperi)%P)/P
        nu = trueAnomoly(ma, e)
        vrad[j] = K*(e*math.cos(w)+math.cos(nu+w)) + gamma
        the = math.tan(nu/2.)*math.sqrt((1.-e)/(1.+e))
        ea = 2.*math.atan(the)
    return vrad

K = 8.399
w = math.radians(78.05)
E = 2453592.7443
P = 16.95347
e = 0.1599
gamma = -21.017

# make a model set of points to generate RVs
bjd_model = np.linspace(min(bjd), max(bjd), 500)
phase = (bjd-E)/P%1
phase_model = (bjd_model-E)/P%1

# get the time of periastron relative to the transit time t0
nutrans = math.pi/2. - w
mtrans = meanAnomoly(nutrans, e)
dttrans = mtrans/2./math.pi * P
tperi = E - dttrans

rv = kepler(bjd, tperi, e, K, w, gamma)
rv_model = kepler(bjd_model, tperi, e, K, w, gamma)

# residuals
res = rv_data - rv

# sort the data points
temp = zip(phase_model, rv_model)
temp = sorted(temp)
phase_model, rv_model = zip(*temp)

plt.plot(phase_model,rv_model,'r-', phase, rv_data, 'bo')
plt.show()

