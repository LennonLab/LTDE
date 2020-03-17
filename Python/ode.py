from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint


def f(y, t, params):
    """
    Your system of differential equations
    """

    N = y[0]
    D = y[1]
    S = y[2]

    try:
        v = params['v'].value
        K = params['K'].value
        c = params['c'].value
        m = params['m'].value
        r = params['r'].value
        d = params['d'].value

    except KeyError:
        v, K, c, m, r, d = params
    # the model equations
    uptake=v*(S/(K+S))*N # substrate uptake capacity of viable cells
    cellsMaintained=uptake/(c*m) # how many cells can meet maintenance costs based on substrate uptqke
    cells2die=d*N # how many cells would die if no maintenance costs were met

    dVdt=cellsMaintained-cells2die

    if (cellsMaintained-cells2die) < 0:
        dDdt = cells2die-cellsMaintained
    else:
        dDdt = 0

    dDdt = dDdt - r*D
    dSdt=r*D*c-uptake     #resource


    return [dVdt, dDdt, dSdt]



def g(t, x0, paras):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(f, x0, t, args=(paras,))
    return x



def residual(paras, t, data):

    """
    compute the residual between actual data and fitted data
    """

    x0 = paras['x10'].value, paras['x20'].value, paras['x30'].value
    model = g(t, x0, paras)

    # you only have data for one of your variables
    x0_model = model[:, 0]
    return (x0_model - data).ravel()


t_measured = np.asarray([10,24,39,52,63,74,88,104,117,131,147,160,174,187,202,216,235,\
            238,258,271,287,299,313,328,349,364,377,392,406,419,434,465,478,\
            496,508,523,537,551,565,580,593,607,620,638,638,663,677,690,702,\
            718,731,745,761,775,788,803,829,840,853,866,880,894,907,921,930])


x0_measured = np.asarray([564000000,260000000,348000000,300000000,372000000,248000000,\
                156000000,136000000,104000000,108000000,68000000,102000000,76000000,52000000,60000000,\
                22000000,13500000,9000000,15000000,17250000,8200000,5200000,4833333,1840000,1000000,\
                1740000,1860000,3620000,6380000,3480000,1180000,1940000,1880000,2460000,1870000,2480000,\
                1410000,1770000,780000,2710000,1200000,4400000,2360000,1280000,700000,920000,3180000,1260000,\
                2340000,880000,1400000,2660000,800000,2000000,1800000,4680000,1060000,298000,2270000,\
                102000,1066000,598000,1010000,2210000,300000])



# initial conditions
x10 = 280000000
x20 = 0
x30 = 0
y0 = [x10, x20, x30]


# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('x10', value=x10, vary=False)
params.add('x20', value=x20, vary=False)
params.add('x30', value=x30, vary=False)
params.add('v', value=220, vary=False)#min=20, max=500)
params.add('K', value=1.4e1, min=2, max=1e10)
params.add('c', value=100,vary=False)# min=10, max=500)
params.add('m', value=5, min=0.1, max=15)
params.add('r', value=0.005,vary=False)#, min=0.0001, max=0.8)
params.add('d', value=0.03,vary=False)#, min=0.0001, max=0.7)


# fit model
result = minimize(residual, params, args=(t_measured, x0_measured), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(np.linspace(0., 1000, 500), y0, result.params)

print(data_fitted)

# plot fitted data
plt.plot(np.linspace(0., 1000, 500), data_fitted[:, 0], '-', linewidth=2, color='red', label='fitted data')
plt.scatter(t_measured, x0_measured)
#plt.legend()
#plt.xlim([0, max(t_measured)])
#plt.ylim([0, 1.1 * max(data_fitted[:, 1])])
plt.yscale('log',basey=10)

# display fitted statistics
report_fit(result)

plt.show()




v = 220
K = 1.4e1
c = 100
m = 5
r = 0.005
d = 0.03

def dP_dt(P, t):
    uptake=v*(P[2]/(K+P[2]))*P[0] # substrate uptake capacity of viable cells
    cellsMaintained=uptake/(c*m) # how many cells can meet maintenance costs based on substrate uptqke
    cells2die=d*P[0] # how many cells would die if no maintenance costs were met

    dVdt=cellsMaintained-cells2die

    if (cellsMaintained-cells2die) < 0:
        dDdt = cells2die-cellsMaintained
    else:
        dDdt = 0

    dDdt = dDdt - r*P[1]
    dSdt=r*P[1]*c-uptake     #resource

    #dDdt=ifelse((cellsMaintained-cells2die)<0,(cells2die-cellsMaintained),0)-r*D               #cells
    return [dVdt, dDdt, dSdt]


#ts = np.linspace(0, 1000, 10000)
#P0 = [int(1e9), 0, 0]
#Ps = odeint(dP_dt, P0, ts)
#N = Ps[:,0]
#S = Ps[:,2]

#plt.plot(ts, N, ".")
#plt.yscale('log',basey=10)

#plt.show()
