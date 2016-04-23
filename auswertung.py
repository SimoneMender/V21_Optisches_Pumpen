import matplotlib as mpl
import numpy as np

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat


# Konstanten
mu0 = 1.25663706e-06
e0 = 1.602176e-19
m0 = 9.109383e-31

# Spule Vertikalkomp. des Erdmagnetfeld
N1 = 20
R1 = 11.735*10**(-2)
I1 = 2.3*0.1

# Spule Horizontalkomp. des Erdmagnetfeldes
N2 = 154
R2 = 15.79*10**(-2)

# Sweep-Spule
N3 = 11
R3 = 16.39*10**(-2)



# Werte

F, H1, S1, H2, S2 =np.genfromtxt('werte.txt', unpack=True)


#Vertikalkomp. des Erdmagnetfeld berechnen

B_erd = mu0 * 8/np.sqrt(125)* (I1*N1)/R1
B_erd=B_erd/10**(-6)
print('B_erd')
print(B_erd)

# Horizontalkomp. des Erdmagnetfeld berechnen

# Isotop 1

I2 = H1 * 0.3
I3 = S1 * 0.1
F = F

B_hor = (mu0 * 8/np.sqrt(125)* (I3*N3)/R3) + (mu0 * 8/np.sqrt(125)* (I2*N2)/R2)

B_hor*=10**(6)

print('dfg')
print(B_hor)

def f(B_hor, a, b):
    return a*B_hor+ b

params, covariance = curve_fit(f, B_hor, F)


errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

a1 = ufloat(params[0], errors[0])
b1 = ufloat(params[1], errors[1])

# Nullstelle ist B-Feld

x1 = -b1/a1
print('B_hor')
print(x1)

# Isotop 2

I4 = H2 * 0.3
I5 = S2 * 0.1

B_hor = mu0 * 8/np.sqrt(125)* (I5*N3)/R3 + mu0 * 8/np.sqrt(125)* (I4*N2)/R2
#print('dfg')

B_hor*=10**6
def f(B_hor, c, d):
    return c*B_hor+ d

params, covariance = curve_fit(f, B_hor, F)


errors = np.sqrt(np.diag(covariance))

print('c =', params[0], '±', errors[0])
print('d =', params[1], '±', errors[1])

c1 = ufloat(params[0], errors[0])
d1 = ufloat(params[1], errors[1])

# Nullstelle ist B-Feld

x2=-d1/c1
print('B_hor')
print(x2)

# Gesamtwert der Horizontalkomp.

x3 = (x1+x2)/2
print('B_hor.ges.')
print(x3)

# Lande-Faktoren berechnen

g1 = (4*np.pi*m0*a1)/(e0)
print('g1')
print(g1*10**6)
g2 = (4*np.pi*m0*c1)/(e0)
print('g2')
print(g2*10**6)

# Kernspin berechnen
