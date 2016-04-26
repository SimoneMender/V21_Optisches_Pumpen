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

#B_hor*=10**(6)


def f(F, a, b):
    return a*F+ b

params, covariance = curve_fit(f, F, B_hor)


errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0], 'in kHz/T')
print('b =', params[1], '±', errors[1], 'in kHz')

a1 = ufloat(params[0], errors[0])
b1 = ufloat(params[1], errors[1])

# Nullstelle ist B-Feld

#x1 = -b1/a1
print('B_hor')
print(b1)

# Graph Isotop 1

plt.plot(F, B_hor, 'bx', label='Werte Isotop 1')

B_plot = np.linspace(0, 1000, 20000)
plt.plot(B_plot, f(B_plot, *params), 'c-', label='Ausgleichsgerade Isotop 1')
plt.xlabel('Frequenz in kHz')
plt.ylabel('B-Feld in T')




# Isotop 2

I4 = H2 * 0.3
I5 = S2 * 0.1

B_hor = mu0 * 8/np.sqrt(125)* (I5*N3)/R3 + mu0 * 8/np.sqrt(125)* (I4*N2)/R2
#print('dfg')

#B_hor*=10**6
def f(F, c, d):
    return c*F+ d

params, covariance = curve_fit(f, F, B_hor)


errors = np.sqrt(np.diag(covariance))

print('c =', params[0], '±', errors[0], 'in kHz/T')
print('d =', params[1], '±', errors[1], 'in kHz')

c1 = ufloat(params[0], errors[0])
d1 = ufloat(params[1], errors[1])

# Nullstelle ist B-Feld

#x2=-d1/c1
print('B_hor')
print(d1)

# Gesamtwert der Horizontalkomp.

x3 = (b1+d1)/2
print('B_hor.ges.')
print(x3)


# Graph Isotop 2

plt.plot(F, B_hor, 'rx', label='Werte Isotop 2')

plt.plot(B_plot, f(B_plot, *params), 'k-', label='Ausgleichsgerade Isotop 2')
plt.xlabel('Frequenz in kHz')
plt.ylabel('B-Feld in T')

plt.legend(loc="best")

plt.tight_layout()
plt.savefig('graph2.pdf')


# Lande-Faktoren berechnen

g1 = (4*np.pi*m0*1000)/(e0*a1)
print('g1')
print(g1)
g2 = (4*np.pi*m0*1000)/(e0*c1)
print('g2')
print(g2)

# Kernspin berechnen

gj = (3.0023*0.5*(0.5+1)+1.0023*(0.5*(0.5+1)))/(2*0.5*(0.5+1))
print('gj')
print(gj)

#I1 = -1 +(gj/(4*g1))+( (1-(gj/(4*g1)))**2 - (1-gj/g1) )**(0.5)
#print('I_1')
#print(I1)

#I2 = -1+ gj/(4*g2)+((1-(gj/(4*g2)))**2 - (1-gj/g2))**(0.5)
#print('I_2')
#print(I2)

# Kernspin Test
muF = 5.05*10**(-27)
muB = 9.27e-24


x1=-(4*g1-gj)/(4*g1)+(((gj-4*g1)/(4*g1))**2-3/4+3/4*(gj)/(g1))**(0.5)
print('I1')
print(x1)

x2=-(4*g2-gj)/(4*g2)+(((gj-4*g2)/(4*g2))**2-3/4+3/4*(gj)/(g2))**(0.5)
print('I2')
print(x2)



# Quadratischer Zeeman-Effekt

U1=g1*muB*B_hor+g1**2*muB**2*B_hor**2*(1-2)/(4.53e-24)
#print(U1)

U2=g2*muB*B_hor+g2**2*muB**2*B_hor**2*(1-2)/(2.01e-24)
#print(U2)
