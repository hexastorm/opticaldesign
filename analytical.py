# Author: Henri Starmans
# Company: Hexastorm
# Date: 3-3-2017

import numpy as np
from scipy.integrate import dblquad
from sympy import symbols, diff, cos, sin, sqrt


# Parameters used:
n=1.53                       # refractive index should equal s.refractive_index(405e-9,4)
facets = 4                   # number of polygon facets
utiltmax = np.radians(90-(180-360/facets)/2) # maximum angle of incidence possible [radians]
utilt=np.radians(28)         # maximum angle of incidence used [radians]
wavelength=0.405             # wavelength [micrometers]
T=35                         # T is 2 times inradius polygon [mm]
speed = np.pi*2*350           # radians per second
print(90-(180-360/facets)/2)

# Lens used 12.5mm Dia. x 90mm FL, VIS-NIR, Inked, Achromatic Lens from Edmund Optics
# NOTE: achromatic lens easier to simulate
#       the f-number can also be calculated for aspheric

D=0.8                # diameter bundle [mm]
efl=90               # effective focal length [mm]
bfl=88.47            # back focal length [mm]
ftol=0.02            # focal length tolerance 
d=12.5               # lens diameter [mm]

# F-number
fnumber=efl/D
# Spot size
# source: https://www.newport.com/n/gaussian-beam-optics
waist=2*wavelength/3.14*fnumber
print("The spot radius is "+str(round(waist,2))+" micrometers.")
raileigh=np.pi*waist**2/(wavelength*1000)
print("The rayleigh range is "+str(round(raileigh,3))+" mm.")
# Longitudinal focus shift
#   Wyant page 41 equation 68
slong=(n-1)/n*T
print("The focus point is at "+str(slong+bfl)+' mm.')
# Transversal focus shift
#   Wyant page 41 equation 70
x=symbols('x') # symbol we need for differentation
expr=T*sin(x)*(1-sqrt((1-sin(x)**2)/(n**2-sin(x)**2)))
disp=expr.evalf(subs={x:utilt})
print("The used transversal focus shift is "+str(disp)+" mm.")
dispmax = expr.evalf(subs={x:utiltmax})
print("The maximum transversal focus shift is "+str(dispmax)+" mm.")
# the transversal focus shift is dependent on the position plane
# we use it too calculate the maximum displacement (chain rule)
#  y=f(x) dx/dt=c   x=ct
#  --> dy/dt=dy/dx(x(t))*c 
sdisp=diff(expr,x)

print("The speed at the edges is "+str(round(sdisp.evalf(subs={x:utilt})/100*speed))+ " m/s.")
# fractional speed
fraction=sdisp.evalf(subs={x:0})/sdisp.evalf(subs={x:utilt})
print("The speed at the center is "+str(round(fraction*100,2))+" %"+" of the speed at the edges.")
# 3rd order Seidel aberrations in [mm]
#  Wyant page 42 equation 72; spherical aberration
sabr=-T/pow(fnumber,4)*((pow(n,2)-1)/(128*pow(n,3)))
#   Wyant page 44 equation 75; coma
#   Note: cosine has been omitted, will be added back when function f is defined
coma=-T*utilt/pow(fnumber,3)*((pow(n,2)-1)/(16*pow(n,3)))
#   Wyant page 45 equation 77 ; astigmatism
#   Note: cosine has been omitted, will be added back when function f is defined
astig=-T*pow(utilt,2)/pow(fnumber,2)*((pow(n,2)-1)/(8*pow(n,3)))
# To calculate the combined RMS we add back the functional forms to the coefficients
# They have been obtained from Wyant, page 17, table 2
# The function f defines the wavefront aberration denoted by Wyant as w
def f(theta,rho):
	return sabr*(rho**4)+astig*(rho**2)*(np.cos(theta)**2)+coma*(rho**3)*np.cos(theta)
# We now try to evaluate equation 62, page 37, Wyant
# First we define two auxiliary functions
def ws(theta,rho):
   return f(theta,rho)**2*rho
def w(theta,rho):
   return f(theta,rho)*rho
# Hereafter we evaluate the integral, i.e. equation 62
var=1/np.pi*dblquad(ws,0,1,lambda rho: 0, lambda rho: 2*np.pi)[0]\
-1/(np.pi**2)*(dblquad(w,0,1,lambda rho: 0, lambda rho: 2*np.pi)[0]**2)
# the RMS value is obtained by the root
rms=np.sqrt(var)
# The lambda RMS is expressed in [mm] and will be converted to wavelength units
lambdarms=rms/(wavelength*1E-3)
print("The lambda OPD RMS is " +str(round(lambdarms,6)))
# The Strehl radius is calculated with the first three terms of its Taylor series
# Wyant page 39 equation 67
rstrehl=1-(2*np.pi*lambdarms)**2+(2*np.pi*lambdarms)**4/2
print("The Strehl radius is "+str(round(rstrehl,2)))
