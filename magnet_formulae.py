# From Magnet Formulas, Eric Dennison
import numpy as np
from scipy.special import ellipk, ellipe, ellipkm1
from numpy import pi, sqrt, linspace
from pylab import plot, xlabel, ylabel, suptitle, legend, show

uo = 4E-7*pi     # Permeability constant - units of H/m
Bo = lambda I, a, u=uo: I*u/2./a    # Central field = f(current, loop radius, perm. constant)
al = lambda r, R: r/R               # Alpha = f(radius of measurement point, radius of loop)
be = lambda z, a: z/a               # Beta = f(axial distance to meas. point, radius of loop)
ga = lambda x, r: x/r               # Gamma = f(axial distance, radius to meas. point)
Q = lambda r, x, a: (1 + al(r,a))**2 + be(x,a)**2   # Q = f(radius, distance to meas. point, loop radius)
k = lambda r, x, a: sqrt(4*al(r,a)/Q(r,x,a))       # k = f(radius, distance to meas. point, loop radius)
K = lambda k: ellipk(k**2.0)          # Elliptic integral, first kind, as a function of k
E = lambda k: ellipe(k**2.0)          # Elliptic integral, second kind, as a function of k

# On-Axis field = f(current and radius of loop, x of measurement point)
def Baxial(I, R, z, u=uo):
    if R == 0:
        if z == 0:
            return np.NaN
        else:
            return 0.0
    else:
        return (u*I*R**2)/2.0/(R**2 + z**2)**(1.5)

# Axial field component = f(current and radius of loop, r and x of meas. point)
def ring_Bz(I, R, z, r):
    if r == 0:
        if z == 0:
            return Bo(I,R)         # central field
        else:
            return Baxial(I,R,z)   # axial field
    else:                          # axial component, any location
        return Bo(I,R)*\
            (E(k(r,z,R))*((1.0-al(r,R)**2-be(z,R)**2)/(Q(r,z,R)-4*al(r,R))) + K(k(r,z,R)))\
            /pi/sqrt(Q(r,z,R))
        
# Radial field component = f(current and radius of loop, r and x of meas. point)
def ring_Br(I, R, z, r):
    if r == 0:
        return 0                   # no radial component on axis!
    else:                          # radial component, any location other than axis.
        return Bo(I,R)*ga(z,r)*\
            (E(k(r,z,R))*((1.0+al(r,R)**2+be(z,R)**2)/(Q(r,z,R)-4*al(r,R))) - K(k(r,z,R)))\
            /pi/sqrt(Q(r,z,R))
    
class CoilGeometry():
    def __init__(self, coil_list): #r, z, current
        self.coil_list = coil_list

    def calculate_field(self, r, z):
        TESLA_TO_GAUSS = 1e4  # Conversion factor from Tesla to Gauss
        Br = 0
        Bz = 0

        for coil_r, coil_z, current in self.coil_list:

            Br += ring_Br(current, coil_r, z - coil_z, r)
            Bz += ring_Bz(current, coil_r, z - coil_z, r)

        return Br * TESLA_TO_GAUSS, Bz * TESLA_TO_GAUSS