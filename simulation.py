%pylab inline
import warnings
import numpy as np
import matplotlib.pyplot as plt
import rayopt as ro

# Lens used 12.5mm Dia. x 90mm FL, VIS-NIR, Inked, Achromatic Lens from Edmund Optics
# LINK: http://www.edmundoptics.com/document/download/391099
filename='zmax_49332ink.zmx'   
with open(filename) as file:
    data=file.read()

# Parameters:
wavelength=405e-9        # wavelength [m]
D=0.8                             # diameter bundle [mm] see, s.scale=0.001 [m]
T=35                              # plate thickness [mm]
utilt=np.radians(0)       # tilt angle plate [radians]

# Radius plate
#  I can't remember why I use tangent, should not matter 
#   as long diameter is large enough
spol=T*np.tan(np.pi/8)

# Create the system
s=ro.zemax.zmx_to_system(data)
s.object.pupil.radius = D/2
# Ensures rays created with function ray_point are in the [-D/2,D/2] range
s.object.pupil.update_radius = False 
s.object.angle = np.radians(0)   # [radians]                     
s.wavelengths = [wavelength]
s.update()  
# changes needed to make the Zemax data compatible with Rayopt
del s[0]
# set physical size of the offset surface, i.e. the left line in the drawing
s[0].radius = 20  # [mm]
# sets the length between the first virtual offset surface and the lens
s[1].distance = 0  # [mm]
# add parallel plate to the system
s.insert(4,ro.elements.Spheroid(distance=10,material='SCHOTT/N-BK7',
                                                   diameter=spol*2,angles=[utilt,0,0])) 
s.insert(5,ro.elements.Spheroid(distance=T/np.cos(utilt),material='basic/air',
                                                   diameter=spol*2,angles=[utilt,0,0]))
#NOTE:  due to rotation the thickness increases to T/np.cos(utilt) 
#             if this is not done the transversal focus shift displacement
#             does not agree with the theoretical model
s.update()
#s.align(s.paraxial.n) # used by jordens for astigmatic focus shift, destroys rotation

# astigmatic focus shift , can also be obtained from print(q) and looking at table
#print("Astigmatic focus shift "+str(abs(q.waist_position.T[0][-1])-abs(q.waist_position.T[1][-1])))+" mm.")

# Geometric trace
g = ro.GeometricTrace(s)
# In my system, I am only interested in one field
# with a field angle equal to zero radians
# Several distribution can be chosen; hexapolar, random, radau
# The radau scheme should be able to give a good result while using not so many rays
fieldangle=0
g.rays_point((0, fieldangle), wavelength=wavelength, nrays=20, 
                     distribution="radau", filter=False, clip=False)
# Geometric focus [used]
g.refocus()
q = ro.GaussianTrace(s)
if utilt==0:  
    fig, ax = plt.subplots()
    s.plot(ax)
    q.plot(ax, color="red", scale=1)
print("The spot radius is "+str(q.spot_radius[-1][0]*1000))
print("The Gaussian waist radius is "+str(round(q.spot_radius[-1][0]*1000,2))+"   micrometers.")
print("The Rayleigh range is "+str(q.rayleigh_range[-1][0])+ " mm.")
# The geometric RMS spotsize is then calculated at the focal point
# i.e. RMS= <(W-<W>)2>1/2
# on default Rayopt specifies the focal point at the last surface 
# as it sets i=surface equal to -1.
# all rays are given the same "weight"
print("RMS geometric spotsize is "+str(g.rms()*1000)+" micrometers.")
# The focus point distance is measured with respect to the lens
print("The focus point distance from the lens is "+str(g.path[-1]-g.path[3])+" mm.")
print("The transversal displacement is "+str(g.y[-1,-1,1])+" mm.")
p, qq, opd = g.opd(resample=False)
print("The lambda OPD RMS is "+str(np.sqrt((opd**2 * g.w).sum()/g.w.sum())))
#
p = ro.ParaxialTrace(s)
print("The Airy radius is "+ str(p.airy_radius[1]*1000)+" micrometers.")
# paraxial focus [not used]
#s.paraxial.refocus()
ro.Analysis(s,refocus_full=False, update=False)
# Gaussian trace
#   plot only works at ultilt is 0 degrees

# Seidel aberrations
#z = ro.PolyTrace(s)
#str(z)
# Retrieve seidel
#print("\n".join(z.print_seidel()))
