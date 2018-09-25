from Py6S import *
import os
print('... test')
s = SixS(os.getcwd() + "\\sixsV1_1.exe")

# Atmospheric data
s.atmos_profile = AtmosProfile.UserWaterAndOzone(3.6819, 0.257)
s.aot550 = 0.0900

# Observation angles
s.altitudes.set_sensor_satellite_level()
s.altitudes.set_target_sea_level()
s.geometry = Geometry.User()
s.geometry.solar_z = 51.9052 # Forth data in angles: sunzen >> File 201505111300_angles.tif
s.geometry.solar_a = 138.0212 # Second data in angles: sunaz
s.geometry.view_z  = 59.2438  # Third data in angles: satzen
s.geometry.view_a  = 103.1633 # First data in angles: sataz


s.wavelength = Wavelength (0.56, 0.71) #In (Schmetz, 2002) VIS06(0.56, 0.71) VIS08(0.74, 0.88)
'''
s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.0938)
#s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(100)

s.run()
print (s.outputs.pixel_reflectance)
print(s.outputs.values)
'''

#s.wavelength = Wavelength (0.74, 0.88) #In (Schmetz, 2002) VIS06(0.56, 0.71) VIS08(0.74, 0.88)
#s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.3360)
s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(60)

s.run()
print(s.outputs.fulltext)
xa = s.outputs.values['coef_xa']
xb = s.outputs.values['coef_xb']
xc = s.outputs.values['coef_xc']
#L = s.outputs.values['pixel_radiance']
L = 41.44
print('xa: %.4f' %xa)
print('xb: %.4f' %xb)
print('xc: %.4f' %xc)
print('L:  %.4f' %L)
y = xa*L-xb
print('y: %.4f' %y)
acr = y/(1+xc*y)
ref = 0.0988
print('acr: %.4f ref: %.4f' %(acr, ref))
#print(s.outputs.pixel_radiance)
#print(s.outputs.environmental_irradiance)
#print(s.outputs.total_gaseous_transmittance)
