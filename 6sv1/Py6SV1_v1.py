from Py6S import *
import os

s = SixS(os.getcwd() + "\\sixsV1_1.exe")
s.atmos_profile = AtmosProfile.UserWaterAndOzone(3.6, 0.9)
# Set the atmosphere profile to be base
s.wavelength = Wavelength (0.56, 0.71) #VIS(0.56, 0.71)
# Set the wavelength to be that of
# Data file Py6S_artigo_ex3.py
s.ground_reflectance = GroundReflectance.HomogeneousMODISBRDF(0.037 , 0.0, 0.133) 
# Set the surf
s.geometry = Geometry.Meteosat()
s.geometry.month = 7
s.geometry.day = 14
s.geometry.year = 2014
s.geometry.gmt_decimal_hour = 7.00
s.geometry.line = 630 # To Goias 630
s.geometry.column = 730 # To Goias 730
s.geometry.solar_z = 30.2610 # Forth data in angles: sunzen >> File 201505111300_angles.tif
s.geometry.solar_a = 83.9880 # Second data in angles: sunaz
s.geometry.view_z = 58.1128  # Third data in angles: satzen
s.geometry.view_a = 104.4975 # First data in angles: sataz
s.run()
print (s.outputs.pixel_radiance)
print (s.outputs.background_radiance)
print (s.outputs.single_scattering_albedo)
print (s.outputs.transmittance_water.downward)
#SixSHelpers.Angles.run_and_plot_360(s, 'view', 'pixel_radiance', colorbarlabel=r"At -sensor Spectral Radiance ($W/m^2\!/\ mu m$)")
#SixSHelpers.Angles.run_and_plot_360(s, 'view', 'pixel_reflectance')
print(s.outputs.values)
