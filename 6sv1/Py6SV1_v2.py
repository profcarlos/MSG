from Py6S import *
import os
print('... test')
s = SixS(os.getcwd() + "\\sixsV1_1.exe")
s.atmos_profile = AtmosProfile.UserWaterAndOzone(3.6, 0.9)
# Set the atmosphere profile to be base
s.wavelength = Wavelength (0.56, 0.71) #In (Schmetz, 2002) VIS06(0.56, 0.71) VIS08(0.74, 0.88)
# Set the wavelength to be that of
# Data file Py6S_artigo_ex3.py
s.ground_reflectance = GroundReflectance.HomogeneousMODISBRDF(0.037 , 0.0, 0.133) 
# Set the surf
s.geometry = Geometry.User()
s.geometry.solar_z = 30.2610 # Forth data in angles: sunzen >> File 201505111300_angles.tif
s.geometry.solar_a = 83.9880 # Second data in angles: sunaz
s.geometry.view_z = 58.1128  # Third data in angles: satzen
s.geometry.view_a = 104.4975 # First data in angles: sataz
s.run()
print (s.outputs.pixel_reflectance)

print(s.outputs.values)
