from Py6S import *
import os

#s = SixS('C:\\Users\\carlos.silveira\\Dropbox\\newPython\\brdf\\sixsV1_1_dell.exe')
s = SixS('C:\\Users\\carlos.silveira\\Dropbox\\newPython\\6SV1\\sixsV1_1_lab.exe')
s.produce_debug_report()

#classmethod UserWaterAndOzone(water, ozone)
#Set 6S to use an atmosphere defined by an amount of water vapour and ozone.
#Arguments:
#•water – The total amount of water in a vertical path through the atmosphere (in g/cm^2)
#•ozone – The total amount of ozone in a vertical path through the atmosphere (in cm-atm)
s.atmosprofile = AtmosProfile.UserWaterAndOzone(0.1, 0.02)
s.aot550 = 0.01

#class Py6S.Wavelength
#Select one or more wavelengths for the 6S simulation.
s.wavelength = Wavelength(0.56, 0.71)
print('wavelength:')
print(s.wavelength)

#class Py6S.AtmosCorr
#Class representing options for selecting atmospheric correction settings for 6S.
#classmethod AtmosCorrBRDFFromReflectance(reflectance)
#Set 6S to perform atmospheric correction using a fully BRDF-represented surface, using a given reflectance value.
s.atmos_corr = AtmosCorr.AtmosCorrBRDFFromReflectance(0.4)
print('atmos_corr: ')
print(s.atmos_corr)

#classmethod HomogeneousMODISBRDF(par1, par2, par3)
#Parameterisation for a surface BRDF based on the MODIS Operational BRDF model.
s.ground_reflectance = GroundReflectance.HomogeneousMODISBRDF(0.4, 0.001, 0.05)
s.ground_altitude = 1.2     # Not change results!
s.ground_pressure = 1000
print('ground_reflectance:')
print(s.ground_reflectance)
#class Geometry.User Stores parameters for a user-defined geometry for 6S.
#Attributes:
#    •solar_z – Solar zenith angle
#    •solar_a – Solar azimuth angle
#    •view_z – View zenith angle •view_a – View azimuth angle
#    •day – The day the image was acquired in (1-31)
#    •month – The month the image was acquired in (0-12)

s.geometry = Geometry.User()
s.geometry.solar_z = 30.2610 # Forth data in angles: sunzen >> File 201505111300_angles.tif
s.geometry.solar_a = 83.9880 # Second data in angles: sunaz
s.geometry.view_z = 58.1128  # Third data in angles: satzen
s.geometry.view_a = 104.4975 # First data in angles: sataz

#class Py6S.AeroProfile Class representing options for Aerosol Profile

s.aeroprofile = AeroProfile.UserProfile(AeroProfile.Continental)
s.aeroprofile.add_layer(5, 0.34) # Add a 5km-thick layer with an AOT of 0.34
s.aeroprofile.add_layer(10, 0.7) # Add a 10km-thick layer with an AOT of 0.7
s.aeroprofile.add_layer(100, 0.01) # Add a 100km-thick layer with an AOT of 0.01

#print(s.aeroprofile.values)


s.run()
print("----")
print( "aot550" + "\t" + str(s.outputs.aot550))
print( "apparent_radiance" + "\t" + str(s.outputs.apparent_radiance))
print( "apparent_reflectance" + "\t" + str(s.outputs.apparent_reflectance))
print( "atmos_corrected_reflectance_brdf" + "\t" + str(s.outputs.atmos_corrected_reflectance_brdf))
print( "atmos_corrected_reflectance_lambertian" + "\t" + str(s.outputs.atmos_corrected_reflectance_lambertian))
print( "atmospheric_intrinsic_radiance" + "\t" + str(s.outputs.atmospheric_intrinsic_radiance))
print( "atmospheric_intrinsic_reflectance" + "\t" + str(s.outputs.atmospheric_intrinsic_reflectance))
print( "azimuthal_angle_difference" + "\t" + str(s.outputs.azimuthal_angle_difference))
print( "background_radiance" + "\t" + str(s.outputs.background_radiance))
print( "background_reflectance" + "\t" + str(s.outputs.background_reflectance))
print( "coef_xa" + "\t" + str(s.outputs.coef_xa))
print( "coef_xb" + "\t" + str(s.outputs.coef_xb))
print( "coef_xc" + "\t" + str(s.outputs.coef_xc))
print( "day" + "\t" + str(s.outputs.day))
print( "diffuse_solar_irradiance" + "\t" + str(s.outputs.diffuse_solar_irradiance))
print( "diffuse_solar_irradiance" + "\t" + str(s.outputs.diffuse_solar_irradiance))
print( "environmental_irradiance" + "\t" + str(s.outputs.environmental_irradiance))
print( "ground_altitude" + "\t" + str(s.outputs.ground_altitude))
print( "ground_pressure" + "\t" + str(s.outputs.ground_pressure))
print( "int_funct_filt" + "\t" + str(s.outputs.int_funct_filt))
print( "int_solar_spectrum" + "\t" + str(s.outputs.int_solar_spectrum))
print( "measured_radiance" + "\t" + str(s.outputs.measured_radiance))
print( "month" + "\t" + str(s.outputs.month))
print( "percent_diffuse_solar_irradiance" + "\t" + str(s.outputs.percent_diffuse_solar_irradiance))
print( "percent_direct_solar_irradiance" + "\t" + str(s.outputs.percent_direct_solar_irradiance))
print( "percent_environmental_irradiance" + "\t" + str(s.outputs.percent_environmental_irradiance))
print( "pixel_radiance" + "\t" + str(s.outputs.pixel_radiance))
print( "pixel_reflectance" + "\t" + str(s.outputs.pixel_reflectance))
print( "scattering_angle" + "\t" + str(s.outputs.scattering_angle))
print( "solar_a" + "\t" + str(s.outputs.solar_a))
print( "solar_z" + "\t" + str(s.outputs.solar_z))
print( "total_gaseous_transmittance" + "\t" + str(s.outputs.total_gaseous_transmittance))
print( "view_a" + "\t" + str(s.outputs.view_a))
print( "view_z" + "\t" + str(s.outputs.view_z))
print( "visibility" + "\t" + str(s.outputs.visibility))
print( "wv_above_aerosol" + "\t" + str(s.outputs.wv_above_aerosol))
print( "wv_mixed_with_aerosol" + "\t" + str(s.outputs.wv_mixed_with_aerosol))
print( "wv_under_aerosol" + "\t" + str(s.outputs.wv_under_aerosol))





