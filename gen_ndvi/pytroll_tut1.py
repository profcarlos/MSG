######################################
# GNC-A Blog Pytroll Example
######################################
 
######################################
# Required Libraries
######################################
 
# Turn on the Debug mode
from mpop.utils import debug_on
debug_on()
# Import the mpop GEO factory
from mpop.satellites import GeostationaryFactory
# Import the mpop area definition module
from mpop.projector import get_area_def
# Import the XRIT decompression utility module
from mipp import xrit
# Import the Image module from the Python Imaging Library
from PIL import Image
# Use pycoast to add shapefiles
from pycoast import ContourWriter
# Import the module to manipulate dates and times
import datetime
# Commom path name manipulations
import os.path
 
######################################
# Choosing the HRIT file
######################################
 
# Define the image time slot (year, month, day, hour, minute)
time_slot = datetime.datetime(2017, 8, 30, 12, 00)
global_data = GeostationaryFactory.create_scene("Meteosat-10", "", "seviri", time_slot)
 
# Define the area accorging to the areas.def file
area = get_area_def("met09globeFull")
 
######################################
# Decompressing HRIT files if required
######################################
 
# Create a list with all SEVIRI bands
bands = ['VIS006', 'VIS008','IR_016','IR_039','WV_062','WV_073','IR_087','IR_097','IR_108','IR_120','IR_134']
 
# Input and Output directories for the compressed and uncompressed data
indir = '/home/pytroll/msg/'
outdir = '/home/pytroll/msg/'
 
# Loop through all bands
for y in range (0, len(bands)):
   print 'Decompressing ' + bands[y]
   # Loop through all segments
   for each band for x in range (1,9):
      # File to uncompress
      msg_compressed_data = indir + 'H-000-MSG3__-MSG3________-' + bands[y] + '___-00000' + str(x) + '___-201708301200-C_'
      # Check if the uncompressed file for this band and segment exists
      msg_uncompressed_data = indir + 'H-000-MSG3__-MSG3________-' + bands[y] + '___-00000' + str(x) + '___-201708301200-__'
      # Uncompress only if the uncompressed file doesn't exist
      if (os.path.isfile(msg_uncompressed_data) == False):
            fn=xrit.decompress(msg_compressed_data,outdir)
 
######################################
# Plotting RGB composites
######################################
 
cw = ContourWriter('/local_disk/data/shapes')
 
######################################
# FOG
######################################
 
print 'Creating the fog RGB'
global_data.load(global_data.image.fog.prerequisites, area_extent=area.area_extent)
img = global_data.image.fog()
img.save("/home/pytroll/output/fog.png")
img = Image.open('/home/pytroll/output/fog.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/fog.png")
 
######################################
# AIRMASS
######################################
 
print 'Creating the airmass RGB'
global_data.load(global_data.image.airmass.prerequisites, area_extent=area.area_extent)
img = global_data.image.airmass()
img.save("/home/pytroll/output/airmass.png")
img = Image.open('/home/pytroll/output/airmass.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/airmass.png")
 
######################################
# GREEN SNOW
######################################
 
print 'Creating the green_snow RGB'
global_data.load(global_data.image.green_snow.prerequisites, area_extent=area.area_extent)
img = global_data.image.green_snow()
img.save("/home/pytroll/output/green_snow.png")
img = Image.open('/home/pytroll/output/green_snow.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/green_snow.png")
 
######################################
# DUST
######################################
 
print 'Creating the dust RGB'
global_data.load(global_data.image.dust.prerequisites, area_extent=area.area_extent)
img = global_data.image.dust()
img.save("/home/pytroll/output/dust.png")
img = Image.open('/home/pytroll/output/dust.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/dust.png")
 
######################################
# ASH
######################################
 
print 'Creating the ash RGB'
global_data.load(global_data.image.ash.prerequisites, area_extent=area.area_extent)
img = global_data.image.ash()
img.save("/home/pytroll/output/ash.png")
img = Image.open('/home/pytroll/output/ash.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/ash.png")
 
######################################
# NATURAL COLORS
######################################
 
print 'Creating the natural RGB'
global_data.load(global_data.image.natural.prerequisites, area_extent=area.area_extent)
img = global_data.image.natural()
img.save("/home/pytroll/output/natural.png")
img = Image.open('/home/pytroll/output/natural.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/natural.png")
 
######################################
# OVERVIEW
######################################
 
print 'Creating the overview RGB'
global_data.load(global_data.image.overview.prerequisites, area_extent=area.area_extent)
img = global_data.image.overview()
img.save("/home/pytroll/output/overview.png")
img = Image.open('/home/pytroll/output/overview.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/overview.png")
 
######################################
# CONVECTION
######################################
 
print 'Creating the convection_co2 RGB'
global_data.load(global_data.image.convection_co2.prerequisites, area_extent=area.area_extent)
img = global_data.image.convection_co2()
img.save("/home/pytroll/output/convection_co2.png")
img = Image.open('/home/pytroll/output/convection_co2.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/convection_co2.png")
 
######################################
# DAY MICROPHYSICS
######################################
 
print 'Creating the Day Microphysics RGB'
global_data.load(global_data.image.day_microphysics.prerequisites, area_extent=area.area_extent)
img = global_data.image.day_microphysics()
img.save("/home/pytroll/output/day_microphysics.png")
img = Image.open('/home/pytroll/output/day_microphysics.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/day_microphysics.png")
 
######################################
# CLOUDTOP
######################################
 
print 'Creating the cloudtop RGB'
global_data.load(global_data.image.cloudtop.prerequisites, area_extent=area.area_extent)
img = global_data.image.cloudtop()
img.save("/home/pytroll/output/cloudtop.png")
img = Image.open('/home/pytroll/output/cloudtop.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/cloudtop.png")
 
######################################
# NIGHT OVERVIEW
######################################
 
print 'Creating the night_overview RGB'
global_data.load(global_data.image.night_overview.prerequisites, area_extent=area.area_extent)
img = global_data.image.night_overview()
img.save("/home/pytroll/output/night_overview.png")
img = Image.open('/home/pytroll/output/night_overview.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/night_overview.png")
 
######################################
# NIGHT FOG
######################################
 
print 'Creating the night_fog RGB'
global_data.load(global_data.image.night_fog.prerequisites, area_extent=area.area_extent)
img = global_data.image.night_fog()
img.save("/home/pytroll/output/night_fog.png")
img = Image.open('/home/pytroll/output/night_fog.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/night_fog.png")
 
######################################
# NIGHT MICROPHYSICS
######################################
 
print 'Creating the night_microphysics RGB'
global_data.load(global_data.image.night_microphysics.prerequisites, area_extent=area.area_extent)
img = global_data.image.night_microphysics()
img.save("/home/pytroll/output/night_microphysics.png")
img = Image.open('/home/pytroll/output/night_microphysics.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/night_microphysics.png")
 
######################################
# SNOW
######################################
 
print 'Creating the snow RGB'
global_data.load(global_data.image.snow.prerequisites, area_extent=area.area_extent)
img = global_data.image.snow()
img.save("/home/pytroll/output/snow.png")
img = Image.open('/home/pytroll/output/snow.png')
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/ne_10m_admin_0_countries.shp",outline="white")
cw.add_shapefile_shapes(img,area,"/home/pytroll/shapefiles/BRA_adm1.shp",outline="white")
img.save("/home/pytroll/output/snow.png")