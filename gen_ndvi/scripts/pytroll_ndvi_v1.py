######################################
# Required libraries
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
# Import the module to manipulate dates and times
import datetime
# Commom path name manipulations
import os.path
# Used to generate NDVI
#import satpy
#from satpy.dataset import combine_attrs

 
######################################
# Choosing a time slot
######################################
 
# Define the image time slot (year, month, day, hour, minute)
time_slot = datetime.datetime(2017, 9, 03, 12, 00)
global_data = GeostationaryFactory.create_scene("Meteosat-10", "", "seviri", time_slot)
 
######################################
# Choosing a region from the areas.def
######################################
 
# Define the area accorging to the areas.def file
area = get_area_def("met09globeFull")
######################################
# Decompressing HRIT files if required
######################################
 
# Create a list with all SEVIRI bands to decompress
bands = ['VIS006', 'VIS008','IR_016','IR_039','WV_062','WV_073','IR_087','IR_097','IR_108','IR_120','IR_134']
 
# Input and Output directories for the compressed and uncompressed data
indir = '/home/GEONETCAST/pytroll/msg/'
outdir = '/home/GEONETCAST/pytroll/msg/output/'
 
# Loop through all bands
for y in range (0, len(bands)):
    print 'Decompressing ' + bands[y]
    # Loop through all segments for each band
    for x in range (1,9):
        # File to uncompress
        msg_compressed_data = indir + 'H-000-MSG3__-MSG3________-' + bands[y] + '___-00000' + str(x) + '___-' + time_slot.strftime('%Y%m%d%H%M') + '-C_'
        # Check if the uncompressed file for this band and segment exists
        msg_uncompressed_data = indir + 'H-000-MSG3__-MSG3________-' + bands[y] + '___-00000' + str(x) + '___-' + time_slot.strftime('%Y%m%d%H%M') + '-__'
        # Uncompress only if the uncompressed file doesn't exist
        if (os.path.isfile(msg_uncompressed_data) == False):
            fn=xrit.decompress(msg_compressed_data,indir)
######################################
# Visualizing a given channel
######################################
 
# Loading the 0.6 um visible channel
global_data.load([0.6],  area_extent=area.area_extent)
global_data.load([0.8],  area_extent=area.area_extent)
global_data.load([10.8], area_extent=area.area_extent)
# Verify the loaded channels
print '..... global data'
print global_data
print '..... calc NDVI'
global_data['ndvi'] = (global_data[0.8] - global_data[0.6]) / (global_data[0.8] + global_data[0.6])
print '..... Save file'
img = global_data.image.overview()
img.save(outdir + '/myoverview.tif')
