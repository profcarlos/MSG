
#-------------------------------------------------------------------------------------------------
Projection facility

Satellite class loader
mpop.satellites is the module englobes all satellite specific modules. In itself, it hold the mighty mpop.satellites.get_satellite_class() method.

class mpop.satellites.GeostationaryFactory

    Factory for geostationary satellite scenes.

    static create_scene(satname, satnumber, instrument, time_slot, area=None, variant='')

        Create a compound satellite scene.

#--------------------------------------------------------------------------------------------------
https://github.com/pytroll/satpy/blob/master/satpy/etc/areas.def

Hi,

Yeah, you can derive your own NDVI and save it geotiff. Is that what you
want?

See http://pytroll.org/quickstart_avhrr.html for channel arithmetics if
needed.

Here is my small sample code using trollimage (master branch) and mpop
in master branch:

from mpop.satellites import GeostationaryFactory
from datetime import datetime
from mpop.projector import get_area_def
import mpop.imageo.geo_image as geo_image

europe = get_area_def("met09globeFull")
tslot = datetime(2014, 7, 1, 12, 00)
glbd = GeostationaryFactory.create_scene("meteosat", "10", "seviri", tslot)
prfx = tslot.strftime('%Y%m%d_%H%M')
glbd.load(area_extent=europe.area_extent)
glbd.area = europe

ndvi = (glbd["VIS008"].data - glbd["VIS006"].data) / \
     (glbd["VIS008"].data + glbd["VIS006"].data)

img = geo_image.GeoImage(
     (ndvi, ndvi, ndvi), glbd.area, tslot, fill_value=(0, 0, 0), mode="RGB")
img.save('ndvi.tif')


Somehow there is a mismatch using mpop pre-master and trollimage master
branch. It fails writing the geotiff. I will talk to Martin about that
next week.

Hope this helps?

-Adam 