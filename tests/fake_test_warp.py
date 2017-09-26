
from pytileproj import warp
from equi7grid.equi7grid import Equi7Grid

gdal_path = r'C:\OSGeo4W64_qgis218\bin'

image = r"D:\Arbeit\fteaching\2016WS\AF_VU\Uebung1_Schneekartierung\LC81920272016080LGN00_north" \
      r"\LC81920272016080LGN00_B4.tif"
outpath = r"D:\Arbeit\fteaching\2016WS\AF_VU\Uebung1_Schneekartierung"


e7g500 = Equi7Grid(500)

def fake_test_warp():
    pass

    warp.warp2tiledgeotiff(e7g500, image, outpath, gdal_path=gdal_path, subgrid_ids=['EU'],
                           resampling_type='near', image_nodata=0.0, tile_nodata=-9999.0)


    pass

if __name__ == '__main__':
    fake_test_warp()



