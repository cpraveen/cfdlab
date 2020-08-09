import numpy as np
from mayavi.mlab import *

def f(x, y):
    sin, cos = np.sin, np.cos
    return sin(x + y) + sin(2 * x - y) + cos(3 * x + 4 * y)

x, y = np.mgrid[-7.:7.05:0.1, -5.:5.05:0.05]

figure(1,fgcolor=(0,0,0),bgcolor=(1,1,1))
s = contour_surf(x, y, f, contours=10, warp_scale=0)
outline()
axes(y_axis_visibility=False)
view(azimuth=0,elevation=0)
title('Contour plot')
#savefig('contour.png',size=(300,300))

figure(2,fgcolor=(0,0,0),bgcolor=(1,1,1))
s = imshow(x, y, f, colormap='jet')
colorbar(orientation='vertical')
outline()
axes(y_axis_visibility=False)
view(azimuth=0,elevation=0)
title('Colour plot')
#savefig('colour.png',size=(300,300))

figure(3,fgcolor=(0,0,0),bgcolor=(1,1,1))
s = surf(x, y, f)
axes()
title('Surface plot')
#savefig('surface.png',size=(300,300))

show()
