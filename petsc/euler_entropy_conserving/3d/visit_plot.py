'''
Run this as follows
    visit -cli -nowin -s ./vis_ot.py
'''
import sys
from visit_utils import *


if len(sys.argv) == 1:
    print 'Specify RHO|P|U|V'
    exit()

variable = sys.argv[1]

OpenDatabase("sol-*.plt database")

v = GetView2D()
v.windowCoords = (0.0, 1.0, 0.0, 1.0)
v.viewportCoords = (0.1, 0.85, 0.07, 0.95)
#v.fullFrameActivationMode = v.On
SetView2D(v)

s = SaveWindowAttributes()
s.format = s.POSTSCRIPT                   # PNG or POSTSCRIPT
s.fileName = "sol_"+variable+"_"
s.width, s.height = 1200,1000
s.screenCapture = 0
#s.family = 1
s.resConstraint = s.NoConstraint
#s.advancedMultiWindowSave = 0
SetSaveWindowAttributes(s)

a = AnnotationAttributes()
a.axes2D.visible = 0
a.axes2D.xAxis.title.visible = 0  # (0,1): Disable/enable display of axis title
a.axes2D.yAxis.title.visible = 0
a.userInfoFlag = 0                # disables username display
a.databaseInfoFlag = 0            # diables database display

AddPlot("Pseudocolor", variable)
SetAnnotationAttributes(a)

p = PseudocolorAttributes()
#p.colorTableName = "hot_desaturated"
#p.minFlag, p.maxFlag = 1, 1
#p.min, p.max = -1.0, 1.0
SetPlotOptions(p)
DrawPlots()

# Set legend
plotName = GetPlotList().GetPlots(0).plotName
legend = GetAnnotationObject(plotName)
legend.numberFormat = "%1.2f"
legend.managePosition = 0
legend.position = (0.85,0.80)
legend.fontFamily = legend.Arial
legend.fontBold = 1
legend.fontItalic = 1
legend.xScale = 1.2
legend.yScale = 2.5

# Save images of all timesteps and add each image filename to a list.
for state in range(TimeSliderGetNStates()):
   SetTimeSliderState(state)
   # Save the image
   SaveWindow()

quit()
