# state file generated using paraview version 6.0.0
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.Set(
    ViewSize=[2840, 1810],
    AxesGrid='Grid Axes 3D Actor',
    OrientationAxesVisibility=0,
    CenterOfRotation=[30.0, 2.0, 5.0],
    CameraPosition=[-41.76523486985973, 25.66681385795997, 57.92621336982935],
    CameraFocalPoint=[17.998571019640842, -8.323171882019706, -1.1966946362036242],
    CameraViewUp=[0.2542478820220702, 0.9269496163351384, -0.2759029236224785],
    CameraFocalDisk=1.0,
    CameraParallelScale=50.30845383820511,
    EnableRayTracing=1,
    BackEnd='OSPRay pathtracer',
    SamplesPerPixel=10,
    OSPRayMaterialLibrary=materialLibrary1,
)

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView3)
layout1.SetSize(2840, 1810)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView3)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
mem_O_solpvd = PVDReader(registrationName='mem_O_sol.pvd', FileName='/Users/ocolomesgene/surfdrive - Oriol Colomes Gene@surfdrive.surf.nl/Research/Tests/LRHS/resTime_resonatingCase/mono_time_lrmm/mem_O_sol.pvd')
mem_O_solpvd.Set(
    CellArrays=['cell'],
    PointArrays=['phi'],
)

# create a new 'Sphere'
sphere1 = Sphere(registrationName='Sphere1')
sphere1.Set(
    Center=[30.0, 2.0, 1.0],
    PhiResolution=50,
)

# create a new 'PVD Reader'
mem_O_solpvd_1 = PVDReader(registrationName='mem_O_sol.pvd', FileName='/Users/ocolomesgene/surfdrive - Oriol Colomes Gene@surfdrive.surf.nl/Research/Tests/LRHS/resTime_resonatingCase/mono_time_lrmm/mem_O_sol.pvd')
mem_O_solpvd_1.Set(
    CellArrays=['cell'],
    PointArrays=['phi', 'q'],
)

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=mem_O_solpvd_1)
calculator1.Set(
    CoordinateResults=1,
    Function='(30)*iHat + (3+1*q)*jHat + 2*kHat',
)

# create a new 'Calculator'
calculator3 = Calculator(registrationName='Calculator3', Input=mem_O_solpvd_1)
calculator3.Set(
    CoordinateResults=1,
    Function='(30)*iHat + (3+1*q)*jHat + 8*kHat',
)

# create a new 'Box'
box1 = Box(registrationName='Box1')
box1.Set(
    XLength=100.0,
    YLength=15.0,
    ZLength=10.0,
    Center=[30.0, -2.5, 5.0],
)

# create a new 'Extract Surface'
extractSurface6 = ExtractSurface(registrationName='ExtractSurface6', Input=mem_O_solpvd_1)

# create a new 'Linear Extrusion'
linearExtrusion6 = LinearExtrusion(registrationName='LinearExtrusion6', Input=extractSurface6)
linearExtrusion6.ScaleFactor = 10.0

# create a new 'PVD Reader'
mem_Ge_solpvd = PVDReader(registrationName='mem_Ge_sol.pvd', FileName='/Users/ocolomesgene/surfdrive - Oriol Colomes Gene@surfdrive.surf.nl/Research/Tests/LRHS/resTime_resonatingCase/mono_time_lrmm/mem_Ge_sol.pvd')
mem_Ge_solpvd.Set(
    CellArrays=['cell'],
    PointArrays=['eta'],
)

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=mem_Ge_solpvd)

# create a new 'Linear Extrusion'
linearExtrusion1 = LinearExtrusion(registrationName='LinearExtrusion1', Input=extractSurface1)
linearExtrusion1.ScaleFactor = 10.0

# create a new 'PVD Reader'
mem_Gk_solpvd = PVDReader(registrationName='mem_Gk_sol.pvd', FileName='/Users/ocolomesgene/surfdrive - Oriol Colomes Gene@surfdrive.surf.nl/Research/Tests/LRHS/resTime_resonatingCase/mono_time_lrmm/mem_Gk_sol.pvd')
mem_Gk_solpvd.Set(
    CellArrays=['cell'],
    PointArrays=['kappa'],
)

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=mem_O_solpvd_1)
calculator2.Set(
    CoordinateResults=1,
    ResultArrayName='Result2',
    Function='(30)*iHat + (3+1*q)*jHat + 5*kHat',
)

# create a new 'Glyph'
glyph2 = Glyph(registrationName='Glyph2', Input=calculator2,
    GlyphType='Sphere')
glyph2.Set(
    OrientationArray=['POINTS', 'No orientation array'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=5.0,
    GlyphTransform='Transform2',
    GlyphMode='All Points',
)

# init the 'Sphere' selected for 'GlyphType'
glyph2.GlyphType.Set(
    Radius=0.2,
    ThetaResolution=15,
    PhiResolution=15,
)

# create a new 'PVD Reader'
mem_Ge_solpvd_1 = PVDReader(registrationName='mem_Ge_sol.pvd', FileName='/Users/ocolomesgene/surfdrive - Oriol Colomes Gene@surfdrive.surf.nl/Research/Tests/LRHS/resTime_resonatingCase/mono_time_lrmm/mem_Ge_sol.pvd')
mem_Ge_solpvd_1.Set(
    CellArrays=['cell'],
    PointArrays=['eta'],
)

# create a new 'Extract Surface'
extractSurface4 = ExtractSurface(registrationName='ExtractSurface4', Input=mem_Ge_solpvd_1)

# create a new 'Extract Surface'
extractSurface2 = ExtractSurface(registrationName='ExtractSurface2', Input=mem_Gk_solpvd)

# create a new 'Linear Extrusion'
linearExtrusion2 = LinearExtrusion(registrationName='LinearExtrusion2', Input=extractSurface2)
linearExtrusion2.ScaleFactor = 10.0

# create a new 'Warp By Scalar'
warpByScalar2 = WarpByScalar(registrationName='WarpByScalar2', Input=linearExtrusion2)
warpByScalar2.Set(
    Scalars=['POINTS', 'kappa'],
    ScaleFactor=10.0,
    Normal=[0.0, 1.0, 0.0],
)

# create a new 'Extract Surface'
extractSurface3 = ExtractSurface(registrationName='ExtractSurface3', Input=mem_O_solpvd)

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=calculator1,
    GlyphType='Sphere')
glyph1.Set(
    OrientationArray=['POINTS', 'No orientation array'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=5.0,
    GlyphTransform='Transform2',
    GlyphMode='All Points',
)

# init the 'Sphere' selected for 'GlyphType'
glyph1.GlyphType.Set(
    Radius=0.2,
    ThetaResolution=15,
    PhiResolution=15,
)

# create a new 'Glyph'
glyph3 = Glyph(registrationName='Glyph3', Input=calculator3,
    GlyphType='Sphere')
glyph3.Set(
    OrientationArray=['POINTS', 'No orientation array'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=5.0,
    GlyphTransform='Transform2',
    GlyphMode='All Points',
)

# init the 'Sphere' selected for 'GlyphType'
glyph3.GlyphType.Set(
    Radius=0.2,
    ThetaResolution=15,
    PhiResolution=15,
)

# create a new 'PVD Reader'
mem_Gk_solpvd_1 = PVDReader(registrationName='mem_Gk_sol.pvd', FileName='/Users/ocolomesgene/surfdrive - Oriol Colomes Gene@surfdrive.surf.nl/Research/Tests/LRHS/resTime_resonatingCase/mono_time_lrmm/mem_Gk_sol.pvd')
mem_Gk_solpvd_1.Set(
    CellArrays=['cell'],
    PointArrays=['kappa'],
)

# create a new 'Linear Extrusion'
linearExtrusion4 = LinearExtrusion(registrationName='LinearExtrusion4', Input=extractSurface4)
linearExtrusion4.ScaleFactor = 10.0

# create a new 'Linear Extrusion'
linearExtrusion7 = LinearExtrusion(registrationName='LinearExtrusion7', Input=linearExtrusion4)
linearExtrusion7.Set(
    ScaleFactor=0.8,
    Vector=[0.0, 1.0, 0.0],
)

# create a new 'Warp By Scalar'
warpByScalar3 = WarpByScalar(registrationName='WarpByScalar3', Input=linearExtrusion7)
warpByScalar3.Set(
    Scalars=['POINTS', 'eta'],
    ScaleFactor=10.0,
    Normal=[0.0, 1.0, 0.0],
)

# create a new 'Warp By Scalar'
warpByScalar1 = WarpByScalar(registrationName='WarpByScalar1', Input=linearExtrusion1)
warpByScalar1.Set(
    Scalars=['POINTS', 'eta'],
    ScaleFactor=10.0,
    Normal=[0.0, 1.0, 0.0],
)

# create a new 'Extract Surface'
extractSurface5 = ExtractSurface(registrationName='ExtractSurface5', Input=mem_Gk_solpvd_1)

# create a new 'Linear Extrusion'
linearExtrusion3 = LinearExtrusion(registrationName='LinearExtrusion3', Input=extractSurface3)
linearExtrusion3.ScaleFactor = 10.0

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=box1)
clip1.ClipType = 'Plane'

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Set(
    Origin=[30.0, 3.0, 5.0],
    Normal=[0.0, 1.0, 0.0],
)

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [30.0, -2.5, 5.0]

# create a new 'Linear Extrusion'
linearExtrusion5 = LinearExtrusion(registrationName='LinearExtrusion5', Input=extractSurface5)
linearExtrusion5.ScaleFactor = 10.0

# create a new 'Warp By Scalar'
warpByScalar4 = WarpByScalar(registrationName='WarpByScalar4', Input=linearExtrusion5)
warpByScalar4.Set(
    Scalars=['POINTS', 'kappa'],
    ScaleFactor=10.0,
    Normal=[0.0, 1.0, 0.0],
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView3'
# ----------------------------------------------------------------

# show data from warpByScalar4
warpByScalar4Display = Show(warpByScalar4, renderView3, 'GeometryRepresentation')

# trace defaults for the display properties.
warpByScalar4Display.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    OSPRayMaterial='PBR_Water',
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
warpByScalar4Display.ScaleTransferFunction.Points = [-0.09867239015793143, 0.0, 0.5, 0.0, 0.07592379303436993, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
warpByScalar4Display.OpacityTransferFunction.Points = [-0.09867239015793143, 0.0, 0.5, 0.0, 0.07592379303436993, 1.0, 0.5, 0.0]

# show data from warpByScalar3
warpByScalar3Display = Show(warpByScalar3, renderView3, 'GeometryRepresentation')

# trace defaults for the display properties.
warpByScalar3Display.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    OSPRayMaterial='PBR_Cork',
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
warpByScalar3Display.ScaleTransferFunction.Points = [-0.046082298946279146, 0.0, 0.5, 0.0, 0.009671677622572018, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
warpByScalar3Display.OpacityTransferFunction.Points = [-0.046082298946279146, 0.0, 0.5, 0.0, 0.009671677622572018, 1.0, 0.5, 0.0]

# show data from clip1
clip1Display = Show(clip1, renderView3, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    Opacity=0.05,
    SelectNormalArray='Normals',
    SelectTCoordArray='TCoords',
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from glyph2
glyph2Display = Show(glyph2, renderView3, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph2Display.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    OSPRayMaterial='Metal_Iron_cast',
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph2Display.ScaleTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph2Display.OpacityTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# show data from glyph1
glyph1Display = Show(glyph1, renderView3, 'GeometryRepresentation')

# get 2D transfer function for 'q'
qTF2D = GetTransferFunction2D('q')

# get color transfer function/color map for 'q'
qLUT = GetColorTransferFunction('q')
qLUT.Set(
    TransferFunction2D=qTF2D,
    RGBPoints=GenerateRGBPoints(
        range_min=-0.49208839845673963,
        range_max=0.030227782204747256,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
glyph1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'q'],
    LookupTable=qLUT,
    SelectNormalArray='Normals',
    OSPRayMaterial='Metal_Iron_cast',
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# show data from glyph3
glyph3Display = Show(glyph3, renderView3, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph3Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'q'],
    LookupTable=qLUT,
    SelectNormalArray='Normals',
    OSPRayMaterial='Metal_Iron_cast',
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph3Display.ScaleTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph3Display.OpacityTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'q'
qPWF = GetOpacityTransferFunction('q')
qPWF.Set(
    Points=[-0.49208839845673963, 0.0, 0.5, 0.0, 0.030227782204747256, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper
timeKeeper1.SuppressedTimeSources = [mem_Ge_solpvd, mem_Gk_solpvd, mem_O_solpvd, mem_Ge_solpvd_1, mem_Gk_solpvd_1, mem_O_solpvd_1]

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=renderView3,
    Cues=timeAnimationCue1,
    AnimationTime=0.0,
    EndTime=78.53981633974607,
    NumberOfFrames=120,
)

# initialize the animation scene

# ----------------------------------------------------------------
# restore active source
SetActiveSource(None)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------