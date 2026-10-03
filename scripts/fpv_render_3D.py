# ParaView (6.0) state for the hinged floating-solar (FPV) array of
# src/tmp_3D_hinged_FPV.jl (Nymphaea Aurora-type, 14 × 14 floaters of 7 m) in
# a short-crested North Sea wind sea.  Same scene set-up as plate_render_3D.py.
#
# Inputs (written by tmp_3D_hinged_FPV.jl in DATA_DIR):
#   solution_κ.pvd   free-surface elevation κ on Γκ (open water + damping zones)
#   solution_η.pvd   array elevation η on Γη, with the per-floater fields `floater`
#                    (floater index) and `checker` (0/1 alternating floaters)
#
# Scene (z up, free surface at z = 0, vertical motions exaggerated by WARP):
#   - water surface warped by κ (PBR_Water)
#   - floaters warped by η, sunk by DRAFT and extruded to a height FLOATER_T,
#     coloured as PV panels with alternating shades per floater
#   - hinge lines as tubes on the deck, cut from the warped array so they
#     follow the relative rotations between floaters
#
# Usage:  File > Load State (Python)  or
#         pvbatch fpv_render_3D.py    (set LRMM_SAVE=1 to write an animation)
# Override the data folder with the LRMM_DATA_DIR environment variable and
# disable ray tracing for a quick preview with LRMM_RAYTRACE=0.

import os

import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 0

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# User settings
# ----------------------------------------------------------------
try:
    _HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:   # 'Load State' in the GUI does not define __file__
    _HERE = os.path.join(os.path.expanduser('~'), 'Progs', 'git_repos', 'CMOE',
                         'HydroElasticFEM_xxx2026.jl', 'scripts')
DATA_DIR = os.environ.get('LRMM_DATA_DIR',
                          os.path.join(_HERE, '..', 'data', 'tmp_3D_hinged_FPV'))

A_FLOATER = 7.0    # floater side length [m] (`a` in tmp_3D_hinged_FPV.jl)
WARP      = 1.0    # vertical exaggeration of η and κ (Hs = 2.5 m: true scale reads well)
FLOATER_T = 0.8    # visual floater height (frame + floats) [m]
DRAFT     = 0.3    # part of the floater below the waterline [m]
HINGE_R   = 0.15   # hinge tube radius [m]
WINDOW    = None   # half-size [m] of a square water window around the array; None = whole tank

# PV-panel colours of the two checkerboard floater families, and hinge colour
PANEL_RGB = ([0.10, 0.16, 0.32], [0.16, 0.25, 0.45])
HINGE_RGB = [0.75, 0.75, 0.78]

# Orbiting camera (one PythonAnimationCue, driven by the normalised animation
# time s ∈ [0, 1]): general view at s = 0 → close-up at s = 0.5 → general view
# at s = 1, while turning N_ORBITS times around the array.
ORBIT_CAMERA = True
N_ORBITS = 1.0     # full turns over the whole animation
AZ0      = 225.0   # start azimuth [deg] (225 = upwave, -y side, looking downwave)
R_FAR    = 500.0   # camera distance in the general view [m] (whole tank in view)
R_NEAR   = 110.0   # camera distance in the close-up [m] (array ≈ 98 m wide)
EL_FAR   = 40.0    # elevation in the general view [deg]
EL_NEAR  = 20.0    # elevation in the close-up [deg]
Z_FOCUS  = 0.0     # array deck
RAY_TRACING = os.environ.get('LRMM_RAYTRACE', '1') != '0'

# HDRI environment (ray tracing only), shared with plate_render_3D.py.
# Override the map with LRMM_HDRI; LRMM_HDRI='' keeps the plain white
# background.  SHOW_SKY = False keeps the HDRI lighting/reflections but hides
# the sky behind the scene (white backplate).
HDRI = os.environ.get('LRMM_HDRI', os.path.join(
    _HERE, '..', 'data', 'tmp_3D', 'animation',
    'kloofendal_48d_partly_cloudy_puresky_4k.hdr'))
SHOW_SKY = True

# ----------------------------------------------------------------
# Readers (the array bounds set the camera target, water window and hinges)
# ----------------------------------------------------------------
kappa_pvd = PVDReader(registrationName='solution_κ.pvd',
                      FileName=os.path.join(DATA_DIR, 'solution_κ.pvd'))
kappa_pvd.Set(PointArrays=['κ'])

eta_pvd = PVDReader(registrationName='solution_η.pvd',
                    FileName=os.path.join(DATA_DIR, 'solution_η.pvd'))
eta_pvd.UpdatePipeline()
xmin, xmax, ymin, ymax, _, _ = eta_pvd.GetDataInformation().GetBounds()
XC, YC = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
NFX = int(round((xmax - xmin) / A_FLOATER))
NFY = int(round((ymax - ymin) / A_FLOATER))
HAS_CHECKER = eta_pvd.PointData.GetArray('checker') is not None

# ----------------------------------------------------------------
# View
# ----------------------------------------------------------------
materialLibrary1 = GetMaterialLibrary()

renderView1 = CreateView('RenderView')
renderView1.Set(
    ViewSize=[2840, 1810],
    OrientationAxesVisibility=0,
    CameraViewUp=[0.0, 0.0, 1.0],
    CameraFocalDisk=1.0,
    UseColorPaletteForBackground=0,
    Background=[1.0, 1.0, 1.0],
)
if RAY_TRACING:
    renderView1.Set(
        EnableRayTracing=1,
        BackEnd='OSPRay pathtracer',
        SamplesPerPixel=10,
        OSPRayMaterialLibrary=materialLibrary1,
    )
    if HDRI and os.path.isfile(HDRI):
        # ParaView 6.0: environment lighting lives in the view Background
        # group (Background = Skybox + Use As Environment Lighting)
        renderView1.Set(
            BackgroundColorMode='Skybox',
            BackgroundTexture=CreateTexture(HDRI),
            UseEnvironmentLighting=1,
            EnvironmentNorth=[0.0, 0.0, 1.0],   # z up: sky above the water
            EnvironmentEast=[1.0, 0.0, 0.0],
        )
        if not SHOW_SKY:
            renderView1.Backgroundmode = 'Both'   # environment light + colour backplate
    elif HDRI:
        print('HDRI not found, using a plain background: ' + HDRI)

layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(2840, 1810)
SetActiveView(renderView1)

# ----------------------------------------------------------------
# Pipelines
# ----------------------------------------------------------------

# Water surface: κ on Γκ (optionally cropped to a window around the array), warped by κ
if WINDOW:
    kappa_clip = Clip(registrationName='WaterWindow', Input=kappa_pvd)
    kappa_clip.ClipType = 'Box'
    kappa_clip.ClipType.Set(
        Position=[XC - WINDOW, YC - WINDOW, -1.0],
        Length=[2 * WINDOW, 2 * WINDOW, 2.0],
    )
    kappa_clip.Invert = 1   # keep the inside of the box
    water_in = kappa_clip
else:
    water_in = kappa_pvd

kappa_surf = ExtractSurface(registrationName='WaterSurface', Input=water_in)
kappa_warp = WarpByScalar(registrationName='WaterWarp', Input=kappa_surf)
kappa_warp.Set(Scalars=['POINTS', 'κ'], ScaleFactor=WARP, Normal=[0.0, 0.0, 1.0])

# Floaters: η on Γη, warped by η, sunk by DRAFT and extruded to FLOATER_T
eta_surf = ExtractSurface(registrationName='ArraySurface', Input=eta_pvd)
eta_warp = WarpByScalar(registrationName='ArrayWarp', Input=eta_surf)
eta_warp.Set(Scalars=['POINTS', 'η'], ScaleFactor=WARP, Normal=[0.0, 0.0, 1.0])
eta_sunk = Transform(registrationName='ArrayDraft', Input=eta_warp)
eta_sunk.Transform.Translate = [0.0, 0.0, -DRAFT]
floaters = LinearExtrusion(registrationName='Floaters', Input=eta_sunk)
floaters.Set(ScaleFactor=FLOATER_T, Vector=[0.0, 0.0, 1.0])

# Hinges: cut the warped array along the interior connection lines x = const
# and y = const, lift the cuts to the deck and draw them as tubes
eta_deck = Transform(registrationName='ArrayDeck', Input=eta_warp)
eta_deck.Transform.Translate = [0.0, 0.0, FLOATER_T - DRAFT + 0.5 * HINGE_R]
hinge_parts = []
for name, normal, nf, origin in (('HingesX', [1.0, 0.0, 0.0], NFX, [xmin, YC, 0.0]),
                                 ('HingesY', [0.0, 1.0, 0.0], NFY, [XC, ymin, 0.0])):
    if nf < 2:
        continue
    cut = Slice(registrationName=name, Input=eta_deck)
    cut.SliceType = 'Plane'
    cut.SliceType.Set(Origin=origin, Normal=normal)
    cut.SliceOffsetValues = [i * A_FLOATER for i in range(1, nf)]
    hinge_parts.append(cut)
hinges = None
if hinge_parts:
    hinge_lines = AppendDatasets(registrationName='HingeLines', Input=hinge_parts)
    hinge_poly = ExtractSurface(registrationName='HingePoly', Input=hinge_lines)
    hinges = Tube(registrationName='Hinges', Input=hinge_poly)
    hinges.Set(Radius=HINGE_R, NumberofSides=12)

# ----------------------------------------------------------------
# Display
# ----------------------------------------------------------------
waterDisplay = Show(kappa_warp, renderView1, 'GeometryRepresentation')
waterDisplay.Set(Representation='Surface', ColorArrayName=[None, ''],
                 DiffuseColor=[0.25, 0.5, 0.75], OSPRayMaterial='PBR_Water')

floaterDisplay = Show(floaters, renderView1, 'GeometryRepresentation')
floaterDisplay.Set(Representation='Surface', Specular=0.6, SpecularPower=60.0)
if HAS_CHECKER:
    ColorBy(floaterDisplay, ('POINTS', 'checker'))
    checkerLUT = GetColorTransferFunction('checker')
    checkerLUT.Set(
        InterpretValuesAsCategories=1,
        Annotations=['0', 'A', '1', 'B'],
        IndexedColors=PANEL_RGB[0] + PANEL_RGB[1],
    )
    floaterDisplay.SetScalarBarVisibility(renderView1, False)
else:
    floaterDisplay.Set(ColorArrayName=[None, ''], DiffuseColor=PANEL_RGB[0])

if hinges is not None:
    hingeDisplay = Show(hinges, renderView1, 'GeometryRepresentation')
    hingeDisplay.Set(Representation='Surface', ColorArrayName=[None, ''],
                     DiffuseColor=HINGE_RGB, Specular=0.8)
    if RAY_TRACING:
        hingeDisplay.OSPRayMaterial = 'Metal_Aluminium_brushed'

# ----------------------------------------------------------------
# Animation: follow the time steps of the pvd files
# ----------------------------------------------------------------
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()
animationScene1.PlayMode = 'Snap To TimeSteps'

# ----------------------------------------------------------------
# Orbiting camera: the same code sets the initial pose and runs on every tick
# ----------------------------------------------------------------
CAMERA_CODE = """
import math
from paraview.simple import GetAnimationScene, GetRenderViews

_F = [{xc}, {yc}, {zf}]

def _pose(s):
    w  = 0.5 * (1.0 - math.cos(2.0 * math.pi * s))             # 0 → 1 → 0
    r  = {r_far} * ({r_near} / {r_far}) ** w                     # geometric zoom
    el = math.radians({el_far} + ({el_near} - {el_far}) * w)
    az = math.radians({az0}) + 2.0 * math.pi * {n_orbits} * s
    return [_F[0] + r * math.cos(el) * math.cos(az),
            _F[1] + r * math.cos(el) * math.sin(az),
            _F[2] + r * math.sin(el)]

def _set_camera(s):
    for v in GetRenderViews():
        v.CameraFocalPoint = _F
        v.CenterOfRotation = _F
        v.CameraPosition = _pose(s)
        v.CameraViewUp = [0.0, 0.0, 1.0]

def _progress():
    sc = GetAnimationScene()
    span = sc.EndTime - sc.StartTime
    return 0.0 if span <= 0 else (sc.AnimationTime - sc.StartTime) / span

def start_cue(self):
    _set_camera(_progress())

def tick(self):
    _set_camera(_progress())

def end_cue(self):
    pass
""".format(xc=XC, yc=YC, zf=Z_FOCUS, r_far=R_FAR, r_near=R_NEAR, el_far=EL_FAR,
           el_near=EL_NEAR, az0=AZ0, n_orbits=N_ORBITS)

_cam_ns = {}
exec(CAMERA_CODE, _cam_ns)
_cam_ns['_set_camera'](0.0)
if ORBIT_CAMERA:
    cameraCue = PythonAnimationCue()
    cameraCue.Script = CAMERA_CODE
    animationScene1.Cues.append(cameraCue)


SetActiveSource(None)

if os.environ.get('LRMM_SAVE'):
    SaveAnimation(os.path.join(DATA_DIR, 'fpv_3D.png'), renderView1,
                  ImageResolution=[1420, 905])
