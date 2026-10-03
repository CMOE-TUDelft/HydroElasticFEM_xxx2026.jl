# ParaView (6.0) state for the 3D floating plate of src/tmp_3D.jl (no
# resonators) in a short-crested sea.  Same scene as lrmm_render_3D.py
# without the resonator balls and springs, for side-by-side comparison.
#
# Inputs (written by tmp_3D.jl in DATA_DIR):
#   solution_κ.pvd   free-surface elevation κ on Γκ (open water + damping zones)
#   solution_η.pvd   plate elevation η on Γη
#
# Scene (z up, free surface at z = 0, vertical motions exaggerated by WARP):
#   - water surface warped by κ (PBR_Water)
#   - plate warped by η, sunk by DRAFT and extruded to a thickness PLATE_T (PBR_Cork)
#
# Usage:  File > Load State (Python)  or
#         pvbatch plate_render_3D.py    (set LRMM_SAVE=1 to write an animation)
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
                          os.path.join(_HERE, '..', 'data', 'tmp_3D'))

WARP     = 2.0    # vertical exaggeration of η and κ
PLATE_T  = 2.0    # visual plate thickness [m]
DRAFT    = 0.5    # part of the plate below the waterline [m] (d₀ = ρb·hb/ρ ≈ 0.51 m)
WINDOW   = None   # half-size [m] of a square water window around the plate; None = whole tank

# Orbiting camera (one PythonAnimationCue, driven by the normalised animation
# time s ∈ [0, 1]): general view at s = 0 → close-up at s = 0.5 → general view
# at s = 1, while turning N_ORBITS times around the plate.  The focal point
# stays on the plate/resonators.
ORBIT_CAMERA = True
N_ORBITS = 1.0     # full turns over the whole animation
AZ0      = 225.0   # start azimuth [deg] (225 = upwave, -y side, looking downwave)
R_FAR    = 650.0   # camera distance in the general view [m] (whole tank in view)
R_NEAR   = 100.0   # camera distance in the close-up [m]
EL_FAR   = 45.0    # elevation in the general view [deg]
EL_NEAR  = 22.0    # elevation in the close-up [deg]
Z_FOCUS  = 0.0    # plate
RAY_TRACING = os.environ.get('LRMM_RAYTRACE', '1') != '0'

# HDRI environment (ray tracing only): the .hdr map is loaded as a Skybox
# background and used as environment lighting, so the PBR water reflects the
# sky.  Override the map with LRMM_HDRI; LRMM_HDRI='' keeps the plain white
# background.  SHOW_SKY = False keeps the HDRI lighting/reflections but hides
# the sky behind the scene (white backplate).
HDRI = os.environ.get('LRMM_HDRI', os.path.join(
    _HERE, '..', 'data', 'tmp_3D', 'animation',
    'kloofendal_48d_partly_cloudy_puresky_4k.hdr'))
SHOW_SKY = True

# ----------------------------------------------------------------
# Readers (the plate bounds set the camera target and water window)
# ----------------------------------------------------------------
kappa_pvd = PVDReader(registrationName='solution_κ.pvd',
                      FileName=os.path.join(DATA_DIR, 'solution_κ.pvd'))
kappa_pvd.Set(PointArrays=['κ'])

eta_pvd = PVDReader(registrationName='solution_η.pvd',
                    FileName=os.path.join(DATA_DIR, 'solution_η.pvd'))
eta_pvd.Set(PointArrays=['η'])
eta_pvd.UpdatePipeline()
xmin, xmax, ymin, ymax, _, _ = eta_pvd.GetDataInformation().GetBounds()
XC, YC = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)

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

# Water surface: κ on Γκ (optionally cropped to a window around the plate), warped by κ
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

# Plate: η on Γη, warped by η, sunk by DRAFT and extruded to a visual thickness
eta_surf = ExtractSurface(registrationName='PlateSurface', Input=eta_pvd)
eta_warp = WarpByScalar(registrationName='PlateWarp', Input=eta_surf)
eta_warp.Set(Scalars=['POINTS', 'η'], ScaleFactor=WARP, Normal=[0.0, 0.0, 1.0])
eta_sunk = Transform(registrationName='PlateDraft', Input=eta_warp)
eta_sunk.Transform.Translate = [0.0, 0.0, -DRAFT]
plate = LinearExtrusion(registrationName='Plate', Input=eta_sunk)
plate.Set(ScaleFactor=PLATE_T, Vector=[0.0, 0.0, 1.0])

# ----------------------------------------------------------------
# Display
# ----------------------------------------------------------------
waterDisplay = Show(kappa_warp, renderView1, 'GeometryRepresentation')
waterDisplay.Set(Representation='Surface', ColorArrayName=[None, ''],
                 DiffuseColor=[0.25, 0.5, 0.75], OSPRayMaterial='PBR_Water')

plateDisplay = Show(plate, renderView1, 'GeometryRepresentation')
plateDisplay.Set(Representation='Surface', ColorArrayName=[None, ''],
                 DiffuseColor=[0.8, 0.65, 0.45], OSPRayMaterial='PBR_Cork')

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
    SaveAnimation(os.path.join(DATA_DIR, 'plate_3D.png'), renderView1,
                  ImageResolution=[1420, 905])
