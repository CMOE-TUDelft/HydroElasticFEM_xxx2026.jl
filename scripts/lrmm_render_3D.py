# ParaView (6.0) state for the 3D locally resonant metaplate of
# src/tmp_3D_resonator.jl: floating plate + array of point resonators in a
# short-crested sea.
#
# Inputs (written by tmp_3D_resonator.jl in DATA_DIR):
#   solution_κ.pvd   free-surface elevation κ on Γκ (open water + damping zones)
#   solution_η.pvd   plate elevation η on Γη, plus two N-component fields that
#                    are constant over the plate: q (resonator displacement)
#                    and eta_r (plate displacement at each attachment point)
#   resonators.csv   x, y, M, K, C, omega_r of each resonator (row i ↔ component i)
#
# Scene (z up, free surface at z = 0, vertical motions exaggerated by WARP):
#   - water surface warped by κ (PBR_Water)
#   - plate warped by η, sunk by DRAFT and extruded to a thickness PLATE_T (PBR_Cork)
#   - one spring (tube) per resonator from the plate top to its mass
#   - one ball per resonator at z = Z_RES + WARP·qᵢ, coloured by qᵢ; the radius
#     maps the natural frequency ωrᵢ linearly from R_MAX (lowest ωr) to R_MIN
#     (highest ωr), i.e. bigger ball = slower resonator
#
# Usage:  File > Load State (Python)  or
#         pvbatch lrmm_render_3D.py    (set LRMM_SAVE=1 to write an animation)
# Override the data folder with the LRMM_DATA_DIR environment variable.

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
                          os.path.join(_HERE, '..', 'data', 'tmp_3D_resonator'))
CSV = os.path.join(DATA_DIR, 'resonators.csv')

WARP     = 2.0    # vertical exaggeration of η, κ and q
PLATE_T  = 2.0    # visual plate thickness [m]
DRAFT    = 0.5    # part of the plate below the waterline [m] (d₀ = ρb·hb/ρ ≈ 0.51 m)
Z_TOP    = PLATE_T - DRAFT      # plate top at rest
SPRING_L = 6.0    # spring rest length above the plate top [m]
Z_RES    = Z_TOP + SPRING_L     # ball centre at rest
R_MAX    = 2.5    # ball radius at the lowest ωr [m] (< half the resonator spacing)
R_MIN    = 1.0    # ball radius at the highest ωr [m]
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
Z_FOCUS  = Z_RES / 2    # between plate top and balls
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

# Resonator centroid → camera target and water window
with open(CSV) as f:
    rows = [list(map(float, l.split(','))) for l in f.readlines()[1:] if l.strip()]
XC = sum(r[0] for r in rows) / len(rows)
YC = sum(r[1] for r in rows) / len(rows)

# ----------------------------------------------------------------
# Programmable-filter scripts: build resonator balls / springs from the
# N-component fields q and eta_r of the plate output and resonators.csv
# ----------------------------------------------------------------
_COMMON = """
import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray

res = np.loadtxt(r'{csv}', delimiter=',', skiprows=1, ndmin=2)
x, y, M = res[:, 0], res[:, 1], res[:, 2]
# natural frequency: omega_r column if present, else sqrt(K/M)
wr = res[:, 5] if res.shape[1] > 5 else np.sqrt(res[:, 3] / M)
N = len(x)
# q and eta_r are constant over the plate: take the first point
q  = np.asarray(inputs[0].PointData['q']).reshape(-1, N)[0]
er = np.asarray(inputs[0].PointData['eta_r']).reshape(-1, N)[0]
out = self.GetPolyDataOutput()

def set_points(xyz):
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(xyz, dtype=float), deep=1))
    out.SetPoints(pts)

def add_array(name, a):
    arr = numpy_to_vtk(np.ascontiguousarray(a, dtype=float), deep=1)
    arr.SetName(name)
    out.GetPointData().AddArray(arr)
""".format(csv=CSV)

BALLS_SCRIPT = _COMMON + """
set_points(np.column_stack([x, y, {z_res} + {warp} * q]))
verts = vtkCellArray()
for i in range(N):
    verts.InsertNextCell(1)
    verts.InsertCellPoint(i)
out.SetVerts(verts)
add_array('q', q)
add_array('M', M)
add_array('omega_r', wr)
span = wr.max() - wr.min()
s = (wr.max() - wr) / span if span > 0 else np.ones(N)
add_array('radius', {r_min} + ({r_max} - {r_min}) * s)
""".format(z_res=Z_RES, warp=WARP, r_max=R_MAX, r_min=R_MIN)

SPRINGS_SCRIPT = _COMMON + """
bottom = np.column_stack([x, y, {z_top} + {warp} * er])
top    = np.column_stack([x, y, {z_res} + {warp} * q])
set_points(np.vstack([bottom, top]))
lines = vtkCellArray()
for i in range(N):
    lines.InsertNextCell(2)
    lines.InsertCellPoint(i)
    lines.InsertCellPoint(N + i)
out.SetLines(lines)
add_array('stretch', np.concatenate([q - er, q - er]))
""".format(z_top=Z_TOP, warp=WARP, z_res=Z_RES)

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
kappa_pvd = PVDReader(registrationName='solution_κ.pvd',
                      FileName=os.path.join(DATA_DIR, 'solution_κ.pvd'))
kappa_pvd.Set(PointArrays=['κ'])

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
eta_pvd = PVDReader(registrationName='solution_η.pvd',
                    FileName=os.path.join(DATA_DIR, 'solution_η.pvd'))
eta_pvd.Set(PointArrays=['η', 'q', 'eta_r'])

eta_surf = ExtractSurface(registrationName='PlateSurface', Input=eta_pvd)
eta_warp = WarpByScalar(registrationName='PlateWarp', Input=eta_surf)
eta_warp.Set(Scalars=['POINTS', 'η'], ScaleFactor=WARP, Normal=[0.0, 0.0, 1.0])
eta_sunk = Transform(registrationName='PlateDraft', Input=eta_warp)
eta_sunk.Transform.Translate = [0.0, 0.0, -DRAFT]
plate = LinearExtrusion(registrationName='Plate', Input=eta_sunk)
plate.Set(ScaleFactor=PLATE_T, Vector=[0.0, 0.0, 1.0])

# Resonator masses: one point per resonator → sphere glyph scaled by radius
balls_pts = ProgrammableFilter(registrationName='ResonatorPoints', Input=eta_pvd)
balls_pts.Set(OutputDataSetType='vtkPolyData', Script=BALLS_SCRIPT)

balls = Glyph(registrationName='ResonatorBalls', Input=balls_pts, GlyphType='Sphere')
balls.Set(
    OrientationArray=['POINTS', 'No orientation array'],
    ScaleArray=['POINTS', 'radius'],
    ScaleFactor=1.0,
    GlyphMode='All Points',
)
balls.GlyphType.Set(Radius=1.0, ThetaResolution=24, PhiResolution=24)

# Springs: plate top at the attachment point → ball centre
springs_lines = ProgrammableFilter(registrationName='ResonatorSpringLines', Input=eta_pvd)
springs_lines.Set(OutputDataSetType='vtkPolyData', Script=SPRINGS_SCRIPT)
springs = Tube(registrationName='ResonatorSprings', Input=springs_lines)
springs.Set(Radius=0.25, NumberofSides=12)

# ----------------------------------------------------------------
# Display
# ----------------------------------------------------------------
waterDisplay = Show(kappa_warp, renderView1, 'GeometryRepresentation')
waterDisplay.Set(Representation='Surface', ColorArrayName=[None, ''],
                 DiffuseColor=[0.25, 0.5, 0.75], OSPRayMaterial='PBR_Water')

plateDisplay = Show(plate, renderView1, 'GeometryRepresentation')
plateDisplay.Set(Representation='Surface', ColorArrayName=[None, ''],
                 DiffuseColor=[0.8, 0.65, 0.45], OSPRayMaterial='PBR_Cork')

springsDisplay = Show(springs, renderView1, 'GeometryRepresentation')
springsDisplay.Set(Representation='Surface', ColorArrayName=[None, ''],
                   DiffuseColor=[0.6, 0.6, 0.6], OSPRayMaterial='Metal_Iron_cast')

qLUT = GetColorTransferFunction('q')
qLUT.ApplyPreset('Cool to Warm', True)
qLUT.AutomaticRescaleRangeMode = 'Clamp and update every timestep'

ballsDisplay = Show(balls, renderView1, 'GeometryRepresentation')
ballsDisplay.Set(Representation='Surface', ColorArrayName=['POINTS', 'q'],
                 LookupTable=qLUT, OSPRayMaterial='Metal_Iron_cast')

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
    SaveAnimation(os.path.join(DATA_DIR, 'lrmm_3D.png'), renderView1,
                  ImageResolution=[1420, 905])
