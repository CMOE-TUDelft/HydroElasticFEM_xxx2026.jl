# ParaView (6.0) state: maximum-over-time elevation along a line through the
# plate, comparing the plate without resonators (src/tmp_3D.jl) and the
# locally resonant metaplate (src/tmp_3D_resonator.jl) in the same
# short-crested sea (same seed).
#
# Inputs:
#   PLATE_DIR/solution_κ.pvd, solution_η.pvd   plate only (κ also has ηin, the
#                                               incident elevation)
#   LRMM_DIR/solution_κ.pvd,  solution_η.pvd    plate + resonators
#
# For each case and field (κ on the free surface Γκ, η on the plate Γη):
#   reader → ExtractTimeSteps (t ≥ T_START) → PlotOverLine (P1 → P2)
#          → TemporalStatistics (max, min) → mask points off the field's
#            surface with NaN (κ has a hole under the plate, η lives only on it)
# so the statistics are exact at every sample point of the line.  The line
# chart breaks the curves at the NaNs: κ is drawn up- and downwave of the
# plate, η over it.
#
# Usage:  File > Load State (Python)  or
#         pvbatch lrmm_max_elevation_lines.py
#         (set LRMM_SAVE=1 to write <OUT_DIR>/max_elevation_line.{png,pdf}
#          and one csv per curve for pgfplots)
# Override the data folders with PLATE_DATA_DIR and LRMM_DATA_DIR.

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
_DATA = os.path.join(_HERE, '..', 'data')
PLATE_DIR = os.environ.get('PLATE_DATA_DIR', os.path.join(_DATA, 'tmp_3D'))
LRMM_DIR  = os.environ.get('LRMM_DATA_DIR', os.path.join(_DATA, 'tmp_3D_resonator'))
OUT_DIR   = os.path.join(_DATA, 'tmp_3D_comparison')

T_START = 50.0    # discard the initial transient [s]

# Sampling line (free surface z = 0).  Default: tank centreline y = 0 over the
# physical region between the inlet (x < 230 m) and outlet (x > 570 m)
# damping zones; the plate spans x ∈ [370, 430] m, y ∈ [-25, 25] m.
P1 = [230.0, 0.0, 0.0]
P2 = [570.0, 0.0, 0.0]
DS = 0.5          # sample spacing [m] (keep (P2-P1)/DS integer so the plate edges are sampled)
X_ARRAY = 'Points_X'   # abscissa: 'Points_X' (line along x), 'Points_Y' or 'arc_length'
X_LABEL = 'x [m]'

SHOW_MIN      = False  # also plot the minimum over time (troughs), dashed
SHOW_INCIDENT = True   # max of the incident elevation ηin (no plate), dotted

# Curve colours (Okabe–Ito): no resonators / resonators / incident
C_PLATE = [0.0, 0.0, 0.0]
C_LRMM  = [0.0, 0.447, 0.698]
C_INC   = [0.5, 0.5, 0.5]

# ----------------------------------------------------------------
# Pipeline builders
# ----------------------------------------------------------------
_L = sum((b - a) ** 2 for a, b in zip(P1, P2)) ** 0.5
RESOLUTION = max(1, int(round(_L / DS)))

MASK_SCRIPT = """
import numpy as np
inp = inputs[0]
output.ShallowCopy(inp.VTKObject)
pd = inp.PointData
mk = 'vtkValidPointMask_minimum'   # 1 only if on the surface at every step
valid = np.asarray(pd[mk]) > 0 if mk in pd.keys() else np.ones(inp.GetNumberOfPoints(), bool)
for src, dst in {pairs}:
    if src in pd.keys():
        a = np.array(pd[src], dtype=float)
        a[~valid] = np.nan
        output.PointData.append(a, dst)
"""


def max_over_line(data_dir, field, arrays, tag):
    """Max/min over t ≥ T_START of `arrays` of solution_<field>.pvd along P1→P2.
    Returns a source with arrays <name>_max / <name>_min (NaN off the surface)."""
    reader = PVDReader(registrationName='%s_%s.pvd' % (tag, field),
                       FileName=os.path.join(data_dir, 'solution_%s.pvd' % field))
    reader.Set(PointArrays=list(arrays))
    reader.UpdatePipelineInformation()
    times = list(reader.TimestepValues)
    i0 = next((i for i, t in enumerate(times) if t >= T_START - 1e-9), len(times) - 1)

    steps = ExtractTimeSteps(registrationName='%s_%s_t>%g' % (tag, field, T_START),
                             Input=reader)
    steps.Set(SelectionMode='Select Time Range',
              TimeStepRange=[i0, len(times) - 1], TimeStepInterval=1)

    line = PlotOverLine(registrationName='%s_%s_line' % (tag, field), Input=steps)
    line.Set(Point1=P1, Point2=P2, Resolution=RESOLUTION, SamplingPattern='Sample Uniformly')

    stats = TemporalStatistics(registrationName='%s_%s_stats' % (tag, field), Input=line)
    stats.Set(ComputeAverage=0, ComputeStandardDeviation=0,
              ComputeMaximum=1, ComputeMinimum=1)

    pairs = [('%s_%s' % (a, s), '%s_%s' % (a, s[:3]))
             for a in arrays for s in ('maximum', 'minimum')]
    masked = ProgrammableFilter(registrationName='%s_%s_max' % (tag, field), Input=stats)
    masked.Set(OutputDataSetType='Same as Input', Script=MASK_SCRIPT.format(pairs=pairs))
    print('%s %s: %d steps in t ∈ [%g, %g] s' % (tag, field, len(times) - i0,
                                                 times[i0], times[-1]))
    return masked


plate_k = max_over_line(PLATE_DIR, 'κ', ['κ', 'ηin'] if SHOW_INCIDENT else ['κ'], 'plate')
plate_e = max_over_line(PLATE_DIR, 'η', ['η'], 'plate')
lrmm_k  = max_over_line(LRMM_DIR, 'κ', ['κ'], 'lrmm')
lrmm_e  = max_over_line(LRMM_DIR, 'η', ['η'], 'lrmm')

# ----------------------------------------------------------------
# Line chart
# ----------------------------------------------------------------
chart = CreateView('XYChartView')
chart.Set(
    ViewSize=[1200, 600],
    ChartTitle='',
    ShowLegend=1,
    LegendLocation='TopRight',
    LeftAxisTitle='max elevation over t > %g s [m]' % T_START,
    BottomAxisTitle=X_LABEL,
    LeftAxisTitleFontSize=18, BottomAxisTitleFontSize=18,
    LeftAxisLabelFontSize=16, BottomAxisLabelFontSize=16,
    LegendFontSize=15,
    BottomAxisUseCustomRange=1,
    BottomAxisRangeMinimum=P1[0] if X_ARRAY == 'Points_X' else 0.0,
    BottomAxisRangeMaximum=P2[0] if X_ARRAY == 'Points_X' else _L,
    LeftAxisGridColor=[0.9, 0.9, 0.9],
    BottomAxisGridColor=[0.9, 0.9, 0.9],
)
if X_ARRAY == 'Points_Y':
    chart.Set(BottomAxisRangeMinimum=P1[1], BottomAxisRangeMaximum=P2[1])

layout1 = CreateLayout(name='Max elevation')
layout1.AssignView(0, chart)
layout1.SetSize(1200, 600)
SetActiveView(chart)

_SOLID, _DASH, _DOT = '1', '2', '3'


def show_curves(src, curves):
    """curves: list of (array, label, rgb, width, style)."""
    d = Show(src, chart, 'XYChartRepresentation')
    d.Set(UseIndexForXAxis=0, XArrayName=X_ARRAY,
          SeriesVisibility=[c[0] for c in curves])
    d.SeriesLabel = [v for c in curves for v in (c[0], c[1])]
    d.SeriesColor = [v for c in curves for v in [c[0]] + ['%g' % x for x in c[2]]]
    d.SeriesLineThickness = [v for c in curves for v in (c[0], str(c[3]))]
    d.SeriesLineStyle = [v for c in curves for v in (c[0], c[4])]
    d.SeriesMarkerStyle = [v for c in curves for v in (c[0], '0')]
    return d


def case_curves(k_name, e_name, name, rgb):
    # thin: free surface κ, thick: plate η (the chart font has no Greek glyphs)
    c = [(k_name + '_max', 'free surface, ' + name, rgb, 2, _SOLID),
         (e_name + '_max', 'plate, ' + name, rgb, 4, _SOLID)]
    if SHOW_MIN:
        c += [(k_name + '_min', 'free surface min, ' + name, rgb, 2, _DASH),
              (e_name + '_min', 'plate min, ' + name, rgb, 4, _DASH)]
    return c


_pc = case_curves('κ', 'η', 'no resonators', C_PLATE)
_lc = case_curves('κ', 'η', 'resonators', C_LRMM)
plate_k_disp = show_curves(plate_k, [c for c in _pc if c[0].startswith('κ')] +
                           ([('ηin_max', 'incident wave', C_INC, 2, _DOT)]
                            if SHOW_INCIDENT else []))
plate_e_disp = show_curves(plate_e, [c for c in _pc if c[0].startswith('η')])
lrmm_k_disp  = show_curves(lrmm_k, [c for c in _lc if c[0].startswith('κ')])
lrmm_e_disp  = show_curves(lrmm_e, [c for c in _lc if c[0].startswith('η')])

SetActiveSource(None)

if os.environ.get('LRMM_SAVE'):
    os.makedirs(OUT_DIR, exist_ok=True)
    base = os.path.join(OUT_DIR, 'max_elevation_line')
    SaveScreenshot(base + '.png', chart, ImageResolution=[1200, 600])
    ExportView(base + '.pdf', view=chart)
    # one csv per curve: x, y, s (arc length) and the _max/_min columns,
    # only the points on the field's surface (ASCII headers for pgfplots)
    import numpy as np
    from vtkmodules.numpy_interface import dataset_adapter as dsa
    for tag, src in [('plate_kappa', plate_k), ('plate_eta', plate_e),
                     ('lrmm_kappa', lrmm_k), ('lrmm_eta', lrmm_e)]:
        d = dsa.WrapDataObject(servermanager.Fetch(src))
        xyz = np.asarray(d.Points)
        cols = {'x': xyz[:, 0], 'y': xyz[:, 1],
                's': np.asarray(d.PointData['arc_length_minimum'])}
        for a in d.PointData.keys():
            if a.endswith(('_max', '_min')):
                cols[a.replace('κ', 'kappa').replace('ηin', 'eta_in')
                      .replace('η', 'eta')] = np.asarray(d.PointData[a])
        data = np.column_stack(list(cols.values()))
        data = data[~np.isnan(data[:, 3])]
        np.savetxt(os.path.join(OUT_DIR, 'max_line_%s.csv' % tag), data,
                   delimiter=',', header=','.join(cols), comments='', fmt='%.6g')
    print('saved ' + base + '.{png,pdf} and csv files in ' + OUT_DIR)
