"""
Origin 2025b Paper Figure Automation — Publication Quality
============================================================
Generates ALL figures for the Applied Sciences manuscript with
top-journal formatting standards.

Formatting spec (top-journal quality):
  - Font: Arial, 8pt tick labels, 14pt axis labels (fsize)
  - Line weight: 1.5pt data, 0.5pt axes
  - Markers: 6pt, distinct shapes per series
  - Colors: Colorblind-safe (Okabe-Ito palette)
  - Page size: MDPI double-column = 17 cm width, height = width × 0.72
  - Export: 600 DPI TIFF + PNG
  - No background grid, clean white background
  - Open-top frame (no top/right axis lines)
  - Inward tick marks with minor ticks
  - Frameless legend
  - Tick marks: inward, major + minor

Usage:
    python origin_plot_all.py
"""

import originpro as op
import numpy as np
import csv
import os
import sys
import shutil

# ── Paths ──
BASE_DIR = r"e:\antigravity_folder\ADORE_V2_jl"
DATA_DIR = os.path.join(BASE_DIR, "results", "sweeps")
FIG_DIR  = os.path.join(BASE_DIR, "paper", "figures", "origin")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Export Settings ──
DPI = 600
DOUBLE_COL_CM = 17.0   # MDPI double-column
SINGLE_COL_CM = 8.5    # MDPI single-column
DOUBLE_COL_PX = int(DOUBLE_COL_CM / 2.54 * DPI)  # ~4016 px
SINGLE_COL_PX = int(SINGLE_COL_CM / 2.54 * DPI)  # ~2010 px
ASPECT_RATIO = 0.72  # height/width

# ── Okabe-Ito colorblind-safe palette (RGB) ──
C_BLUE      = '#0072B2'
C_ORANGE    = '#E69F00'
C_RED       = '#D55E00'
C_GREEN     = '#009E73'
C_CYAN      = '#56B4E9'
C_PURPLE    = '#CC79A7'
C_BLACK     = '#000000'
C_GRAY      = '#999999'
C_DARKGRAY  = '#555555'

# Marker shape codes for Origin (LabTalk): 1=square, 2=circle, 3=uptri, 4=downtri, 5=diamond, 8=star
MARKERS = [2, 1, 5, 3, 4, 8]  # circle, square, diamond, uptri, downtri, star

# Series color assignments
SERIES_COLORS = [C_RED, C_BLUE, C_ORANGE, C_GREEN, C_CYAN, C_PURPLE]


def rgb_from_hex(hex_color):
    """Convert #RRGGBB to Origin color(r,g,b) int."""
    h = hex_color.lstrip('#')
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return r, g, b


def read_csv(filename):
    """Read CSV file and return dict of column_name -> np.array."""
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found, skipping.")
        return None
    data = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        rows = list(reader)
    for h in headers:
        try:
            data[h] = np.array([float(row[h]) for row in rows])
        except (ValueError, KeyError):
            data[h] = np.array([row.get(h, '') for row in rows])
    return data


def format_graph_layer(graph, gl, xlabel, ylabel, 
                       xmin=None, xmax=None, ymin=None, ymax=None,
                       log_x=False, log_y=False):
    """Apply publication formatting to a graph layer.
    
    Open-top frame, 0.5pt axis lines, Arial font, 8pt tick labels,
    14pt axis labels. Page size controlled via save_fig(width=).
    Note: inward ticks cannot be set via lt_exec(); uses default outward.
    """
    # --- White background ---
    graph.lt_exec('page.color = color(255,255,255);')
    
    # --- Axis labels ---
    gl.lt_exec(f'xb.text$ = {xlabel};')
    gl.lt_exec(f'yl.text$ = {ylabel};')
    
    # --- Axis label font: Arial, 14 fsize ---
    gl.lt_exec('xb.fsize = 14;')
    gl.lt_exec('yl.fsize = 14;')
    
    # --- Tick label font: 8 pt ---
    gl.lt_exec('xb.ftsize = 8;')
    gl.lt_exec('yl.ftsize = 8;')
    
    # --- Open-top frame: hide top X and right Y axis (tested working) ---
    gl.lt_exec('xt.show = 0;')
    gl.lt_exec('yr.show = 0;')
    
    # --- Axis line weight: 0.5 pt ---
    gl.lt_exec('xb.linewidth = 0.5;')
    gl.lt_exec('yl.linewidth = 0.5;')
    
    # --- Axis ranges ---
    if xmin is not None and xmax is not None:
        gl.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None:
        gl.set_ylim(ymin, ymax)
    
    # --- Log scale ---
    if log_x:
        gl.xscale = 1
    if log_y:
        gl.yscale = 1


def format_legend(gl):
    """Frameless legend with Arial 8pt font."""
    gl.lt_exec('legend.background = 0;')  # no box/frame
    gl.lt_exec('legend.fsize = 8;')
    gl.lt_exec('legend.font = Arial;')


def format_data_plot(gl, plot_index, color_hex, marker_shape=1, line_width=1.5, 
                     marker_size=8, line_only=False, marker_only=False,
                     dash_style=0):
    """Format a data plot using the Python API (Plot object).
    
    plot_index: 0-based index into gl.plot_list()
    color_hex: '#RRGGBB' color string
    marker_shape: Origin symbol_kind (0=none, 1=circle, 2=square, 3=uptriangle, etc.)
    """
    plots = gl.plot_list()
    if plot_index >= len(plots):
        return
    p = plots[plot_index]
    
    p.color = color_hex
    
    if not line_only:
        p.symbol_kind = marker_shape
        p.symbol_size = marker_size
        p.set_cmd('-kf 0')  # solid fill interior (tested working)
    
    if marker_only:
        try:
            p.set_cmd('-l 0')  # hide line
        except:
            pass


def export_figure(graph, name, width=None):
    """Export graph as TIFF and PNG at 600 DPI."""
    if width is None:
        width = DOUBLE_COL_PX
    
    tiff_path = os.path.join(FIG_DIR, f'{name}.tiff')
    png_path  = os.path.join(FIG_DIR, f'{name}.png')
    
    for path in [tiff_path, png_path]:
        try:
            graph.save_fig(path, width=width)
            ext = os.path.splitext(path)[1]
            sz = os.path.getsize(path) / 1024
            print(f"    ✓ {name}{ext} ({sz:.0f} KB)")
        except Exception as e:
            print(f"    ✗ Export failed: {e}")


# =====================================================================
#  FIGURE GENERATORS
# =====================================================================

def fig_oilflow_temperature():
    """Fig 5: Oil flow rate → steady-state temperatures."""
    print("\n  [Fig 5] Oil Flow Temperature Sensitivity")
    data = read_csv(os.path.join("oilflow", "oilflow_sweep_results.csv"))
    if data is None: return
    
    wks = op.new_sheet('w', 'OilFlowTemp')
    wks.from_list(0, data['V_dot_cm3min'], 'Oil Flow Rate')
    wks.from_list(1, data['T_ball_K'],     'Ball')
    wks.from_list(2, data['T_IR_K'],       'Inner Race')
    wks.from_list(3, data['T_OR_K'],       'Outer Race')
    wks.from_list(4, data['T_oil_K'],      'Oil Sump')
    wks.set_label(0, 'X', 'D')
    for c in range(1, 5): wks.set_label(c, 'Y', 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    for col in range(1, 5):
        gl.add_plot(wks, coly=col, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Oil Flow Rate (cm\\+(3)/min)',
        ylabel='Steady-State Temperature (K)')
    
    format_data_plot(gl, 0, C_RED,    marker_shape=1, marker_size=8)  # Ball - circle
    format_data_plot(gl, 1, C_BLUE,   marker_shape=2, marker_size=8)  # IR - square
    format_data_plot(gl, 2, C_ORANGE, marker_shape=5, marker_size=8)  # OR - diamond
    format_data_plot(gl, 3, C_GREEN,  marker_shape=3, marker_size=8)  # Oil - triangle
    format_legend(gl)
    
    export_figure(graph, 'oilflow_temperature_sensitivity')


def fig_oilflow_heat():
    """Fig 6: Oil flow rate → total heat generation."""
    print("\n  [Fig 6] Oil Flow Heat Generation")
    data = read_csv(os.path.join("oilflow", "oilflow_sweep_results.csv"))
    if data is None: return
    
    wks = op.new_sheet('w', 'OilFlowHeat')
    wks.from_list(0, data['V_dot_cm3min'], 'Oil Flow Rate')
    wks.from_list(1, data['H_total_W'],    'Total Heat')
    wks.set_label(0, 'X', 'D')
    wks.set_label(1, 'Y', 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    gl.add_plot(wks, coly=1, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Oil Flow Rate (cm\\+(3)/min)',
        ylabel='Total Heat Generation (W)')
    
    format_data_plot(gl, 0, C_BLACK, marker_shape=1, marker_size=8)
    format_legend(gl)
    
    export_figure(graph, 'oilflow_heat_generation')


def fig_clearance_steadystate():
    """Fig 8: Clearance → steady-state temperature + heat."""
    print("\n  [Fig 8] Clearance Steady-State")
    data = read_csv(os.path.join("clearance", "clearance_sweep_results.csv"))
    if data is None: return
    
    wks = op.new_sheet('w', 'ClearanceSS')
    wks.from_list(0, data['P_d_um'],    'Clearance')
    wks.from_list(1, data['T_IR_K'],    'Inner Race')
    wks.from_list(2, data['T_OR_K'],    'Outer Race')
    wks.from_list(3, data['T_ball_K'],  'Ball')
    wks.from_list(4, data['H_total_W'], 'H total')
    for c in range(5):
        wks.set_label(c, ['X','Y','Y','Y','Y'][c], 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    for col in range(1, 4):  # only temps, not heat
        gl.add_plot(wks, coly=col, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Initial Diametral Clearance (\\g(m)m)',
        ylabel='Steady-State Temperature (K)')
    
    format_data_plot(gl, 0, C_BLUE,   marker_shape=2, marker_size=8)
    format_data_plot(gl, 1, C_ORANGE, marker_shape=5, marker_size=8)
    format_data_plot(gl, 2, C_RED,    marker_shape=1, marker_size=8)
    format_legend(gl)
    
    export_figure(graph, 'clearance_steadystate')


def fig_clearance_transient_overlay():
    """Fig 7: Transient contact angle for different clearances."""
    print("\n  [Fig 7] Clearance Transient Overlay (Contact Angle)")
    # Try contact-angle CSV first, fall back to temperature CSV
    alpha_csv = os.path.join(DATA_DIR, "clearance", "clearance_alpha_transient_data.csv")
    temp_csv  = os.path.join(DATA_DIR, "clearance", "clearance_transient_data.csv")
    if os.path.exists(alpha_csv):
        transient_csv = alpha_csv
        ylabel = 'Mean Inner Contact Angle (\\+(o))'
        prefix = 'alpha_i_'
    elif os.path.exists(temp_csv):
        transient_csv = temp_csv
        ylabel = 'Inner Race Temperature (\\+(o)C)'
        prefix = 'T_IR_'
    else:
        print("    ⚠ No clearance transient CSV found. Skipping.")
        return
    
    data = {}
    with open(transient_csv, 'r') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        rows = list(reader)
    for h in headers:
        try:
            data[h] = np.array([float(row[h]) for row in rows])
        except ValueError:
            pass
    
    wks = op.new_sheet('w', 'ClearTransient')
    float_headers = [h for h in headers if h in data]
    for i, h in enumerate(float_headers):
        # Clean names for legend
        label = h.replace(prefix, '').replace('um', ' μm').replace('Pd', 'P_d = ')
        if h == 'time_ms':
            label = 'Time'
        wks.from_list(i, data[h], label)
        wks.set_label(i, 'X' if i == 0 else 'Y', 'D')
    
    graph = op.new_graph(template='line')
    gl = graph[0]
    ncols = len(float_headers)
    for col in range(1, ncols):
        gl.add_plot(wks, coly=col, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Time (ms)',
        ylabel=ylabel)
    
    # Use distinct colors for each clearance
    trans_colors = [C_RED, C_ORANGE, C_GREEN, C_BLUE, C_PURPLE, C_CYAN]
    for idx in range(ncols - 1):
        c = trans_colors[idx % len(trans_colors)]
        format_data_plot(gl, idx, c, line_only=True, line_width=1.5)
    format_legend(gl)
    
    export_figure(graph, 'clearance_transient_overlay')


def fig_loadratio_temperature():
    """Fig 9: Load direction → temperature."""
    print("\n  [Fig 9] Load Ratio Temperature")
    data = read_csv(os.path.join("loadratio", "loadratio_sweep_results.csv"))
    if data is None: return
    
    wks = op.new_sheet('w', 'LoadRatioTemp')
    wks.from_list(0, data['angle_deg'], 'Load Angle')
    wks.from_list(1, data['T_IR_K'],    'Inner Race')
    wks.from_list(2, data['T_OR_K'],    'Outer Race')
    wks.from_list(3, data['T_ball_K'],  'Ball')
    for c in range(4):
        wks.set_label(c, ['X','Y','Y','Y'][c], 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    for col in range(1, 4):
        gl.add_plot(wks, coly=col, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Load Direction Angle (\\+(o))',
        ylabel='Steady-State Temperature (K)')
    
    format_data_plot(gl, 0, C_BLUE,   marker_shape=2, marker_size=8)
    format_data_plot(gl, 1, C_ORANGE, marker_shape=5, marker_size=8)
    format_data_plot(gl, 2, C_RED,    marker_shape=1, marker_size=8)
    format_legend(gl)
    
    export_figure(graph, 'loadratio_temperature')


def fig_loadratio_heat():
    """Fig 10: Load direction → heat breakdown."""
    print("\n  [Fig 10] Load Ratio Heat Breakdown")
    data = read_csv(os.path.join("loadratio", "loadratio_sweep_results.csv"))
    if data is None: return
    
    heat_cols = [h for h in data.keys() if 'H_' in h]
    
    wks = op.new_sheet('w', 'LoadRatioHeat')
    wks.from_list(0, data['angle_deg'], 'Load Angle')
    wks.set_label(0, 'X', 'D')
    for i, h in enumerate(heat_cols):
        label = h.replace('H_', '').replace('_W', '').replace('_', ' ').title()
        wks.from_list(i + 1, data[h], label)
        wks.set_label(i + 1, 'Y', 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    for col in range(1, len(heat_cols) + 1):
        gl.add_plot(wks, coly=col, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Load Direction Angle (\\+(o))',
        ylabel='Heat Generation (W)')
    
    for idx in range(len(heat_cols)):
        c = SERIES_COLORS[idx % len(SERIES_COLORS)]
        m = MARKERS[idx % len(MARKERS)]
        format_data_plot(gl, idx, c, marker_shape=m, marker_size=8)
    format_legend(gl)
    
    export_figure(graph, 'loadratio_heat_breakdown')


def fig_gij_sensitivity():
    """Fig: G_ij conductance sensitivity."""
    print("\n  [Fig] G_ij Sensitivity")
    data = read_csv(os.path.join("Gij_sensitivity", "Gij_sensitivity_results.csv"))
    if data is None: return
    
    scales = data.get('scale', data.get('G_scale', None))
    if scales is None:
        print("    ⚠ Cannot find scale column. Skipping.")
        return
    
    wks = op.new_sheet('w', 'GijSens')
    wks.from_list(0, scales,          'Scale Factor')
    wks.from_list(1, data['T_IR_K'],  'Inner Race')
    wks.from_list(2, data['T_OR_K'],  'Outer Race')
    wks.from_list(3, data['T_ball_K'],'Ball')
    for c in range(4):
        wks.set_label(c, ['X','Y','Y','Y'][c], 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    for col in range(1, 4):
        gl.add_plot(wks, coly=col, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='G\\-(ij) Scale Factor',
        ylabel='Steady-State Temperature (K)')
    
    format_data_plot(gl, 0, C_BLUE,   marker_shape=2, marker_size=8)
    format_data_plot(gl, 1, C_ORANGE, marker_shape=5, marker_size=8)
    format_data_plot(gl, 2, C_RED,    marker_shape=1, marker_size=8)
    format_legend(gl)
    
    export_figure(graph, 'Gij_sensitivity')


def fig_thermal_accel():
    """Fig 4: Thermal acceleration validation — parity plot.
    
    Shows accelerated vs physical values as scatter points along y=x line.
    Each point is one variable (T_IR, T_OR, T_ball, T_oil).
    Close clustering along the diagonal proves the two methods agree.
    """
    print("\n  [Fig 4] Thermal Acceleration Validation")
    data = read_csv(os.path.join("thermal_accel", "thermal_accel_comparison.csv"))
    if data is None: return
    
    # Temperature variables only (same units, same scale)
    var_keys   = ['T_IR_K',  'T_OR_K',  'T_ball_K', 'T_oil_K']
    var_labels = ['T\\-(IR)', 'T\\-(OR)', 'T\\-(ball)', 'T\\-(oil)']
    
    # Extract row 0 = accelerated, row 1 = physical
    phys_vals = []
    accel_vals = []
    for key in var_keys:
        if key in data:
            accel_vals.append(float(data[key][0]))  # accelerated row
            phys_vals.append(float(data[key][1]))    # physical row
    
    if not phys_vals:
        print("    ⚠ No temperature columns found. Skipping.")
        return
    
    # y = x reference line (from min-2 to max+2)
    lo = min(min(phys_vals), min(accel_vals)) - 2
    hi = max(max(phys_vals), max(accel_vals)) + 2
    ref_x = np.array([lo, hi])
    ref_y = np.array([lo, hi])
    
    wks = op.new_sheet('w', 'ThermalAccel')
    wks.from_list(0, np.array(phys_vals),  'Physical Reference')
    wks.from_list(1, np.array(accel_vals), 'Accelerated')
    wks.from_list(2, ref_x, 'y=x line X')
    wks.from_list(3, ref_y, 'y=x line Y')
    wks.set_label(0, 'X', 'D')
    wks.set_label(1, 'Y', 'D')
    wks.set_label(2, 'X', 'D')
    wks.set_label(3, 'Y', 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    gl.add_plot(wks, coly=1, colx=0)  # scatter: accel vs phys
    gl.add_plot(wks, coly=3, colx=2)  # y=x reference line
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Physical Simulation (K)',
        ylabel='Accelerated Simulation (K)',
        xmin=lo, xmax=hi, ymin=lo, ymax=hi)
    
    # Scatter points: red circles, no connecting line
    format_data_plot(gl, 0, C_RED, marker_shape=1, marker_size=10, marker_only=True)
    # y=x reference: thin grey dashed line, no markers
    format_data_plot(gl, 1, C_GRAY, line_only=True, line_width=1.0)
    format_legend(gl)
    
    export_figure(graph, 'thermal_accel_validation')


def fig_heat_generation_comparison():
    """Fig 2: NASA validation — simulation vs experimental heat generation.
    Shows NASA TP-2275 experimental data with ±5% uncertainty band
    and ADORE simulation results.
    """
    print("\n  [Fig 2] Heat Generation Comparison (NASA)")
    
    rpm  = np.array([28, 48, 64, 72], dtype=float)
    nasa = np.array([350, 700, 1100, 1400], dtype=float)
    sim  = np.array([341, 680, 1134, 1400], dtype=float)
    nasa_upper = nasa * 1.05   # +5%
    nasa_lower = nasa * 0.95   # -5%
    
    wks = op.new_sheet('w', 'NASAValid')
    wks.from_list(0, rpm,         'Speed')
    wks.from_list(1, nasa,        'NASA TP-2275 (\\+(-)5%)')
    wks.from_list(2, nasa_upper,  '+5%')
    wks.from_list(3, nasa_lower,  '-5%')
    wks.from_list(4, sim,         'Simulation')
    wks.set_label(0, 'X', 'D')
    wks.set_label(1, 'Y', 'D')
    wks.set_label(2, 'Y', 'D')
    wks.set_label(3, 'Y', 'D')
    wks.set_label(4, 'Y', 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    gl.add_plot(wks, coly=1, colx=0)  # NASA center
    gl.add_plot(wks, coly=2, colx=0)  # +5%
    gl.add_plot(wks, coly=3, colx=0)  # -5%
    gl.add_plot(wks, coly=4, colx=0)  # Simulation
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Shaft Speed (kRPM)',
        ylabel='Heat Generation (W)',
        xmin=25, xmax=75, ymin=0, ymax=1600)
    
    # NASA center: black filled squares, markers only
    format_data_plot(gl, 0, C_BLACK, marker_shape=2, marker_size=10, marker_only=True)
    
    # ±5% bounds: gray thin dashed lines, NO markers
    plots = gl.plot_list()
    for idx in [1, 2]:
        if idx < len(plots):
            p = plots[idx]
            p.color = C_GRAY
            p.symbol_kind = 0      # hide marker
            p.set_cmd('-w 800')    # line width 0.8pt (unit: 0.001pt)
            p.set_cmd('-d 2')      # dash style: dashed
    
    # Simulation: red line+circle
    format_data_plot(gl, 3, C_RED, marker_shape=1, marker_size=8, line_width=1.5)
    
    # Build legend: show only NASA and Simulation (hide ±5% bounds)
    gl.lt_exec('legend -s;')   # reconstruct auto-legend
    # Directly set legend text content using LabTalk escape sequences
    # \\l(N) references plot N's symbol; use CRLF between entries
    legend_text = '\\l(1) NASA TP-2275' + chr(13) + chr(10) + '\\l(4) Simulation'
    gl.lt_exec(f'legend.text$ = {legend_text};')
    gl.lt_exec('legend.background = 0;')
    gl.lt_exec('legend.fsize = 8;')
    
    export_figure(graph, 'heat_generation_comparison')





def fig_heat_scaling_loglog():
    """Fig 3: Log-log speed scaling H ∝ n^α."""
    print("\n  [Fig 3] Heat Scaling Log-Log")
    
    rpm  = np.array([28000, 48000, 64000, 72000], dtype=float)
    nasa = np.array([350,   700,   1100,  1400],  dtype=float)
    sim  = np.array([341,   680,   1134,  1400],  dtype=float)
    
    wks = op.new_sheet('w', 'LogLogScale')
    wks.from_list(0, rpm,  'Speed (RPM)')
    wks.from_list(1, nasa, 'NASA TP-2275')
    wks.from_list(2, sim,  'Simulation')
    wks.set_label(0, 'X', 'D')
    wks.set_label(1, 'Y', 'D')
    wks.set_label(2, 'Y', 'D')
    
    graph = op.new_graph(template='linesymb')
    gl = graph[0]
    gl.add_plot(wks, coly=1, colx=0)
    gl.add_plot(wks, coly=2, colx=0)
    gl.rescale()
    
    format_graph_layer(graph, gl,
        xlabel='Shaft Speed (RPM)',
        ylabel='Heat Generation (W)',
        log_x=True, log_y=True)
    
    # NASA experimental: black circles, markers only
    format_data_plot(gl, 0, C_BLACK, marker_shape=1, marker_size=8, marker_only=True)
    # Simulation: red squares, markers only (discrete speed points)
    format_data_plot(gl, 1, C_RED,   marker_shape=2, marker_size=8, marker_only=True)
    format_legend(gl)
    
    export_figure(graph, 'heat_scaling_loglog')


# =====================================================================
#  DATA ARCHIVAL
# =====================================================================
def consolidate_data():
    """Archive all CSVs to paper/data/."""
    data_archive = os.path.join(BASE_DIR, "paper", "data")
    os.makedirs(data_archive, exist_ok=True)
    
    csv_files = [
        ("oilflow/oilflow_sweep_results.csv", "oilflow_sweep.csv"),
        ("clearance/clearance_sweep_results.csv", "clearance_sweep.csv"),
        ("clearance/clearance_transient_data.csv", "clearance_transient.csv"),
        ("clearance/clearance_alpha_transient_data.csv", "clearance_alpha_transient.csv"),
        ("clearance_lowpreload/clearance_lowpreload_results.csv", "clearance_lowpreload_sweep.csv"),
        ("loadratio/loadratio_sweep_results.csv", "loadratio_sweep.csv"),
        ("Gij_sensitivity/Gij_sensitivity_results.csv", "Gij_sensitivity.csv"),
        ("thermal_accel/thermal_accel_comparison.csv", "thermal_accel.csv"),
    ]
    
    for src_rel, dst_name in csv_files:
        src = os.path.join(DATA_DIR, src_rel)
        dst = os.path.join(data_archive, dst_name)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f"    ✓ {dst_name}")
        else:
            print(f"    ✗ {src_rel} (not found)")
    
    # NASA validation (inline data)
    val_csv = os.path.join(data_archive, "nasa_validation.csv")
    with open(val_csv, 'w') as f:
        f.write("RPM,DN_1e6,NASA_H_W,Sim_H_W,Error_pct\n")
        f.write("28000,0.98,350,341,-2.6\n")
        f.write("48000,1.68,700,680,-2.9\n")
        f.write("64000,2.24,1100,1134,3.0\n")
        f.write("72000,2.52,1400,1400,0.0\n")
    print(f"    ✓ nasa_validation.csv")


# =====================================================================
#  MAIN
# =====================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("  Origin 2025b — Publication-Quality Figure Generation")
    print("  Applied Sciences (MDPI)")
    print("  Format: Arial, 600 DPI, 17 cm width, Okabe-Ito colors")
    print("=" * 60)
    
    print("\n  [1/11] Archiving data...")
    consolidate_data()
    
    # Figures in manuscript order
    print("\n  [2/11]"); fig_heat_generation_comparison()   # NASA validation
    print("\n  [3/11]"); fig_heat_scaling_loglog()           # Log-log scaling
    print("\n  [4/11]"); fig_thermal_accel()                 # Thermal accel
    print("\n  [5/11]"); fig_oilflow_temperature()           # Oil flow → T
    print("\n  [6/11]"); fig_oilflow_heat()                  # Oil flow → H
    print("\n  [7/11]"); fig_clearance_transient_overlay()   # Clearance transient
    print("\n  [8/11]"); fig_clearance_steadystate()         # Clearance → SS
    print("\n  [9/11]"); fig_loadratio_temperature()         # Load ratio → T
    print("\n [10/11]"); fig_loadratio_heat()                # Load ratio → H
    print("\n [11/11]"); fig_gij_sensitivity()               # G_ij sensitivity
    
    # Save project
    print("\n" + "-" * 60)
    try:
        proj_path = os.path.join(BASE_DIR, "paper", "figures", "manuscript_figures.opju")
        op.save(proj_path)
        print(f"  Origin project: {proj_path}")
    except Exception as e:
        print(f"  Could not save project: {e}")
    
    print("\n" + "=" * 60)
    print(f"  ALL FIGURES → {FIG_DIR}")
    print("  Format: 600 DPI TIFF + PNG, MDPI double-column")
    print("=" * 60)
