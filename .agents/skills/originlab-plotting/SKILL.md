---
name: OriginLab Publication Plotting
description: Automate publication-quality figure generation using the `originpro` Python API with Origin 2025b. Covers worksheet setup, graph creation, formatting, color palettes, export, and known pitfalls.
---

# OriginLab Publication Plotting Skill

Generate publication-quality figures in Origin 2025b via the `originpro` Python API.
This skill encodes all proven patterns and hard-won pitfalls from the ADORE project.

## Prerequisites

| Item | Detail |
|------|--------|
| **Origin** | Origin 2025b, installed at `d:\Program Files\OriginLab\Origin2025b\` |
| **Python** | `C:\Users\yuxinlu_qh\miniconda3\python.exe` |
| **Package** | `originpro` (ships with Origin, auto-discovered) |
| **Kill first** | Always `Stop-Process -Name "Origin*" -Force` before running |

## Publication Formatting Spec

```
Font:       Arial, 8pt ticks, 14pt axis labels (fsize units)
Line width: 1.5pt data lines, 0.5pt axes
Markers:    10pt default, filled (symbol_interior=1)
Colors:     Okabe-Ito colorblind-safe palette (see below)
Page size:  17 cm (MDPI double-column) width
Export:     600 DPI, TIFF + PNG dual export
Background: White, no grid
```

## Okabe-Ito Colorblind-Safe Palette

```python
C_BLUE      = '#0072B2'
C_ORANGE    = '#E69F00'
C_RED       = '#D55E00'
C_GREEN     = '#009E73'
C_CYAN      = '#56B4E9'
C_PURPLE    = '#CC79A7'
C_BLACK     = '#000000'
C_GRAY      = '#999999'
```

## Core Workflow

### Step 1: Read CSV Data

```python
import csv, numpy as np, os

def read_csv(filename):
    """Read CSV → dict of column_name → np.array."""
    data = {}
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    for h in reader.fieldnames:
        try:
            data[h] = np.array([float(row[h]) for row in rows])
        except (ValueError, KeyError):
            data[h] = np.array([row.get(h, '') for row in rows])
    return data
```

### Step 2: Create Worksheet

```python
import originpro as op

wks = op.new_sheet('w', 'SheetName')

# from_list(col_index, np.array, long_name)
wks.from_list(0, x_data,  'X Label')
wks.from_list(1, y1_data, 'Series 1')
wks.from_list(2, y2_data, 'Series 2')

# Mark column designations: 'X' or 'Y'
wks.set_label(0, 'X', 'D')
wks.set_label(1, 'Y', 'D')
wks.set_label(2, 'Y', 'D')
```

> [!IMPORTANT]
> The **long name** (3rd arg of `from_list`) appears in the legend.
> Always set column designation with `set_label(col, 'X'|'Y', 'D')`.

### Step 3: Create Graph & Add Plots

```python
# Templates: 'linesymb' (markers+lines), 'line' (lines only)
graph = op.new_graph(template='linesymb')
gl = graph[0]  # first graph layer

for col in range(1, num_y_cols + 1):
    gl.add_plot(wks, coly=col, colx=0)

gl.rescale()  # auto-fit axes to data
```

> [!WARNING]
> **Never use the `'column'` template** for comparison plots — it creates stacked bars
> that are nearly impossible to format correctly. Use `'linesymb'` with `marker_only=True` instead.

### Step 4: Format Graph Layer

```python
def format_graph_layer(graph, gl, xlabel, ylabel,
                       xmin=None, xmax=None, ymin=None, ymax=None,
                       log_x=False, log_y=False):
    # White background
    graph.lt_exec('page.color = color(255,255,255);')

    # Axis labels (LabTalk — proven stable)
    gl.lt_exec(f'xb.text$ = {xlabel};')
    gl.lt_exec(f'yl.text$ = {ylabel};')

    # Axis label font size
    gl.lt_exec('xb.fsize = 14;')
    gl.lt_exec('yl.fsize = 14;')

    # Open-top frame: hide top-X and right-Y axes (proven working)
    gl.lt_exec('xt.show = 0;')
    gl.lt_exec('yr.show = 0;')

    # Frameless legend (proven working)
    gl.lt_exec('legend.background = 0;')

    # Axis ranges (Python API)
    if xmin is not None and xmax is not None:
        gl.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None:
        gl.set_ylim(ymin, ymax)

    # Log scale
    if log_x: gl.xscale = 1
    if log_y: gl.yscale = 1
```

> [!WARNING]
> **Inward tick marks** cannot be set via `lt_exec()`. The documented LabTalk
> property `layer.x.ticks = 5` does not take effect through `gl.lt_exec()`,
> `graph.lt_exec()`, or `op.lt_exec()`. Default outward ticks are acceptable
> for publication.

### Step 5: Format Individual Data Plots

```python
def format_data_plot(gl, plot_index, color_hex,
                     marker_shape=1, line_width=1.5,
                     marker_size=10, line_only=False, marker_only=False):
    plots = gl.plot_list()
    if plot_index >= len(plots):
        return
    p = plots[plot_index]

    p.color = color_hex  # '#RRGGBB' string

    if not line_only:
        p.symbol_kind = marker_shape
        p.symbol_size = marker_size
        # IMPORTANT: p.symbol_interior = 1 does NOT work reliably.
        # Use set_cmd instead for filled markers:
        p.set_cmd('-kf 0')       # 0 = solid fill interior
        p.set_cmd(f'-csf {color_hex}')  # fill color matches line color

    if marker_only:
        try:
            p.set_cmd('-l 0')  # hide connecting line
        except:
            pass
```

**Marker shape codes** (`symbol_kind` values):

| Code | Shape |
|------|-------|
| 0 | None |
| 1 | Circle |
| 2 | Square |
| 3 | Up triangle |
| 4 | Down triangle |
| 5 | Diamond |
| 8 | Star |

### Step 6: Export

```python
DPI = 600
WIDTH_PX = int(17.0 / 2.54 * DPI)  # ≈ 4016 px for 17 cm

def export_figure(graph, name, width=WIDTH_PX):
    fig_dir = r"path\to\figures"
    for ext in ['.tiff', '.png']:
        graph.save_fig(os.path.join(fig_dir, f'{name}{ext}'), width=width)
```

## LabTalk Special Characters for Axis Labels

| Syntax | Renders as | Example |
|--------|-----------|---------|
| `\\+(text)` | Superscript | `cm\\+(3)` → cm³ |
| `\\-(text)` | Subscript | `G\\-(ij)` → G_ij |
| `\\g(m)` | Greek μ | `\\g(m)m` → μm |
| `\\+(o)` | Degree symbol | `\\+(o)C` → °C |

## Figure Patterns

### Multi-Series Line + Marker Plot (most common)

```python
graph = op.new_graph(template='linesymb')
gl = graph[0]
for col in range(1, 4):
    gl.add_plot(wks, coly=col, colx=0)
gl.rescale()

format_graph_layer(graph, gl, xlabel='...', ylabel='...')
format_data_plot(gl, 0, C_BLUE,   marker_shape=2, marker_size=10)  # square
format_data_plot(gl, 1, C_ORANGE, marker_shape=5, marker_size=10)  # diamond
format_data_plot(gl, 2, C_RED,    marker_shape=1, marker_size=10)  # circle
```

### Parity / Validation Scatter Plot

Good for showing agreement between two methods (e.g., accelerated vs physical).

```python
# Data: 4 scatter points + y=x reference line
wks.from_list(0, phys_vals,  'Physical')
wks.from_list(1, accel_vals, 'Accelerated')
wks.from_list(2, ref_x,     'y=x X')
wks.from_list(3, ref_y,     'y=x Y')

graph = op.new_graph(template='linesymb')
gl = graph[0]
gl.add_plot(wks, coly=1, colx=0)  # scatter
gl.add_plot(wks, coly=3, colx=2)  # y=x line

format_data_plot(gl, 0, C_RED, marker_shape=1, marker_size=12, marker_only=True)
format_data_plot(gl, 1, C_GRAY, line_only=True, line_width=1.0)
```

### Multi-Line Transient Overlay (no markers)

```python
graph = op.new_graph(template='line')  # 'line' not 'linesymb'
gl = graph[0]
for col in range(1, ncols):
    gl.add_plot(wks, coly=col, colx=0)
gl.rescale()

colors = [C_RED, C_ORANGE, C_GREEN, C_BLUE, C_PURPLE, C_CYAN]
for idx in range(ncols - 1):
    format_data_plot(gl, idx, colors[idx % len(colors)],
                     line_only=True, line_width=1.5)
```

### Log-Log Scaling Plot

```python
format_graph_layer(graph, gl, xlabel='...', ylabel='...',
                   log_x=True, log_y=True)
```

## Known Pitfalls & Fixes

> [!CAUTION]
> **Critical: Kill Origin before running.** If a stale Origin session is open,
> the script attaches to it and inherited state causes silent formatting failures.
> Always run: `Stop-Process -Name "Origin*" -Force; Start-Sleep -Seconds 2`

1. **Function defined but not called** — If a figure function exists but isn't
   listed in `main()`, the figure on disk is from an older run with no formatting.
   Always verify every `fig_*()` function is called.

2. **`column` template** — Creates stacked/grouped bars that are very hard to
   format via the Python API. Prefer `linesymb` with `marker_only=True` for
   comparison plots.

3. **Axis labels reset by `rescale()`** — Origin can auto-set axis labels from
   worksheet column names after `rescale()`. Always call `format_graph_layer()`
   **after** `gl.rescale()`.

4. **Colors not applied** — If all series appear grey/black, check:
   - Is the function actually being called in `main()`?
   - Is `gl.plot_list()` returning the expected number of plots?
   - Is the template correct (`linesymb` not `line` if markers are needed)?

5. **LabTalk crashes** — Avoid `set -d` range syntax, `wks.col` property chains,
   and `layer.plot` manipulation via LabTalk. Stick to the Python API for
   per-plot formatting.

6. **Worksheet column long names** — These appear in the legend. Set meaningful
   names via `wks.from_list(col, data, 'Pretty Label')`.

## Reference Script

The canonical implementation is:
```
e:\antigravity_folder\ADORE_V2_jl\paper\origin_plot_all.py
```

This script generates 10 publication figures for the Applied Sciences manuscript,
covering: NASA validation, log-log scaling, thermal acceleration parity,
oil-flow sensitivity, clearance sweeps, load-ratio analysis, and G_ij sensitivity.
