"""Test: layer -plot addplot LabTalk X-Function for error bars."""
import originpro as op
import numpy as np
import os

rpm = np.array([28, 48, 64, 72], dtype=float)
nasa = np.array([350, 700, 1100, 1400], dtype=float)
nasa_err = nasa * 0.05
sim = np.array([341, 680, 1134, 1400], dtype=float)

wks = op.new_sheet('w', 'NASAFinal')
wks.from_list(0, rpm, 'Speed')
wks.from_list(1, nasa, 'NASA TP-2275')
wks.from_list(2, nasa_err, 'Err')
wks.from_list(3, sim, 'Simulation')
wks.cols_axis('xyey')

# Get names
op.lt_exec('string wbn$ = page.name$;')
wb = op.get_lt_str('wbn$')
op.lt_exec('string shn$ = layer.name$;')
sh = op.get_lt_str('shn$')
print(f"Workbook: {wb}, Sheet: {sh}")

# Create graph
graph = op.new_graph(template='linesymb')
gl = graph[0]
op.lt_exec('string grn$ = page.name$;')
gr = op.get_lt_str('grn$')
print(f"Graph: {gr}")

# Add NASA via Python API first (no error bars)
gl.add_plot(wks, coly=1, colx=0)

# Now try LabTalk range to associate error data to existing plot
# The key is: range syntax for error data assignment
# In Origin, after adding a plot, we can assign error data using:
# set plotN -ye rangeRef
# where plotN can be accessed as %C when selected

# Activate graph
op.lt_exec(f'win -a {gr};')
# Select layer 1 
op.lt_exec('layer -s 1;')    # select first data plot (1-indexed in layer -s)

# Now %C should refer to the first data plot
# Try assigning yErr data using full range reference
op.lt_exec(f'set %C -ye [{wb}]{sh}!C;')
print("Attempted set -ye with range")

# Alternative: try set command with column number
op.lt_exec(f'set %C -ye [{wb}]{sh}!col(3);')
print("Attempted set -ye with col(3)")

# Alternative: try the -e switch directly
op.lt_exec(f'set %C -e 1;')  # enable error bars
print("Attempted set -e 1")

# Add simulation
gl.add_plot(wks, coly=3, colx=0)
gl.rescale()

plots = gl.plot_list()
print(f"Plots: {len(plots)}")

graph.save_fig(os.path.join(r'e:\antigravity_folder\ADORE_V2_jl\paper\figures\origin', 'test_errbar8.png'), width=2000)
print("Exported test_errbar8.png")

op.exit()
