"""
Origin 2025b Python Control - Paper Figure Generation
=====================================================
Uses the `originpro` package to control Origin externally.

IMPORTANT: Origin must be running before executing this script,
OR this script will launch it automatically via COM.

Usage:
    python origin_plot.py
"""

import originpro as op
import numpy as np
import os
import sys

# === Configuration ===
ORIGIN_PATH = r"d:\Program Files\OriginLab\Origin2025b"
FIGURE_DIR = r"e:\antigravity_folder\ADORE_V2_jl\paper\figures"

def setup_origin():
    """Initialize Origin connection."""
    # originpro connects via COM automation
    # Origin must be installed and registered
    if op.oext is None:
        print("ERROR: Cannot connect to Origin. Make sure Origin 2025b is installed.")
        sys.exit(1)
    print(f"Connected to Origin: {op.get_lt_str('system.version$')}")
    return True


def create_demo_plot():
    """Create a simple demo plot to verify the connection works."""
    # Create a new workbook
    wks = op.new_sheet('w', 'DemoData')
    
    # Generate demo data
    x = np.linspace(0, 10, 50)
    y1 = np.sin(x)
    y2 = np.cos(x)
    
    # Set data into worksheet
    wks.from_list(0, x, 'X')
    wks.from_list(1, y1, 'sin(x)')
    wks.from_list(2, y2, 'cos(x)')
    
    # Set column designations
    wks.set_label(0, 'X')
    wks.set_label(1, 'Y')
    wks.set_label(2, 'Y')
    
    # Create a graph
    graph = op.new_graph(template='line')
    gl = graph[0]  # first layer
    
    # Plot data
    plot1 = gl.add_plot(wks, coly=1, colx=0)
    plot2 = gl.add_plot(wks, coly=2, colx=0)
    
    # Rescale axes
    gl.rescale()
    
    # Export
    out_path = os.path.join(FIGURE_DIR, 'origin_demo.png')
    graph.save_fig(out_path, width=800)
    print(f"Demo figure saved to: {out_path}")
    
    return True


def create_validation_plot(csv_path, output_name):
    """
    Generic function to create a publication-quality plot from CSV data.
    
    Parameters:
        csv_path: path to CSV file with columns
        output_name: output filename (without extension)
    """
    wks = op.new_sheet('w', output_name)
    wks.from_file(csv_path)
    
    graph = op.new_graph(template='line')
    gl = graph[0]
    
    # Add all Y columns
    ncols = wks.cols
    for i in range(1, ncols):
        gl.add_plot(wks, coly=i, colx=0)
    
    gl.rescale()
    
    # Publication formatting
    gl.set_xlabel('X Label')
    gl.set_ylabel('Y Label')
    
    # Export at 600 DPI for journal
    out_path = os.path.join(FIGURE_DIR, f'{output_name}.tiff')
    graph.save_fig(out_path, width=2400)  # 2400px at 600dpi = 4 inches = ~10cm
    print(f"Figure saved: {out_path}")


# === Main ===
if __name__ == '__main__':
    print("=" * 50)
    print("Origin-Python Integration Test")
    print("=" * 50)
    
    try:
        setup_origin()
        create_demo_plot()
        print("\n✓ Origin-Python connection verified!")
        print("You can now use this script as a template for paper figures.")
    except Exception as e:
        print(f"\nERROR: {e}")
        print("\nTroubleshooting:")
        print("1. Make sure Origin 2025b is running")
        print("2. Or try: origin_path\\Origin64.exe /automation")
        print(f"   {ORIGIN_PATH}\\Origin64.exe /automation")
    finally:
        # Don't exit Origin when done (keep it open for inspection)
        pass
