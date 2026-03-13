import csv
rows = list(csv.DictReader(open(r'e:\antigravity_folder\ADORE_V2_jl\results\sweeps\loadratio\loadratio_sweep_results.csv')))
hk = [k for k in rows[0].keys() if 'H_' in k]
print('Heat columns:', hk)
print()
for r in rows:
    vals = ', '.join(f'{k}={float(r[k]):>8.3f}' for k in hk)
    print(f'angle={r["angle_deg"]:>5s}: {vals}')
