#!/usr/bin/env python3
"""Siemens MAMMOMAT Revelation against Hologic Selenia Dimensions on the L2 phantom: acquisition data from the DICOM
headers and the numbers of the tuned recombination, one line per exposure and a summary per system and LE kV.

Reads:  the Siemens exposures (28kV, 30kV, 32kV folders) and the Excel file of l2_tuned_recombination.py;
        the Hologic exposures and the Excel file of l2_tuned_recombination_multisystem.py (in the hologic folder).
Writes: system_comparison.xlsx in the data folder (sheets per_exposure, summary, notes).
Usage:  python3 compare_systems.py [data folder, default /Volumes/SanDisk/new_L2_measurements]
Every exposure is taken on its own; the summary gives median and range over the exposures of a group."""
import os, sys, re, importlib.util, numpy as np, pandas as pd, pydicom

DATA = sys.argv[1] if len(sys.argv) > 1 else '/Volumes/SanDisk/new_L2_measurements'
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ms', os.path.join(HERE, 'l2_tuned_recombination_multisystem.py')); MS = importlib.util.module_from_spec(spec); spec.loader.exec_module(MS)

rows = []
# --- Siemens: headers of every triplet, method numbers from the Excel file of the Siemens script
TS = MS.scan([os.path.join(DATA, k) for k in ('28kV', '30kV', '32kV') if os.path.isdir(os.path.join(DATA, k))], ('.ima',), ('4AFC',))
XS = pd.read_excel(os.path.join(DATA, 'L2_tuned_recombination_numbers.xlsx'), sheet_name=None)
tun = XS['tuned_images'].set_index('id'); chk = XS['siemens_weighting_check'].set_index('id')
for t in TS:
    hd = [pydicom.dcmread(t[k], stop_before_pixels=True) for k in ('L', 'H', 'R')]
    level = re.search(r'(LOW|NORMAL|HIGH)MAS', t['session']); r = dict(system='Siemens', id=t['id'], session=t['session'], dose_level_folder=level.group(1).lower() if level else '')
    r.update(MS.acquisition(*hd))
    if t['id'] in tun.index:
        u = tun.loc[t['id']]
        r.update(header_W=u['siemens_header_W'], system_effective_weight=u['siemens_effective_weight'], w_star=u['w_star'], system_left_in_share_percent=u['siemens_left_in_share_percent'],
                 tuned_weight=u['tuned_weight'], tuned_weight_measured=u['tuned_weight_measured'], tuned_left_in_share_percent=u['left_in_share_percent'], regions_used=u['regions_used'],
                 residual_0p2_system_mm2=u['residual_0p2_siemens_mm2'], residual_0p2_tuned_mm2=u['residual_0p2_tuned_mm2'])
    if t['id'] in chk.index:
        c = chk.loc[t['id']]
        r.update(check_interior_correlation=c['check_interior_correlation'], check_interior_difference_rms=c['check_interior_difference_rms'], interior_sd_R=c['check_interior_image_sd'],
                 check_share_reproduced_percent=c['check_share_reproduced_percent'], check_edge_15mm_correlation=c['check_edge_15mm_correlation'])
    rows.append(r)
# --- Hologic: everything from the summary of the multi-system script (it already holds the header data)
XH = pd.read_excel(os.path.join(DATA, 'hologic', MS.EXCEL), sheet_name=None)
ht = XH['tuned_images']; ht = ht[ht['strength'] == 1.0].set_index('id'); hc = XH['own_weighting_check'].set_index('id')
for i, u in ht.iterrows():
    r = dict(system='Hologic', id=i, session=u['session'], dose_level_folder='AEC (AutoFilter)')
    r.update({k: u[k] for k in MS.ACQ_KEYS if k in u})
    r.update(header_W='', system_effective_weight=u['system_effective_weight'], w_star=u['w_star'], system_left_in_share_percent=u['system_left_in_share_percent'],
             tuned_weight=u['tuned_weight'], tuned_weight_measured=u['tuned_weight_measured'], tuned_left_in_share_percent=u['left_in_share_percent'], regions_used=u['regions_used'],
             residual_0p2_system_mm2=u['residual_0p2_system_mm2'], residual_0p2_tuned_mm2=u['residual_0p2_tuned_mm2'],
             LE_raw_median_interior=u['LE_raw_median_interior'], HE_raw_median_interior=u['HE_raw_median_interior'])
    if i in hc.index:
        c = hc.loc[i]
        r.update(check_interior_correlation=c['check_interior_correlation'], check_interior_difference_rms=c['check_interior_difference_rms'], interior_sd_R=c['check_interior_image_sd'],
                 check_share_reproduced_percent=c['check_share_reproduced_percent'], check_edge_15mm_correlation=c['check_edge_15mm_correlation'])
    rows.append(r)
P = pd.DataFrame(rows)
for k in ('residual_0p2_system_mm2', 'residual_0p2_tuned_mm2', 'kV_LE', 'mAs_LE', 'mAs_HE', 'AGD_mGy_total'):
    P[k] = pd.to_numeric(P[k], errors='coerce')
P['mAs_ratio_HE_LE'] = (P['mAs_HE'] / P['mAs_LE']).round(3)
P['AGD_share_HE_percent'] = (100 * pd.to_numeric(P['AGD_mGy_HE'], errors='coerce') / P['AGD_mGy_total']).round(1)
P['system_weight_over_w_star'] = (P['system_effective_weight'] / P['w_star']).round(4)

# --- summary per system, LE kV (and Siemens dose level): median [lowest-highest]
NUM = ['mAs_LE', 'mAs_HE', 'mAs_ratio_HE_LE', 'exposure_ms_LE', 'exposure_ms_HE', 'HVL_mmAl_LE', 'HVL_mmAl_HE', 'AGD_mGy_LE', 'AGD_mGy_HE', 'AGD_mGy_total', 'AGD_share_HE_percent',
       'ESD_mGy_LE', 'ESD_mGy_HE', 'thickness_mm', 'force_N', 'header_W', 'system_effective_weight', 'w_star', 'system_left_in_share_percent', 'tuned_weight', 'tuned_left_in_share_percent',
       'residual_0p2_system_mm2', 'residual_0p2_tuned_mm2', 'interior_sd_R', 'check_interior_correlation', 'check_share_reproduced_percent', 'check_edge_15mm_correlation']
TXT = ['model', 'software', 'target_LE', 'filter_LE', 'filter_mm_LE', 'target_HE', 'filter_HE', 'filter_mm_HE', 'kV_HE', 'detector', 'pixel_mm', 'SID_mm', 'magnification', 'grid', 'bits_raw',
       'bits_recombined', 'AEC']


def cell(v):
    v = pd.to_numeric(v, errors='coerce').dropna()
    if not len(v):
        return ''
    f = (lambda x: '%.3g' % x) if v.abs().max() < 1e-2 else (lambda x: ('%.4f' % x).rstrip('0').rstrip('.'))
    return f(v.median()) if v.min() == v.max() else '%s [%s-%s]' % (f(v.median()), f(v.min()), f(v.max()))


S = []
P['group_kV'] = P['kV_LE'].round().astype(int)
for (sy, kv, lv), g in P.groupby(['system', 'group_kV', 'dose_level_folder'], sort=True):
    r = dict(system=sy, LE_kV=kv, dose_level=lv, exposures=len(g))
    r.update({k: ', '.join(sorted(set(str(x) for x in g[k] if str(x) not in ('', 'nan')))) for k in TXT if k in g})
    r.update({k: cell(g[k]) for k in NUM if k in g})
    S.append(r)
S = pd.DataFrame(S)
notes = pd.DataFrame([
    ('per_exposure', 'one line per exposure: acquisition data from the DICOM headers and the numbers of the tuned recombination (strength 1) and of the check at the system\'s own weighting'),
    ('summary', 'per system, LE tube voltage and (Siemens) dose-level folder: median [lowest-highest] over the exposures'),
    ('AGD_mGy', 'Organ Dose (0040,0316) of the header, converted from dGy: the manufacturer\'s own estimate of the average glandular dose, not recalculated'),
    ('ESD_mGy', 'Entrance Dose in mGy (0040,8302) of the header'),
    ('HVL_mmAl', 'Half Value Layer (0040,0314) of the header; Siemens writes none'),
    ('header_W', 'weight in the header of the recombined image; Siemens only (private tag); Hologic writes none, its weight is known only as measured'),
    ('system_effective_weight', 'weight on the LE image contained in the system\'s recombined image, measured from it at 0.06-0.15 cy/mm'),
    ('w_star', 'weight that cancels the L2 sphere background, measured from the raw LE and HE images at 0.06-0.15 cy/mm'),
    ('system_left_in_share_percent', 'share of the sphere contrast that the system\'s own recombination leaves, 1 - weight / w_star'),
    ('tuned_left_in_share_percent', 'share left by the tuned recombination at strength 1 (breast target measured on Siemens; used unchanged for Hologic)'),
    ('residual_0p2 ... mm2', 'power of the log-signal residual at 0.2-0.3 cy/mm (noise included), under the system\'s weighting and tuned'),
    ('interior_sd_R', 'standard deviation of the system\'s recombined image more than 15 mm inside the outline, in its own grey values (the grey scales of the two systems differ: gain about 7000 (Siemens) and 3350 (Hologic) grey values per log unit)'),
    ('check_...', 'our recombination of the raw images at the system\'s own weighting against the system\'s image: interior correlation, share of its power reproduced at 0.06-1 cy/mm, correlation in the outer 15 mm'),
    ('LE/HE_raw_median_interior', 'median raw value in the interior (Hologic only here): a dark offset, if present, matters more at low values; not corrected on either system'),
], columns=['sheet / quantity', 'meaning'])
out = os.path.join(DATA, 'system_comparison.xlsx')
with pd.ExcelWriter(out) as w:
    cols = ['system', 'id', 'session', 'dose_level_folder'] + [c for c in P.columns if c not in ('system', 'id', 'session', 'dose_level_folder', 'group_kV')]
    P[cols].to_excel(w, sheet_name='per_exposure', index=False); S.to_excel(w, sheet_name='summary', index=False); notes.to_excel(w, sheet_name='notes', index=False)
print('written', out, '| exposures: Siemens %d, Hologic %d' % ((P.system == 'Siemens').sum(), (P.system == 'Hologic').sum()))
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_colwidth', 40)
print(S.set_index(['system', 'LE_kV', 'dose_level']).T.to_string())
