#!/usr/bin/env python3
"""Figures of the method document, all from the real runs:
  - the Excel file written by l2_tuned_recombination.py on the SanDisk (112 exposures, strength 1.0, and the check),
  - tissue_target.json (per-view spectra of the nine breast views; numbers only),
  - one example exposure processed here with the functions of l2_tuned_recombination.py itself,
  - the Excel file written by l2_tuned_recombination_multisystem.py in the hologic folder (15 Hologic exposures).
Usage: python3 make_figures.py [data folder, default /Volumes/SanDisk/new_L2_measurements]"""
import os, sys, json, glob, shutil, importlib.util, warnings, numpy as np, pandas as pd, pydicom
warnings.filterwarnings('ignore')
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, 'figures'); TUN = os.path.dirname(HERE)
DATA = sys.argv[1] if len(sys.argv) > 1 else '/Volumes/SanDisk/new_L2_measurements'
spec = importlib.util.spec_from_file_location('l2', os.path.join(TUN, 'l2_tuned_recombination.py')); M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
plt.rcParams.update({'font.size': 8.5, 'axes.titlesize': 9, 'legend.fontsize': 7.5, 'axes.grid': True, 'grid.alpha': 0.25, 'savefig.bbox': 'tight', 'figure.dpi': 100})
FC = M.FC; EDGES = np.array(M.EDGES); KV = (28, 30, 32); COL = {28: '#1f77b4', 30: '#2ca02c', 32: '#d62728'}
S = pd.read_excel(os.path.join(DATA, 'L2_tuned_recombination_numbers.xlsx'), sheet_name=None); TI, CU, CH, CC = S['tuned_images'], S['tuned_curves'], S['siemens_weighting_check'], S['check_curves']
kv_of = dict(zip(TI['id'], TI['kV'])); W_of = dict(zip(TI['id'], TI['siemens_header_W'])); CU['kV'] = CU['id'].map(kv_of); CU['W'] = CU['id'].map(W_of); CC['kV'] = CC['id'].map(kv_of)
CU['band'] = CU.groupby('id').cumcount(); CC['band'] = CC.groupby('id').cumcount()
def band_stat(df, col, kv):
    a = df[df['kV'] == kv].pivot(index='id', columns='band', values=col).to_numpy(float); return np.nanmedian(a, 0), np.nanmin(a, 0), np.nanmax(a, 0)
def logx(ax, label=True):
    ax.set_xscale('log'); ax.set_xlim(0.06, 5.9)
    if label: ax.set_xlabel('spatial frequency (cy/mm)')
out = {}
# ------------------------------------------------------------------ Siemens' weighting as measured from R
fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
for kv in KV:
    med, lo, hi = band_stat(CU, 'siemens_effective_weight', kv); r2 = band_stat(CU, 'siemens_fit_R2', kv)[0]; ok = r2 >= 0.2; W = float(TI[TI['kV'] == kv]['siemens_header_W'].median())
    ax[0, 0].fill_between(FC[ok], lo[ok], hi[ok], color=COL[kv], alpha=0.2, lw=0); ax[0, 0].plot(FC[ok], med[ok], 'o-', ms=3, color=COL[kv], label='%d kV: measured $w_\\mathrm{eff}(f)$' % kv); ax[0, 0].axhline(W, color=COL[kv], ls=':', lw=1)
    ax[0, 1].plot(FC[ok], (med / W)[ok], 'o-', ms=3, color=COL[kv], label='%d kV' % kv)
    a = band_stat(CU, 'siemens_gain_alpha', kv); ax[1, 0].fill_between(FC[ok], a[1][ok], a[2][ok], color=COL[kv], alpha=0.2, lw=0); ax[1, 0].plot(FC[ok], a[0][ok], 'o-', ms=3, color=COL[kv], label='%d kV' % kv)
    ax[1, 1].plot(FC, r2, 'o-', ms=3, color=COL[kv], label='%d kV' % kv)
    out['weff_%d' % kv] = dict(W=W, weff=np.round(med, 4).tolist(), alpha=np.round(a[0], 0).tolist(), R2=np.round(r2, 3).tolist())
tt = json.load(open(os.path.join(TUN, 'tissue_target.json'))); KB = np.array(tt['K'])
ax[0, 1].plot(FC[:9], KB[:9], 'k--', lw=1.2, label='$K(f)$ used for the breast target')
ax[0, 0].set_title("(a) effective weight $w_\\mathrm{eff}(f)$; dotted: header $W$"); ax[0, 0].set_ylabel('weight'); ax[0, 1].set_title('(b) relative to the header value, $w_\\mathrm{eff}(f)/W$')
ax[1, 0].set_title("(c) gain $\\alpha(f)$ (grey values per unit log signal)"); ax[1, 1].set_title('(d) fit quality $R^2$'); ax[1, 1].axhline(0.2, color='k', ls=':', lw=0.8)
for a_ in ax.ravel(): logx(a_); a_.legend()
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'siemens_weighting.pdf')); plt.close(fig)
# ------------------------------------------------------------------ cancelling weight and coherence
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
for kv in KV:
    med, lo, hi = band_stat(CU, 'w_star', kv); ax[0].fill_between(FC, lo, hi, color=COL[kv], alpha=0.2, lw=0); ax[0].plot(FC, med, 'o-', ms=3, color=COL[kv], label='%d kV' % kv)
    c = band_stat(CU, 'coherence', kv); ax[1].fill_between(FC, c[1], c[2], color=COL[kv], alpha=0.2, lw=0); ax[1].plot(FC, c[0], 'o-', ms=3, color=COL[kv], label='%d kV' % kv); out['ws_%d' % kv] = dict(ws=np.round(med, 4).tolist(), coh=np.round(c[0], 4).tolist())
ax[0].set_title('(a) cancelling weight $w^*(f)=S_{HL}/S_{LL}$'); ax[0].set_ylabel('weight'); ax[1].set_title('(b) coherence $\\gamma^2(f)$ of $\\Lambda_H$ and $\\Lambda_L$')
for a_ in ax: logx(a_); a_.legend()
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'wstar_coherence.pdf')); plt.close(fig)
# ------------------------------------------------------------------ breast target
NS = M.noise_shape(0.085); views = tt['views']; PD = []; TT = []; NO = []; UU = []
for v in views:
    LL, HH, HL = (np.array(v[k]) for k in ('LL', 'HH', 'HL')); WK = v['W'] * KB; pd_ = HH - 2 * WK * HL + WK ** 2 * LL; no = v['HH_raw'][-1] * NS + WK ** 2 * LL[-1] * NS
    PD.append(pd_); NO.append(no); TT.append(np.clip(pd_ - no, 0, None))
    UU.append((HL / LL - WK) ** 2 * LL)          # part of the residual that is tissue structure left by the weight: S_LL (W K - w*)^2
PD, TT, NO, UU = np.array(PD), np.array(TT), np.array(NO), np.array(UU); LOB = FC < 0.3; USH = UU[:, LOB].sum(1) / TT[:, LOB].sum(1)     # its share of the target, 0.06-0.3 cy/mm, per view
out['breast'] = dict(T=np.median(TT, 0).tolist(), PD=np.median(PD, 0).tolist(), noise=np.median(NO, 0).tolist(), T_lo=TT.min(0).tolist(), T_hi=TT.max(0).tolist(), n_roi=[v['n'] for v in views], W=[v['W'] for v in views], kv=[v['kv'] for v in views], thk=[v['thk'] for v in views], K=KB.tolist(), NS=NS.tolist(),
                     check_T=float(np.abs(np.median(TT, 0) - np.array(M.BREAST['T'])).max()), check_PD=float(np.abs(np.median(PD, 0) / np.array(M.BREAST['PD']) - 1).max()))
out['breast'].update(U=np.median(UU, 0).tolist(), U_share_per_view=USH.tolist(), U_share_min_median_max=[float(USH.min()), float(np.median(USH)), float(USH.max())],
                     U_share_of_band_medians=(np.median(UU, 0)[:8] / np.median(TT, 0)[:8]).tolist(), ws_coarse=[float((np.array(v['HL']) / np.array(v['LL']))[:2].mean()) for v in views])
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.9))
for p in PD: ax[0].plot(FC, p, color='0.6', lw=0.7)
ax[0].plot(FC, np.median(PD, 0), 'o-', ms=3, color='#1f77b4', label='median of 9 views'); ax[0].plot(FC, np.median(NO, 0), 's--', ms=3, color='#ff7f0e', label='quantum noise (median)'); ax[0].set_title('(a) residual $P_D(f)$, each view (grey)')
for p in TT: ax[1].plot(FC[:8], np.where(p[:8] > 0, p[:8], np.nan), color='0.6', lw=0.7)
ax[1].plot(FC[:8], np.median(TT, 0)[:8], 'o-', ms=3, color='#1f77b4', label='target $T(f)$: median'); ax[1].plot(FC[:8], np.median(UU, 0)[:8], 's--', ms=3, color='#d62728', label='$U(f)$: left by the weight'); ax[1].axvspan(0.06, 0.3, color='#2ca02c', alpha=0.12, lw=0, label='bands used for $\\rho_0$'); ax[1].set_title('(b) structured residual $P_D-$noise')
ax[2].plot(FC, NS, 'o-', ms=3, color='k'); ax[2].set_title('(c) noise shape $g(f)$ (last band = 1)'); ax[2].set_ylim(0, None)
for a_ in ax[:2]: a_.set_yscale('log'); a_.set_ylabel('mm$^2$'); a_.set_ylim(a_.get_ylim()[0], a_.get_ylim()[1] * 12); a_.legend(loc='upper right')      # room for the legend above the curves
for a_ in ax: logx(a_)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'breast_target.pdf')); plt.close(fig)
# ------------------------------------------------------------------ example exposure with the script's own functions
EXDIR = next(d for d in (os.path.join(DATA, '28kV', 'Measurements', 'CEM_NEWL2_MASAR_NORMALMAS_28KV_02'), os.path.join(DATA, '28kV', 'CEM_NEWL2_MASAR_NORMALMAS_28KV_02')) if os.path.isdir(d))   # folder layout before / after 19 Nov 2026
T = M.scan([EXDIR], '.ima', ('4AFC',)); t = T[0]; out['example'] = dict(id=t['id'], kv=t['kv'], mAsL=t['mAsL'], mAsH=t['mAsH'], W=t['W'])
dsL, dsH, dsR = (pydicom.dcmread(t[k]) for k in ('L', 'H', 'R')); L, H, R = (d.pixel_array.astype(np.float64) for d in (dsL, dsH, dsR)); px = float(dsL.ImagerPixelSpacing[0]); vmax = 2 ** int(dsR.BitsStored) - 1
lL, lH = -np.log(np.clip(L, 1, None)), -np.log(np.clip(H, 1, None)); fp = M.footprint(H); dist = ndi.distance_transform_edt(fp) * px; inner = dist > M.MARGIN_MM; valid = (R > 0) & (R < vmax)
rois = M.Rois(px); pos = rois.positions(inner & valid); nshape = M.noise_shape(px); St = rois.stats(lL, lH, R, pos); keep = M.select(St); m = M.measure(St, keep, rois.norm, nshape); d1 = M.design(m, 1.0)
mid = (FC >= 0.1) & (FC < 1.0); ws = St['HL'][keep].sum(0) / St['LL'][keep].sum(0); le = np.log(np.clip((St['HH'] - 2 * ws * St['HL'] + ws ** 2 * St['LL'])[:, mid].sum(1), 1e-30, None)); med = np.median(le[keep]); mad = 1.4826 * np.median(np.abs(le[keep] - med))
out['example'].update(rois=len(pos), kept=int(keep.sum()), left_out=int((~keep).sum()), rho0=d1['rho0'], ws=np.round(m['ws'], 4).tolist(), weff=np.round(m['w_eff'], 4).tolist(), wt=np.round(d1['wt'], 4).tolist(), C=m['C'].tolist(), E=m['E'].tolist(), LL=m['LL'].tolist(), HH=m['HH'].tolist(), HL=m['HL'].tolist(), alpha=m['alpha'].tolist(), beta=m['beta'].tolist(), R2=m['R2'].tolist(),
                      open_level=float(np.percentile(H, 99.5)), footprint_area_cm2=float(fp.sum() * px * px / 100), inner_area_cm2=float((inner & valid).sum() * px * px / 100), LE_median=float(np.median(L[inner])), HE_median=float(np.median(H[inner])), R_median=float(np.median(R[inner & valid])), R_sd=float(R[inner & valid].std()), threshold_sd=3.5)
fig, ax = plt.subplots(1, 3, figsize=(7.4, 3.3), gridspec_kw=dict(width_ratios=[1, 1, 1.15]))
k = 4; img = M.binned(R, k); rm = np.median(R[inner & valid]); ext = (0, R.shape[1] * px, R.shape[0] * px, 0)
for a_, ttl in ((ax[0], '(a) outline and interior'), (ax[1], '(b) ROIs used and left out')):
    a_.imshow(img, cmap='gray', vmin=rm - 400, vmax=rm + 400, extent=ext); a_.contour(np.arange(fp.shape[1])[::k] * px, np.arange(fp.shape[0])[::k] * px, fp[::k, ::k].astype(float), [0.5], colors='#ffd400', linewidths=0.8)
    a_.contour(np.arange(fp.shape[1])[::k] * px, np.arange(fp.shape[0])[::k] * px, inner[::k, ::k].astype(float), [0.5], colors='#00d0ff', linewidths=0.8); a_.set_title(ttl); a_.set_xlabel('mm'); a_.grid(False)
cy = np.array([p[0] + M.N / 2 for p in pos]) * px; cx = np.array([p[1] + M.N / 2 for p in pos]) * px
ax[1].plot(cx[keep], cy[keep], '.', ms=2.5, color='#33dd33'); ax[1].plot(cx[~keep], cy[~keep], 'x', ms=4, color='#ff3030', mew=1)
ax[2].hist(le[keep], bins=40, color='#33aa33', alpha=0.8, label='used'); ax[2].hist(le[~keep], bins=20, color='#ff3030', alpha=0.8, label='left out'); ax[2].axvline(med + 3.5 * mad, color='k', ls='--', lw=1, label='median + 3.5 robust sd')
ax[2].set_title('(c) residual power at $w^*$'); ax[2].set_xlabel('ln(power)'); ax[2].legend(); ax[2].set_ylabel('regions')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'regions.png'), dpi=200); plt.close(fig)
# ------------------------------------------------------------------ design of the tuned weight (example)
fig, ax = plt.subplots(1, 3, figsize=(7.4, 3.9)); Tb = np.array(M.BREAST['T']); lo = FC < M.FADE[0]
ax[0].plot(FC, m['HH'], 'o-', ms=3, color='0.4', label='$S_{HH}$'); ax[0].plot(FC, m['C'], 'o-', ms=3, color='k', label='$C=S_{HL}^2/S_{LL}$ (cancellable)'); ax[0].plot(FC, np.where(m['E'] > 0, m['E'], np.nan), 's-', ms=3, color='#ff7f0e', label='$E$: non-noise residual at $w^*$')
ax[0].plot(FC[:8], Tb[:8], 'o-', ms=3, color='#1f77b4', label='breast target $T$'); ax[0].axvspan(0.06, 0.3, color='#2ca02c', alpha=0.12, lw=0); ax[0].set_yscale('log'); ax[0].set_ylabel('mm$^2$'); ax[0].set_title('(a) ingredients (example exposure)')
for s_, ls in ((0.5, ':'), (1.0, '-'), (2.0, '--')):
    ds_ = M.design(m, s_); ax[1].plot(FC, 100 * ds_['rho'], 'o' + ls, ms=3, color='#2ca02c', label='strength %g: $\\rho_0$ = %.1f %%' % (s_, 100 * ds_['rho0']))
ax[1].plot(FC[m['R2'] >= 0.2], 100 * (1 - m['w_eff'] / m['ws'])[m['R2'] >= 0.2], 'o-', ms=3, color='#d62728', label='Siemens: $1-w_\\mathrm{eff}/w^*$'); ax[1].set_ylim(-1, 40); ax[1].set_ylabel('%'); ax[1].set_title('(b) share of the sphere contrast left in')
ok = m['R2'] >= 0.2; ax[2].plot(FC[ok], m['w_eff'][ok], 'o-', ms=3, color='#d62728', label="Siemens $w_\\mathrm{eff}$"); ax[2].plot(FC, m['ws'], 'o-', ms=3, color='k', label='$w^*$'); ax[2].plot(FC, d1['wt'], 'o-', ms=3, color='#2ca02c', label='tuned $w_t$, strength 1')
ax[2].set_title('(c) weights on $\\Lambda_L$'); ax[2].set_ylabel('weight')
for a_ in ax: logx(a_); a_.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'design.pdf')); plt.close(fig)
# ------------------------------------------------------------------ the correction filter c(f), one normal-dose exposure per kV
f = np.geomspace(0.004, 5.88, 500); g = lambda x: x ** 2 / (x ** 2 + M.F0 ** 2); fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9)); ex = {}
for kv in KV:
    i0 = sorted(TI[(TI['kV'] == kv) & TI['session'].str.contains('NORMAL')]['id'])[0]; c_ = CU[CU['id'] == i0].sort_values('band'); al = c_['siemens_gain_alpha'].to_numpy(); be = -al * c_['siemens_effective_weight'].to_numpy(); ok = c_['siemens_fit_R2'].to_numpy() >= 0.2
    cb = np.where(ok, al * c_['tuned_weight'].to_numpy() + be, 0.0); filt = np.interp(np.log10(np.clip(f, FC[0], FC[-1])), np.log10(FC), cb); filt[f < FC[0]] = cb[0] * g(f[f < FC[0]]) / g(FC[0]); ex[kv] = i0
    ax[0].plot(f, filt, color=COL[kv], label='%d kV' % kv); ax[0].plot(FC, cb, 'o', ms=3, color=COL[kv]); out['c_%d' % kv] = dict(id=i0, c=np.round(cb, 1).tolist())
ax[0].set_xscale('log'); ax[0].set_xlabel('spatial frequency (cy/mm)'); ax[0].set_ylabel('grey values per unit log signal'); ax[0].set_title('(a) correction filter $c(f)$ applied to $\\Lambda_L$'); ax[0].legend(); ax[0].axvline(FC[0], color='k', ls=':', lw=0.8)
ax[1].plot(f, g(f), 'k'); ax[1].set_xscale('log'); ax[1].set_xlim(0.004, 1); ax[1].set_xlabel('spatial frequency (cy/mm)'); ax[1].set_title('(b) large-scale roll-off $f^2/(f^2+f_0^2)$, $f_0$ = 0.05 cy/mm'); ax[1].axvline(FC[0], color='k', ls=':', lw=0.8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'filter.pdf')); plt.close(fig)
# ------------------------------------------------------------------ residual spectra against the breasts, per kV (all exposures)
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.9), sharey=True)
for a_, kv in zip(ax, KV):
    a_.fill_between(FC, M.BREAST['PD_LO'], M.BREAST['PD_HI'], color='#1f77b4', alpha=0.18, lw=0, label='real breasts (9 views)'); a_.plot(FC, M.BREAST['PD'], color='#1f77b4', lw=1, label='real breasts, median')
    for col, c_, lab in (('residual_siemens_mm2', '#d62728', "L2, Siemens' weighting"), ('residual_full_cancellation_mm2', 'k', 'L2, full cancellation'), ('residual_tuned_mm2', '#2ca02c', 'L2, tuned (strength 1)')):
        md, l_, h_ = band_stat(CU, col, kv); a_.fill_between(FC, l_, h_, color=c_, alpha=0.2, lw=0); a_.plot(FC, md, 'o-', ms=2.5, color=c_, label=lab)
    a_.set_yscale('log'); logx(a_); a_.set_title('%d kV (%d exposures)' % (kv, (TI['kV'] == kv).sum()))
ax[0].set_ylabel('residual power (mm$^2$)'); ax[0].legend(fontsize=6.5)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'residual_spectra.pdf')); plt.close(fig)
# ------------------------------------------------------------------ summary over all exposures
TI['level'] = TI['session'].str.extract(r'(NORMAL|LOW|HIGH)MAS')[0]; fig, ax = plt.subplots(1, 3, figsize=(7.4, 3.5)); x0 = 0
for kv in KV:
    for lv in ('LOW', 'NORMAL', 'HIGH'):
        s_ = TI[(TI['kV'] == kv) & (TI['level'] == lv)]; x = x0 + np.linspace(-0.3, 0.3, len(s_))
        ax[0].plot(x, s_['siemens_left_in_share_percent'], '.', ms=4, color='#d62728'); ax[0].plot(x, s_['left_in_share_percent'], '.', ms=4, color='#2ca02c')
        ax[1].plot(x, s_['siemens_effective_weight'], '.', ms=4, color='#d62728'); ax[1].plot(x, s_['w_star'], '.', ms=4, color='k'); ax[1].plot(x, s_['tuned_weight'], '.', ms=4, color='#2ca02c'); x0 += 1
lab = ['%s\n%s' % (lv[0], kv if lv == 'NORMAL' else '') for kv in KV for lv in ('LOW', 'NORMAL', 'HIGH')]
for a_ in ax[:2]: a_.set_xticks(range(9)); a_.set_xticklabels(lab, fontsize=7); a_.set_xlabel('exposure level (L, N, H) and LE kV')
ax[0].set_title('(a) sphere contrast left in (%)'); ax[0].plot([], [], '.', color='#d62728', label='Siemens'); ax[0].plot([], [], '.', color='#2ca02c', label='tuned, strength 1'); ax[0].legend(loc='center right'); ax[0].set_ylim(0, 32)
ax[1].set_title('(b) weights at 0.06-0.15 cy/mm'); ax[1].plot([], [], '.', color='#d62728', label="Siemens $w_\\mathrm{eff}$"); ax[1].plot([], [], '.', color='k', label='$w^*$'); ax[1].plot([], [], '.', color='#2ca02c', label='tuned $w_t$'); ax[1].legend(loc='upper left', fontsize=6.5)
ax[2].plot(TI['tuned_weight'], TI['tuned_weight_measured'] - TI['tuned_weight'], '.', ms=4, color='#2ca02c'); ax[2].axhline(0, color='k', lw=0.8); ax[2].set_title('(c) measured minus designed $w_t$'); ax[2].set_xlabel('designed $w_t$')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'summary.pdf')); plt.close(fig)
out['summary'] = TI.groupby(['kV', 'level']).agg(n=('id', 'count'), mAsL=('mAs_LE', 'median'), mAsH=('mAs_HE', 'median'), W=('siemens_header_W', 'median'), weff=('siemens_effective_weight', 'median'), ws=('w_star', 'median'), wt=('tuned_weight', 'median'), left=('left_in_share_percent', 'median'), sleft=('siemens_left_in_share_percent', 'median'), used=('regions_used', 'median'), leftout=('regions_left_out', 'median')).round(4).reset_index().to_dict('records')
out['tuned_dev'] = dict(max=float((TI['tuned_weight_measured'] - TI['tuned_weight']).abs().max()), mean=float((TI['tuned_weight_measured'] - TI['tuned_weight']).mean()))
out['res02'] = {int(kv): dict(siemens=float(TI[TI['kV'] == kv]['residual_0p2_siemens_mm2'].astype(float).median()), tuned=float(TI[TI['kV'] == kv]['residual_0p2_tuned_mm2'].astype(float).median())) for kv in KV}; out['res02']['breast'] = float(M.BREAST['PD'][3])
# ------------------------------------------------------------------ check at Siemens' weighting
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.9))
for kv in KV:
    md, l_, h_ = band_stat(CC, 'share_reproduced_percent', kv); ax[0].fill_between(FC, l_, h_, color=COL[kv], alpha=0.2, lw=0); ax[0].plot(FC, md, 'o-', ms=3, color=COL[kv], label='%d kV' % kv); out['share_%d' % kv] = np.round(md, 1).tolist()
    r2 = band_stat(CC, 'fit_R2_in_regions', kv)[0] >= 0.2; a1 = band_stat(CC, 'siemens_effective_weight', kv)[0]; a2 = band_stat(CC, 'ours_effective_weight', kv)[0]; ax[1].plot(FC[r2], a1[r2], 'o-', ms=3, color=COL[kv], label="%d kV, Siemens' $R$" % kv); ax[1].plot(FC[r2], a2[r2], 'x', ms=6, color=COL[kv], label='%d kV, ours' % kv)
    for col, ls in (('power_siemens', '-'), ('power_ours', '--'), ('power_difference', ':')): ax[2].plot(FC, band_stat(CC, col, kv)[0], ls, color=COL[kv])
ax[0].set_title('(a) share reproduced (%)'); ax[0].set_ylim(-3, 103); ax[1].set_title('(b) effective weight'); ax[2].set_title('(c) power: Siemens (solid), ours (dashed),\ndifference (dotted)'); ax[2].set_yscale('log'); ax[2].set_ylabel('grey values$^2$ mm$^2$')
for a_ in ax: logx(a_)
ax[0].legend(); ax[1].legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'check.pdf')); plt.close(fig)
CH['kV2'] = CH['kV']; out['check'] = CH.groupby('kV').agg(n=('id', 'count'), corr=('check_interior_correlation', 'median'), corr_lo=('check_interior_correlation', 'min'), corr_hi=('check_interior_correlation', 'max'), rms=('check_interior_difference_rms', 'median'), sd=('check_interior_image_sd', 'median'), share=('check_share_reproduced_percent', 'median'), share_lo=('check_share_reproduced_percent', 'min'), share_hi=('check_share_reproduced_percent', 'max'), edge=('check_edge_15mm_correlation', 'median'), edge_lo=('check_edge_15mm_correlation', 'min'), edge_hi=('check_edge_15mm_correlation', 'max')).round(3).reset_index().to_dict('records')
# ------------------------------------------------------------------ the comparison images of the real run for the example exposure
for sub, name in (('28kV_strength_1.0', 'comparison_tuned.png'), ('28kV_siemens_weighting', 'comparison_check.png')):
    src = glob.glob(os.path.join(DATA, sub, t['session'], t['id'] + '_comparison.png')); shutil.copyfile(src[0], os.path.join(FIG, name)) if src else print('missing', sub)
from PIL import Image
for name in ('comparison_tuned.png', 'comparison_check.png'):
    im = Image.open(os.path.join(FIG, name)); im.convert('RGB').save(os.path.join(FIG, name.replace('.png', '.jpg')), quality=88); print(name, im.size)
# ------------------------------------------------------------------ Hologic: the same method (l2_tuned_recombination_multisystem.py), its own run
HX = os.path.join(DATA, 'hologic', 'L2_tuned_recombination_multisystem_numbers.xlsx')
if os.path.exists(HX):
    SH = pd.read_excel(HX, sheet_name=None); HT, HU, HC, HK = SH['tuned_images'], SH['tuned_curves'], SH['own_weighting_check'], SH['check_curves']
    HT = HT[HT['strength'] == 1].copy(); HU = HU[HU['strength'] == 1].copy(); hkv = dict(zip(HT['id'], HT['kV']))
    for df in (HU, HK):
        df['kV'] = df['id'].map(hkv); df['band'] = df.groupby('id').cumcount()
    FH = HU[HU['id'] == HU['id'].iloc[0]]['band_centre_cy_mm'].to_numpy(float); HKV = (29, 30); HCOL = {29: '#9467bd', 30: '#e377c2'}
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.8))
    for kv in HKV:
        we = band_stat(HU, 'system_effective_weight', kv)[0]; ws = band_stat(HU, 'w_star', kv)[0]; ok = band_stat(HU, 'system_fit_R2', kv)[0] >= 0.2
        ax[0, 0].plot(FH[ok], (we / we[:2].mean())[ok], 'o-', ms=3, color=HCOL[kv], label='Hologic %d kV' % kv)
        ax[0, 1].plot(FH[ok], 100 * (1 - we / ws)[ok], 'o-', ms=3, color=HCOL[kv], label='Hologic %d kV, own' % kv)
        ax[0, 1].plot(FH, 100 * band_stat(HU, 'left_in_share_rho', kv)[0], 'o:', ms=3, color=HCOL[kv], label='Hologic %d kV, tuned' % kv)
        md, l_, h_ = band_stat(HK, 'share_reproduced_percent', kv); ax[1, 0].fill_between(FH, l_, h_, color=HCOL[kv], alpha=0.2, lw=0); ax[1, 0].plot(FH, md, 'o-', ms=3, color=HCOL[kv], label='Hologic %d kV' % kv)
    for kv in (28, 30):
        we = band_stat(CU, 'siemens_effective_weight', kv)[0]; ws = band_stat(CU, 'w_star', kv)[0]; ok = band_stat(CU, 'siemens_fit_R2', kv)[0] >= 0.2
        ax[0, 0].plot(FC[ok], (we / we[:2].mean())[ok], 's--', ms=3, color=COL[kv], label='Siemens %d kV' % kv)
        ax[0, 1].plot(FC[ok], 100 * (1 - we / ws)[ok], 's--', ms=3, color=COL[kv], label='Siemens %d kV, own' % kv)
        ax[1, 0].plot(FC, band_stat(CC, 'share_reproduced_percent', kv)[0], 's--', ms=3, color=COL[kv], label='Siemens %d kV' % kv)
    gH = HU[HU['id'] == HU['id'].iloc[0]]['noise_shape_g'].to_numpy(float)
    ax[1, 1].plot(FH, gH, 'o-', ms=3, color='#9467bd', label='Hologic: measured in the open field'); ax[1, 1].plot(FC, M.noise_shape(0.085), 's--', ms=3, color='k', label='Siemens: from the published MTF')
    ax[0, 0].set_title('(a) effective weight, relative to 0.06-0.15 cy/mm'); ax[0, 0].set_ylim(0, 1.15)
    ax[0, 1].set_title('(b) share of the sphere contrast left in (%)'); ax[0, 1].set_ylim(-3, 100)
    ax[1, 0].set_title("(c) check: share of the system's image reproduced (%)"); ax[1, 0].set_ylim(-3, 103)
    ax[1, 1].set_title('(d) noise shape $g$ (last band = 1)'); ax[1, 1].set_ylim(0, None)
    for a_ in ax.ravel(): logx(a_, label=a_ in ax[1]); a_.set_xlim(0.06, 7.2); a_.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hologic.pdf')); plt.close(fig)
    num = lambda df, c: pd.to_numeric(df[c], errors='coerce')
    def rng(df, c, nd=4):
        v = num(df, c); return [round(float(v.median()), nd), round(float(v.min()), nd), round(float(v.max()), nd)]
    out['hologic'] = dict(per_kV={}, check={}, n=int(len(HT)))
    for kv in HKV:
        g = HT[HT['kV'] == kv]; c = HC[HC['kV'] == kv]
        out['hologic']['per_kV'][kv] = dict(n=int(len(g)), weff=rng(g, 'system_effective_weight'), ws=rng(g, 'w_star'), wt=rng(g, 'tuned_weight'), own_share=rng(g, 'system_left_in_share_percent', 2),
                                             tuned_share=rng(g, 'left_in_share_percent', 2), rois=[int(num(g, 'regions_used').median()), int(num(g, 'regions_left_out').median())],
                                             res02_system=float(num(g, 'residual_0p2_system_mm2').median()), res02_tuned=float(num(g, 'residual_0p2_tuned_mm2').median()),
                                             E=float(num(g, 'E_sum_0p06_0p3_mm2').median()), C=float(num(g, 'C_sum_0p06_0p3_mm2').median()),
                                             weight_by_band=np.round(band_stat(HU, 'system_effective_weight', kv)[0], 4).tolist(), alpha_by_band=np.round(band_stat(HU, 'system_gain_alpha', kv)[0], 0).tolist(),
                                             R2_by_band=np.round(band_stat(HU, 'system_fit_R2', kv)[0], 3).tolist())
        out['hologic']['check'][kv] = dict(n=int(len(c)), corr=rng(c, 'check_interior_correlation'), rms=round(float(num(c, 'check_interior_difference_rms').median()), 1),
                                           sd=round(float(num(c, 'check_interior_image_sd').median()), 1), share=rng(c, 'check_share_reproduced_percent', 1), edge=rng(c, 'check_edge_15mm_correlation', 3),
                                           weight_diff_max=round(float((num(c, 'check_ours_effective_weight') - num(c, 'system_effective_weight')).abs().max()), 4))
    out['hologic']['tuned_dev_max'] = round(float((num(HT, 'tuned_weight_measured') - num(HT, 'tuned_weight')).abs().max()), 4)
    out['hologic']['noise_shape'] = gH.tolist(); out['hologic']['band_centres'] = FH.tolist()
    src = os.path.join(DATA, 'hologic', 'Hologic_29kV_strength_1.0', 'Hologic_CEM_Oct_22_2026', 'MG14_29kV_comparison.png')
    if os.path.exists(src):
        from PIL import Image
        Image.open(src).convert('RGB').save(os.path.join(FIG, 'comparison_hologic.jpg'), quality=88)
else:
    print('no Hologic results found:', HX)
json.dump(out, open(os.path.join(HERE, 'numbers.json'), 'w'), indent=1); print('figures:', sorted(os.listdir(FIG)))
print('breast table check: largest |T - BREAST[T]| %.2e mm2; PD median within %.2f %%' % (out['breast']['check_T'], 100 * out['breast']['check_PD']))
