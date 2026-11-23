#!/usr/bin/env bash
# Full chain of the GATE 10 Revelation model for one case = phantom model + background + measured reference set:
#   PHANTOM=L2   ./run_all.sh   real L2 container with the new (July 2026) plate, PMMA spheres in paraffin, compared
#                               with ONE July pair (28 kV high mAs, unit 10159; PAIR, default HIGHMAS_28KV2_00),
#                               simulated with that pair's exact mAs; another pair = another run with PAIR=<id>
#   PHANTOM=L2 BACKGROUND=PA12_MCT ./run_all.sh   the REPORT §10.4 background (PA12 + CaCO3 spheres in MCT oil) in the
#                               real L2; at the end its measurement predicted from the current L2 (redesign_compare.py)
#   PHANTOM=L2 BACKGROUND=gland_adipose ./run_all.sh   tissue reference (glandular spheres in adipose tissue, not
#                               buildable): the target the PA12_MCT prediction is judged against
#   PHANTOM=L2v2 ./run_all.sh   suggested redesign (default; background PA12_MCT, reference pair 0035)
# Needs: pip install opengate==10.1.1 (brings Geant4 data on first run, SpekPy, uproot); numpy, scipy, matplotlib, pydicom.
# The phantom geometry is already generated (phantom/ for L2v2, phantom_L2/ for L2); to regenerate it:
#   python phantom_meshes.py --phantom L2   and   python sphere_packing.py --phantom L2   (L2v2: needs CadQuery, ~45 min)
# July set: the exposure values (per pair, no averages) come from the DICOM headers, in data/measured_july28high.json;
#   python measured_set.py --set july28high [--root "<folder 'New L2 Measurements'>"]   rewrites it from the images.
# The comparison at the end needs the July images: in the project folder (session folders, as now), on the SanDisk, or
# under ROOT.
# Look at the model first:  python revelation_gate.py --phantom L2 --visu vrml_file_only  (then view_vrml.py)
set -euo pipefail
cd "$(dirname "$0")"

export REVELATION_PHANTOM=${PHANTOM:-L2v2}         # L2v2 or L2
[ -n "${BACKGROUND:-}" ] && export REVELATION_BACKGROUND=$BACKGROUND   # default per phantom: PA12_MCT (L2v2), PMMA_paraffin (L2)
[ -n "${SET:-}" ] && export REVELATION_SET=$SET    # default per phantom: 0035 (L2v2), july28high (L2)
[ -n "${PAIR:-}" ] && export REVELATION_PAIR=$PAIR # the one measured pair; default per set (July: HIGHMAS_28KV2_00)
ROOT=${ROOT:-}                                     # folder 'New L2 Measurements' (July images), if not on the SanDisk
THREADS=${THREADS:-8}
RUNS=${RUNS:-5}           # runs per beam (seeds 1..RUNS); scatter statistics add up across runs
N_LE=${N_LE:-2e7}         # photons per LE run, total over all threads (5 runs: 1e8 -> LE scatter noise ~1.5 %)
N_HE=${N_HE:-4e6}         # photons per HE run (5 runs: 2e7 -> HE scatter noise < 1 %)
# Measured on the workstation (8 threads, 22 Oct 2026): LE 36 000 photons/s, HE 12 000 photons/s, plus ~1 min set-up
# per run -> defaults about 50 min LE + 35 min HE per case. Phase space ~3 MB (LE) / ~19 MB (HE) per 1e6 photons.

CASE=$(python -c "import config; print(config.case_name())")
SET_NAME=$(python -c "import config; print(config.MEASURED['set'])")
PAIR_NAME=$(python -c "import config; print(config.MEASURED['pair'])")
echo "case $CASE (phantom $REVELATION_PHANTOM), measured set $SET_NAME, pair $PAIR_NAME"
ROOT_ARG=(); [ -n "$ROOT" ] && ROOT_ARG=(--root "$ROOT")

python -c "import config, sys; sys.exit(0 if config.MEASURED['set_ready'] else 1)" || {
  echo "measured set $SET_NAME: exposure values not read yet; reading them from the DICOM headers"
  python measured_set.py --set "$SET_NAME" ${ROOT_ARG[@]+"${ROOT_ARG[@]}"}; }
python spectra.py                       # X-ray spectra (SpekPy), extra Al fitted to the pair's LE entrance dose
python selftest.py                      # grid, path lengths, normalisation, packings
python projector.py --beams LE HE       # noise-free primary images (phantom, open field)

for beam in LE HE; do
  n=$([ "$beam" = LE ] && echo "$N_LE" || echo "$N_HE")
  for s in $(seq 1 "$RUNS"); do
    tag="${beam}_phantom_${CASE}_s${s}"
    python revelation_gate.py --beam "$beam" --n "$n" --threads "$THREADS" --seed "$s" --tag "$tag"
    python postprocess.py "output/$tag" --name "${beam}_${CASE}" --accumulate --delete-phsp
  done
done

# images, w*, Siemens-weight residual, texture, iodine (noise-free and with quantum noise of the measured exposure)
python analysis.py --scatter-le "LE_${CASE}" --scatter-he "HE_${CASE}"
python analysis.py --scatter-le "LE_${CASE}" --scatter-he "HE_${CASE}" --noise

# measurement-like stored-value twin of the measured pair (its exact LE and HE mAs), and the comparison with it
python measlike.py --scatter-le "LE_${CASE}" --scatter-he "HE_${CASE}"
python compare_measured.py ${ROOT_ARG[@]+"${ROOT_ARG[@]}"} || echo "comparison skipped: measured images of set $SET_NAME not reachable (see ROOT)"
# the investigation's own July analysis (july2026/pipe.py, agg.py, discagg.py), unchanged, on the measured pair and on
# its simulated twin, plus the rim profiles; report in output/pipeline/<set>_<case>/report/<pair>/ (July set only)
if [ "$SET_NAME" = "july28high" ]; then
  python pipeline_compare.py all --workers 1 || echo "pipeline comparison failed (see above)"
fi
# another background in the real L2: its measurement predicted from the measured current L2 and the two twins (needs
# the pipeline results of L2_PMMA_paraffin for the same pair); report in output/redesign/<set>_<pair>_PMMA_paraffin_vs_<background>/
if [ "$SET_NAME" = "july28high" ] && [ "$REVELATION_PHANTOM" = "L2" ] && [ "$CASE" != "L2_PMMA_paraffin" ]; then
  python redesign_compare.py --redesign "${CASE#L2_}" || echo "redesign prediction skipped (see above)"
  if [ "$CASE" = "L2_gland_adipose" ]; then   # the tissue reference adds its column to the PA12_MCT prediction
    python redesign_compare.py --redesign PA12_MCT || echo "PA12_MCT prediction not updated (see above)"
  fi
fi

# Optional checks (see README, validation):
#  - projector vs Monte Carlo primary: keep one phase space (omit --delete-phsp) and compare
#    output/mc/mcprimary_<name>.npy with output/projector/<case>/primary_<beam>.npy on 1 mm bins;
#  - other cases: PHANTOM=L2v2 BACKGROUND=PMMA_paraffin ./run_all.sh, ... (outputs are kept apart per case:
#    output/projector/<case>/, output/mc/<beam>_<case>, output/analysis/<case>, output/measlike/<case>/<set>/);
#  - scatter from the open field (paddle/support) that --target phantom leaves out:
#    python revelation_gate.py --beam LE --n 2e8 --target field --tag LE_field_s1
#  - background-only and iodine-only variants: --no-spheres / --no-discs (and projector.py with the same flags)
