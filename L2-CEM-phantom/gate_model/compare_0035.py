#!/usr/bin/env python3
"""
Kept for the earlier commands: the comparison with the measured pair 0035 is now done by compare_measured.py, which
also compares with the July 2026 set. This is the same as

    python compare_measured.py --set 0035 [--phantom L2v2] [--background PMMA_paraffin] [--bin 1]

(output in output/compare/0035_<case>/). Note that 0035 shows the old L2 container with the old plate; the matching
simulation for the real L2 with the new plate is  python compare_measured.py --phantom L2  (July set).
"""
import sys
from compare_measured import main

if __name__ == "__main__":
    if "--set" not in sys.argv:
        sys.argv += ["--set", "0035"]
    if "--background" not in sys.argv:
        sys.argv += ["--background", "PMMA_paraffin"]       # the measured phantom's materials (as before)
    main()
