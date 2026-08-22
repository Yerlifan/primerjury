# -*- coding: utf-8 -*-
"""THIS FILE IS NO LONGER IN USE.  The full version: eski/measure_panel.py

What it did was moved to OPTION (4) in the verification/full_chain.py menu:
    panel_measurement.py  ->  python3 -m screening --mod panel-olc --tam-derinlik

Why it changed (the detail: eski/NEDEN_BURADALAR.md):
  * the pairs are read from the PANEL TSV now, not from a hand kept ciftler.tsv
  * mm<=1 and mm<=3 are measured in a single run and the criterion is written on
    every row
  * it survives interruption (checkpoints) and the membership source can be traced
  * this file's most valuable feature, measuring the old and the new engine bin by
    bin side by side, was moved to option (4) UNCHANGED

The one entry point: verification/full_chain.py in the directory root

"""
# -------------------------------------------------------------------------
# measure_panel.py, DISABLED. It stays only to stop old call routes politely; it
#                makes no measurement.
#
# INPUT  : none. The full version is in eski/measure_panel.py.
# OUTPUT : it prints its own docstring to the screen and ends with sys.exit and an
#          error message.
# CALLED BY: no menu key runs this file. What it did was moved to
#          verification/full_chain.py key 4 (panel_measurement.py, --mod panel-olc
#          --tam-derinlik).
#
# The file was left rather than deleted: an old shortcut or note may still point at
# this path, and stopping plainly and pointing at the right key is better than
# silently running the wrong thing.
# -------------------------------------------------------------------------
import sys
print(__doc__)
sys.exit(u'This script is disabled. Use verification/full_chain.py -> option (4).')
