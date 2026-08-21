# -*- coding: utf-8 -*-
"""BU DOSYA ARTIK KULLANILMIYOR.  Tam hali: eski/measure_panel.py

Yaptigi is verification/full_chain.py menusundeki SECENEK (4)'e tasindi:
    panel_measurement.py  ->  python3 -m screening --mod panel-olc --tam-derinlik

Neden degisti (ayrinti: eski/NEDEN_BURADALAR.md):
  * ciftler artik elle tutulan ciftler.tsv'den degil, PANEL TSV'sinden okunuyor
  * mm<=1 ve mm<=3 tek kosuda olculuyor, olcut her satirda yaziyor
  * kesintiye dayanikli (checkpoint), uyelik kaynagi izlenebilir
  * bu dosyanin en degerli ozelligi - eski/yeni motoru kutu kutu yan yana
    olcmesi - secenek (4)'e AYNEN tasindi

Tek giris noktasi: klasor kokundeki verification/full_chain.py
"""
# ---------------------------------------------------------------------------
# measure_panel.py — DEVRE DISI. Yalnizca eski cagri yollarini kibarca durdurmak
#                icin duruyor; hicbir olcum yapmaz.
#
# GIRDI  : yok. Tam hali eski/measure_panel.py icindedir.
# CIKTI  : kendi docstring'ini ekrana basar ve sys.exit ile hata mesaji vererek
#          sonlanir.
# CAGRAN : hicbir menu tusu bu dosyayi calistirmaz. Yaptigi is
#          verification/full_chain.py tusu 4'e (panel_measurement.py, --mod panel-olc
#          --tam-derinlik) tasindi.
#
# Dosya silinmek yerine birakildi: eski bir kisayol ya da not hala bu yolu
# gosteriyor olabilir ve sessizce yanlis bir sey calistirmaktansa acikca
# durup dogru tusa yonlendirmek dogrudur.
# ---------------------------------------------------------------------------
import sys
print(__doc__)
sys.exit('Bu betik devre disi. verification/full_chain.py -> secenek (4) kullanin.')
