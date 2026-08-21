#!/usr/bin/env bash
# =====================================================================
# CALISTIR.sh
# Toplantı kararlarının tamamını uçtan uca çalıştırır ve her adımı
# tarih saatle log'lar. Yarıda kesilirse aynı komutla kaldığı yerden
# devam eder, çünkü 08 ve 09 checkpoint tutuyor.
#
#   bash CALISTIR.sh
#
# Ortam değişkeniyle yol değiştirilebilir:
#   PT=/baska/yol bash CALISTIR.sh
#
# Adımlar:
#   0  bağımlılık denetimi
#   1  baskın alel konsensüsü (12)   ham okumalardan, belirsizliksiz
#   2  toplu tasarım (08)            hedef başına primer adayı
#   3  özgüllük ve ham okuma doğrulaması (09)
#   4  dış veritabanı taraması (14)
#   5  referans tabanlı tasarım (15)
#   6  Excel teslimatı (13)
#   7  öz denetim: regresyon testi (17) + teslim denetimi (18)
# =====================================================================
set -uo pipefail

PT="${PT:-/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama}"
KONS_IUPAC="${KONS_IUPAC:-$PT/referans_konsensus/self/konsensus}"
BASKIN="${BASKIN:-$PT/referans_konsensus/baskin}"
KONS="${KONS:-$BASKIN/konsensus}"
ADAY="${ADAY:-$PT/primer_adaylari}"
FINAL="${FINAL:-$PT/primer_final}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ANA_LOG="$PT/CALISTIR.log"
MAX_OKUMA="${MAX_OKUMA:-20000}"
TOP="${TOP:-10}"

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
say(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$ANA_LOG"; }

say "================================================================"
say "BASLANGIC"
say "PT          : $PT"
say "IUPAC kons  : $KONS_IUPAC"
say "baskin kons : $KONS"
say "ADAY        : $ADAY"
say "FINAL       : $FINAL"

if [ ! -d "$KONS_IUPAC" ]; then
  say "HATA: konsensus klasoru yok: $KONS_IUPAC"
  say "Once 06 ve 07 --mode self calistirilmali."
  exit 1
fi

# --- 0. bagimlilik denetimi -------------------------------------------
say "----------------------------------------------------------------"
say "ADIM 0/7  bagimlilik denetimi"
python3 - <<'PYX' 2>&1 | tee -a "$ANA_LOG"
eksik = []
for m in ("primer3", "Bio", "mappy", "openpyxl"):
    try:
        __import__(m)
        print("   %-10s TAMAM" % m)
    except ImportError:
        eksik.append(m)
        print("   %-10s EKSIK" % m)
if eksik:
    kur = {"Bio": "biopython", "primer3": "primer3-py"}
    paket = " ".join(kur.get(m, m) for m in eksik)
    # Debian ve Ubuntu tabanli dagitimlarda sistem Python'u PEP 668 ile
    # korumali; duz 'pip install' externally-managed-environment hatasi verir.
    print("   Kurmak icin: pip install --break-system-packages " + paket)
    print("   Alternatif (sistem Python'una dokunmadan):")
    print("     python3 -m venv ~/primer_venv && source ~/primer_venv/bin/activate")
    print("     pip install " + paket + " primer3-py biopython openpyxl mappy")
    if "mappy" in eksik:
        print("   mappy olmadan ayirt edilemez kutu olcumu ve capraz bulasma")
        print("   olcumu YAPILAMAZ; kosu durmaz ama ozgulluk guvencesi zayiflar.")
PYX

# --- 1. baskin alel konsensusu ---------------------------------------
# samtools consensus -A ciktisi degisken pozisyonlari IUPAC koduyla yaziyor.
# Baglanma kurali kodu temsil ettigi butun bazlara uyar saydigi icin, rakip
# konsensusu belirsizlik acisindan zenginse hangi primeri secerseniz secin
# rakipte de bagliyor gorunuyor. Bu adim ayni ham okumalardan ikinci ve
# belirsizliksiz bir konsensus uretir.
say "----------------------------------------------------------------"
say "ADIM 1/7  baskin alel konsensusu (12)"
T0=$(date +%s)
MEVCUT=$(ls "$KONS"/*_konsensus.fasta 2>/dev/null | wc -l)
if [ "$MEVCUT" -gt 0 ] && [ "${YENIDEN_KONS:-0}" != "1" ]; then
  say "  zaten var ($MEVCUT dosya), atlandi. Yeniden uretmek icin YENIDEN_KONS=1"
else
  python3 "$HERE/12_baskin_alel_konsensus.py" \
      --kons "$KONS_IUPAC" \
      --fastq "$PT/fastq files" \
      --out "$BASKIN" 2>&1 | tee -a "$ANA_LOG"
fi
say "ADIM 1 bitti, sure=$(( ($(date +%s)-T0)/60 )) dakika"
N=$(ls "$KONS"/*_konsensus.fasta 2>/dev/null | wc -l)
say "baskin alel konsensus dosyasi: $N"
[ "$N" -gt 0 ] || { say "HATA: baskin konsensus uretilemedi"; exit 1; }

# --- 2. toplu tasarim -------------------------------------------------
say "----------------------------------------------------------------"
say "ADIM 2/7  toplu tasarim (08)"
T0=$(date +%s)
python3 "$HERE/08_toplu_tasarim.py" \
    --kons "$KONS" \
    --hedefler "$HERE/hedefler.tsv" \
    --out "$ADAY" 2>&1 | tee -a "$ANA_LOG"
RC=${PIPESTATUS[0]}
say "ADIM 2 bitti, cikis=$RC, sure=$(( ($(date +%s)-T0)/60 )) dakika"
[ "$RC" -eq 0 ] || { say "ADIM 2 hatali bitti, duruldu"; exit "$RC"; }

A=$(ls "$ADAY"/*__*.tsv 2>/dev/null | wc -l)
say "uretilen aday dosyasi: $A"

# --- 3. ozgulluk ve ham okuma dogrulamasi -----------------------------
say "----------------------------------------------------------------"
say "ADIM 3/7  ozgulluk ve ham okuma dogrulamasi (09)"
T0=$(date +%s)
python3 "$HERE/09_ozgulluk.py" \
    --adaylar "$ADAY" \
    --pt "$PT" \
    --kons "$KONS" \
    --hedefler "$HERE/hedefler.tsv" \
    --out "$FINAL" \
    --top "$TOP" --max-okuma "$MAX_OKUMA" 2>&1 | tee -a "$ANA_LOG"
RC=${PIPESTATUS[0]}
say "ADIM 3 bitti, cikis=$RC, sure=$(( ($(date +%s)-T0)/60 )) dakika"

# --- 4. dis veritabani taramasi ---------------------------------------
say "----------------------------------------------------------------"
say "ADIM 4/7  dis veritabani taramasi (14)"
T0=$(date +%s)
if command -v blastn >/dev/null 2>&1; then
  python3 "$HERE/14_dis_veritabani.py" \
      --final "$FINAL" --db "$PT/REFERANS_DB" \
      --out "$FINAL/dis_veritabani.tsv" 2>&1 | tee -a "$ANA_LOG"
else
  say "  blastn bulunamadi, dis veritabani taramasi atlandi."
  say "  Kurmak icin: sudo apt-get install -y ncbi-blast+"
fi
say "ADIM 4 bitti, sure=$(( ($(date +%s)-T0)/60 )) dakika"

# --- 5. referans tabanli tasarim (kapsanamayan hedefler) --------------
# Numunedeki kutular ayrilamadigi icin cift bulunamayan hedefler icin
# referans veritabanindan tasarim yapilir. Bu ciftler numuneyle
# DOGRULANMAZ; yalnizca numunede destek olup olmadigi olculur ve cikti
# ayri bir sayfada, ayri etiketle sunulur.
say "----------------------------------------------------------------"
say "ADIM 5/7  referans tabanli tasarim (15)"
T0=$(date +%s)
REFC="$PT/primer_referans"
if [ -f "$HERE/hedefler_referans.tsv" ]; then
  python3 "$HERE/15_referans_tasarim.py" \
      --db "$PT/REFERANS_DB" --pt "$PT" \
      --hedefler-ref "$HERE/hedefler_referans.tsv" \
      --out "$REFC" --max-okuma "$MAX_OKUMA" 2>&1 | tee -a "$ANA_LOG"
else
  say "  hedefler_referans.tsv yok, adim atlandi"
fi
say "ADIM 5 bitti, sure=$(( ($(date +%s)-T0)/60 )) dakika"

# --- 6. Excel teslimati -----------------------------------------------
say "----------------------------------------------------------------"
say "ADIM 6/7  Excel teslimati (13)"
REFARG=""
[ -s "$REFC/primer_referans.tsv" ] && REFARG="--referans $REFC/primer_referans.tsv"
python3 "$HERE/13_teslim_excel.py" \
    --aday "$ADAY" --final "$FINAL" \
    --bol "$ADAY/kume_setleri" \
    --adlar "$HERE/taxid_adlari.tsv" \
    --hedefler "$HERE/hedefler.tsv" --kons "$KONS" $REFARG \
    --out "$PT/MicRhoBooster_Primer_Tasarimi.xlsx" 2>&1 | tee -a "$ANA_LOG"

# --- 7. oz denetim ----------------------------------------------------
# Iki ayri denetim calisir:
#   17  kod kurallarini bagimsiz referans uygulamalarla karsilastirir
#   18  teslim edilen tabloyu sifirdan yeniden olcer (tasarim kodunu
#       ice aktarmaz), yazili Tm/GC/urun boyu degerlerini de dogrular
# 18 KRITIK bulgu bulursa cikis kodu 1 doner ve burada acikca raporlanir.
say "----------------------------------------------------------------"
say "ADIM 7/7  oz denetim (17 regresyon + 18 teslim denetimi)"
T0=$(date +%s)
python3 "$HERE/17_regresyon_testi.py" \
    --gercek-veri --aday "$ADAY" --kons "$KONS" 2>&1 | tee -a "$ANA_LOG"
RC17=${PIPESTATUS[0]}
if [ -s "$FINAL/primer_final.tsv" ]; then
  python3 "$HERE/18_teslim_denetimi.py" \
      --final "$FINAL" --kons "$KONS" --hedefler "$HERE/hedefler.tsv" \
      --out "$FINAL/teslim_denetimi.tsv" 2>&1 | tee -a "$ANA_LOG"
  RC18=${PIPESTATUS[0]}
else
  RC18=2
  say "  primer_final.tsv yok, teslim denetimi calistirilamadi"
fi
say "ADIM 7 bitti, 17=$RC17, 18=$RC18, sure=$(( ($(date +%s)-T0)/60 )) dakika"
[ "$RC17" -eq 0 ] || say "DIKKAT: regresyon testinde kalan test var, log'a bakin"
[ "$RC18" -eq 0 ] || say "DIKKAT: teslim denetiminde KRITIK bulgu var, log'a bakin"

# --- ozet -------------------------------------------------------------
say "----------------------------------------------------------------"
if [ -s "$FINAL/primer_final.tsv" ]; then
  G=$(awk -F'\t' 'NR>1 && $4=="GECTI"' "$FINAL/primer_final.tsv" | wc -l)
  H=$(awk -F'\t' 'NR>1 && $4=="GECTI" {print $2}' "$FINAL/primer_final.tsv" | sort -u | wc -l)
  say "SONUC: butun kurallardan gecen aday = $G, kapsanan hedef = $H"
  say "Dosya: $FINAL/primer_final.tsv"
  say "Excel: $PT/MicRhoBooster_Primer_Tasarimi.xlsx"
else
  say "UYARI: primer_final.tsv olusmadi"
fi
if [ -s "$ADAY/ayirt_edilemez.tsv" ]; then
  AE=$(( $(wc -l < "$ADAY/ayirt_edilemez.tsv") - 1 ))
  say "AYIRT EDILEMEZ TAKSON CIFTI: $AE  (dosya: $ADAY/ayirt_edilemez.tsv)"
  say "  Bu ciftler dizi duzeyinde ayrilmiyor; rakip listesinden cikarildilar."
fi
say "Loglar:"
say "  $ANA_LOG"
say "  $ADAY/toplu_tasarim.log"
say "  $FINAL/ozgulluk.log"
say "Checkpointler (yarida kesilirse ayni komut kaldigi yerden devam eder):"
say "  $ADAY/checkpoint.json"
say "  $FINAL/checkpoint.json"
say "BITTI"
