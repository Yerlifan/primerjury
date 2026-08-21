#!/usr/bin/env bash
# =====================================================================
# AGIR_ISLER.sh
# Bulut kutusunda güvenle çalışmayan işleri sizin WSL'inizde sırayla
# çalıştırır. Her adım kendi başına da çalıştırılabilir; sırayla
# gitmesinin tek sebebi çıktıların birbirini beslemesidir.
#
#   bash AGIR_ISLER.sh                # A, B, C, D, E adımları
#   bash AGIR_ISLER.sh --kraken       # F adımını da ekler (dakikalar)
#   SECILEN_ESIK=0.02 bash AGIR_ISLER.sh --kraken   # eşiği doğrudan uygula
#   bash AGIR_ISLER.sh --yalniz C     # yalnız tek adım
#
# Adımlar:
#   A  kutu kurtarma (20)        B-1_2233851 konsensüsünü okumalardan kurar
#   B  hedef kimliği (22)        her hedefin adı ile dizinin gösterdiğini
#                                yan yana koyar
#   C  mfeprimer (19)            dış özgüllüğün ikinci bağımsız ölçümü
#   D  topluluk trendi (21)      bolluk kitabını cins düzeyi esaslı kurar
#   E  öz denetim (17, 18)       regresyon ve teslim denetimi
#   F  Kraken2 güven eşiği        çıktı dosyalarından, veritabanı gerekmez
#                                (önce eşik taraması, sonra seçilen eşik)
#   G  Excel teslimatı (13)      ölçülen kimlik sütunlarıyla, B'den sonra
#   H  geniş dış VT taraması     SILVA, UNITE, ROD, PR2; SAATLER sürer,
#                                yalnız "--yalniz H" ile çalışır
#
# Her adım kendi log'unu $PT/agir_log altına yazar; hiçbiri diğerinin
# çıktısını silmez. Bir adım hata verirse betik DURMAZ, hatayı işaretler
# ve devam eder; sonda hangi adımların başarısız olduğu listelenir.
# =====================================================================
set -uo pipefail

PT="${PT:-/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama}"
HERE="$(cd "$(dirname "$0")" && pwd)"
KONS="${KONS:-$PT/referans_konsensus/baskin/konsensus}"
ADAY="${ADAY:-$PT/primer_adaylari}"
FINAL="${FINAL:-$PT/primer_final}"
LOGD="$PT/agir_log"
MFE="${MFE:-$PT/tools/mfeprimer}"
K2DB="${K2DB:-$HOME/k2db}"
HAM_FASTQ="${HAM_FASTQ:-}"
CONF="${CONF:-0.1}"
IS="${IS:-4}"

KRAKEN=0; YALNIZ=""
while [ $# -gt 0 ]; do
  case "$1" in
    --kraken) KRAKEN=1; shift;;
    --yalniz) YALNIZ="$2"; shift 2;;
    --pt) PT="$2"; shift 2;;
    *) echo "bilinmeyen secenek: $1" >&2; exit 2;;
  esac
done

mkdir -p "$LOGD"
ANA="$LOGD/AGIR_ISLER.log"
ts(){ date '+%Y-%m-%d %H:%M:%S'; }
say(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$ANA"; }
BASARISIZ=()

calistir() {   # $1 = adim harfi, $2 = aciklama, gerisi komut
  local h="$1" ac="$2"; shift 2
  if [ -n "$YALNIZ" ] && [ "$YALNIZ" != "$h" ]; then return 0; fi
  say "----------------------------------------------------------------"
  say "ADIM $h  $ac"
  local t0=$(date +%s)
  "$@" 2>&1 | tee -a "$LOGD/adim_$h.log" | tail -n 40
  local rc=${PIPESTATUS[0]}
  say "ADIM $h bitti, cikis=$rc, sure=$(( ($(date +%s)-t0)/60 )) dakika"
  if [ "$rc" -ne 0 ]; then
    BASARISIZ+=("$h ($ac), cikis=$rc")
    say "  DIKKAT: bu adim hata verdi, tam log: $LOGD/adim_$h.log"
  fi
}

say "================================================================"
say "AGIR ISLER BASLANGIC"
say "PT     : $PT"
say "KONS   : $KONS"
say "loglar : $LOGD"

# --- bagimlilik denetimi ---------------------------------------------
say "bagimlilik denetimi"
EKSIK=()
for k in blastn python3; do command -v "$k" >/dev/null 2>&1 || EKSIK+=("$k"); done
[ -x "$MFE" ] || { chmod +x "$MFE" 2>/dev/null || true; }
[ -x "$MFE" ] || EKSIK+=("mfeprimer (calistirilabilir degil: $MFE)")
python3 - <<'PYX' 2>&1 | tee -a "$ANA"
for m in ("primer3", "Bio", "mappy", "openpyxl"):
    try:
        __import__(m); print("   %-10s TAMAM" % m)
    except ImportError:
        print("   %-10s EKSIK  -> pip install --break-system-packages %s"
              % (m, {"Bio": "biopython", "primer3": "primer3-py"}.get(m, m)))
PYX
if [ ${#EKSIK[@]} -gt 0 ]; then
  say "EKSIK: ${EKSIK[*]}"
  say "  blastn icin: sudo apt-get install -y ncbi-blast+"
fi

# --- A. kutu kurtarma -------------------------------------------------
# B-1_2233851 self referansi %19 IUPAC oldugu icin minimap2 o diziden
# yalnizca 2 minimizer cikarabiliyor ve konsensus sifir uzunlukta cikiyor.
# Okumalar yerinde (5914 tane); kalip okumalarin kendisinden kurulur.
KURT_FQ="$PT/fastq files/B-1/B-1-reads_2233851.fastq"
if [ -f "$KURT_FQ" ] && [ ! -s "$KONS/B-1_2233851_baskin_konsensus.fasta" ]; then
  calistir A "kutu kurtarma (20): B-1_2233851" \
    python3 "$HERE/recover_bins.py" --fastq "$KURT_FQ" \
      --etiket B-1_2233851 --out "$KONS"
else
  say "ADIM A atlandi (fastq yok ya da konsensus zaten var)"
fi

# --- B. hedef kimligi -------------------------------------------------
calistir B "hedef kimligi (22): ad ile dizinin gosterdigi karsilastirmasi" \
  python3 "$HERE/target_identity.py" --kons "$KONS" --db "$PT/REFERANS_DB" \
    --hedefler "$HERE/hedefler.tsv" --adlar "$HERE/taxid_adlari.tsv" \
    --is-parcacigi "$IS" --out "$FINAL/hedef_kimlik.tsv"

# --- C. mfeprimer -----------------------------------------------------
# Bulut kutusunda bacteria.16S adiminda bellek yetmedi. Burada 16 GB var,
# ama yine de veritabani basina zaman asimi konuyor ve asilan veritabani
# "olculemedi" olarak isaretleniyor, sessizce temiz sayilmiyor.
calistir C "mfeprimer (19): dis ozgullugun ikinci olcumu" \
  python3 "$HERE/mfeprimer_layer.py" --final "$FINAL" --db "$PT/REFERANS_DB" \
    --mfe "$MFE" --cpu "$IS" --zaman-asimi 7200 \
    --out "$FINAL/mfeprimer.tsv"

# --- D. topluluk trendi -----------------------------------------------
KIMLIK=()
for f in "$PT"/t_kimlik/kimlik_*.tsv "$HERE"/kimlik_*.tsv \
         "$HERE"/kimlik/kimlik_*.tsv; do
  [ -f "$f" ] && KIMLIK+=("$f")
done
if [ ${#KIMLIK[@]} -eq 0 ]; then
  say "UYARI: kimlik_*.tsv bulunamadi. Tur duzeyi guvenilirlik isareti"
  say "  yalnizca ayirt_edilemez.tsv'den turetilecek, kutu kimligi olcumu"
  say "  devre disi kalacak. Beklenen konum: $HERE/kimlik/"
fi
# Once rutbe kapsamasi olculur: bolluk hangi rutbede okunabilir. Bracken
# CALISTIRILMAZ; gerekcesi 25'in basinda ve bolluk_rutbe_kaniti.md'de.
RUTBEARG=""
if [ -d "$PT/kraken_c${SECILEN_ESIK:-0.02}" ]; then
  calistir D0 "rutbe kapsamasi (25): bolluk hangi rutbede okunabilir" \
    python3 "$HERE/abundance_rank.py" \
      --kraken "$PT/kraken_c${SECILEN_ESIK:-0.02}" \
      --out "$PT/bolluk_rutbe"
  [ -s "$PT/bolluk_rutbe/ozet.tsv" ] && RUTBEARG="--rutbe $PT/bolluk_rutbe"
else
  say "  NOT: kraken_c${SECILEN_ESIK:-0.02} yok, rutbe kapsamasi atlandi."
  say "       Once: SECILEN_ESIK=0.02 bash AGIR_ISLER.sh --yalniz F"
fi
calistir D "topluluk trendi (21): rutbe farkindali bolluk kitabi" \
  python3 "$HERE/community_trends.py" --bracken "$PT/bracken results" \
    --ayirt "$ADAY/ayirt_edilemez.tsv" \
    ${KIMLIK[@]+--kimlik "${KIMLIK[@]}"} \
    --adlar "$HERE/taxid_adlari.tsv" $RUTBEARG \
    --out "$PT/Microbooster_Topluluk_Trend_Analizi.xlsx"

# --- E. oz denetim ----------------------------------------------------
calistir E "oz denetim (17): regresyon takimi" \
  python3 "$HERE/regression_test.py" --gercek-veri --aday "$ADAY" --kons "$KONS"
# 18 KRITIK bulgu bulunca cikis kodu 1 doner. Bu bir COKME DEGIL, bir
# BULGUDUR: ham primer_final.tsv'de alan karisimi satirlari duruyor ve
# Excel bunlari zaten disariya aliyor. Ikisini ayirmazsak her kosuda
# "basarisiz adim" yazar ve gercek cokmeler gozden kacar.
say "----------------------------------------------------------------"
say "ADIM E2  oz denetim (18): teslim denetimi"
T0E2=$(date +%s)
python3 "$HERE/check_deliverables.py" --final "$FINAL" --kons "$KONS" \
  --hedefler "$HERE/hedefler.tsv" --out "$FINAL/teslim_denetimi.tsv" \
  2>&1 | tee -a "$LOGD/adim_E2.log" | tail -n 40
RCE2=${PIPESTATUS[0]}
say "ADIM E2 bitti, cikis=$RCE2, sure=$(( ($(date +%s)-T0E2)/60 )) dakika"
if [ "$RCE2" = 1 ]; then
  say "  E2: KRITIK bulgu var (cokme degil). Bulgular:"
  awk -F'\t' 'NR>1 && $1=="KRITIK"{say[$4]++} END{for(k in say) printf "     %-24s %d\n", k, say[k]}' \
    "$FINAL/teslim_denetimi.tsv" 2>/dev/null | tee -a "$ANA"
  say "  Bu bulgular calisma kitabindan zaten cikariliyor; ham tabloda kaliyor."
elif [ "$RCE2" -ne 0 ]; then
  BASARISIZ+=("E2 (teslim denetimi cokmesi), cikis=$RCE2")
  say "  DIKKAT: bu adim COKTU, tam log: $LOGD/adim_E2.log"
fi

# --- F. Kraken2 guven esigi ------------------------------------------
# ONEMLI: yeniden siniflandirmaya GEREK YOK. Kraken2'nin --output dosyalari
# her okumanin k-mer LCA dizisini zaten tasiyor; guven puani tam olarak
# o diziden hesaplanir. 106 GB veritabani da, ham barkod fastq dosyalari
# da gerekmez. Olculdu: vurus taksonlarinin %99,84'u rapor dosyalarindan
# kurulan agacta yer aliyor.
#
# Esik EZBERDEN SECILMEZ. Kisa okuma verisi icin sik onerilen 0,1 degeri
# bu ONT verisinde okumalarin %69'unu siniflandirilmamis birakiyor, cunku
# ONT okumalarinin k-mer'lerinin cogu veritabaninda karsilik bulmuyor ve
# bu k-mer'ler puanin paydasina giriyor. Once tarama calisir, esik
# tablodan secilir.
if [ "$KRAKEN" = 1 ] || [ "$YALNIZ" = "F" ]; then
  if [ ! -d "$PT/kraken results" ]; then
    say "ADIM F atlandi: 'kraken results' klasoru yok"
  else
    calistir F "Kraken2 guven esigi taramasi (24)" \
      python3 "$HERE/reassign_confidence.py" --kraken "$PT/kraken results" \
        --tarama 0,0.002,0.005,0.01,0.02,0.05,0.1 --tarama-okuma 20000 \
        --out "$PT/kraken_guven"
    say "  Tarama tablosu: $PT/kraken_guven/esik_taramasi.tsv"
    say "  Esigi secip su komutu calistirin:"
    say "    python3 $HERE/reassign_confidence.py \\"
    say "        --kraken '$PT/kraken results' --confidence <ESIK> \\"
    say "        --out '$PT/kraken_c<ESIK>'"
    if [ -n "${SECILEN_ESIK:-}" ]; then
      calistir F2 "Kraken2 guven esigi uygulamasi (24), esik=$SECILEN_ESIK" \
        python3 "$HERE/reassign_confidence.py" --kraken "$PT/kraken results" \
          --confidence "$SECILEN_ESIK" --out "$PT/kraken_c$SECILEN_ESIK"
    fi
  fi
fi

# --- G. Excel teslimatini yenile --------------------------------------
# B adiminda uretilen hedef_kimlik.tsv, calisma kitabina "olculen kimlik"
# sutunu olarak girer. Bu yuzden Excel B'den SONRA yeniden uretilir.
if [ -z "$YALNIZ" ] || [ "$YALNIZ" = "G" ]; then
  REFC="$PT/primer_referans"
  REFARG=""
  [ -s "$REFC/primer_referans.tsv" ] && REFARG="--referans $REFC/primer_referans.tsv"
  calistir G "Excel teslimati (13), olculen kimlik sutunlariyla" \
    python3 "$HERE/export_excel.py" \
      --aday "$ADAY" --final "$FINAL" --bol "$ADAY/kume_setleri" \
      --adlar "$HERE/taxid_adlari.tsv" --hedefler "$HERE/hedefler.tsv" \
      --kons "$KONS" --kimlik "$FINAL/hedef_kimlik.tsv" $REFARG \
      --out "$PT/PrimerJury_Primer_Tasarimi.xlsx"
fi

# --- H. genis dis veritabani taramasi ---------------------------------
# Dar kume (NCBI RefSeq 16S ve ITS) yalnizca tip susu ve temsilci dizileri
# icerir; kulturlenmemis cevresel soylar orada YOKTUR. Bu numunede
# bakteri hedeflerinin ucunun RefSeq'te %90 uzeri akrabasi bile
# bulunamadi, yani topluluğun buyuk kismi dar kumede temsil edilmiyor.
# SILVA SSU NR99, UNITE, ROD ve PR2 cevresel dizileri de tasir.
# UZUN SURER (saatler); bu yuzden ayri adim ve elle baslatilir.
if [ "$YALNIZ" = "H" ]; then
  calistir H "genis dis veritabani taramasi (14 --genis)" \
    python3 "$HERE/external_databases.py" --final "$FINAL" \
      --db "$PT/REFERANS_DB" --genis --is-parcacigi "$IS" --kons "$KONS" \
      --hedefler "$HERE/hedefler.tsv" --adlar "$HERE/taxid_adlari.tsv" \
      --kimlik "$FINAL/hedef_kimlik.tsv" \
      --zaman-asimi 21600 --out "$FINAL/dis_veritabani_genis.tsv"
  calistir H2 "genis kume, ikinci olcum (19 --genis)" \
    python3 "$HERE/mfeprimer_layer.py" --final "$FINAL" \
      --db "$PT/REFERANS_DB" --mfe "$MFE" --genis --cpu "$IS" \
      --zaman-asimi 21600 --blast "$FINAL/dis_veritabani_genis.tsv" \
      --out "$FINAL/mfeprimer_genis.tsv"
fi

# --- ozet -------------------------------------------------------------
say "----------------------------------------------------------------"
say "OZET"
for f in "$FINAL/hedef_kimlik.tsv" "$FINAL/mfeprimer.tsv" \
         "$FINAL/teslim_denetimi.tsv" \
         "$PT/Microbooster_Topluluk_Trend_Analizi.xlsx" \
         "$KONS/B-1_2233851_baskin_konsensus.fasta"; do
  if [ -s "$f" ]; then
    say "  VAR    $(basename "$f")  ($(stat -c%s "$f") bayt)"
  else
    say "  YOK    $(basename "$f")"
  fi
done
if [ ${#BASARISIZ[@]} -gt 0 ]; then
  say "BASARISIZ ADIMLAR:"
  for x in "${BASARISIZ[@]}"; do say "   $x"; done
  say "Tam loglar: $LOGD"
  exit 1
fi
say "butun adimlar tamamlandi"
say "Loglar: $LOGD"
