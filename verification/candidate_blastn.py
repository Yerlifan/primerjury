# -*- coding: utf-8 -*-
"""
candidate_blastn.py - kimlik asamasi icin blastn tabanli ADAY BULUCU.

Kisa liste sorununun tam cozumu icin en guclu aday. Bu dosya PUANLAMAYA
DOKUNMAZ: kimlik yuzdesi yine kimlik_dogrulama.hizala ile, bizim tanimimizla
hesaplanir. Degisen tek sey, hangi kayitlarin hizalanacagina KIMIN karar
verdigidir.
"""
# ---------------------------------------------------------------------------
# candidate_blastn.py
#
# GIRDI  : kutu konsensusu (sorgu) ve REFERANS_DB altindaki BLAST indeksli
#          veritabani (.nin/.nhr/.nsq dosyalari zaten mevcut)
# CIKTI  : [(baslik, dizi), ...] - blastn'in anlamli isabet buldugu kayitlar
# CAGRAN : verification/identity_verification.py ve all_bin_identities.py,
#          yalniz ADAY_BULUCU=blastn verildiginde
#
# NEDEN BU SECENEK EN GUCLU - VE NEDEN DAHA ONCE DUSUNULMEDI
# Bu projede kendi tohum + hizalama hattimizi yazdik. Ama blastn tam olarak bu
# is icin tasarlanmis, otuz yildir sinaniyor ve siralamayi BIT SKORU ile
# yapiyor: eslesme uzunlugu, kimlik ve bosluk cezasi birlikte hesaba giriyor,
# ayrica e-degeri veritabani boyutuna gore duzeltiliyor. Bizim tohum sayimiz
# ise kimligin KOTU BIR VEKILIDIR - uzun ve korunmus bir kayit, alakasiz oldugu
# halde cok tohum toplar.
#
# Daha da onemlisi: blastn BU KLASORDE ZATEN KULLANILIYOR. Kuresel ozgulluk
# katmani onunla kosuldu ve REFERANS_DB altinda dokuz veritabaninin BLAST
# indeksi (.nin) hazir duruyor. Yani hazir, indeksli ve sinanmis bir araci
# yanimizda tutarken kendi hattimizi yazmisiz. Bu, projenin kendi kuralina
# aykiridir: once eldeki olculmus araca bakilir.
#
# KISA LISTE SORUNU NASIL COZULUR
# Bizim yolda adaylar tohum sayisina gore siralanip ilk N alinir; olculdu ki
# kazanan bir sorguda 4171. siradan geldi. blastn'de "ilk N" diye bir on eleme
# yoktur; e-deger esigini gecen HER kayit doner. Kesme noktasi, keyfi bir sira
# sayisi degil ISTATISTIKSEL BIR ESIKTIR ve raporlanabilir.
#
# NEREDE KULLANILMAZ
# Primer baglanma aramasi ve in-siliko PCR. blastn de kisa sorgularda tohum
# ayarina baglidir ve kayipsizlik GARANTISI vermez; oralarda guvercin yuvasi
# motoru kalir.
#
# GUVENLIK: blastn yoksa var_mi() False doner ve cagiran taraf mevcut kisa
# liste yoluna devam eder. Zincir kirilmaz.
# ---------------------------------------------------------------------------

from __future__ import print_function
import os
import subprocess
import sys
import tempfile

_DENENDI = False
_VAR = False
_SEBEP = u''
_SURUM = u''


def var_mi():
    """blastn calisir durumda mi. Bir kez olcer."""
    global _DENENDI, _VAR, _SEBEP, _SURUM
    if _DENENDI:
        return _VAR
    _DENENDI = True
    try:
        p = subprocess.Popen(['blastn', '-version'], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        c = p.communicate(timeout=30)[0].decode('utf-8', 'replace')
        if p.returncode != 0:
            _SEBEP = u'blastn -version sifir disi kod dondurdu'
            return False
        _SURUM = c.strip().splitlines()[0] if c.strip() else u'bilinmiyor'
        _VAR = True
        return True
    except OSError:
        _SEBEP = (u'blastn PATH uzerinde yok. Kurulum: '
                  u'micromamba install -n mikro -c bioconda blast')
        return False
    except Exception as e:
        _SEBEP = u'blastn calistirilamadi: %s: %s' % (type(e).__name__, e)
        return False


def sebep():
    return _SEBEP


def surum():
    var_mi()
    return _SURUM or u'yok'


def indeks_var_mi(fasta_yolu):
    """Veritabaninin BLAST indeksi kurulu mu.

    Indeks yoksa makeblastdb ile kurulabilir ama bu DAKIKALAR surer ve
    kendiliginden yapilmaz: kullanicinin haberi olmadan gigabaytlik indeks
    uretmek dogru olmaz. Bunun yerine False donulur ve cagiran taraf eski
    yola duser.
    """
    return os.path.exists(fasta_yolu + '.nin') or os.path.exists(fasta_yolu + '.nal')


def adaylar(q, fasta_yolu, e_deger=1e-5, en_fazla=5000, iplik=3):
    """blastn ile aday kayitlari bulur.

    Doner: [(baslik, dizi), ...] ya da None (blastn/indeks yoksa).
    None ile bos liste AYRI SEYLERDIR:
      None -> arac ya da indeks yok, ESKI YOLA DUS
      []   -> blastn kostu ve hicbir kayit esigi gecmedi
    Bu ayrim karistirilirsa bir veritabani sessizce atlanir.
    """
    if not var_mi() or not indeks_var_mi(fasta_yolu):
        return None
    fd, sorgu = tempfile.mkstemp(suffix='.fa', prefix='blq_')
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write('>sorgu\n%s\n' % q)
        komut = ['blastn', '-query', sorgu, '-db', fasta_yolu,
                 '-outfmt', '6 sseqid bitscore evalue',
                 '-evalue', str(e_deger),
                 '-max_target_seqs', str(en_fazla),
                 '-num_threads', str(iplik),
                 # dc-megablast: uzak akrabalari da bulur. megablast yakin
                 # eslesmeler icin hizlidir ama %90 altini kacirir; bizim
                 # kutularin bir kismi referanslara %85 civarinda benziyor.
                 '-task', 'dc-megablast']
        p = subprocess.Popen(komut, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cikti, hata = p.communicate()
        if p.returncode != 0:
            sys.stderr.write(u'blastn ERROR (%s): %s\n'
                             % (p.returncode, hata.decode('utf-8', 'replace')[:300]))
            return None
        idler = []
        gorulen = set()
        for satir in cikti.decode('utf-8', 'replace').splitlines():
            a = satir.split('\t')
            if a and a[0] and a[0] not in gorulen:
                gorulen.add(a[0])
                idler.append(a[0])
        if not idler:
            return []
        return _dizileri_cek(fasta_yolu, idler)
    finally:
        try:
            os.unlink(sorgu)
        except OSError:
            pass


def _dizileri_cek(fasta_yolu, idler):
    """blastdbcmd ile dizileri toplu ceker. Tek tek cekmek cok yavas olurdu."""
    fd, liste = tempfile.mkstemp(suffix='.txt', prefix='blid_')
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write('\n'.join(idler) + '\n')
        p = subprocess.Popen(
            ['blastdbcmd', '-db', fasta_yolu, '-entry_batch', liste,
             '-outfmt', '%f'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cikti, hata = p.communicate()
        if p.returncode != 0:
            sys.stderr.write(u'blastdbcmd ERROR: %s\n'
                             % hata.decode('utf-8', 'replace')[:200])
            return None
        out = []
        bas, parca = None, []
        for satir in cikti.decode('utf-8', 'replace').splitlines():
            if satir.startswith('>'):
                if bas is not None:
                    out.append((bas, ''.join(parca).upper()))
                bas, parca = satir[1:].strip(), []
            else:
                parca.append(satir.strip())
        if bas is not None:
            out.append((bas, ''.join(parca).upper()))
        return out
    finally:
        try:
            os.unlink(liste)
        except OSError:
            pass


def secili_mi():
    """ADAY_BULUCU=blastn verilmis ve blastn calisiyor mu."""
    return (os.environ.get('ADAY_BULUCU', '').strip().lower() == 'blastn'
            and var_mi())


if __name__ == '__main__':
    print(u'blastn calisiyor mu : %s' % (u'EVET' if var_mi() else u'HAYIR'))
    if not var_mi():
        print(u'sebep               : %s' % _SEBEP)
        sys.exit(1)
    print(u'version             : %s' % surum())
    kok = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rd = os.path.join(kok, 'REFERANS_DB')
    if os.path.isdir(rd):
        for f in sorted(os.listdir(rd)):
            if f.endswith(('.fasta', '.fna')):
                y = os.path.join(rd, f)
                print(u'  %-40s index: %s'
                      % (f, u'VAR' if indeks_var_mi(y) else u'yok'))
