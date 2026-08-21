"""Havuz uzerinde hizli cift (F+R) urun taramasi; ispcr ile birebir olcut."""
# ---------------------------------------------------------------------------
# cift.py — binlerce aday primer ciftini ham okumalara karsi hizli taramak icin
#           kullanilan "havuz" yolu; olcut ispcr.amplify ile ayni, hiz farkli.
#
# GİRDİ  : tarayici.Havuz nesnesi (butun okumalar tek buyuk dizide birlestirilmis
#          ve SEED=10 icin indekslenmis olarak) + F ve R primer dizileri.
# ÇIKTI  : dosyaya YAZMAZ. urunler() -> (m, boy). m = okuma sayisi uzunlugunda
#          bool maske (hangi okumada urun var), boy = urun boyu -> okuma sayisi.
# ÇAĞRAN : KAPSAMLI_ARAMA/motor.py bu dosyayi ada gore bulup yukler (zorunlu
#          degil; yoksa cift=None kalir) ve kendini_sina.py yuklendigini
#          dogrular. Ayni klasordeki son_sec.py, nd_v2/v2b/v3/v4.py, esle.py,
#          arke_esle.py, arke_skor.py, clo_ara.py, clo_ara2.py,
#          dogrula_cift.py "from cift import urunler" der. Menude motor.py
#          uzerinden olcum tuslarinda yuklenir (P, K, D, T, U, H, 1-9).
#
# NEDEN HAVUZ: bir aday ciftini tek tek her okumada aramak, aday sayisi
# binlere cikinca kabul edilemez. Havuz butun okumalari tek bir buyuk dizide
# birlestirir, aralarina 'N' dolgusu koyar ve bir kez indeksler; sonra her
# aday tek numpy taramasiyla olculur. Olcut degismez, yalnizca ayni hesabin
# duzeni degisir.
#
# NEDEN HAM OKUMALAR, REFERANS DIZI DEGIL: referans/konsensus tek bir
# ozetlenmis dizidir ve numunedeki gercek varyantlari, tur ici farklari,
# nanopore hata desenini icermez. Panelin karar sayilari bu yuzden numunenin
# ham okumalari uzerinde uretilir.
# ---------------------------------------------------------------------------
import sys,pickle,numpy as np
# Oturuma bagli mutlak yollar - bu betikler bir WSL calisma dizininde yazildi.
# motor.py ice aktarmadan once bu yollari gecerli hale getirir.
sys.path.insert(0,'/tmp/wk2/SON_ETAP_betikleri'); sys.path.insert(0,'/tmp/wk2/is')
import ispcr

# ---------------------------------------------------------------------------
# urunler — havuzdaki hangi okumalarda F+R ciftinin urun verdigini bulur.
#
# OLCUT (ispcr.amplify ile birebir):
#   - F ileri yonde baglanir (hv.bul; 3' son 2 baz TAM, uyumsuzluk <= mm)
#   - R'nin kendisi degil TERS TUMLEYENI aranir (hv.bul5). Motor yalniz arti
#     ipligi taradigi icin, rc(R)'nin arti iplikte bulunmasi R'nin gercekte
#     karsi iplige bagliandigi anlamina gelir. bul5, R'nin 3' son 2 bazina
#     karsilik gelen rc(R)'nin ILK 2 bazinin tam tutmasini sart kosar.
#   - Urun boyu penceresi lo..hi. Varsayilan (70, 250) qPCR kisitindan gelir:
#     yapilandirmada URUN_IDEAL 60-150, URUN_KABUL 150-250 bp'dir.
#   - seg >= f + len(F): R'nin baglanma yeri F'nin SAGINDA ve F ile ORTUSMEDEN
#     olmali; iki primer birbirine bakmali. Aksi halde ayni bolgeye binmis iki
#     baglanma sahte urun uretirdi.
# ---------------------------------------------------------------------------
def urunler(hv,F,R,lo=70,hi=250,mm=1):
    fs=hv.bul(F,mm)
    # F hicbir yerde baglanmiyorsa urun de olamaz - bos maske ile erken cikis.
    if fs.size==0: return np.zeros(hv.n,bool),{}
    rs=hv.bul5(ispcr.rc(R),mm)
    if rs.size==0: return np.zeros(hv.n,bool),{}
    # sid: havuzdaki her global konumun HANGI okumaya ait oldugu. Okumalar
    # arasindaki 'N' dolgu bolgelerinde sid = -1'dir.
    fid=hv.sid[fs]; rid=hv.sid[rs]
    # sid < 0 olan vuruslar dolgu bolgesine dusmustur; atilir. Dolgu hem iki
    # okumanin ucunun yanlislikla birlesip sahte bir motif olusturmasini, hem
    # de bir urunun iki ayri okumaya yayilmasini engeller.
    okf=fid>=0; fs=fs[okf]; fid=fid[okf]
    okr=rid>=0; rs=rs[okr]; rid=rid[okr]
    m=np.zeros(hv.n,bool); boy={}
    # R vuruslari okuma kimligine gore siralanip her okumanin araligi
    # searchsorted ile onceden cikarilir. Boylece asagidaki dongude bir F
    # vurusu icin YALNIZ AYNI OKUMADAKI R vuruslarina bakilir - urunun tek bir
    # okuma icinde kalmasi bu sekilde garanti edilir. Ayrica her F icin butun
    # R listesini taramak yerine O(1) dilim alinir; hiz buradan gelir.
    order=np.argsort(rid,kind='stable'); rs=rs[order]; rid=rid[order]
    starts=np.searchsorted(rid,np.arange(hv.n),'left')
    ends=np.searchsorted(rid,np.arange(hv.n),'right')
    for i,f in zip(fid,fs):
        # Okuma zaten urun vermis olarak isaretliyse tekrar bakilmaz: bu bir
        # VARLIK olcumudur (kac okumada urun var), derinlik olcumu degil.
        if m[i]: continue
        seg=rs[starts[i]:ends[i]]
        if seg.size==0: continue
        # Urun boyu: F'nin basindan rc(R)'nin sonuna; primerler urune dahildir.
        sz=seg+len(R)-f
        sel=(sz>=lo)&(sz<=hi)&(seg>=f+len(F))
        if sel.any():
            m[i]=True
            # Boy dagilimina okuma basina TEK deger yazilir (ilk uyan urun).
            v=int(sz[sel][0]); boy[v]=boy.get(v,0)+1
    return m,boy
