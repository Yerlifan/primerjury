# ---------------------------------------------------------------------------
# geometry_core.py — paneldeki 21 hedefin 42 primerinin qPCR GEOMETRI denetimi; her primer
#          icin Tm/GC/hairpin/dimer hesaplanip kural ihlalleri listelenir.
#
# GİRDİ  : dosya okumaz. Denetlenecek 21 hedef x (F, R) listesi P sabitinde
#          dosyanin kendi icine gomuludur. Disaridan yalniz primer3 gerekir.
# ÇIKTI  : calisilan dizine geo.json (her primer icin olculen ham degerler ve
#          ihlal listesi) yazar; ayrica ekrana satir satir GECTI / IHLAL tablosu
#          basar.
# ÇAĞRAN : screening/self_test.py bu dosyayi AYRI BIR SUREC olarak
#          gecici bir dizinde calistirir, urettigi geo.json'u okur ve
#          screening/geometry.py'nin ayni 42 primer icin verdigi
#          gc/tm/hp/hd/son5/uc degerleriyle BIREBIR karsilastirir. Tek bir
#          deger bile tutmazsa kendini sinama duser ve kosu baslamaz. Bu
#          sinama her olcum tusunun basinda kosar: full_chain.py asamalari
#          1, 2, 3, 4, 5, 6, 7, 9 ve dogrudan tus 8. Dosya menuden DOGRUDAN
#          cagrilmaz; referans/altin standart olarak durur.
#
# ÖLÇÜLDÜ - HANGI KOPYA CALISIYOR: self_test.py, yapilandirma.BETIK_YOLLARI
# sirasini gezip ILK buldugu geometry_core.py'yi alir. Sira engine,
# engine, engine'dir ve engine'nde geometry_core.py
# YOKTUR; dolayisiyla fiilen engine/geometry_core.py calistirilir. Bu dosya
# ile o dosya BAYT BAYT AYNIDIR (ayni md5), bu yuzden sonuc degismez - ama
# ikisinden biri degistirilirse sinama sessizce digerini olcmeye devam eder.
# ---------------------------------------------------------------------------
import primer3, json, sys
# Panelin 21 hedefi ve her hedefin (ileri, geri) primer cifti. Denetim bu sabit
# liste uzerinde kosar; boylece sonuc her calistirmada aynidir ve
# screening/geometry.py ile karsilastirilabilir bir altin standart olur.
P = [
 ("Metanojen_universal","GTGGAGCTTGCGGTTTAATTG","CAGGATGCTTCACAGTACGAAC"),
 ("Methanothrix_cinsi","GAGAGGTACTTCAGGGGTAGG","CTAGCTTTCGTCCCTTGCC"),
 ("Petrimonas_cinsi","AAGTCGCGTGAAGGATGAAG","AAAATTTCACCGCCGACTTAAC"),
 ("Metanomikrobiyales_hidrojenotrof","TGGGACCGCCTCTGCTAAAG","CATTGTAGCCCGCGTGTAGC"),
 ("Mantar_universal (F1)","GGTTACCCGCTGAACTTAAGC","CGCTTCACTCGCCGTTAC"),
 ("Methanosarcina_cinsi","TCGCTAGGTGTCAGGCATG","GCGATTCAGGCAAGGTCTTC"),
 ("Metilotrofik_metanojen","CAATCCTGAAACCCGTCCATAG","ATATTCACCGCCTGATGTTGAC"),
 ("Nitrosocosmicus_AOA","ACTCTGAGTGATTTCCGTTAAGG","TGCTTTAGGCCCAATAAACGTC"),
 ("Proteiniphilum_cinsi","GGTTCCTTGAGTGTGGATGAGG","CTTGAGCGTCAGTTATGGCTTAG"),
 ("Proteolitik_Cloacibacillus","AGCTAGTAGGTTGGGTAACGG","GATTTCTTCACCCACGCGG"),
 ("Bacteroidales_kumesi","AGCTAGGATTTGGTTGCTGTGG","CCTCAGCGTCAATATTGGCTTTG"),
 ("Mantar_universal (F2)","GTGCATGGCCGTTCTTAGTTG","CAAACTTCCATCGGCTTGAGC"),
 ("Microascaceae_askomikot","ATCAATAAGCGGAGGAAAAGAAACC","CCTCTTCAAATTACAACTCGGACTG"),
 ("Sakarolitik_Sphaerochaeta","ATCTGGCCATGTACTGACGC","CTGGTGCACATCGTTTACTGTG"),
 ("Asetoklastik_metanojenler","CCGGGAGAGGTGAGAGGTAC","CGGGTATCTAATCCGGTTCGTG"),
 ("Bakteri_universal","ACAAGCGGTGGAGCATGTG","ACGACAGCCATGCAGCAC"),
 ("Petriella_musispora","GGAGTCGTCCTAATATGCGAGTG","CAAATCCATCCGAGAACATCAGG"),
 ("Proteolitik_Cloacimonas","TTAAAGGCAGCGGCTCACC","GAACCCGACACCTAGTGATTATCG"),
 ("Arke_universal","CTGCGGTTTAATTGGATTCAACGC","GAACTGACGACGGCCATGC"),
 ("Methanothrix_soehngenii_turu","AATGTAGCAATACATGGCGAACTG","TTCCAGCAATCGAGACCTATCG"),
 ("Methanosarcina_mazei_turu","GCCCTTGGGACCGGCATAA","TCGCTGGCTAGTAGGTACATTACA"),
]
# Termodinamik kosullar. Tm bir dizinin degil, dizi + TAMPON'un ozelligidir:
# tuz ve dNTP derisimi degisince Tm de degisir. Bu degerler kullanilan qPCR
# karisimini (QuantiNova SYBR Green) temsil eder; config.py'de ayni
# degerler P3 sabitinde "geometry_core.py ile birebir" notuyla tekrarlanmistir. Bu
# sayilar gercek reaksiyondan saparsa asagidaki 58-62 Tm penceresi anlamsizlasir.
#   mv_conc   = 50   mM tek degerlikli katyon (Na+/K+)
#   dv_conc   = 1.5  mM iki degerlikli katyon (Mg++)
#   dntp_conc = 0.6  mM dNTP
#   dna_conc  = 50   nM oligo derisimi
KW=dict(mv_conc=50, dv_conc=1.5, dntp_conc=0.6, dna_conc=50)
# Erime sicakligi (nearest-neighbor, primer3).
def tm(p): return round(primer3.calc_tm(p, **KW),2)
# GC yuzdesi.
def gc(p): return round(100.0*sum(p.count(c) for c in 'GC')/len(p),1)
# Hairpin (sac tokasi): primerin KENDI uzerine katlanmasinin erime sicakligi.
def hp(p): return round(primer3.calc_hairpin(p, **KW).tm,1)
# Homodimer: primerin KENDI kopyasiyla eslesmesinin erime sicakligi.
def hd(p): return round(primer3.calc_homodimer(p, **KW).tm,1)
# Heterodimer: F ile R'nin BIRBIRIYLE eslesmesinin erime sicakligi.
def het(a,b): return round(primer3.calc_heterodimer(a,b, **KW).tm,1)
# ---------------------------------------------------------------------------
# viol — tek bir primerin ihlal ettigi kurallarin listesi. Bos liste = GECTI.
# Her kuralin neden var oldugu ve qPCR'deki karsiligi asagida satir satir.
# ---------------------------------------------------------------------------
def viol(p):
    v=[]
    # UZUNLUK 18-25 baz. 18'in altinda primer genomda/toplulukta rastgele de
    # bulunabilecek kadar kisalir, ozgulluk duser; 25'in ustunde Tm asiri
    # yukselir, katlanma ve yanlis baglanma sansi artar, sentez maliyeti buyur.
    if not (18<=len(p)<=25): v.append("uz %d"%len(p))
    g=gc(p)
    # GC %40-60. G/C ciftinin uc hidrojen bagi vardir, A/T'nin iki. Cok dusuk GC
    # -> dupleks 60 C'de tutunamaz, verim duser; cok yuksek GC -> primer hedef
    # disi GC zengini bolgelere de yapisir ve sablon 95 C'de tam acilmayabilir.
    if not (40<=g<=60): v.append("GC %%%.1f"%g)
    t=tm(p)
    # Tm PENCERESI 58-62 C. Panelin tamami TEK bir protokolde, ayni cihaz
    # programiyla (60 C annealing hedefi) kosacaktir. Bir primerin Tm'i
    # pencerenin altina duserse o kuyuda baglanma zayiflar (yalanci negatif),
    # ustune cikarsa 60 C'de hedef disi bolgelere de baglanir (yalanci pozitif).
    # Pencere dar tutulur ki 21 hedefin hepsi ayni kosuda calisabilsin.
    if not (58<=t<=62): v.append("Tm %.2f"%t)
    # HAIRPIN esigi 45 C. Primer kendi uzerine katlanip sac tokasi yaparsa, o
    # yapi annealing sicakliginda (60 C civari) hala ayakta oldugunda primer
    # sablona baglanmak yerine kendisiyle mesgul olur -> verim duser, Cq gecikir.
    # Esik annealing sicakliginin belirgin altina cekilmistir; 45 C'nin ustunde
    # eriyen bir yapi guvenli sayilmaz. SYBR Green'de bu olcut ELEYICIDIR.
    if hp(p)>=45: v.append("hairpin Tm %.1f"%hp(p))
    # HOMODIMER esigi 45 C. Primerin kendi kopyasiyla eslesmesi primer-dimer
    # uretir. Primer-dimer hem primeri ve dNTP'yi tuketir hem de SYBR Green
    # cift ipligin TAMAMINI boyadigi icin hedef yokken bile SINYAL VERIR;
    # erime egrisinde dusuk Tm'li ikinci bir tepe olarak gorunur. Prob temelli
    # bir kimyada uyari sayilabilecek bu olcut SYBR'de eleyicidir.
    if hd(p)>=45: v.append("homodimer Tm %.1f"%hd(p))
    # 3' SON BAZ G ya da C olmali ("GC clamp"). Polimeraz primeri 3' ucundan
    # uzatir; o uctaki uc hidrojen bagli G/C, primerin uzatmanin basladigi yerde
    # saglam oturmasini saglar. 3' ucu A/T olan primer uctan "solur" ve verim
    # degisken olur.
    if p[-1] not in 'GC': v.append("3' uc %s (G/C degil)"%p[-1])
    n5=sum(1 for c in p[-5:] if c in 'GC')
    # 3' SON 5 BAZDA EN COK 3 G/C. Bir onceki kuralin karsi agirligi. 3' ucu
    # G/C ile fazla doldurmak ("GC clamp"i abartmak) primerin YALNIZCA son
    # birkac bazi tutan hedef disi bolgelerde bile saglam yapismasina yol acar;
    # uc oturdugu icin polimeraz uzatmaya baslar ve ozgul olmayan urun cikar.
    # Yani: ucta yeterli tutus olsun (bir onceki kural), ama fazlasi olmasin.
    if n5>3: v.append("3' son 5 bazda %d G/C"%n5)
    # DEJENERE BAZ yasak. ACGT disi bir harf (N, R, Y...) primerin tek bir oligo
    # degil oligo KARISIMI olarak sentezlenmesi demektir: her bir varyantin
    # etkin derisimi duser, Tm tek bir sayi olmaktan cikar ve yukaridaki
    # termodinamik hesaplarin hepsi anlamini yitirir.
    if any(c not in 'ACGT' for c in p): v.append("dejenere baz")
    return v

# ---------------------------------------------------------------------------
# BURADAN ASAGISI GOSTERI BOLUMUDUR - yalniz betik DOGRUDAN calistirildiginda
# kosar. 2026-08-11 duzeltmesi: bu bolum eskiden ICE AKTARIMDA da kosuyordu.
# Sonucu:
#   1) geometry_core.py'yi import eden her betik, ekrana yukaridaki P listesinin (2
#      Agustos panelinin) 63 satirlik geometri dokumunu basiyordu. O listede
#      Petrimonas ve Bacteroidales gibi ciftlerin ESKI dizileri var; konsolu
#      ya da gunlugu okuyan biri onlari GUNCEL saniyordu.
#   2) import, calisilan dizine geo.json YAZIYORDU. Ice aktarmanin dosya
#      yazmasi beklenmeyen bir yan etkidir.
# self_test.py bu dosyayi AYRI BIR SUREC olarak (cwd=gecici dizin)
# calistirip geo.json okuyor; orada __name__ == "__main__" oldugu icin
# davranis DEGISMEZ. Yalniz "import geo" sessizlesti.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rows=[]
    for name,f,r in P:
        for lbl,p in (("Ileri",f),("Geri",r)):
            rows.append(dict(hedef=name,primer=lbl,dizi=p,uz=len(p),gc=gc(p),tm=tm(p),
                             hp=hp(p),hd=hd(p),uc=p[-1],gc5=sum(1 for c in p[-5:] if c in 'GC'),
                             ihlal=viol(p)))
        # ---- CIFT DUZEYI KURALLAR (tek primere bakarak gorunmezler) ----
        d=abs(tm(f)-tm(r)); h=het(f,r); pv=[]
        # dTm < 1.5 C. F ve R AYNI kuyuda, AYNI annealing sicakliginda calisir.
        # Tm'leri ayrisirsa dusuk Tm'li primer o sicaklikta zayif baglanir; bir
        # iplik digerinden hizli cogalir (asimetrik amplifikasyon), verim ve Cq
        # tekrarlanabilirligi bozulur.
        if d>=1.5: pv.append("dTm %.2f"%d)
        # HETERODIMER esigi 45 C. F ve R birbirine yapisirsa olusan primer-dimer,
        # homodimerle ayni zarari verir ama daha siktir: iki farkli dizi bulusur.
        # SYBR Green'de dogrudan sahte sinyal kaynagidir.
        if h>=45: pv.append("heterodimer Tm %.1f"%h)
        rows.append(dict(hedef=name,primer="CIFT",dizi="dTm=%.2f  het=%.1f"%(d,h),uz="",gc="",tm="",hp="",hd="",uc="",gc5="",ihlal=pv))
    # geo.json CALISILAN DIZINE yazilir (mutlak yol yok). self_test.py bu dosyayi
    # gecici bir dizinde cwd=td vererek calistirir ve ciktiyi oradan okur; bu yuzden
    # goreli yol kasitlidir, proje agacini kirletmez.
    json.dump(rows,open('geo.json','w'))
    # Ekran ciktisi. CIFT satirinda hem cift duzeyi ihlaller hem de o hedefin iki
    # primerinin ihlalleri TEK satirda toplanir: karar hedef bazinda verilir, bir
    # hedefin tek primeri bile duserse o hedef siparise giremez.
    for r in rows:
        if r['primer']=='CIFT':
            pr = [x['ihlal'] for x in rows if x['hedef']==r['hedef'] and x['primer']!='CIFT']
            allv = r['ihlal'] + [i for s in pr for i in s]
            print("%-34s CIFT  %-24s -> %s" % (r['hedef'], r['dizi'], "IHLAL: "+"; ".join(allv) if allv else "GECTI"))
        else:
            print("%-34s %-5s %-26s uz%3s GC%5s Tm%6s hp%5s hd%6s 3'%s g5=%s  %s" % (
                r['hedef'],r['primer'],r['dizi'],r['uz'],r['gc'],r['tm'],r['hp'],r['hd'],r['uc'],r['gc5'],
                ("IHLAL: "+"; ".join(r['ihlal'])) if r['ihlal'] else "GECTI"))
