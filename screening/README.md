# screening — paketin ici

Bu klasor bir Python paketidir; **dogrudan calistirilmaz**.

Tek giris noktasi bir ust klasordeki:

```
screening.bat
```

Kullanim kilavuzu (Turkce, ayrintili):

```
KAPSAMLI_ARAMA_NASIL_CALISTIRILIR.md
```

Okuma motoru hatasinin anlatimi kilavuzun **8. bolumundedir**.
Bu dosyanin onceki, daha uzun hali: `eski/README.md`.

## Paketin parcalari

| dosya | is |
|---|---|
| `__main__.py` | ana akis, butun modlar |
| `motor.py` | projenin mevcut olcum kodunu ice aktarir (yeniden yazmaz) |
| `okuma_motoru.py` | DUZELTILMIS ham okuma motoru (kayipsiz tohumlama) |
| `kaba_kuvvet.py` | tohumsuz bagimsiz referans - dogrulama icin |
| `test_motor.py` / `bagimsiz_dogrulama.py` | motor esitlik testleri |
| `geometri.py` | toplanti geometri kurallari (parametrelestirilmis) |
| `hedefler.py` | panel / uyelik / kutu okuma |
| `hedef_uyelik.tsv` | **uye-rakip tanimi - elle duzenlenebilir** |
| `uretec.py` | pencere, cift, ARMS, 144 hucreli izgara |
| `numune.py` | ham okuma in-silico PCR + Wilson |
| `referans.py` / `kuresel_tarama.py` | referans kapsam / kuresel ozgulluk |
| `panel_olcum.py` | secenek (4) panel yeniden olcum |
| `uyelik_denetimi.py` | secenek (5) uyelik duyarlilik analizi |
| `konsensus_uret.py` | secenek (6) konsensus yeniden uretim |
| `hepsi.py` | secenek (9) her seyi sirayla kos |
| `rapor.py` / `kontrol.py` | rapor uretimi / checkpoint |
| `panel_guncelle.py` | duzeltmeyi panel xlsx'ine isler - ELLE calistirilir, menude yok |
| `ciftler.tsv` | diger oturumun uye tanimi; secenek (5) bunu karsilastirmaya alir |
| `eski/` | yerini yenisi alan dosyalar (silinmedi) |
