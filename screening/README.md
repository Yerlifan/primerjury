# screening — paketin ici

Bu klasor bir Python paketidir; **dogrudan calistirilmaz**.

Tek giris noktasi bir ust klasordeki:

```
verification/full_chain.py
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
| `engine_gateway.py` | projenin mevcut olcum kodunu ice aktarir (yeniden yazmaz) |
| `read_engine.py` | DUZELTILMIS ham okuma motoru (kayipsiz tohumlama) |
| `brute_force.py` | tohumsuz bagimsiz referans - dogrulama icin |
| `engine_test.py` / `independent_check.py` | motor esitlik testleri |
| `geometry.py` | toplanti geometri kurallari (parametrelestirilmis) |
| `targets.py` | panel / uyelik / kutu okuma |
| `hedef_uyelik.tsv` | **uye-rakip tanimi - elle duzenlenebilir** |
| `generator.py` | pencere, cift, ARMS, 144 hucreli izgara |
| `sample.py` | ham okuma in-silico PCR + Wilson |
| `reference.py` / `global_scan.py` | referans kapsam / kuresel ozgulluk |
| `panel_measurement.py` | secenek (4) panel yeniden olcum |
| `membership_check.py` | secenek (5) uyelik duyarlilik analizi |
| `build_consensus.py` | secenek (6) konsensus yeniden uretim |
| `run_all.py` | secenek (9) her seyi sirayla kos |
| `report.py` / `checks.py` | rapor uretimi / checkpoint |
| `update_panel.py` | duzeltmeyi panel xlsx'ine isler - ELLE calistirilir, menude yok |
| `ciftler.tsv` | diger oturumun uye tanimi; secenek (5) bunu karsilastirmaya alir |
| `eski/` | yerini yenisi alan dosyalar (silinmedi) |
