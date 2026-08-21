# MicRhoBooster Primer Tasarımı: Ayrıntılı Çalışma Kaydı

Son güncelleme: 2026-08-01
Proje klasörü: `C:\Users\yerli\Masaüstü\PrimerTasarlama`
Betikler: `steps/`

Bu belge, yapılan bütün işleri, ölçülen bütün sayıları, bulunan bütün hataları
ve verilen bütün kararları kaydeder. Amaç, hiçbir sonucun kaynağının
belirsiz kalmaması ve her iddianın hangi ölçümden geldiğinin izlenebilmesidir.

---

## 1. Görev ve kabul ölçütleri

Anaerobik çürütücü (MicRhoBooster) örneklerinden alınan ONT nanopor rDNA
amplikon verisi Kraken2 ile sınıflandırıldı. Görev, toplantı kararlarında
belirlenen taksonlar için PCR primer çiftleri tasarlamak.

### 1.1 Oligo kuralları

- Yalnız A/C/G/T; dejenere baz yok
- Uzunluk 18-25
- GC %40-60 (sert sınır %35-65)
- 3' uç G ya da C olmalı
- Son beş bazda en fazla üç G/C
- En fazla dört ardışık aynı baz

### 1.2 Termodinamik

- Tm 58-62 °C (sert sınır 57-63)
- Saç tokası ΔG ≥ -3000
- Kendi kendine dimer ≥ -6000
- Hetero-dimer ≥ -6000
- Tm **iki bağımsız kütüphaneyle** ölçülür (primer3 ve Biopython), sistematik
  fark veriden ölçülür

### 1.3 Ürün

- 70-250 bp (en fazla 300)
- 90-150 en iyi puanı alır
- Çift içi Tm farkı 1,5 °C altında
- Ürünün başı ileri primere, sonu geri primerin ters tümleyenine **birebir**
  eşit olmalı; makineyle denetlenir

### 1.4 Bağlanma kuralı (PCR fiziği)

- Son iki baz birebir eşleşmeli
- Son beş bazda en fazla bir uyumsuzluk
- Toplamda en fazla üç uyumsuzluk
- 5' çıkıntı serbest
- Primerler karşıt zincirlerde ve 3' uçları birbirine bakacak

### 1.5 Özgüllük

- Ürün hedeflenen **her** üyede oluşmalı
- Rakiplerin **hiçbirinde** oluşmamalı
- En az bir primer rakiplerde hiç bağlanma yeri bulmamalı (yetim primer)
- Çiftler ham okumalar üzerinde yeniden sınanır
- Rakip oranları **Wilson alt sınırıyla** değerlendirilir

### 1.6 Genel kural

Hiçbir karar tek bir kod yoluna bırakılmaz. Her adım iki bağımsız ölçüm alır;
ayrışma olursa reddedilir ve kaydedilir. Primer BLAST zorunludur.

### 1.7 Toplantı hedefleri ve istenen düzey

`hedefler.tsv` içinde kodlu. Yirmi hedef, dört karar başlığı.

**Karar 1, tür özgül (6 hedef):** Methanosarcina mazei, Methanothrix
soehngenii, Methanosarcina barkeri, Podospora pseudopauciseta, Dictyostelium
discoideum, Trichoderma asperellum

**Karar 2, cins özgül (4 hedef):** Bacteroides, Alistipes, Proteiniphilum,
Petrimonas

**Karar 3, işlevsel grup (7 hedef):** Hidrojenotrofik metanojenler,
Metilotrofik metanojen, Asetoklastik metanojenler, Sakarolitik bakteriler,
Proteolitik sintrofik bakteriler, Nitrosocosmicus AOA, Trichoderma cinsi

**Karar 4, evrensel (3 hedef):** Bakteri, Arke, Mantar

**Tür özgüllüğünde hoşgörü:** 1-2 çapraz **tür** kabul edilir. Ölçü çapraz tür
sayısıdır, o türlerde oluşan ürün sayısı değil.

### 1.8 Çalışma kısıtları

- Ağır işler kullanıcının WSL'inde çalışır, bulut kutusunda değil
- Her şey tarih ve saatle loglanır
- Kesinti için kontrol noktaları tutulur
- Her şey çapraz kontrol edilir
- Önbellekten sağlama toplamı doğrulanmadan okuma yapılmaz
- Veriye bağlı değerler koda gömülmez
- Teslimden önce her zaman öz denetim

---

## 2. Boru hattı

### 2.1 Ana zincir

| Betik | İş |
|---|---|
| `check_environment.sh` | bağımlılık denetimi |
| `reclassify_kraken2.sh` | Kraken2 sınıflandırma |
| `analyze_ambiguous_bases.sh` | N oranı analizi |
| `generate_primer_candidates.py` | oligo adayları, kompozisyon ve termodinamik |
| `design_group_primers.py` | çift oluşturma, bağlanma ve ürün kuralı motoru |
| `split_clusters.py` | kutu bölme |
| `anchored_reference_consensus.sh` | referans çapalı konsensüs |
| `freeze_reference.sh` | referans sabitleme |
| `batch_design.py` | hedef başına toplu tasarım |
| `specificity.py` | ham okumalarda özgüllük, Wilson alt sınırı |
| `indistinguishable_targets.py` | ayırt edilemez kutuların bulunması |
| `check_bin_identity.py` | kutu kimlik denetimi |
| `dominant_allele_consensus.py` | baskın alel konsensüsü |
| `export_excel.py` | Excel teslimatı |

### 2.2 Doğrulama ve dış ölçüm

| Betik | İş |
|---|---|
| `external_databases.py` | dış veritabanı taraması (blastn), kapsam denetimi, takson ayrımı |
| `design_from_reference.py` | referans veritabanından tasarım |
| `check_primer_geometry.py` | geometri denetimi |
| `regression_test.py` | regresyon takımı, **137 test** |
| `check_deliverables.py` | bağımsız teslim denetimi |
| `mfeprimer_layer.py` | ikinci bağımsız özgüllük ölçümü |
| `recover_bins.py` | referanssız konsensüs kurtarma |
| `community_trends.py` | topluluk trend çalışma kitabı |
| `target_identity.py` | hedef adı ile verinin karşılaştırılması |
| `reassign_confidence.py` | Kraken2 güven eşiği, çevrimdışı |
| `abundance_rank.py` | bolluğun okunabildiği rütbe |
| `reference_identity.py` | referans primerlerin numunede ne çoğalttığı |
| `check_taxonomic_level.py` | tür/cins özgüllüğünün doğrudan sınanması |
| `alan_denetimi.py` | alan (arke/bakteri/mantar) tutarlılığı, ortak modül |

### 2.3 Sürücüler

| Betik | İş |
|---|---|
| `CALISTIR.sh` | ana zincir |
| `AGIR_ISLER.sh` | ağır adımlar (A-H), WSL'de |
| `SENKRON.sh` | SHA256 manifesto, iki taraflı doğrulama |

---

## 3. Bulunan hatalar ve düzeltmeleri

Bu bölüm çalışmanın en önemli kısmıdır. Her madde ölçülmüş bir hatayı,
sayısını ve düzeltmesini taşır. Bunların çoğu benim hatamdı ve öyle
işaretlenmiştir.

### 3.1 Ürün boyu, geri primerin uzunluğunu unutuyordu

Ürün boyu iki 3' ucu arasındaki mesafe olarak ölçülüyordu. Doğrusu ileri
primerin 5' ucundan geri primerin 5' ucuna kadarki mesafedir. Fark tipik
olarak 44 baz; gerçek 70-94 bp'lik hedef dışı ürünler alt sınırın altında
kalıp hiç sayılmıyordu. Düzeltildi ve regresyon testine bağlandı: test artık
gerçek bir kalıp kurup primerleri o kalıptan kesiyor, beklenen değer kodun
kendi ölçüsünden değil dizinin kendisinden geliyor.

### 3.2 Denetim betiğim yanlış DNA derişimi kullanıyordu

`check_deliverables.py` Tm'i yeniden ölçerken `--dna 250` nM kullanıyordu,
oysa 03 ve 04 varsayılanı 50 nM. Yeniden ölçülen Tm sistematik olarak kayıyor
ve **197 sahte KRİTİK bulgu** üretiyordu. Varsayılan 50 nM'ye çekildi.

### 3.3 Alan karışımı: bakteriyel hedef, mantar lokusunda

`Sakarolitik_bakteriler` hedefinin kutuları B=6, F=2 dağılımında. Doğru sınıf
olan B'de hiçbir çift bulunamayınca, geriye yalnız yanlış lokustaki F2 çifti
kalıyordu ve tablo hedefi kapsanmış gösteriyordu. Kural denetimi bunu
yakalamaz, çünkü ortada kural ihlali yok. `alan_denetimi.py` ortak modülü
yazıldı; alan bilgisi elle yazılmış bir tablodan değil, konsensüs dosya
adlarından ve `hedefler.tsv` taxid listesinden çıkarılıyor. 08 çapraz alan
tasarımı hiç üretmiyor, 13 bunları çalışma kitabından çıkarıyor, 18 her koşuda
KRİTİK olarak işaretliyor.

### 3.4 09, dışlanan takson dosyasını konuma göre okuyordu

`dislanan_takson.tsv` sütunları konuma göre okunuyordu; dosya biçimi
değişseydi dışlama **sessizce** devre dışı kalacaktı. Başlığa göre okumaya
çevrildi ve tanınmayan başlıkta yüksek sesle duruyor. Beş test eklendi.

### 3.5 Mantar kimliği: Trichoderma değil, Petriella

Kraken2 mantar kutularını Trichoderma / Metarhizium / Podospora diye
etiketliyordu. Ben de bir noktada "aslında Trichoderma" dedim. **Yanlıştı.**
Baskın organizma Microascaceae familyasından, Petriella'ya en yakın:
ITS bölgesinde RefSeq'te %98,2, UNITE'ta %98,8-98,9.

Ara bir çelişki de çözüldü: UNITE'a tam konsensüsle sorulduğunda Hypocreales
çıkıyordu (bitscore 5280'e 872). Bu bir uzunluk yapaylığıydı; UNITE'ın tam
operon kayıtları 3700 bp'lik hizalamalar üretiyor ve hizalamaya korunmuş
18S/28S hâkim oluyor. Yalnız ITS penceresiyle sorulduğunda ilk 15 vuruşun
15'i Microascaceae çıktı. Kanıt: `primer_final/mantar_kimlik_kaniti.md`.

### 3.6 Kutu kurtarmada yön yapaylığı ve birleştirilmiş oran

`recover_bins.py` iki yarı arasında %72 ıraksama bildiriyordu; sebep
tohumların karşıt zincirlerden gelmesiydi. Hizalama tabanlı, yön farkında
karşılaştırmaya çevrildi. Ardından ikinci bir hata çıktı: kapsama eksiği
ikame farkının içine katlanıyordu. Gerçek ikame oranı 0,0021 iken
birleştirilmiş sayı 0,0145 çıkıyordu. Üç sayı ayrı verilir oldu: ikame oranı,
indel oranı, kapsama.

Blastochloris kutusu (B-1_2233851) kurtarıldı: kendi referansı %19 IUPAC
olduğu için minimap2 yalnız 2 minimizer bulabiliyordu; okumaların kendisinden
1449 bp tam kapsamlı konsensüs kuruldu.

### 3.7 Bayat dosya, sağlama toplamı olmadan okundu

`device_stage_files` trend çalışma kitabının **eski** bir kopyasını getirdi
(sağlama toplamı e58d7414, cihazdaki 62b027a7). Bundan sonra her getirilen
dosya cihazdakiyle sağlama toplamı karşılaştırılarak doğrulanır oldu.
Ayrıca bir kez kendi aynamdaki kopyayı değiştirip eşliği bozdum; düzeltilmiş
sürüm iki tarafa birden yazıldı.

### 3.8 Rütbe yüzdeleri %368 topluyordu

`abundance_rank.py` iç içe aynı rütbedeki düğümleri iki kez sayıyordu.
Kraken2 raporunda gerçek şube "P", alt şube "P1", onun altı "P2" diye kodlanır
ve hepsi aynı ana rütbeye katlanır; üst düğümün klanı alt düğümünkini zaten
içerdiği için ikisini birden saymak çift sayımdır. Düzeltmeden önce
barcode10'un yüzdeleri toplamı **%368,54** çıkıyordu. Ata tabanlı ayrıklıkla
düzeltildi; yirmi barkodun hepsi %99,85-100,00 aralığına oturdu.
Bunu kendi aritmetik denetimimle teslimden önce yakaladım.

### 3.9 "Aynı organizma" uyarısı çoğunluk yerine çokluğa bakıyordu

Uyarı, Asetoklastik metanojenler için 19 kutunun 6'sında (çokluk, çoğunluk
değil) tetikleniyordu. Gerçek çoğunluk şartına çevrildi (`d * 2 > n`) ve
destek oranları her zaman gösterilir oldu.

### 3.10 26, yorum satırlarını başlık sanıyordu

`hedefler_referans.tsv` `#` ile başlayan satırlarla açılıyor; `DictReader`
bunları başlık sanıp KeyError veriyordu. Filtrelendi.

### 3.11 26, cinse bakıp "uyuşuyor" diyordu

M. barkeri ürünleri yalnız cins düzeyinde eşleştiği hâlde "uyuşuyor"
yazılıyordu. Tür ve cins uyumu ayrı ayrı raporlanır oldu.

### 3.12 Kimlik sıralaması uzunluğa bakıyordu

`target_identity.py` vuruşları önce uzunluğa göre sıralıyordu; 524 bp'de
%94,47 eşleşme, 504 bp'de %98,21 eşleşmenin önüne geçiyor ve kimlik yanlış
çıkıyordu. Bitscore'a çevrildi.

### 3.13 Kraken2 güven eşiği ezberden seçilemez

Kısa okuma verisi için sık önerilen 0,1 değeri bu ONT verisinde okumaların
**%69'unu** sınıflandırılmamış bırakıyor, çünkü ONT okumalarının k-mer'lerinin
çoğu veritabanında karşılık bulmuyor ve bu k-mer'ler puanın paydasına giriyor.
Eşik bir taramayla veriden seçildi: 0,02.

`reassign_confidence.py` yeniden sınıflandırmaya gerek bırakmıyor; Kraken2'nin
`--output` dosyaları her okumanın k-mer LCA dizisini zaten taşıyor ve güven
puanı tam olarak oradan hesaplanıyor. 106 GB veritabanı da, ham fastq da
gerekmiyor. Ölçüldü: vuruş taksonlarının **%99,84'ü** rapor dosyalarından
kurulan ağaçta yer alıyor. Öz test: güven 0'da tam olarak 0 okuma yer
değiştiriyor.

### 3.14 Bracken önerimi kendi ölçümümle reddedildi

Bracken çalıştırmayı önerdim, sonra güven düzeltmesinden sonra okumaların
yalnız **%17,3'ünün** cins düzeyine ulaştığını ölçüp kendi önerimi geri
çektim. Bracken veritabanı bütünlüğünü varsayar; bu veri o varsayımı
karşılamıyor. Gerekçe: `bolluk_rutbe_kaniti.md`.

### 3.15 ROD ökaryot veritabanıdır, arke/bakteriye atanmıştı

**Bu, oturumun en ciddi bulgularından biri.**

`ROD_v1.2_operon_variants.fasta` başlıklarındaki soyağacı sayıldı:
60320 kaydın **60320'si Eukaryota**, 0 Bacteria, 0 Archaea, içinde 9753 mantar.
ROD, Krabberød ve arkadaşlarının Ribosomal Operon Database'idir: genom
derlemelerinden çıkarılmış tam uzunlukta **ökaryot** rDNA operonları.

Buna rağmen `SINIF_DB_GENIS` içinde ROD A1/A2/B sınıflarına (arke ve bakteri)
atanmış, F1/F2'ye atanmamıştı. Sonucu iki yönlüydü:

- **Sahte temizlik:** 71 arke/bakteri çifti için "ROD'da hedef dışı ürün yok"
  yazıldı. Bu özgüllük kanıtı değildir; o veritabanında o alandan tek dizi
  bile yok. Aynı taramada SILVA SSU 71 kaydın 71'inde ürün buldu (41324 ürün),
  yani primerler ölçülemez durumda değildi, ölçülen şey boştu.
- **Kaçırılan ölçüm:** ROD'un 9753 mantar tam operonu mantar primerleri için
  hiç taranmadı.

`fungi.18SrRNA.fna` (4037 dizi, indeksli) ise hiçbir sınıfın listesinde yoktu.
Bu da benim atlamamdı.

**Düzeltme:** sınıf başına elle veritabanı seçimi kaldırıldı. Her sınıf eldeki
on rDNA veritabanının hepsini görüyor, liste tek bir ortak kaynaktan
türetiliyor. Gerekçe: gerçek PCR'da primer yalnız kendi alanının rDNA'sıyla
değil ortamdaki tüm DNA ile karşılaşır; bir arke primerinin bakteri 23S'inde
yanlış bağlanması tam olarak aranan hata türüdür.

### 3.16 İndeks iddiam eksikti

`ref_all.fna` ve `ref_all2.fna` için "indeks yok" demiştim. Doğrusu: BLAST
indeksi yok, ama **mfeprimer indeksi var** (eski `.primerqc` biçimi, artı
`.fai` ve `.json`). Yalnız `.nin` aradığım için eksik söyledim. Sekiz BLAST
bileşeninin (`.ndb .nhr .nin .njs .not .nsq .ntf .nto`) hepsi ve çok ciltli
`.00.*` kalıbı ayrı ayrı tarandı.

Ayrıca o indeksleri kurmaya gerek de yok. `.fai` kimlik kümeleri
karşılaştırıldı:

```
archaea.16S     1160
bacteria.16S   26877
fungi.ITS      20394
fungi.18SrRNA   4037
fungi.28SrRNA  12890
toplam         65358   <-- ref_all2.fna da tam 65358
ref_all2'de olup beşte olmayan : 0
beşte olup ref_all2'de olmayan : 0
```

`ref_all2.fna` bu beş dosyanın birebir birleşimi, `ref_all.fna` ise 28S ve
18S'siz hâli (65358 - 12890 - 4037 = 48431). İkisi de zaten taranan
dosyaların alt kümesi. Önceki `makeblastdb` önerim geri alındı.

`SILVA_138.2_LSUParc.fasta` bilerek dışarıda: Parc kümesi kısmi ve düşük
kaliteli kayıtları da içerir, `LSURef_NR99` aynı kapsamın %99 tekrarsız
küratörlü temsilcisidir. SSU tarafında da Parc değil NR99 kullanılıyor.

Yan bulgu: `SILVA_SSURef_NR99.fasta` (824 MB) ve `SILVA_LSURef_NR99.fasta`
(297 MB) hem `REFERANS_DB` hem `_to_delete/REFERANS_DB` içinde duruyor ve
`cmp` ile birebir aynı oldukları doğrulandı. 1,1 GB gereksiz kopya.

Ayrıntı: `veritabani_kapsami_kaniti.md`.

### 3.17 Boş ölçüm, temiz ölçüm gibi görünüyordu (kapsam denetimi)

ROD hatasının asıl dersi buydu. `external_databases.py`'ye kapsam denetimi
eklendi: her (sınıf, veritabanı) ikilisi için, taramadan önce o sınıfın kendi
konsensüs dizileri veritabanına megablast ile aranıyor ve en uzun hizalama
ölçülüyor. Eşik uydurulmuyor, veriden geliyor: aranan ürün en fazla
`prod_max` baz uzunluğunda olduğu için, veritabanındaki en uzun benzer bölge
`prod_max`'tan kısaysa o ürünün orada oluşması zaten mekanik olarak
imkânsızdır.

Sonuç: 970 kaydın **474'ü KAPSAM_YOK**. Yani "hedef dışı ürün yok" diye
okunabilecek 674 satırın yalnız 201'i gerçekten bir şey ölçüyordu.

### 3.18 "Hedef dışı ürün" sayısı, hedefin kendisini de sayıyordu

Vuruş başlıklarını veritabanlarından çözdüm. En yüksek sayıları veren
kayıtların bir kısmı **hedefin kendisiydi**:

| Kayıt | Ham sayı | Vuruşlar gerçekte |
|---|---|---|
| Asetoklastik × archaea.16S | 306 | *Methanothrix soehngenii* GP6, Opfikon |
| Asetoklastik × SILVA | 1143 | hepsi Methanosaetaceae;Methanothrix |
| Methanothrix soehngenii × SILVA | 308 | hepsi Methanothrix |
| **Nitrosocosmicus × SILVA** | **1119** | **Halobacteria, Methanoperedenaceae, Cenarchaeum** |
| **Petrimonas × SILVA** | **707** | **Clostridium, Bacteroides** |
| **Trichoderma × fungi.28S** | **1150** | **Calonectria, Acremonium, Trichothecium** |

İlk üçü hata değil, primerin işini yapmasıdır. Ham sayıya göre sıralamak
kullanıcıyı **yanlış primerleri düzeltmeye** gönderiyordu.

Düzeltme: her ürün, oluştuğu referansın taksonuna göre `urun_kendi_taksonda`,
`urun_yabanci_taksonda`, `urun_takson_bilinmiyor` diye ayrılıyor. Takson
adları elle yazılmıyor; `hedefler.tsv` + `taxid_adlari.tsv` (beyan edilen ad)
ve `hedef_kimlik.tsv` (ölçülen ad) olmak üzere iki kaynaktan geliyor.
Evrensel hedefler (`*A`, `*B`, `*F`) ayrıca işaretleniyor; onlarda çok taksonu
birden çoğaltmak amaç olduğu için hiçbir ürün yabancı sayılmıyor.

### 3.19 "Yabancı takson" da tek başına yeterli değildi

İkinci turdan sonra ölçüldü: hedeflerin bir kısmı tek takson değil **işlevsel
grup**.

| Hedef | "Yabancı" vuruşlar | Gerçekte |
|---|---|---|
| Hidrojenotrofik metanojenler | *Methanobacterium alcaliphilum*, *Methanosphaera stadtmanae* | ikisi de hidrojenotrofik metanojen |
| Nitrosocosmicus AOA | *Nitrosotalea*, *Nitrosopumilus* | ikisi de amonyak oksitleyen arke |
| Zoopagomycota | *Piptocephalis moniliformis* | kendisi Zoopagomycota |
| Asetoklastik | *Methanolobus*, *Methanosalsum*, *Methanohalophilus* | aynı aile |

Buna karşılık Petrimonas → *Flavobacterium*, *Phocaeicola vulgatus* ve
Trichoderma → *Calonectria*, *Acremonium* gerçekten uzak.

Düzeltme: uzaklık elle yazılmış bir "işlevsel grup" tablosuyla değil
**veriden** ölçülüyor. SILVA, ROD, UNITE ve PR2 başlıkları soyağacını taşıyor;
hedefin kendi soyağacı, kendi taksonunda ürün veren referanslardan
çıkarılıyor; her yabancı vuruşun bu soyağacıyla paylaştığı derinlik
ölçülüyor. Aile düzeyine kadar paylaşan `yakin`, paylaşmayan `uzak`.

### 3.20 Uzaklık eşiği kendi vuruş sayısına göre kayıyordu

Referans soyağacını önce **ortak ön ekten** alıyordum. O zaman tek vuruşta
derinlik 7, iki çeşitli vuruşta 4 çıkıyor ve aynı yabancı vuruş bir kayıtta
uzak, ötekinde yakın sayılabiliyordu. **Baskın** soyağacına çevrildi, eşik
sabitlendi. Testle bağlandı.

### 3.21 Sıfır ürün veren çift "temiz" görünüyordu (çift düzeyinde)

`Sakarolitik_bakteriler` F2'nin üç çifti sekiz kapsanan veritabanının
hiçbirinde ürün vermiyordu, **kendi taksonunda bile** (`kendi=0`). Özet bunları
"hiçbir yabancı taksonda ürün vermeyen" diye sayıyordu. ROD hatasının çift
düzeyindeki aynısı. Artık `kendi=0` olan çiftler "temiz" sayılmıyor,
`OLCUM GECERSIZ` diye ayrı raporlanıyor.

### 3.22 27, referans çiftlerini sessizce düşürüyordu

`primer_referans.tsv`'de ad `Methanosarcina_barkeri_referans`; `_referans`
ekini soyunca `Methanosarcina_barkeri` çıkıyor, oysa `hedefler.tsv`'deki ad
`Methanosarcina_barkeri_turu`. Eşleşme tutmayınca hedefin **tek** primer takımı
(de novo hiç çifti yok) sessizce düştü ve hedef `CIFT_YOK` göründü. Proteiniphilum
ile Podospora'da adlar birebir tuttuğu için hata tek hedefte göze çarpmadı.
Önce birebir, sonra ön ek eşleşmesi deneniyor; hiçbiri tutmazsa yüksek sesle
"ÖLÇÜM DIŞI KALDI" yazılıyor.

### 3.23 "sp" tür adı sayılıyordu

UNITE'ta `s__Trichoderma_sp` biçimindeki kayıtlar tür adı sayılıyor ve panelde
ayrı bir "tür" gibi duruyordu: yalnız Trichoderma için 16910, Marasmius için
1712, Podospora için 1326 kayıt. Tür adı belli olmayan bir kayıtta ürün
oluşması, tür ayrımının başarısız olduğunu göstermez. `sp`, `spp`, `cf`, `aff`
artık tür adı sayılmıyor. RefSeq'teki `Trichoderma sp.` noktalı olduğu için
zaten takılmıyordu; kaçan yalnız alt çizgili biçimlerdi.

### 3.24 Cins düzeyi ölçümü, cins dışını gizliyordu

Eski sayım yalnız "panelin kaç türünde ürün var" diyordu. Proteiniphilum'un
beş çiftinden ikisi *Fermentimonas caenicola*'yı da çoğaltıyor, ki başka bir
cins; eski etiket ikisini de `CINS_ICI_3_3` yazıyor, yani cins özgüllüğünü
ihlal eden çift en geniş kapsamlı çift gibi görünüyordu. Cins kararı
`CINS_OZGUL` / `CINS_AYRIMI_YOK` / `CINS_ICINDE_URUN_YOK` diye ayrıldı. Ölçülen
kimlik burada hedef sayılmıyor: Proteiniphilum istendiyse Fermentimonas çapraz
çoğaltmadır, kutulardaki organizma o olsa bile.

### 3.25 RefSeq tip suşu kayıtları panele hiç girmiyordu

RefSeq ITS ve 28S kayıtları `NR_172285.1 Petriella musispora CBS 745.69 ITS
region; from TYPE material` biçimindedir. Önceki sürümde başlıkta noktalı
virgül varsa RefSeq dalı hiç çalışmıyor, dolayısıyla **tip materyali**
kayıtlarının tür adı çıkarılamıyor ve panele alınmıyorlardı. Panel tam da bu
kayıtlara en çok ihtiyaç duyulan yerdi: tür ayrımının altın standardı tip suşu
dizileridir. Yalnız Petriella için 35 kayıt böyle düşüyordu. SILVA'nın
soyağaçlı başlıkları bu daldan zaten geçemez, ayrı korumaya gerek yok.

### 3.26 Referans tasarımında rakip kümesi, doğrulama panelinden dardı

**Yeni tasarımların kök sebebi.** `design_from_reference.py` rakipleri tek
veritabanından ve ad başına altı kayıtla topluyordu. Podospora tasarımı 29
rakip diziyle yapılmış (`fungi.ITS.fna`'da toplam 14 Podospora kaydı var);
27'nin doğrulama paneli ise UNITE dahil 242 kayıt ve 50 tür. Tasarım
*P. anserina* ile *P. comata*'yı dışladığını sanırken geniş panelde ikisini de
çoğaltıyordu. Yani başarısızlığın sebebi tasarım motoru değil, tasarım anında
**görülmeyen** rakiplerdi.

Sentetik olarak yeniden üretildi: dar rakip kümesiyle **3548 çift** geçiyor,
geniş kümeyle **0**. Eski davranışta ikisi de 3548 veriyordu.

Düzeltme: `veritabani` sütunu virgülle birden çok dosya alıyor; rakip listesi
elle yazılanla sınırlı değil, hedefin cinsindeki öteki türler **veriden**
bulunup ekleniyor; tür adı tanımı 27'den alınıyor, iki ayrı kopya olmasın diye.

### 3.27 Cins adı alt dize olarak aranıyordu

`sec()` fonksiyonu `ic='Bacteroides'` yazıldığında *Parabacteroides* ve
*Acetobacteroides* kayıtlarını da hedef **üyesi** sayıyordu.
`bacteria.16S.fna`'da sözcük sınırlı Bacteroides 86 kayıt verirken, alt dize
araması bunlara 17 Parabacteroides ile 1 Acetobacteroides kaydını da katardı;
primer o zaman başka cinslerde de ürün vermek zorunda kalır ve cins
özgüllüğü daha tasarım anında kaybedilirdi. Sözcük sınırı eşleşmesine
çevrildi.

---

## 4. Bilimsel bulgular

### 4.1 Baskın mantar organizması

Microascaceae, *Petriella*'ya en yakın; Kraken2'nin Trichoderma / Metarhizium /
Podospora etiketleri yanlış. İki bağımsız veritabanı ayırt edici ITS
bölgesinde uyuşuyor. Kanıt: `mantar_kimlik_kaniti.md`.

### 4.2 Bolluk yalnız verinin desteklediği rütbede okunabilir

- Arke: cins düzeyi, okumaların %86-97'si
- Bakteri: şube düzeyi, %8,8
- Mantar: sınıf/şube düzeyi, %0,1-1,5

Eşik yapaylığı olmadığı güven-0 kontrol koşusuyla doğrulandı; okuma uzunluğu
ve okuma sayısı karıştırıcıları açıkça dışlandı. Kanıt:
`bolluk_rutbe_kaniti.md`.

### 4.3 Referans tasarımlı primerlerin numunede çoğalttıkları

- **Proteiniphilum doğrulandı:** 5 çiftin 5'i → *Proteiniphilum acetatigenes*,
  %98,3-100
- **M. barkeri tür düzeyinde tutmuyor:** ürünler *M. thermophila* %97,89 ve
  *M. flavescens* %98,53. Cins düzeyinde etiketlenmeli
- **Podospora yok:** 0 ürün okuması, üçüncü bağımsız doğrulama

Kanıt: `primer_referans/referans_ve_mfeprimer_bulgulari.md`.

### 4.4 Bir gerilim ve çözümü

27, *M. barkeri* primerlerinin referans dizilere karşı tür ayrımı yaptığını
söylüyor (16 türlük panelde 1 çift çaprazsız, 4 çift eşik içi). 26 ise bu
primerlerin **bu numunede** çoğalttığı ürünlerin *M. thermophila* %97,9 ve
*M. flavescens* %98,5 çıktığını ölçtü. Çelişki değil: numunedeki organizma
adlandırılmış hiçbir türle özdeş değil, ikisine de yaklaşık %2 uzak. Primer
tasarımı sağlam, numunedeki organizma *M. barkeri* değil. Raporda bu ayrım
açıkça yazılmalı.

### 4.5 mfeprimer ile blastn ayrışması

31 ayrışan kayıt / 27 çift, hepsi mfeprimer-buluyor / blastn-bulmuyor
yönünde. Öncelik: Sakarolitik F2 (12760, zaten dışlandı), Alistipes B
(545, 485, 19), Zoopagomycota F1 (137/135/135/104), Metilotrofik A1 (60, 11),
Nitrosocosmicus (44, 44, 2).

Ölçülmüş bir mfeprimer davranışı: `--misEnd 3` ile çalıştırıldığında, 3'
uçtaki son bazı değiştirilmiş bir primer bozulmamış primerle **aynı** sayıda
amplikon veriyor (323'e 323). mfeprimer'ın sessiz başarısızlığı da yakalanıyor
(`"no valid db"` çıktısı ile dönüş kodu 0).

---

## 5. Ölçüm sonuçları

### 5.1 Geniş dış veritabanı taraması (14, üçüncü tur)

970 kayıt. 674'ünde hiç ürün yok, bunların kapsamı doğrulanmış olanı 201.
Kapsanan veritabanıyla ölçülen özgül çift: 67. Kendi taksonunda ürün vererek
ölçülebilir olanı: 47. Uzak taksonda hiç ürün vermeyen: 2. Kendi taksonunda
bile ürün vermeyen (ölçüm geçersiz): 20.

**Gerçekten temiz iki çift:**

```
Nitrosocosmicus_AOA  A1   kendi=26  yakin=10  uzak=0
   F ACCGCGTGTCACTATCGC
   R TCGATAGTACCAATTAGGCACCAC

Microascaceae_askomikot  F2   kendi=1  yakin=0  uzak=0
   F ATCAATAAGCGGAGGAAAAGAAACC
   R CCTCTTCAAATTACAACTCGGACTG
```

İkincisinde `kendi=1`, yani dayanak tek kayıt; zayıf kanıt.

**Hedef başına en iyi ölçülebilir çift (uzak takson ürünü):**

| Hedef | En iyi uzak | Yakın | Kendi |
|---|---|---|---|
| Nitrosocosmicus_AOA A1 | 0 | 10 | 26 |
| Microascaceae_askomikot F2 | 0 | 0 | 1 |
| Metilotrofik_metanojen A1 | 1 | 2 | 5 |
| Nitrosocosmicus_AOA A2 | 2 | 97 | 43 |
| Methanothrix_soehngenii A1 | 8 | 0 | 72 |
| Hidrojenotrofik A2 | 8 | 0 | 44 |
| Methanothrix_soehngenii A2 | 10 | 0 | 116 |
| Hidrojenotrofik A1 | 11 | 0 | 34 |
| Methanosarcina_mazei A2 | 42 | 10 | 58 |
| Petrimonas_cinsi B | 43 | 26 | 19 |
| Zoopagomycota F1 | 50 | 0 | 6 |
| Asetoklastik A1 | 227 | 0 | 749 |
| Trichoderma_cinsi F2 | 2192 | 85 | 215 |

**Hedef içi yayılım kritik:** Nitrosocosmicus A1'in beş çiftinden biri 0, bir
başkası 1098 veriyor. Çift seçimi sonucu tamamen değiştiriyor; Excel'deki
"önerilen çiftler" listesi bu ölçüte göre yeniden kurulmalı.

**Hiç ölçülemeyen (kendi taksonunda ürün yok):** Alistipes_cinsi B 5 çift,
Bacteroidaceae_ailesi B 5 çift, Sakarolitik_bakteriler F2 5 çift,
Zoopagomycota F1 4 çift, Trichoderma_cinsi F2 1 çift.

### 5.2 Düzey denetimi (27)

Çapraz tür hoşgörüsü 2.

**Tür özgüllüğü:**

| Hedef | Sonuç |
|---|---|
| *Methanothrix soehngenii* | 10 çiftin 10'u TUR_OZGUL, hiç çapraz yok. **Karşılıyor** |
| *Methanosarcina barkeri* | 1 çift TUR_OZGUL, 4 çift eşik içi (*M. baltica* vb.). **Karşılıyor** |
| *Methanosarcina mazei* | 1 çift eşik içi, tek çapraz *M. soligelidi*. **Karşılıyor** |
| *Dictyostelium discoideum* | 5 çiftin 1'i eşik içi, çapraz *Marasmius rhyssophyllus*; hedef türde çoğalttığı gerçekten *D. discoideum*. **Karşılıyor** |
| *Trichoderma asperellum* | 2 çift 3-4 çapraz (*Petriella guttulata*, *setifera*), 1 çift hedef türde ürün yok. Hedef türde çoğalttığı *T. asperellum* değil *Petriella musispora*. **Karşılamıyor** |
| *Podospora pseudopauciseta* | 5 çiftin 5'i hedef türde hiç ürün vermiyor, kardeş türleri çoğaltıyor. **Karşılamıyor** |

**Cins özgüllüğü:**

| Hedef | Sonuç |
|---|---|
| **Proteiniphilum** | 5 çiftin 3'ü CINS_OZGUL; 2'si *Fermentimonas caenicola*'yı da çoğaltıyor. **Karşılıyor** |
| **Petrimonas** | 5 çiftin 5'i panel içinde CINS_OZGUL (*P. mucosa*, *P. sulfuriphila*). Ancak geniş taramada panel dışında uzak taksonlar var. **Kısmen** |
| Bacteroides | 65 türün 0'ı. **Karşılamıyor** |
| Alistipes | 22 türün 0'ı. **Karşılamıyor** |

**On hedefin altısı istenen düzeyi karşılıyor.**

Bir düzeltme: Dictyostelium için önce "panel yok" demiştim. Yanlıştı. Üç
mantar veritabanında sıfır kaydı olduğu doğru, ama PR2'de (ökaryot SSU) 112
kayıt var ve çift bunları çoğaltıyor. "Mantar veritabanlarında yok"
gözlemini "hiçbir yerde yok" diye genelleştirmiştim.

### 5.3 Toplantı hedeflerinin genel durumu

Yirmi hedefin **on beşi** kullanılabilir çiftle kapsandı. Kapsanmayanlar:

- **Proteolitik sintrofik bakteriler:** hiç çift yok, ne de novo ne referans.
  17664 elemenin 17623'ünü tek başına *Cloacibacillus porcorum* yapıyor.
- **Sakarolitik bakteriler:** 5 çift geçmiş görünüyor ama hepsi F2 sınıfında,
  yani mantar lokusunda. 18 her koşuda KRİTİK işaretliyor, Excel'e girmiyorlar.
- **Podospora pseudopauciseta:** organizma numunede yok, konusuz kalıyor.

*Methanosarcina barkeri* ile *Proteiniphilum* yalnız referans setiyle
kapsanıyor.

---

## 6. Primer BLAST durumu

**Yapıldı**, rDNA tarafında NCBI'nin varsayılanından geniş:

- `external_databases.py`: `blastn -task blastn-short`, on rDNA veritabanı, beş
  sınıf, 970 kayıt, kapsam denetimi, kendi/yabancı/yakın/uzak ayrımı
- `mfeprimer_layer.py`: mfeprimer `spec`, ikinci bağımsız yöntem
- `check_taxonomic_level.py`: kardeş tür panellerine karşı blastn
- `22` ve `26`: kimlik blastn'leri

Primer-BLAST'ın algoritmik olarak yaptığı iş (primerleri ayrı ayrı ara, ters
yönde ve ürün boyu aralığında buluşan vuruş çiftlerini bul) yapıldı.

**Açık kalan iki yüzey:**

1. Elimizdeki veritabanlarının hepsi rDNA. Bir primerin genomun rDNA dışı bir
   yerinde bağlanıp bağlanmadığına bakılmadı; NCBI Primer-BLAST'ın `nt`
   araması bunu kapsar, koşulmadı.
2. `14` ve `19` yalnız `primer_final.tsv` okuyor. Referanstan tasarlanan
   çiftler dış veritabanı taramasına hiç girmedi. Oysa Proteiniphilum'un
   teslim edilebilir çiftleri tam olarak onlar.

---

## 7. Regresyon takımı

`regression_test.py`, **137 test, 137/137 geçiyor.** Her testin beklenen
sonucu toplantı kararlarından ya da bilinen matematikten türetilir; kodun
kendi yardımcı fonksiyonlarına güvenilmez.

| Grup | Konu |
|---|---|
| 1 | ters tümleyen |
| 2 | kompozisyon kuralları |
| 3 | IUPAC zincir kuralı |
| 4 | bağlanma kuralı, kaba kuvvetle karşılaştırma (400 deneme) |
| 5 | ürün geometrisi, gerçek kalıpla |
| 6-9 | Wilson, konsensüs, ayırt edilemezlik, hizalama arka ucu |
| 10 | ham okuma taraması ile tasarım kuralı tutarlılığı |
| 11 | dış veritabanı ürün boyu |
| 12 | 08 ile 09 arasında dosya biçimi sözleşmesi |
| 13 | dış veritabanı küme eşlemesi ve kapsam denetimi |
| 14 | kendi taksonu ile yabancı takson ayrımı |
| 15 | yabancı vuruşun uzaklığı (soyağacı derinliği) |
| 16 | karar düzeyi denetimi, tür ve cins özgüllüğü |
| 17 | referans tasarımında rakip kümesi |

---

## 8. Senkronizasyon

`SENKRON.sh` SHA256 manifestosu üretir; `[BETIK]` ve `[CIKTI]` bölümleri var,
`--dogrula` ile iki taraf karşılaştırılır. Kasıtlı bir kayma yaratılarak
yakaladığı doğrulandı. 42 betik dosyası iki tarafta bayt bayt aynı.

Bulut kutusundan cihaza yazarken Windows yolu kullanılmalı
(`C:\Users\yerli\Masaüstü\PrimerTasarlama\...`); `/sessions/...` biçimi
reddediliyor. `device_bash` dosya silemiyor; silinecekler `_to_delete/`
altına taşınıyor.

---

## 9. Yeni tasarım hazırlığı (2026-08-01)

Karar: numunede bulunmayan organizmalar için **hem** karardaki adla referanstan
tasarım **hem** gerçek organizmaya göre tasarım yapılacak; başarısız üç grubun
hepsi denenecek.

`hedefler_referans.tsv` dokuz satıra çıkarıldı:

| Ad | Sınıf | Veritabanı | Hedef |
|---|---|---|---|
| Methanosarcina_barkeri_referans | A1 | archaea.16S | *M. barkeri* |
| Proteiniphilum_cinsi_referans | B | bacteria.16S | *P. saccharofermentans*, *P. acetatigenes* |
| Podospora_pseudopauciseta_referans | F1 | fungi.ITS + UNITE | *P. pseudopauciseta* |
| Bacteroidaceae_ailesi_referans | B | bacteria.16S | Bacteroides cinsi |
| Alistipes_cinsi_referans | B | bacteria.16S | Alistipes cinsi |
| Microascaceae_askomikot_referans | F2 | fungi.ITS + UNITE | *T. asperellum* |
| Petriella_musispora_referans | F2 | fungi.ITS + 28S + UNITE | *P. musispora* |
| Proteolitik_Cloacibacillus_referans | B | bacteria.16S | *C. porcorum*, *C. evryensis* |
| Sakarolitik_Sphaerochaeta_referans | B | bacteria.16S | Sphaerochaeta cinsi |

Elle yazılan kardeş tür adları listeden **çıkarıldı**: artık veriden
bulunuyorlar ve elle liste eksik kalabileceği için yanlış güven veriyordu.
Elle kalanlar yalnız başka cinslerden rakipler.

*Petriella musispora* için taxid kesin bulunamadı ve **uydurulmadı**.
`hedefler.tsv`'ye isteğe bağlı bir `hedef_tur` (7.) sütunu eklendi; tür adı
oraya doğrudan yazılıyor. Sütun taxid adını **eklemiyor, değiştiriyor**: eklese
hedef tür kümesine Trichoderma da girerdi ve Trichoderma çoğaltan bir çift
"hedef türde ürün var" sayılırdı.

Üç yeni satır `karar 5` diye numaralandırıldı ve not sütununa "ÖLÇÜMDEN
TÜRETİLDİ, onaya sunulur" yazıldı; toplantıda konuşulmadılar, kararmış gibi
görünmesinler.

---

## 10. Kalan işler

1. Rapor (docx) ve Excel'in bütün yeni ölçümlerle yeniden kurulması: kapsam
   denetimi, kendi/yabancı/yakın/uzak ayrımı, düzey denetimi, rütbe kapsaması,
   kimlik bulguları
2. Yeni referans tasarımının çalıştırılması ve 27 ile doğrulanması
3. Referans setinin 14 ve 19 taramalarına dahil edilmesi
4. NCBI Primer-BLAST ile `nt` taraması (rDNA dışı hedef dışı bağlanma)
5. Sakarolitik için 08 yamalı tam yeniden koşu (~2,5 saat), ham tablodaki
   alan karışımı satırlarının kalkması
6. `_to_delete/` klasörünün elle silinmesi (1,1 GB gereksiz kopya)
7. Islak laboratuvar doğrulaması

---

## 11. Çalıştırma komutları

```bash
PT="/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama"
cd "$PT/steps"

# Regresyon takımı
python3 regression_test.py

# Senkron doğrulama
bash SENKRON.sh --dogrula

# Geniş dış veritabanı taraması (12 dakika)
bash AGIR_ISLER.sh --yalniz H

# Düzey denetimi (5 dakika)
python3 check_taxonomic_level.py \
  --hedefler hedefler.tsv --adlar taxid_adlari.tsv \
  --final "$PT/primer_final" --referans "$PT/primer_referans/primer_referans.tsv" \
  --db "$PT/REFERANS_DB" --kimlik "$PT/primer_final/hedef_kimlik.tsv" \
  --is-parcacigi 4 --out "$PT/primer_final/duzey_denetimi.tsv"

# Yeni referans tasarımı (yarım saate kadar)
python3 design_from_reference.py --db "$PT/REFERANS_DB" --pt "$PT" \
  --hedefler-ref hedefler_referans.tsv --out "$PT/primer_referans"
```

---

## 12. Kanıt belgeleri

| Dosya | İçerik |
|---|---|
| `mantar_kimlik_kaniti.md` | Petriella bulgusu, iki veritabanı, ITS penceresi |
| `bolluk_rutbe_kaniti.md` | rütbe kapsaması, Bracken'in neden çalıştırılmadığı |
| `referans_ve_mfeprimer_bulgulari.md` | referans setin numunede çoğalttıkları, mfeprimer ayrışmaları |
| `veritabani_kapsami_kaniti.md` | indeks dökümü, ROD hatası, kapsam denetimi |
| `CALISMA_KAYDI.md` | bu belge |

## 13. Arşivlenen ara çıktılar

| Dosya | Neden arşivlendi |
|---|---|
| `primer_final/_eski/dis_veritabani_genis_ROD_hatali_20260801.tsv` | ROD yanlış sınıfa atanmışken |
| `primer_final/_eski/dis_veritabani_genis_taksonsuz_20260801.tsv` | kendi/yabancı ayrımı yokken |
| `primer_final/_eski/dis_veritabani_genis_uzakliksiz_20260801.tsv` | yakın/uzak ayrımı yokken |
| `primer_final/_eski/duzey_denetimi_ilk_tur_20260801.tsv` | referans eşleşme hatası ve "sp" hatası varken |
| `primer_final/_eski/duzey_denetimi_cinssiz_20260801.tsv` | cins kararı eklenmeden önce |
| `primer_referans/_eski/primer_referans_dar_rakip_20260801.tsv` | dar rakip kümesiyle tasarlanmışken |
