# WSL betikleri, yeniden sınıflandırma ve N analizi

Kararınız uyarınca A ve C paralel yürütülüyor. A sizin WSL'de çalıştıracağınız
Kraken2 kolu, C benim hazırladığım N analizi kolu. C, sınıflandırmayı beklemiyor
çünkü mevcut konsensüs ve fastq dosyalarıyla çalışıyor.

## Kraken2, 106 GB veritabanı ve 16 GB RAM

`01` betiği veritabanı boyutunu kullanılabilir RAM ile karşılaştırıp
`--memory-mapping` bayrağını kendiliğinden açıyor, yani 16 GB RAM ile 106 GB
veritabanı çalıştırılabilir. Ali de zaten pluspf16'yı bellek eşlemesiyle
kullanmış. Daha büyük bir veritabanınız varsa `--db` ile onu verin, betik aynı
kararı onun boyutuna göre yeniden verir.

Tek bilinmeyen hız. Bu yüzden betikte `--only-benchmark` kipi var: en küçük
dosyayı işleyip gerçek süreyi ölçüyor ve toplam için sayı basıyor. Bu bir engel
değil, planlama verisi; ölçümü görüp tam koşuyu `--yes` ile başlatırsınız.
Süreyi kısaltmak isterseniz `--threads` ile çekirdek sayısını yükseltmek ve
veritabanını en hızlı diske almak en doğrudan iki müdahale.

## Ham okuma dosyaları

`PrimerTasarlama\fastq files` klasöründeki dosyalar ham okuma değil, her
örnekten en bol beş taksonun çıkarılmış hali; bu hocanın istediği biçim.
Yeniden sınıflandırmayı ham dosyalar üzerinden yapacağız, `01` betiğine
`--in` ile ham fastq klasörünü verin. Alt kümelerle yeniden sınıflandırma
yapılmayacak, çünkü o küme zaten pluspf16'nın kararlarıyla seçilmiş ve
Zoopagomycota ile Microascaceae aramasını kısıtlar.

## Ali'nin konsensüs betiğinde sessiz bir hata var

Bunu sandbox'a samtools 1.19.2 kurup sentetik bir BAM ile doğruladım, tahmin
değil. `consensus2.sh` betiği `samtools consensus` komutunu iki kez çağırıyor:

```
samtools consensus -a --use-qual              # yorumda "degenerate consensus"
samtools consensus -a --use-qual -c 0.9       # yorumda "strict"
```

Sizin samtools 1.22.1 yardım çıktısına göre öntanımlı mod `bayesian`, ve
`-c/--call-fract` ile `-q/--use-qual` yalnızca `simple` modun seçenekleri. Test
sonucu şu:

| Çağrı | Yüzde 60 / 40 dağılan pozisyonda çıktı |
|---|---|
| `-a --use-qual` | **N** |
| `-a --use-qual -c 0.9` | **N**, birinciyle bayt bazında aynı |
| `-a -A` | **M**, yani A veya C |
| `-a -m simple --use-qual -c 0.9` | N |
| `-a -m simple -A -c 0.9` | M |

Yani `-c 0.9` bayesian modda hiçbir etki yapmıyor, betikteki iki çağrının çıktısı
birebir aynı. Daha önemlisi `-A/--ambig` verilmediği için samtools hiçbir koşulda
IUPAC kodu basamıyor; gerçekten iki allelli bir pozisyonda bile N yazıyor.

Sonuç: projedeki `_consensus_strict.fasta` dosyalarındaki 1.897 N, iki bambaşka
nedeni ayırt edilemez biçimde birleştiriyor. Bir kısmı düşük okuma derinliği, bir
kısmı gerçek suş içi değişkenlik. `02` betiği bu ikisini ayırıyor. Klasördeki
"dejenere konsensüs" hiç dejenere olmamış, o yüzden karşılaştırma için yeniden
üretilmesi gerekiyor.

## Sıra

### 1. Ortam kontrolü, hiçbir şeyi değiştirmez

```bash
cd "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama/WSL_betikleri"
bash 00_ortam_kontrol.sh 2>&1 | tee ortam_raporu.txt
```

Çekirdek ve RAM ölçer, araç sürümlerini listeler, `kraken2`, `bracken`,
`samtools consensus`, `minimap2`, `blastn` yardım çıktılarını
`yardim_ciktilari/` klasörüne yazar, diskte `hash.k2d` arayarak Kraken2
veritabanlarını bulur, her biri için boyutu ve Bracken `kmer_distrib`
dosyalarının varlığını raporlar, ham barkod fastq dosyalarını arar ve
`kraken results\*_output` dosyalarından medyan okuma uzunluğunu ölçer.

### 2. Kraken2 hız ölçümü

```bash
bash 01_kraken2_yeniden_siniflandir.sh \
     --db  ~/k2db \
     --in  /ham/fastq/klasoru \
     --out "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama/kraken_yeni" \
     --only-benchmark
```

Betik 16 GB RAM ile 106 GB veritabanını görünce `--memory-mapping` bayrağını
kendiliğinden açar ve bunu loga yazar. En küçük dosyayı işler, süreyi ölçer,
tahmini basar ve durur. Tahmini gördükten sonra `--only-benchmark` yerine
`--yes` ile tam koşuya geçilir.

Betiğin güvenlik önlemleri:

- Kullanacağı her bayrağı önce `kraken2 --help` çıktısında arar, bulamazsa durur.
- Bellek eşlemesi kararını RAM ile veritabanı boyutundan türetir, elle yazmaz.
  `--force-mmap` ve `--no-mmap` ile geçersiz kılınabilir.
- İş parçacığı sayısını `nproc` eksi iki seçer, elle verilen değer çekirdek
  sayısını aşarsa uyarır.
- Taban ad çakışmasını baştan tespit eder ve durur. `barcode03.fastq` ile
  `barcode03.fastq.gz` aynı çıktı adına ineceği için biri sessizce düşerdi;
  bu projede `fastq files\F2-1` klasöründeki kopyalar tam bu türden bir sessiz
  çift sayıma yol açmıştı.
- Var olan çıktının üzerine yazmaz, atlar.

### 3. N analizi, sınıflandırmayı beklemez

```bash
bash 02_N_analizi.sh \
     --pt  "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama" \
     --out "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama/N_analizi"
```

Her hedef için `minimap2 -ax map-ont` ile okumaları kendi konsensüsüne hizalar,
sonra iki bağımsız ölçüm alır: `samtools consensus -a -A` ile IUPAC kodlarına
izin veren konsensüs, ve `samtools consensus -a -A -f PILEUP` ile pozisyon başına
derinlik, IUPAC çağrı ve baz dizisi. Sınıflandırma PILEUP pozisyon numaralarına
dayanır, FASTA uzunluğuna değil; ikisi arasında ayrılık varsa sayısı loglanır.

Her N pozisyonu dört sınıftan birine giriyor:

| Sınıf | Ölçüt | Sonuç |
|---|---|---|
| `iki_allelli` | derinlik yeterli, IUPAC iki bazlı kod, ikinci bazın oranı eşiğin üstünde | gerçek suş içi değişkenlik, primer ayağı için kalıcı yasak |
| `dusuk_derinlik` | derinlik eşiğin altında | okuma desteği yok, 85 maskesi kapsamında yasak |
| `kurtarilabilir` | derinlik yeterli, tek baz açık çoğunlukta | dizi değiştirilmiyor, yalnızca not ediliyor |
| `belirsiz` | yukarıdakilerin dışı | yasak, ayrıca incelenecek |

Eşikler veriden türetiliyor. Derinlik eşiği hedef başına medyan derinliğin
medyanının yüzde 10'u, tabanı 5; ikinci baz oranı eşiği samtools'un kendi
öntanımlı `--het-fract` değeri olan 0,15. İkisi de `--min-depth` ve
`--het-fract` ile geçersiz kılınabilir ve kullanılan değer ile türetme biçimi
çıktıya yazılıyor.

Nanopore veriniz R10.4 ise `--config r10.4_sup` ya da `--config r10.4_dup`
eklemek kalite kalibrasyonunu düzeltir. Basecalling modelinizi bilmediğim için
öntanımlı olarak eklemedim; betik istenen değerin bu samtools sürümünde tanımlı
olup olmadığını kontrol ediyor.

Çıktılar: `N_pozisyonlari.tsv` (her N pozisyonu ayrıntılı), `hedef_ozeti.tsv`
(hedef başına özet), `maske/*.bed` (primer yerleşimi için yasak bölgeler),
`ambig/*.fa` ve `pileup/*.txt` (ham ölçümler).

### 4. Bracken yeniden tahmin

Henüz kesinleşmedi, iki bilgi gerekiyor: yeni veritabanında hangi
`database<K>mers.kmer_distrib` dosyalarının bulunduğu ve ölçülen medyan okuma
uzunluğu. Ali `database300mers.kmer_distrib` kullanmış, ancak `kraken results`
içindeki okuma uzunlukları 4 kb dolayında görünüyor, yani uygun uzunlukta bir
dağılım dosyasının `bracken-build` ile üretilmesi gerekebilir. `00` betiği her
ikisini de ölçüyor.

Ali'nin `bracken_species.sh` betiğinde iki hata var, düzeltilmiş sürümü
üreteceğim: yolların başındaki bölü işareti eksik ve `/"$VAR"` ile yamanmış,
ayrıca `basename ... _report.txt` kalıbı `.report` uzantılı dosyalarla
eşleşmiyor.

## Temizlik yapıldı

`_to_delete\` klasörüne 26 GB taşındı. Bağlı klasörlerde silme yetkisi
olmadığı için taşıma yöntemi kullanıldı; yer boşalması için o klasörü Windows
tarafından silmeniz gerekiyor. İçerik:

- `REFERANS_DB\` altında yarım kalmış UNITE mfeprimer indeks çalışma alanı
  (24 GB) ve onun 412 MB'lık logu. Logun son 1 MB'ı
  `REFERANS_DB\UNITE_ITS.fasta.log.son1MB.txt` olarak yerinde bırakıldı.
- md5 ile doğrulanmış SILVA kopyaları: `SILVA_SSURef_NR99.fasta` ve `.gz`,
  `SILVA_LSURef_NR99.fasta` ve `.gz`. `SILVA_138.2_*` sürümleri korundu.
- `fastq_F2-1_yanlis_kopyalar\` altında F2-1 klasörüne yanlışlıkla kopyalanmış
  beş F1-4 fastq dosyası.
- `diger\` altında `ARACLAR\linux-x64` (Bun çalışma zamanı, primer aracı değil),
  bracken A2 `report copy.txt`, `A1_1_methanosarcina_mazei.fasta`.

## mfeprimer

Windows `mfeprimer.exe` dosyasını kullanmayın. O 2025-09-29 derlemesi, yani 4.2
kuşağı, ve `REFERANS_DB` içindeki `v5-single-strand-64` indeks biçimini
okuyamaz. `ARACLAR\mfeprimer` (Linux ELF, modül sürümü
`v0.0.0-20260717054035`) 4.4 kuşağından ve o indeksleri üreten ikili. Temiz bir
etiketli sürüm isterseniz:

```bash
wget https://github.com/quwubin/MFEprimer-3.0/releases/download/v4.4.0/mfeprimer-4.4.0-linux-amd64.gz
gunzip mfeprimer-4.4.0-linux-amd64.gz && chmod +x mfeprimer-4.4.0-linux-amd64
./mfeprimer-4.4.0-linux-amd64 --version
```

v4.4.0'ın getirdiği `--bind-amp-only` bayrağı özgüllük süzgecinde işimize
yarayacak: rakip veritabanlarında genom çapındaki bütün k-mer eşleşmelerine
termodinamik hesap yapmak yerine yalnızca öngörülen amplikon içindeki
bağlanmaları raporluyor. 16 GB RAM ile büyük referanslarda bu fark önemli.

## Belirsiz nükleotid durumu, özet

Doksan dokuz konsensüs dosyasının tamamı kaynak diskten tarandı. Alfabe yalnızca
A (56.064), C (66.639), G (59.428), T (60.139) ve N (1.897); R, Y, S, W, K, M, B,
D, H, V kodlarından hiçbiri yok, boşluk da yok. Toplam 244.167 bp içinde binde
7,8 oranında N var. Doksan sekiz dosyanın 19'u tamamen temiz.

Dağılım bunun ağırlıklı olarak kapsama sorunu olduğunu gösteriyor. En kötü dosya
`A2-3-reads_118126`, 4381 bp içinde 137 N ve bunlar 118 ayrı blok halinde, en
uzun blok 3 bp, ve o taksonun fastq dosyası yalnızca 376 KB. Tek yapısal istisna
`F1-1-reads_2093779`, orada tek parça 43 bp'lik bir N bloğu var.

N artı üç baz güvenlik payı maskelendikten sonra 98 dosyanın tamamında en az bir
25 bp temiz pencere kalıyor, yalnızca iki dosyada (`B-4-reads_1642646` ve
`B-4-reads_1642647`) en uzun temiz blok 145 bp'de kalıyor. Sizin kurallarınıza
göre yalnızca iki primer ayak izinin temiz olması gerektiği, ürünün iç kısmındaki
N'in kabul edilebilir olduğu için gerçek kısıt bundan çok daha gevşek.

---

# Ek: 03_primer_aday_uret.py, toplantı kararlarının uygulanması

Bu betik toplantı kararlarındaki oligo, termodinamik, ürün ve bölge kurallarını
doğrudan uyguluyor. Sınıflandırmayı ve N analizini beklemez; maske dosyası
verilmezse tüm konsensüs üzerinde çalışır, verilirse maskeli pozisyonlar primer
ayak izine hiç girmez.

```bash
python3 03_primer_aday_uret.py \
    --consensus "consensus sequences/A1-1/A1-1-reads_2209_consensus_strict.fasta" \
    --mask      "N_analizi/maske/A1-1_2209_maske.bed" \
    --out       "primer_adaylari/A1-1_2209.tsv"
```

## İki bağımsız Tm ölçümü, gerçek veriyle sınandı

Kural şuydu: sıcaklık iki bağımsız kütüphaneyle hesaplanır ve aralarındaki
sabit kaymadan tolerans kadar sapan oligo elenir. Betik kaymayı elle almıyor,
o hedefteki bütün adayların farkının medyanından ölçüyor. Gerçek
*Methanosarcina mazei* konsensüsünde (A1-1, taxid 2209, 1445 bp) 5.782 aday
üzerinde ölçüm:

| Ölçüm | Değer |
|---|---|
| primer3 eksi Biopython, medyan kayma | eksi 1,50 derece |
| kaymanın standart sapması | 0,44 derece |
| tolerans dışında kalan oligo | 0 |

Kayma sistematik ve dar, yani kuralın dayandığı varsayım verinizde doğrulanıyor.
Tolerans öntanımlı 2,0 derece, `--tm-cross-tol` ile değiştirilebilir.

## Uygulanan kurallar

Oligo: uzunluk 18 ile 25, GC yüzde 40 ile 60 tercih ve 35 ile 65 sert sınır,
3' uç G ya da C, son beş bazda en fazla üç G veya C, aynı bazın en fazla dört
tekrarı. Termodinamik: Tm 58 ile 62 tercih ve 57 ile 63 sert sınır, hairpin dG
en az eksi 3000, self-dimer en az eksi 6000, çift içi hetero-dimer en az eksi
6000. Ürün: 70 ile 250 ve 300'e kadar, puanlamada 90 ile 150 en iyi, ürün GC 40
ile 60 tercih, iki primerin Tm farkı 1,5 dereceden az. Her eşik komut satırından
değiştirilebilir ve kullanılan değer çıktının başına yazılır.

Ürünün makine doğrulaması her adayda yapılıyor: kalıptan kesilen parçanın başı
ileri primere, sonu geri primerin ters tümleyenine birebir eşit mi diye
karşılaştırılıyor. Aynı kalıptan kesildiği için bu aşamada sıfır başarısızlık
bekleniyor; kontrolün asıl işlevi rakip konsensüslere ve ham okumalara karşı
çalıştırıldığında ortaya çıkacak.

## Dejenere baz politikası

Öntanımlı `--degeneracy-budget 0`, yani toplantı kararınıza uygun olarak
oligoda dejenere baz yok. Ancak Karar 3'teki işlev grubu primerleri ve Karar
4'teki anaerobik universal primer için tek bir ACGT dizisinin bütün üyeleri
yakalaması çoğu zaman mümkün olmuyor. Bu iki karar için sınırlı dejenerelik
açılabilir:

```bash
--degeneracy-budget 2 --degeneracy-fold-max 4
```

Betik bu durumda en fazla iki dejenere pozisyona ve toplam dört kata izin
veriyor, ve dejenere bazın son beş bazda bulunmasını her koşulda yasaklıyor
çünkü uzama oradan başlıyor. Tür ve cins özgül primerlerde (Karar 1 ve 2)
bütçenin sıfır kalmasını öneriyorum, orada özgüllük her şeyden önemli.

## Gerçek hedeflerde ölçülen sonuç

| Hedef | Uzunluk | N | Maskeli oran | Geçerli çift |
|---|---|---|---|---|
| A1-1 taxid 2209, *M. mazei* | 1445 bp | 0 | yüzde 0 | 56.454 |
| A2-3 taxid 118126, en kötü N yoğunluğu | 4381 bp | 137 | yüzde 3,13 | 101.245 |
| B-4 taxid 1642646, en kısa temiz blok | 1493 bp | 65 | yüzde 4,35 | 20.197 |

En kötü durumdaki hedeflerde bile on binlerce geçerli çift kalıyor, yani N'ler
tasarımı hiçbir hedefte tıkamıyor. Darboğaz özgüllük aşamasında olacak. Çıktı
ceza puanına göre sıralanıyor ve öntanımlı olarak ilk 5.000 satıra kırpılıyor,
`--max-pairs` ile değiştirilebilir.

Bir hedefin işlenmesi 1,5 kb için yaklaşık 14 saniye, 4,4 kb için yaklaşık bir
dakika sürüyor; 98 hedef için toplam yarım saatin altında.

---

# Ek: 04_grup_primer.py, Karar 2, 3 ve 4 için çok üyeli motor

Karar 2 (cins özgül), Karar 3 (işlev grubu) ve Karar 4 (universal) aynı motoru
kullanır; aradaki tek fark hedef kümesinin büyüklüğü ve dejenerelik bütçesidir.
Hizalama kullanılmaz. Adaylar bir çapa konsensüsten üretilir, sonra her aday
toplantı kararındaki bağlanma kuralıyla her üyeye ve her rakibe karşı taranır.

Kabul ölçütleri toplantı kararındaki gibi: hedeflenen her üyede ürün oluşmalı,
rakiplerin hiçbirinde oluşmamalı, ve ayrım sağlam olmalı, yani primerlerden en
az biri rakiplerde hiç bağlanma yeri bulamamalı. İki primerin de zayıfça
bağlanıp yalnızca birlikte yetersiz kalmasıyla oluşan temizlik reddedilir.

## Konsensüs dosyalarında yön normalizasyonu yok, motor bunu düzeltiyor

`consensus2.sh` her takson için bir çekirdek okuma seçip konsensüsü o okumanın
yönünde kuruyor. Sonuç olarak konsensüsler birbirine göre yön normalizasyonu
yapılmamış durumda. Ölçtüm: asetoklastik metanojen kümesindeki beş üyeden biri,
`A1-4-reads_3078083`, diğerlerine göre ters saklanmış. Kırk korunmuş probdan
26'sı eksi zincire bağlanıyor, artı zincire hiçbiri bağlanmıyor.

Bunun etkisi sessiz ve tam yıkıcı: ileri primer o üyede artı zincir yerine eksi
zincire bağlandığı için o üyede hiç ürün oluşmuyor, dolayısıyla bütün çiftler
"bir üyede ürün yok" diye eleniyor. Normalizasyon eklenmeden önce asetoklastik
grup için sıfır geçerli çift çıkıyordu; eklendikten sonra 4.755 çift çıktı.

Motor her koşuda çapadan 40 korunmuş prob üretip her diziyi oyluyor ve ters
olanları ters tümleyene çeviriyor, hangi dizinin çevrildiğini log'a yazıyor.
Kaynak dosyalara dokunulmuyor.

## Doğrulama: motor bilinen bir universal primeri kendi başına buldu

Motorun bulduğu korunmuş oligolardan biri `GGTTACCTTGTTACGACTTA`. Bu, klasik
1492R universal 16S primerinin çekirdek bölgesi. Motor bu diziyi hiçbir referans
listesine bakmadan, yalnızca beş üyenin konsensüsündeki korunmuşluktan buldu.
Bağlanma taraması ayrıca kaba kuvvet karşılaştırmasıyla sınandı ve üç örnek
oligoda her üyede aynı uyumsuzluk sayısını verdi.

## Karar 3 örneği: asetoklastik metanojenler

Hedef kümesi taxid 2209, 2223, 2208, 3078083, 1434102. Rakipler taxid 394967,
83984, 224719, 1406512, yani hidrojenotrofik ve metilotrofik metanojenler.

```
2.041 oligo x 9 dizi tarandı
her hedef üyeye bağlanan oligo   : 897
rakiplerde hiç bağlanmayan oligo : 188
elenen, bir üyede ürün yok       : 77.746
elenen, rakipte ürün oluşuyor    : 6.044
elenen, yetim primer yok         : 1.064
elenen, hetero-dimer dG          : 186
geçerli çift                     : 4.755
```

En iyi aday:

```
F  ATCTCCGGGCTCTTGCTCTC   Tm 61,4
R  TGGGTCTGCGGCCTATCAG    Tm 61,4
   ürün beş üyenin tamamında 88 bp, ileri primer rakiplerde hiç bağlanmıyor
```

Beş üyenin hepsinde aynı ürün uzunluğu çıkması bu bölgenin gerçekten korunmuş
olduğunu gösteriyor. Bir hedef grubu üç saniyede işleniyor.

---

# Ek: Bağımsız çapraz denetim ve düzeltmeler

İki bağımsız denetçi `03` ve `04` betiklerini kural kural sınadı, sentetik ve
gerçek veriyle test yazdı. Bulunan kusurların tamamı düzeltildi ve düzeltme
sonrası koşu tekrarlandı.

## 04'te bulunan ağır kusurlar

| Kusur | Etkisi | Düzeltme |
|---|---|---|
| `product_len` yalnızca "F artı zincirde, R eksi zincirde" konfigürasyonunu tarıyordu | Kalıp çift sarmal olduğu için ters konfigürasyon da gerçek üründür. Rakip yanlışlıkla "temiz", üye yanlışlıkla "ürün yok" sayılıyordu | İki konfigürasyon da taranıyor |
| Yön oylaması berabere kalınca (artı=eksi=0) sessizce "ters değil" diyordu; gerçek veride 21 dizide oluyor, 3'ü gerçekten ters | Ters bir rakip özgüllük denetiminden tamamen muaf kalıyordu | İkinci bağımsız ölçüt eklendi: bütün canlılarda korunmuş SSU motifleri. İkisi de karar veremezse dizi adı yazılarak uyarılıyor |
| `find_bindings` kalıptaki N'i seed taramasında göremiyordu | 98 konsensüste 3' uç hizalama pozisyonlarının yüzde 3'ü kör, en kötü dosyada yüzde 21,8. Yanlış "yetim primer" ilanları üretiyordu | Seed varyantları artık N de üretiyor |
| 5' sarkma serbest kuralı hiç uygulanmamıştı | Her dizinin iki ucunda yaklaşık 24 bp kör bant; tam da kapsamanın düştüğü yer | 5' sarkma serbest, örtüşen kısım en az `--min-overlap` baz |
| Rakipteki ürün yalnızca hedef uzunluk penceresindeyse reddediliyordu | 370 bp'lik bir rakip bandı "temiz" sayılıyordu | Rakipte artık her bant reddediliyor, `--competitor-prod-max 0` sınırsız |
| Aynı etikete düşen iki üye sessizce tek diziye çöküyordu | Çıktı iki üyeyi kapsadığını iddia ediyordu | Etiket çakışmasında betik duruyor |
| Rakip verilmediğinde her oligo "yetim" ilan ediliyordu | Asılsız özgüllük güvencesi | `rakip_verilmedi` yazılıyor ve uyarı basılıyor |
| Maske araması hem fazla geniş hem sessizce boş | Başka taksonların koordinatları çapaya bindiriliyor, ya da maske sessizce yok sayılıyordu | Grup ve taxid birlikte eşleştiriliyor, bulunamazsa duruyor |
| `--degeneracy-budget 2` betiği çökertiyordu | N meşru dejenere baz sayılıyor, primer3 patlıyordu | N `DEGEN_FOLD`'dan çıkarıldı, N içeren oligo hiç üretilmiyor |
| "En kötü uyumsuzluk" `min` ile hesaplanıyordu | Ad ile hesap çelişiyordu | Üye içinde en iyi, üyeler arasında en kötü olarak netleştirildi |

## 03'te bulunan kusurlar

`--degeneracy-budget > 0` çökertiyordu (aynı N sorunu, düzeltildi). `--prod-max`
ölü parametreydi, artık yumuşak üst sınır olarak cezalandırılıyor. F ve R ayak
izlerinin çakışmasına karşı koruma yoktu, eklendi. Var olmayan maske yolu
sessizce maskelemeyi kapatıyordu, artık hata veriyor. Çok kayıtlı FASTA sessizce
birleştiriliyordu, artık duruyor. BED'in kontig sütunu yok sayılıyordu,
`--mask-contig` eklendi. İki kütüphane toleransı 2,0 derecede sabitti ve gözlenen
dağılımın 5,5 standart sapması olduğu için hiçbir oligoyu elemiyordu; artık
tolerans veriden türetiliyor (`--tm-cross-k` çarpı standart sapma). Çift Tm farkı
kuralı 1,5'i kabul ediyordu, artık kesin olarak altında olmalı. `--gc-clamp-last 0`
kuralı kapatmak yerine tüm oligoya uyguluyordu. Aynı lokusun 1-2 bp kaymış
kopyaları çıktıyı dolduruyordu, `--min-locus-spacing` eklendi.

Ürün makine doğrulaması bu aşamada tanım gereği geçiyor, çünkü ürün primerlerin
türetildiği kalıbın aynı kopyasından kesiliyor. Betik bunu artık açıkça yazıyor
ve sütun adı `evet_ayni_kalip` oldu. Kontrolün asıl işlevi 05 özgüllük
aşamasında ortaya çıkacak.

## Yön denetimi: sorun sandığımdan çok daha yaygın

Bütün canlılarda korunmuş beş SSU motifinin genişletilmiş varyantlarıyla 98
konsensüsün tamamı bağımsız olarak denetlendi:

| Grup | Düz | Ters | Motif yok |
|---|---|---|---|
| A1 | 3 | 17 | 0 |
| A2 | 5 | 15 | 0 |
| B | 6 | 14 | 0 |
| F1 | 0 | 0 | 19 |
| F2 | 7 | 12 | 0 |
| **Toplam** | **21** | **58** | **19** |

Motifi saptanabilen 79 dizinin 58'i ters yönde saklanmış. Yön rastgele, çünkü
`consensus2.sh` her takson için bir çekirdek okuma seçip konsensüsü o okumanın
yönünde kuruyor ve nanopore okumaları iki yönde de geliyor. F1 grubunda SSU
motifi bulunmuyor, çünkü o amplikon ITS ve 18S bölgesi; onların yönü çapa
problarıyla belirleniyor.

Bunun anlamı: yön normalizasyonu küçük bir düzeltme değil, dosyaların yüzde
59'unu etkileyen zorunlu bir adım. Normalizasyon olmadan yapılan her diziler
arası karşılaştırma yanlıştır.

## Alan başına universal sonuçları

| Alan | Üye | Rakip | Korunmuş oligo | Geçerli çift | Ürün |
|---|---|---|---|---|---|
| Bakteri (B) | 20 | 40 | 137 | 141 | 129-140 bp |
| Arke (A1 artı A2) | 40 | 20 | 112 | 268 | 76-87 bp |
| Mantar uzun operon (F2) | 19 | 15 | 975 | 6.919 | 145-150 bp |
| Mantar ITS/18S (F1) | 19 | 15 | 28 | **0** | yok |

En iyi adaylar:

```
Bakteri universal
F  ACGACAGCCATGCAGCAC     Tm 61,0
R  ACAAGCGGTGGAGCATGTG    Tm 61,0     ürün 129-140 bp, geri primer yetim

Arke universal
F  CGGCGTTGAGTCCAATTAAAC  Tm 58,1
R  CGCAAGGCTGAAACTTAAAGG  Tm 58,1     ürün 79 bp, ileri primer yetim

Mantar universal (uzun operon)
F  TAAGAACGGCCATGCACCAC   Tm 61,0
R  AATTGACGGAAGGGCACCAC   Tm 60,9     ürün 145-147 bp, geri primer yetim
```

F1 grubunda çift bulunamamasının tek sebebi `F1-1-reads_1159556`
(*Ustilaginoidea virens*): 28 korunmuş oligonun ürün verebileceği bütün
kombinasyonları tek başına engelliyor (93 elemenin 93'ü). Bu grup ITS bölgesi
olduğu için türler arası ayrışma zaten yüksek; ya o üye kümeden çıkarılmalı ya da
F1 için tek bir universal yerine iki ayrı primer seti tasarlanmalı.

Arke tarafında `A1-1-reads_2209` en kısıtlayıcı üye (1.163 elemenin 1.163'ü);
bakteri tarafında `B-1-reads_1129264` (*Sphaerochaeta associata*, bir spiroket)
1.878 elemeden sorumlu. İkisi de beklenen davranış: en uzak akraba, korunmuş
pencereyi daraltan üyedir.

## Hâlâ yapılması gereken

Bu adayların kabul ölçütü sizin 98 konsensüsünüzde çalışması değil. Universal
primerlerin gerçek kapsamı SILVA 138.2 SSU'ya karşı ölçülmeli; `05` özgüllük
betiği bunu yapacak. Aynı betik `mfeprimer` ile `--bind-amp-only` ve `blastn`
ile REFERANS_DB taramasını, rakip oranının Wilson alt sınırıyla
değerlendirilmesini ve ham okumalarda doğrulamayı da içerecek.

---

# Ek: F1 için iki set, ve bu sırada çıkan daha büyük sorun

Önce adlandırma: konsensüs klasörlerindeki `F1-*` Ali'nin adlandırmasıdır ve
KISA amplikonu (1155-1649 bp) tutar. Ekran görüntünüze göre kısa fungus primeri
F2'dir ve barcode09-12'ye karşılık gelir. Aşağıda "F1" derken konsensüs
klasörünün adını kullanıyorum, yani kısa fungus amplikonunu.

## Bölme tasarlanabilirliğe göre yapıldı, benzerliğe göre değil

`05_kume_bol.py` betiği kümeyi 04 ile deniyor, geçerli çift çıkmazsa 04'un
raporladığı en çok engelleyen üyeyi çıkarıp tekrar deniyor, ilk geçerli çift
çıktığında kalanları bir set sayıp çıkarılanları yeni kümeye alıyor.

K-mer benzerliğine göre kümelemeyi denedim ve bıraktım: F1 grubunda ikili
Jaccard benzerliği 0,028 ile 0,035 arasında kalıyor, yani ölçüm gürültüden
ibaret. ITS bölgesinde bu beklenen bir durum.

## Sonuç

| Set | Üye | Geçerli çift | En iyi aday |
|---|---|---|---|
| SET1 | 16 | 12 | F `GTACTTGTTCGCTATCGGTCTC` (Tm 58,6) R `CGAGTCGGGTTGTTTGGG` (Tm 58,4), ürün 90-91 bp |
| SET2 | 3 | 844 | F `TGTGCGTTCAAAGATTTGATGATTC` (Tm 59,1) R `TTGGTTCTCGCAACGATGAAG` (Tm 59,2), ürün 90 bp |

SET2 üyeleri: `F1-1-reads_1159556`, `F1-1-reads_2093779`, `F1-1-reads_2093780`.
Hepsi F1-1 barkodundan. Bölme taksonomik değil barkod temelli çıktı ve bu
şüpheli olduğu için araştırdım.

## Asıl sorun: aynı taksonun konsensüsleri yıllar arasında farklı bölgeleri kapsıyor

Aynı taksonun farklı barkodlardaki konsensüslerini 12-mer benzerliğiyle
karşılaştırdım:

| Takson | Karşılaştırma | Benzerlik | Yorum |
|---|---|---|---|
| 44689 *Dictyostelium* | F1-2 vs F1-3 | 0,987 | aynı bölge |
| 44689 | F1-1 vs F1-2 | 0,065 | farklı bölge |
| 44689 | F1-2 vs F1-4 | 0,105 | farklı bölge |
| 2093779 *Podospora* | F1-1 vs F1-4 | 0,750 (ters) | aynı bölge, ters yön |
| 2093779 | F1-1 vs F1-3 | 0,056 | farklı bölge |
| 2093780 *Podospora* | F1-1 vs F1-4 | 0,049 | farklı bölge |
| 2545709 *S. osmophilus* | F1-2 vs F1-3 | 0,061 | farklı bölge |
| 101201 *Trichoderma* | F1-2 vs F1-4 | 0,922 | aynı bölge |

Yani tek bir taksonun dört yıla ait konsensüsü, rDNA'nın üç ayrı parçasını
kapsayabiliyor. Bu suş farkı değil, kapsanan pencerenin farkı.

Sebebi `consensus2.sh`: her takson için okumalardan bir çekirdek okuma seçiyor
ve konsensüsü onun etrafında kuruyor. Okuma uzunlukları ve kapsama değiştiğinde
çekirdek okuma operonun başka bir parçasına düşüyor, konsensüs de oraya
oturuyor. Yani F1 konsensüsleri tanımlı bir amplikon değil, keyfî bir çekirdek
okumanın etrafındaki yerel birleştirmedir.

Sonucu iki tane. Birincisi, F1 için yaptığım iki set aslında taksonomik bir
bölme değil BÖLGE bölmesidir; SET2 üç F1-1 dosyasını topluyor çünkü onlar başka
bir pencereye oturuyor. İkincisi ve daha önemlisi, F1 grubunda yıllar arası
karşılaştırma bu haliyle geçersizdir.

Bu sorun yalnızca F1'e özgü. A1 grubunda beş üyenin tamamı aynı 88 bp ürünü
veriyor, bakteri grubunda 20 üyede 129-140 bp, F2 grubunda 19 üyede 145-150 bp.
Onlar tutarlı.

## Önerilen çözüm

F1 konsensüsleri keyfî çekirdek okuma yerine ORTAK BİR REFERANSA göre yeniden
kurulmalı. Elinizde `UNITE_ITS.fasta` ve `fungi.18SrRNA.fna` BLAST indeksli
duruyor. Yöntem: her taksonun okumalarını o taksonun UNITE ya da 18S referans
dizisine `minimap2` ile hizala, konsensüsü referans koordinat sisteminde çağır.
Böylece dört yılın konsensüsü aynı pencereye oturur, hem yıllar arası
karşılaştırma hem tek bir universal primer mümkün olur.

Bu yapılana kadar yukarıdaki iki set kullanılabilir, ancak SET2'nin biyolojik
bir grup olmadığını, yalnızca farklı bir rDNA penceresine oturan üç dosya
olduğunu akılda tutmak gerekir.

---

# Ek: Klasör adları düzeltildi ve bir çıkarımım çürütüldü

## Klasör ile barkod eşlemesi bağımsız olarak kanıtlandı

Klasör adlarına hiç güvenmeden, her konsensüs ve fastq klasörünün takson
kümesini o barkodun Bracken çıktısındaki ilk beş taksonla eşleştirdim.
`sequence_extraction.sh` zaten ilk beşi çıkardığı için bu eşleşme tam olmalı ve
öyle oldu: 20 klasörün 20'sinde de kesişim 5/5 ya da 4/4, ikinci en iyi aday
her zaman geride. Konsensüs ve fastq eşlemeleri birbirinin aynısı.

| Klasör | Barkod | Klasör | Barkod |
|---|---|---|---|
| A1-1..A1-4 | barcode01-04 | F1-1..F1-4 | barcode13-16 |
| A2-1..A2-4 | barcode05-08 | F2-1..F2-4 | barcode09-12 |
| B-1..B-4 | barcode17-20 | | |

Bu, ekran görüntüsündeki grup ile barkod eşlemesiyle birebir aynı. Yani
konsensüs ve fastq klasörleri baştan beri doğruydu.

## Ne yeniden adlandırıldı

Yalnızca `bracken results\genus`, `bracken results\species` ve `kraken results`
klasörlerindeki F1 ve F2 alt klasörleri takas edildi; onlar barkodlarla ters
eşleşiyordu. Şimdi hepsi uyumlu:

```
bracken results/genus/F1   -> barcode13 14 15 16
bracken results/genus/F2   -> barcode09 10 11 12
bracken results/species/F1 -> barcode13 14 15 16
bracken results/species/F2 -> barcode09 10 11 12
kraken results/F1          -> barcode13 14 15 16
kraken results/F2          -> barcode09 10 11 12
```

Manifest ve geri alma betiği `_yeniden_adlandirma\` klasöründe.

## Çürütülen çıkarım: Uzun ve Kısa etiketleri

Daha önce konsensüs uzunluklarına bakıp "ekran görüntüsü doğru, Ali'nin klasör
adları ters" demiştim. Yanlıştı. Kraken çıktılarından okuma uzunluklarını
doğrudan ölçünce şu çıktı:

| Barkod | Medyan okuma | Grup | Ekran görüntüsü diyor ki | Ölçüm diyor ki |
|---|---|---|---|---|
| barcode01-04 | 1419-1428 bp | A1 | Kısa | Kısa, doğru |
| barcode05-08 | 4300-4324 bp | A2 | Uzun | Uzun, doğru |
| barcode09-12 | 3697-3700 bp | F2 | Kısa | **Uzun** |
| barcode13-16 | 1204-1479 bp | F1 | Uzun | **Kısa** |
| barcode17-20 | 1478-1483 bp | B | Bakteri paneli | 1,5 kb |

Ekran görüntüsündeki grup ile barkod eşlemesi doğru, ama fungus grupları için
Uzun ve Kısa açıklamaları ters yazılmış. Arke gruplarında açıklama doğru.

Excel'deki grup açıklamaları buna göre düzeltildi: F1 artık Fungus Kısa Primer,
F2 artık Fungus Uzun Primer. Barkod, yıl ve renk atamaları değişmedi çünkü
onlar zaten doğruydu; yalnızca açıklama sütunu düzeltildi. Doğrulama ve Log
sayfasına ölçüm değerleri yazıldı.

Bu düzeltmenin daha önceki F1 bulgularına etkisi yok, yalnızca adlandırma
değişti: "F1 konsensüsleri farklı bölgeleri kapsıyor" bulgusu barcode13-16 yani
kısa fungus amplikonu için geçerli.

# Ek: 06_referans_capali_konsensus.sh

Her taksonun konsensüsünü keyfî bir çekirdek okuma yerine ortak bir referansın
koordinat sisteminde kurar.

```bash
bash 06_referans_capali_konsensus.sh \
     --pt  "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama" \
     --out "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama/referans_konsensus" \
     --groups F1,F2
```

Takson başına akış: okumalardan bir örneklem alınır, uygun BLAST
veritabanına karşı `blastn` ile en yüksek toplam bit skorunu veren referans
seçilir, referans `blastdbcmd` ya da `samtools faidx` ile çıkarılır, SILVA gibi
RNA alfabeli kaynaklarda U harfleri T'ye çevrilir, bütün okumalar `minimap2` ile
o referansa hizalanır, sonra `samtools consensus` iki bağımsız ölçüm üretir:
`-a -A` ile IUPAC kodlu konsensüs, `-a -A -f PILEUP` ile pozisyon başına
derinlik ve baz dizisi.

Veritabanı seçimi gruba göre: A grupları `archaea.16S.fna`, B grubu
`bacteria.16S.fna`, F grupları `fungi.ITS.fna` artı `fungi.18SrRNA.fna` artı
`fungi.28SrRNA.fna` arasından en iyi skoru veren.

Betik sonunda hizalılık denetimi yapıyor: aynı taksonun farklı barkodlarda aynı
referansa oturup oturmadığını kontrol edip oturmayanları listeliyor. Oturmayan
varsa o taksonda referansı elle sabitlemek gerekir, aksi halde yıllar yine aynı
koordinat sistemine gelmez.

Çıktı referans koordinatında olduğu için aynı taksonun bütün yılları hizalı olur
ve yön normalizasyonu da gereksizleşir; 04 yine de kontrol eder.

Sıra: `06` ile konsensüsleri yeniden kur, `02_N_analizi.sh` ile maske üret,
sonra `04_grup_primer.py` ve `05_kume_bol.py` ile tasarıma dön.
