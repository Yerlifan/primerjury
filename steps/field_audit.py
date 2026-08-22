#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
field_audit.py
ALAN (domain) TUTARLILIK DENETIMI, tek kaynak.

Sorun: bir hedefin taxid listesi, o hedefin ait olmadigi bir lokus
kitapliginda da gecebiliyor. Ornek: Sakarolitik bakteriler hedefi
bakteriyeldir (Bacteroides, Proteiniphilum, Sphaerochaeta,
Acetobacteroides) ve kutulari B sinifindadir; ama taxid 1760811
(Acetobacteroides) mantar ITS/28S kitapliginda da iki kutuda goruluyor.
Dogru sinif olan B'de hicbir cift bulunamayinca, geriye yalniz yanlis
lokustaki tek kutudan tasarlanmis F2 cifti kaliyor. Tablo o hedefi
kapsanmis gosteriyor, oysa kapsanmiyor. Kural denetimi bunu yakalamaz,
cunku ortada kural ihlali yoktur.

Cozum: hedefin uye kutularinin hangi alanlarda (A arke, B bakteri,
F mantar) bulundugu VERIDEN cikarilir. Hedef birden cok alanda kutu
tasiyorsa, azinlik alandaki tasarim isaretlenir.

Bu dosya elle yazilmis bir alan tablosu icermez. Alan bilgisi yalnizca
konsensus dosya adlarindan ve targets.tsv'deki taxid listesinden gelir.
"""
import os
import re

ETIKET = re.compile(r"((?:A1|A2|B|F1|F2))-\d+_(\d+)$")


def taxid_alanlari(kons_klasoru):
    """Konsensus dosya adlarindan taxid -> {alan: kutu_sayisi} cikarir."""
    d = {}
    if not kons_klasoru or not os.path.isdir(kons_klasoru):
        return d
    for ad in sorted(os.listdir(kons_klasoru)):
        if not ad.endswith(".fasta"):
            continue
        et = re.sub(r"_(baskin|ref|self)?_?konsensus\.fasta$", "", ad)
        m = ETIKET.match(et)
        if not m:
            continue
        sinif, taxid = m.group(1), m.group(2)
        d.setdefault(taxid, {})
        d[taxid][sinif[0]] = d[taxid].get(sinif[0], 0) + 1
    return d


def hedef_taxidleri(hedefler_tsv):
    """targets.tsv -> {hedef_adi: [taxid, ...]}"""
    d = {}
    if not hedefler_tsv or not os.path.exists(hedefler_tsv):
        return d
    with open(hedefler_tsv, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] == "karar":
                continue
            d[p[1]] = [t.strip() for t in p[3].split(",") if t.strip()]
    return d


def alan_dagilimi(hedef, sinif, kons_klasoru=None, hedefler_tsv=None,
                  taxid_alan=None, hedef_taxid=None):
    """Doner: (uyumsuz_mu, dagilim_sozlugu, baskin_alan)

    uyumsuz_mu True ise bu (hedef, sinif) ikilisi hedefin BASKIN alaninda
    degildir ve teslim edilmemelidir. Hedef tek alanda kutu tasiyorsa ya da
    veri eksikse (False, {}, None) doner: karar verilemiyorsa suclamayiz.
    """
    if taxid_alan is None:
        taxid_alan = taxid_alanlari(kons_klasoru)
    if hedef_taxid is None:
        hedef_taxid = hedef_taxidleri(hedefler_tsv)
    tl = hedef_taxid.get(hedef)
    if not tl or not taxid_alan:
        return False, {}, None
    dagilim = {}
    for t in tl:
        for alan, n in taxid_alan.get(t, {}).items():
            dagilim[alan] = dagilim.get(alan, 0) + n
    if len(dagilim) < 2:
        return False, dagilim, (max(dagilim, key=lambda x: dagilim[x])
                                if dagilim else None)
    baskin = max(dagilim, key=lambda x: dagilim[x])
    bu = (sinif or "")[:1]
    return (bool(bu) and bu != baskin), dagilim, baskin


def aciklama(dagilim, sinif, baskin):
    ad = {"A": "archaea", "B": "bacteria", "F": "fungi"}
    return ("the target's member bins are %s, and this pair is in the %s field "
            "while the dominant field is %s. The design was made from the reads "
            "of a locus library the target does not belong to, so it does not "
            "represent the target."
            % (", ".join("%s=%d" % (ad.get(k, k), v)
                         for k, v in sorted(dagilim.items())),
               ad.get((sinif or "")[:1], sinif), ad.get(baskin, baskin)))
