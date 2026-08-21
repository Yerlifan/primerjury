# -------------------------------------------------------------------------
# geometry_core.py - the qPCR GEOMETRY audit of the 42 primers of the panel's 21
#          targets; Tm, GC, hairpin and dimer are computed for every primer and the
#          rule violations are listed.
#
# INPUT  : it reads no file. The list of 21 targets x (F, R) to audit is embedded in
#          the file itself, in the P constant. From outside, only primer3 is needed.
# OUTPUT : it writes geo.json into the working directory (the raw measured values
#          and the violation list for every primer), and prints a line by line
#          PASSED / VIOLATION table to the screen.
# CALLED BY: screening/self_test.py runs this file AS A SEPARATE PROCESS in a
#          temporary directory, reads the geo.json it produces, and compares it
#          EXACTLY against the gc/tm/hp/hd/son5/uc values screening/geometry.py
#          gives for the same 42 primers. If even one value does not match, the
#          self test fails and the run does not start. That test runs at the head of
#          every measuring key: full_chain.py stages 1, 2, 3, 4, 5, 6, 7, 9 and key 8
#          directly. The file is not called DIRECTLY from the menu; it stands as the
#          reference, the gold standard.
#
# WHICH COPY RUNS: self_test.py walks the config.BETIK_YOLLARI order and takes the
# FIRST geometry_core.py it finds. In the source study that list held three
# directories and two of them carried a copy of this file, so editing one could leave
# the test silently measuring the other. In this repository BETIK_YOLLARI is a single
# entry, engine/, and there is exactly one geometry_core.py; the ambiguity is gone.
# Keep it that way.
# -------------------------------------------------------------------------
import primer3, json, sys
# The panel's 21 targets and each target's (forward, reverse) primer pair. The audit
# runs over this fixed list, so the result is the same on every run and it makes a
# gold standard comparable against screening/geometry.py.
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
# The thermodynamic conditions. A Tm is a property not of a sequence but of the
# sequence PLUS THE BUFFER: change the salt or dNTP concentration and the Tm changes
# too. These values represent the qPCR mix in use (QuantiNova SYBR Green); the same
# values are repeated in config.py in the P3 constant with the note "identical to
# geometry_core.py". If these numbers drift from the real reaction, the 58-62 Tm
# window below stops meaning anything.
#   mv_conc   = 50   mM monovalent cation (Na+/K+)
#   dv_conc   = 1.5  mM divalent cation (Mg++)
#   dntp_conc = 0.6  mM dNTP
#   dna_conc  = 50   nM oligo concentration
KW=dict(mv_conc=50, dv_conc=1.5, dntp_conc=0.6, dna_conc=50)
# Erime sicakligi (nearest-neighbor, primer3).
def tm(p): return round(primer3.calc_tm(p, **KW),2)
# GC yuzdesi.
def gc(p): return round(100.0*sum(p.count(c) for c in 'GC')/len(p),1)
# Hairpin: the melting temperature of the primer folding back ON ITSELF.
def hp(p): return round(primer3.calc_hairpin(p, **KW).tm,1)
# Homodimer: the melting temperature of the primer pairing with ITS OWN copy.
def hd(p): return round(primer3.calc_homodimer(p, **KW).tm,1)
# Heterodimer: the melting temperature of F and R pairing WITH ONE ANOTHER.
def het(a,b): return round(primer3.calc_heterodimer(a,b, **KW).tm,1)
# -------------------------------------------------------------------------
# viol - the list of rules a single primer violates. An empty list = PASSED.
# Why each rule exists, and what it corresponds to in qPCR, is written out below.
# -------------------------------------------------------------------------
def viol(p):
    v=[]
    # UZUNLUK 18-25 baz. 18'in altinda primer genomda/toplulukta rastgele de
    # bulunabilecek kadar kisalir, ozgulluk duser; 25'in ustunde Tm asiri
    # yukselir, katlanma ve yanlis baglanma sansi artar, sentez maliyeti buyur.
    if not (18<=len(p)<=25): v.append("uz %d"%len(p))
    g=gc(p)
    # GC 40-60%. A G/C pair has three hydrogen bonds, an A/T two. Too low a GC and the
    # duplex cannot hold at 60 C, so the yield drops; too high a GC and the primer also
    # sticks to off-target GC rich regions, while the template may not open fully at 95 C.
    if not (40<=g<=60): v.append("GC %%%.1f"%g)
    t=tm(p)
    # THE Tm WINDOW 58-62 C. The whole panel will run in ONE protocol, under the same
    # instrument program (a 60 C annealing target). If a primer's Tm falls below the
    # window, binding weakens in that well (a false negative); above it, the primer binds
    # off-target regions at 60 C too (a false positive). The window is kept narrow so
    # that all 21 targets can run in the same conditions.
    if not (58<=t<=62): v.append("Tm %.2f"%t)
    # THE HAIRPIN THRESHOLD 45 C. If a primer folds back on itself into a hairpin and
    # that structure is still standing at the annealing temperature (around 60 C), the
    # primer is busy with itself instead of binding the template; the yield drops and
    # the Cq is delayed. The threshold is set well below the annealing temperature; a
    # structure melting above 45 C does not count as safe. Under SYBR Green this
    # criterion is ELIMINATING.
    if hp(p)>=45: v.append("hairpin Tm %.1f"%hp(p))
    # THE HOMODIMER THRESHOLD 45 C. A primer pairing with its own copy produces a
    # primer-dimer. A primer-dimer consumes both the primer and the dNTPs, and because
    # SYBR Green stains ALL double stranded DNA it GIVES A SIGNAL even with no target
    # present; it appears on the melt curve as a second peak at a low Tm. In a probe
    # based chemistry this criterion might count as a warning; under SYBR it eliminates.
    if hd(p)>=45: v.append("homodimer Tm %.1f"%hd(p))
    # THE LAST BASE AT THE 3' END must be G or C (a "GC clamp"). The polymerase extends
    # the primer from its 3' end, and a G/C with its three hydrogen bonds at that end
    # keeps the primer seated firmly where extension starts. A primer with A/T at the 3'
    # end "breathes" at the tip and the yield becomes variable.
    if p[-1] not in 'GC': v.append("3' uc %s (G/C degil)"%p[-1])
    n5=sum(1 for c in p[-5:] if c in 'GC')
    # AT MOST 3 G/C IN THE LAST 5 BASES AT THE 3' END. The counterweight to the previous
    # rule. Filling the 3' end with too much G/C (overdoing the "GC clamp") makes the
    # primer stick firmly even to off-target regions that hold only its last few bases;
    # because the tip is seated, the polymerase starts extending and a non-specific
    # product comes out. So: enough grip at the tip (the previous rule), but no more.
    if n5>3: v.append("3' son 5 bazda %d G/C"%n5)
    # A DEGENERATE BASE is forbidden. A letter outside ACGT (N, R, Y and so on) means
    # the primer is synthesised not as one oligo but as a MIXTURE of oligos: the
    # effective concentration of each variant drops, the Tm stops being a single number,
    # and every thermodynamic calculation above loses its meaning.
    if any(c not in 'ACGT' for c in p): v.append("dejenere baz")
    return v

# -------------------------------------------------------------------------
# EVERYTHING BELOW THIS IS THE DEMONSTRATION SECTION - it runs only when the script
# is run DIRECTLY. The 2026-08-11 fix: this section used to run ON IMPORT as well.
# The consequences:
#   1) Every script importing geometry_core.py printed the 63 line geometry dump of
#      the P list above (the panel of 2 August) to the screen. That list holds the
#      OLD sequences of pairs such as Petrimonas and Bacteroidales, and anyone
#      reading the console or the log took them for CURRENT ones.
#   2) The import WROTE geo.json into the working directory. An import writing a
#      file is an unexpected side effect.
# self_test.py runs this file AS A SEPARATE PROCESS (cwd=a temporary directory) and
# reads geo.json; there __name__ == "__main__", so the behaviour IS UNCHANGED. Only
# "import geo" went quiet.
# -------------------------------------------------------------------------
if __name__ == "__main__":
    rows=[]
    for name,f,r in P:
        for lbl,p in (("Ileri",f),("Geri",r)):
            rows.append(dict(hedef=name,primer=lbl,dizi=p,uz=len(p),gc=gc(p),tm=tm(p),
                             hp=hp(p),hd=hd(p),uc=p[-1],gc5=sum(1 for c in p[-5:] if c in 'GC'),
                             ihlal=viol(p)))
        # ---- CIFT DUZEYI KURALLAR (tek primere bakarak gorunmezler) ----
        d=abs(tm(f)-tm(r)); h=het(f,r); pv=[]
        # dTm < 1.5 C. F and R work in the SAME well at the SAME annealing temperature. If
        # their Tm values diverge, the primer with the lower Tm binds weakly at that
        # temperature; one strand amplifies faster than the other (asymmetric
        # amplification), and the yield and the reproducibility of the Cq both suffer.
        if d>=1.5: pv.append("dTm %.2f"%d)
        # THE HETERODIMER THRESHOLD 45 C. If F and R stick to one another, the primer-dimer
        # that forms does the same damage as a homodimer but more often: two different
        # sequences meet. Under SYBR Green it is a direct source of false signal.
        if h>=45: pv.append("heterodimer Tm %.1f"%h)
        rows.append(dict(hedef=name,primer="CIFT",dizi="dTm=%.2f  het=%.1f"%(d,h),uz="",gc="",tm="",hp="",hd="",uc="",gc5="",ihlal=pv))
    # geo.json is written INTO THE WORKING DIRECTORY (there is no absolute path).
    # self_test.py runs this file in a temporary directory with cwd=td and reads the
    # output from there, which is why the relative path is deliberate; it does not
    # pollute the project tree.
    json.dump(rows,open('geo.json','w'))
    # The screen output. On a PAIR line, both the pair level violations and the
    # violations of that target's two primers are gathered on ONE line: the decision is
    # made per target, and if even one primer of a target fails, that target cannot go
    # to order.
    for r in rows:
        if r['primer']=='CIFT':
            pr = [x['ihlal'] for x in rows if x['hedef']==r['hedef'] and x['primer']!='CIFT']
            allv = r['ihlal'] + [i for s in pr for i in s]
            print("%-34s CIFT  %-24s -> %s" % (r['hedef'], r['dizi'], "IHLAL: "+"; ".join(allv) if allv else "GECTI"))
        else:
            print("%-34s %-5s %-26s uz%3s GC%5s Tm%6s hp%5s hd%6s 3'%s g5=%s  %s" % (
                r['hedef'],r['primer'],r['dizi'],r['uz'],r['gc'],r['tm'],r['hp'],r['hd'],r['uc'],r['gc5'],
                ("IHLAL: "+"; ".join(r['ihlal'])) if r['ihlal'] else "GECTI"))
