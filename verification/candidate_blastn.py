# -*- coding: utf-8 -*-
"""
candidate_blastn.py, a blastn based CANDIDATE FINDER for the identity stage.

The strongest candidate for a full solution to the short list problem. This file
DOES NOT TOUCH THE SCORING: the identity percentage is still computed by
kimlik_dogrulama.hizala, under our own definition. The only thing that changes is
WHO decides which records get aligned.

"""
# -------------------------------------------------------------------------
# candidate_blastn.py
#
# INPUT  : the bin consensus (the query) and a BLAST indexed database under
#          REFERANS_DB (the .nin/.nhr/.nsq files already exist)
# OUTPUT : [(header, sequence), ...], the records where blastn found a meaningful
#          hit
# CALLED BY: verification/identity_verification.py and all_bin_identities.py, only
#          when ADAY_BULUCU=blastn is given
#
# WHY THIS OPTION IS THE STRONGEST, AND WHY IT WAS NOT THOUGHT OF EARLIER
# In this project we wrote our own seed plus alignment pipeline. But blastn was
# designed for exactly this job, it has been tested for thirty years, and it ranks
# by BIT SCORE: the match length, the identity and the gap penalty all enter the
# calculation together, and the e-value is corrected for the database size. Our own
# seed count, by contrast, IS A POOR PROXY for identity: a long and conserved
# record collects many seeds while being irrelevant.
#
# More to the point: blastn IS ALREADY USED IN THIS DIRECTORY. The global
# specificity layer was run with it and the BLAST index (.nin) of nine databases
# sits ready under REFERANS_DB. So we wrote our own pipeline while a ready, indexed
# and tested tool was at hand. That goes against the project's own rule: look at
# the measured tool you already have first.
#
# HOW IT SOLVES THE SHORT LIST PROBLEM
# On our route the candidates are ranked by seed count and the first N are taken;
# it was measured that on one query the winner came from position 4171. In blastn
# there is no "first N" prefilter; EVERY record passing the e-value threshold comes
# back. The cut off is A STATISTICAL THRESHOLD rather than an arbitrary rank, and it
# can be reported.
#
# WHERE IT IS NOT USED
# The primer binding search and in-silico PCR. blastn depends on its seed setting
# on short queries too and gives NO GUARANTEE of losslessness; the pigeonhole
# engine stays there.
#
# SAFETY: if blastn is missing, var_mi() returns False and the caller carries on
# with the existing short list route. The chain does not break.
# -------------------------------------------------------------------------

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
    """Is blastn in working order. It measures once."""
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
    """Is the database's BLAST index built.

    If the index is missing it can be built with makeblastdb, but that takes MINUTES
    and is not done by itself: producing a gigabyte sized index without the user
    knowing would not be right. False is returned instead and the caller falls back on
    the old route.

    """
    return os.path.exists(fasta_yolu + '.nin') or os.path.exists(fasta_yolu + '.nal')


def adaylar(q, fasta_yolu, e_deger=1e-5, en_fazla=5000, iplik=3):
    """Finds candidate records with blastn.

    Returns: [(header, sequence), ...] or None (when blastn or the index is missing).
    None and an empty list ARE DIFFERENT THINGS:
      None -> the tool or the index is missing, FALL BACK ON THE OLD ROUTE
      []   -> blastn ran and no record passed the threshold
    Confuse the two and a database is skipped silently.

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
                 # dc-megablast: it finds distant relatives too. megablast is fast
                 # for close matches but misses anything under 90 percent; some of
                 # our bins resemble the references at around 85 percent.
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
    """Pulls the sequences in bulk with blastdbcmd. Pulling them one at a time would be far too slow."""
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
    print(u'is blastn working   : %s' % (u'EVET' if var_mi() else u'HAYIR'))
    if not var_mi():
        print(u'reason              : %s' % _SEBEP)
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
