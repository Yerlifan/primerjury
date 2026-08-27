# -*- coding: utf-8 -*-
"""A CHECKPOINT LEDGER: it stores expensive measurements and rescues interrupted work.

WHY IT EXISTS
  The alignments in this pipeline take hours. If a run is cut short (a power cut,
  memory, a window closed by accident) starting again from the beginning is not
  acceptable. But a badly built cache is worse still: it hands back an old
  measurement as though it belonged to the new parameters, and nobody notices.

THREE PRINCIPLES
  1. THE RAW MEASUREMENT is stored, THE VERDICT is not.
     The hits blastn returned are kept; the decision "this species is such and
     such" is derived again on every run. So when a threshold or a naming rule
     changes, the result is brought up to date without repeating the scan, in
     seconds rather than hours.
  2. THE SIGNATURE carries the conditions the measurement was taken under.
     The database list, the blast parameters, a digest of the query file. When the
     signature changes that entry is invalid and is measured again. Measuring
     again rather than silently returning the old value is expensive, but correct.
  3. WRITING IS ATOMIC.
     It writes to a temporary file first and then moves it into place. If it is
     interrupted while writing, the old ledger stays intact; there is no such
     thing as half a JSON file.

TO USE IT
    from checkpoint import Defter
    d = Defter(os.path.join(OUT, 'checkpoint', 'three_way.json'),
               imza={'db': GRUP_DB, 'threshold': TUR_ESIGI, 'version': 'v1'})
    if d.var(anahtar):
        isabetler = d.al(anahtar)
    else:
        isabetler = pahali_olcum()
        d.yaz(anahtar, isabetler)      # goes to disk AFTER every measurement
    d.kapat()

  `anahtar` has to describe the measurement unambiguously (the bin plus the
  database, for instance).
"""
from __future__ import print_function
import hashlib
import io
import json
import os
import tempfile


def _cogul(n, tekil, cogul):
    """One measurement, two measurements. A report a person reads should not
    stumble over its own grammar."""
    return '%d %s' % (n, tekil if n == 1 else cogul)


def _imza_ozeti(imza):
    """Reduce the signature to an ordered, repeatable digest. sort_keys is
    essential: the digest must not change when the dictionary order changes,
    otherwise the ledger is counted invalid on every run and the cache is of no
    use at all."""
    ham = json.dumps(imza, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.md5(ham.encode('utf-8')).hexdigest()[:16]


class Defter(object):
    def __init__(self, yol, imza, yaz_sikligi=1, sessiz=False):
        self.yol = yol
        self.imza = _imza_ozeti(imza)
        self.yaz_sikligi = max(1, int(yaz_sikligi))
        self.sessiz = sessiz
        self.kayit = {}
        self.bekleyen = 0
        self.kurtarilan = 0
        self.yeni = 0
        klasor = os.path.dirname(os.path.abspath(yol))
        if klasor and not os.path.isdir(klasor):
            os.makedirs(klasor)
        self._oku()

    def _oku(self):
        if not os.path.exists(self.yol):
            return
        try:
            with io.open(self.yol, encoding='utf-8') as fh:
                d = json.load(fh)
        except (ValueError, IOError, OSError) as e:
            # A corrupt ledger is not swallowed in silence: what happened is said
            # out loud, and then the measurement is taken from scratch. Swallowing
            # it silently makes somebody rescan for hours and then ask "why did the
            # cache not work".
            self._not(u'the checkpoint could not be read (%s), it will be measured '
                      u'from scratch: %s' % (type(e).__name__, self.yol))
            return
        if d.get('imza') != self.imza:
            _n = len(d.get('kayit') or {})
            self._not(u'the checkpoint signature has changed (the parameters differ), '
                      u'%s counted invalid' % _cogul(_n, 'entry was', 'entries were'))
            return
        self.kayit = d.get('kayit') or {}
        self.kurtarilan = len(self.kayit)
        if self.kurtarilan:
            self._not(u'%s rescued from the checkpoint: %s'
                      % (_cogul(self.kurtarilan, 'measurement was',
                                'measurements were'),
                         os.path.basename(self.yol)))

    def _not(self, mesaj):
        if not self.sessiz:
            print(u'  [checkpoint] %s' % mesaj)

    def var(self, anahtar):
        return anahtar in self.kayit

    def al(self, anahtar, varsayilan=None):
        return self.kayit.get(anahtar, varsayilan)

    def yaz(self, anahtar, deger):
        self.kayit[anahtar] = deger
        self.yeni += 1
        self.bekleyen += 1
        if self.bekleyen >= self.yaz_sikligi:
            self.kaydet()

    def kaydet(self):
        """An atomic write: a temporary file plus os.replace. The temporary file is
        opened in the same directory, because os.replace is only atomic within one
        file system (a system temporary directory can be another partition)."""
        if not self.bekleyen:
            return
        klasor = os.path.dirname(os.path.abspath(self.yol))
        gec = None
        try:
            fd, gec = tempfile.mkstemp(dir=klasor, suffix='.tmp')
            with io.open(fd, 'w', encoding='utf-8', newline='\n') as fh:
                fh.write(json.dumps({'imza': self.imza, 'kayit': self.kayit},
                                    ensure_ascii=False, sort_keys=True))
            os.replace(gec, self.yol)
            gec = None
            self.bekleyen = 0
        except (IOError, OSError) as e:
            self._not(u'the checkpoint COULD NOT BE WRITTEN (%s), the measurements are '
                      u'in memory and will be lost if the run is cut short'
                      % type(e).__name__)
        finally:
            if gec and os.path.exists(gec):
                try:
                    os.remove(gec)
                except OSError:
                    pass

    def kapat(self):
        self.kaydet()
        self._not(u'%s rescued, %s taken'
                  % (_cogul(self.kurtarilan, 'measurement', 'measurements'),
                     _cogul(self.yeni, 'new measurement', 'new measurements')))

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.kapat()
        return False
