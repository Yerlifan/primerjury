# -*- coding: utf-8 -*-
"""SPLIT DEMULTIPLEXED BAM FILES INTO ONE FASTQ PER BARCODE.

WHAT IT IS FOR
--------------
dorado (basecaller ... --kit-name ...) writes the barcode of every read into
the BC:Z tag of its BAM output. This step turns a directory of such BAMs into
<out>/barcodeNN.fastq (+ unclassified.fastq), which is what bin_reads.py and
the rest of the chain expect. The quality string and the read id are carried
over unchanged.

SAFETY
------
Nothing is overwritten. If an earlier run left fastq files behind without the
SPLIT.done marker they are renamed with a .partial suffix rather than deleted.
A BAM whose ".done" marker is missing is treated as unfinished and the step
stops, because dorado writes to a pipe and a truncated BAM still opens
(samtools quickcheck cannot tell: the BGZF EOF block is absent either way).
Pass --no-markers when the BAMs did not come from a marker-writing run.

RUN
---
    python3 steps/split_barcodes.py --bam-dir basecalled --out sequences/barcodes [--force] [--no-markers]

OUTPUT
------
    <out>/barcodeNN.fastq, <out>/unclassified.fastq, <out>/BARCODE_COUNTS.tsv, <out>/SPLIT.done
"""
from __future__ import print_function
import argparse
import collections
import io
import os
import shutil
import subprocess
import sys
import time


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--bam-dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--force', action='store_true', help='split again even if SPLIT.done exists')
    ap.add_argument('--no-markers', action='store_true', help='do not require <bam>.done markers')
    a = ap.parse_args()
    if not shutil.which('samtools'):
        raise SystemExit('STOPPED: samtools not found')
    bam_dir, out = os.path.abspath(a.bam_dir), os.path.abspath(a.out)
    done = os.path.join(out, 'SPLIT.done')
    if os.path.exists(done) and not a.force:
        print(u'  already split (%s). Use --force to redo.' % done)
        return 0
    bams = sorted(f for f in os.listdir(bam_dir) if f.endswith('.bam'))
    if not bams:
        raise SystemExit('STOPPED: no .bam in %s' % bam_dir)
    if not a.no_markers:
        missing = [f for f in bams if not os.path.exists(os.path.join(bam_dir, f + '.done'))]
        if missing:
            raise SystemExit('STOPPED: %d BAM(s) carry no .done marker (possibly truncated): %s'
                             % (len(missing), ', '.join(missing[:3])))
    if not os.path.isdir(out):
        os.makedirs(out)
    for f in os.listdir(out):
        if f.endswith('.fastq'):
            os.rename(os.path.join(out, f), os.path.join(out, f + '.partial_%d' % int(time.time())))
    writers = {}
    reads = collections.Counter()
    bases = collections.Counter()
    t0 = time.time()
    for i, f in enumerate(bams):
        p = subprocess.Popen(['samtools', 'view', os.path.join(bam_dir, f)], stdout=subprocess.PIPE)
        n = 0
        for raw in p.stdout:
            s = raw.decode('utf-8', 'replace').rstrip('\n').split('\t')
            if len(s) < 11:
                continue
            bc = 'unclassified'
            for field in s[11:]:
                if field.startswith('BC:Z:'):
                    v = field[5:]
                    bc = ('barcode' + v.rsplit('barcode', 1)[1][:2]) if 'barcode' in v else v
                    break
            g = writers.get(bc)
            if g is None:
                g = writers[bc] = io.open(os.path.join(out, bc + '.fastq'), 'w', encoding='ascii')
            g.write(u'@%s\n%s\n+\n%s\n' % (s[0], s[9], s[10]))
            reads[bc] += 1
            bases[bc] += len(s[9])
            n += 1
        p.wait()
        if p.returncode != 0:
            raise SystemExit('STOPPED: samtools view %s -> %d' % (f, p.returncode))
        print(u'  [%2d/%d] %-52s %7d reads  (%.0f s)' % (i + 1, len(bams), f[:52], n, time.time() - t0))
        sys.stdout.flush()
    for g in writers.values():
        g.close()
    rows = [u'barcode\treads\tbases']
    for bc in sorted(reads):
        rows.append(u'%s\t%d\t%d' % (bc, reads[bc], bases[bc]))
    io.open(os.path.join(out, 'BARCODE_COUNTS.tsv'), 'w', encoding='utf-8',
            newline='\n').write(u'\n'.join(rows) + u'\n')
    io.open(done, 'w', encoding='utf-8').write(u'date=%s\nbam=%d\nreads=%d\n'
                                              % (time.strftime('%Y-%m-%d %H:%M:%S'),
                                                 len(bams), sum(reads.values())))
    print(u'  %d reads, %d barcodes; %s' % (sum(reads.values()), len(reads), out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
