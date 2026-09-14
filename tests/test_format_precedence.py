# -*- coding: utf-8 -*-
"""The format-operator precedence trap.

    u'%s' % v.startswith('YES')

reads like "format v, then test the prefix", but Python calls v.startswith first and formats
its result: the expression is the text 'True' or 'False', and a non-empty text is always
true. A colour rule written this way never fails: every "NO" row turns green, every warning
cell turns red, a count built from it comes out as zero. The right form is

    (u'%s' % v).startswith('YES')

This test scans every Python file in the repository for the wrong form.

RUN
    python3 tests/test_format_precedence.py
"""
from __future__ import print_function
import io
import os
import re
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATTERN = re.compile(r"""u?['"]%s['"]\s*%\s*[A-Za-z_]\w*(?:\[[^\]]*\])*\.(startswith|endswith)\(""")
SKIP = {'.git', '__pycache__', 'REFERENCE_DB'}


def findings(root=KOK):
    out = []
    for d, dirs, files in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP]
        for f in files:
            if not f.endswith('.py'):
                continue
            path = os.path.join(d, f)
            if os.path.abspath(path) == os.path.abspath(__file__):
                continue          # this file quotes the wrong form on purpose
            for i, line in enumerate(io.open(path, encoding='utf-8', errors='replace'), 1):
                code = line.split('#', 1)[0]
                if PATTERN.search(code):
                    out.append((os.path.relpath(path, root), i, line.strip()))
    return out


def test_pattern_itself():
    bad = "x = 1 if u'%s' % v.startswith(u'YES') else 2"
    good = "x = 1 if (u'%s' % v).startswith(u'YES') else 2"
    ok = PATTERN.search(bad) is not None and PATTERN.search(good) is None
    print('  pattern: %s' % ('ok' if ok else 'FAIL'))
    return ok


def test_repository():
    f = findings()
    for rel, i, line in f:
        print('  %s:%d  %s' % (rel, i, line[:100]))
    print('  repository: %s (%d finding)' % ('ok' if not f else 'FAIL', len(f)))
    return not f


if __name__ == '__main__':
    r = [test_pattern_itself(), test_repository()]
    print('format precedence: %s' % ('PASS' if all(r) else 'FAIL'))
    sys.exit(0 if all(r) else 1)
