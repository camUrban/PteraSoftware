#!/usr/bin/env python

'''
Capture function
----------------------------------

'''

import sys
from contextlib import contextmanager

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO


@contextmanager
def capture_sys_output():
    capture_out, capture_err = StringIO(), StringIO()
    current_out, current_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = capture_out, capture_err
        yield capture_out, capture_err
    finally:
        sys.stdout, sys.stderr = current_out, current_err
