
import argparse
import sys

from rst2ansi import rst2ansi


def b(s):
    # Useful for very coarse version differentiation.
    PY2 = sys.version_info[0] == 2
    PY3 = sys.version_info[0] == 3
    if PY3:
        return s.encode("utf-8")
    else:
        return s


class RSTHelpFormatter(argparse.RawTextHelpFormatter):
    """
    Custom formatter class that is capable of interpreting ReST.
    """
    def format_help(self):
        ret = rst2ansi(b(super(RSTHelpFormatter, self).format_help()) +
                       b('\n'))
        return ret.encode(sys.stdout.encoding,
                          'replace').decode(sys.stdout.encoding)

    def format_usage(self):
        ret = rst2ansi(b(super(RSTHelpFormatter, self).format_usage()) +
                       b('\n'))
        return ret.encode(sys.stdout.encoding,
                          'replace').decode(sys.stdout.encoding)
