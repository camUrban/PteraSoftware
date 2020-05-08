__version__ = '0.6.4'

try:
    from mando.core import Program
except ImportError as e:  # pragma: no cover
    # unfortunately the only workaround for Python2.6, argparse and setup.py
    e.version = __version__
    raise e

main = Program()
command = main.command
arg = main.arg
parse = main.parse
execute = main.execute
