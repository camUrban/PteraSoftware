import pytest
from mando import Program

from . import capture

program = Program('example.py', '1.0.10')


@program.command(doctype='numpy')
def simple_numpy_docstring(arg1, arg2='string'):
    '''One line summary.

    Extended description.

    Parameters
    ----------
    arg1 : int
        Description of `arg1`
    arg2 : str
        Description of `arg2`
    Returns
    -------
    str
        Description of return value.
    '''
    return int(arg1) * arg2


GENERIC_COMMAND_CASES = [
    ('simple_numpy_docstring 2 --arg2=test', 'testtest'),
]


@pytest.mark.parametrize('args,result', GENERIC_COMMAND_CASES)
def test_generic_command(args, result):
    args = args.split()
    assert result == program.execute(args)
    assert program.parse(args)[0].__name__ == program._current_command


NUMPY_DOCSTRING_HELP_CASES = [
    ('simple_numpy_docstring --help 2 --arg2=test', '''usage: example.py simple_numpy_docstring [-h] [--arg2 ARG2] arg1

Extended description.

positional arguments:
  arg1         Description of `arg1`

optional arguments:
  -h, --help   show this help message and exit
  --arg2 ARG2  Description of `arg2`
'''),
]


@pytest.mark.parametrize('args,result', NUMPY_DOCSTRING_HELP_CASES)
def test_numpy_docstring_help(args, result):
    args = args.split()
    with pytest.raises(SystemExit):
        with capture.capture_sys_output() as (stdout, stderr):
            program.execute(args)
    assert result == stdout.getvalue()
