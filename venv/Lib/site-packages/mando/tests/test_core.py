import pytest
from mando import Program


program = Program('example.py', '1.0.10')


def NoopCompleter(prefix, **kwd):
    return []


program.option(
    "-f", "--foo", dest='foo', default='bar', completer=NoopCompleter,
    help="Real programmers don't comment their code. \
          If it was hard to write, it should be hard to read."
)

program.add_subprog('sub')
program.sub.option(
    "-i", "--inc", dest='inc', type=int, default=0,
    help="Some help text."
)


@program.command
def getopt(name):
    '''
        :param name: Name of option to return.
    '''
    # also allows for: Script.foo
    return getattr(program, name)


@program.sub.command
def powOfSub(b, e):
    '''
        :param b: Base.
        :param e: Exponent.
    '''
    return int(b) ** int(e) + program.inc


@program.sub.command('powOfSub2')
def powOfSub2_impl(b, e):
    '''
        :param b: Base.
        :param e: Exponent.
    '''
    return int(b) ** int(e) - program.inc


@program.command
def goo(pos, verbose=False, bar=None):
    pass


@program.command
def vara(pos, foo, spam=24, *vars):
    '''
    :param vars: Yeah, you got it right, the variable arguments.
    '''
    pass


@program.command
def another(baw, owl=42, json=False, tomawk=None):
    '''This yet another example showcasing the power of Mando!

    :param baw: That's the positional argument, obviously.
    :param -o, --owl: Yeah, I know, this is too much.
    :param -j, --json: In case you want to pipe it through something.
    :param -t, --tomawk: Well, in this case -t isn't for time.'''
    pass


@program.command('alias')
def analiased(a, b=4):
    pass


@program.command
def power(x, y=2):
    return int(x) ** y


@program.command('more-power')
def more_power(x, y=2):
    '''This one really shows off complete power.

    :param x <int>: Well, the base.
    :param -y <int>: You got it, the exponent.'''

    return x ** y


@program.command('more-powerful')
@program.arg('x', type=int, completer=NoopCompleter)
@program.arg('y', '-y', '--epsilon', type=int)
def more_power_2(x, y=2):
    return x ** y


@program.command
@program.arg('x', type=int)
@program.arg('y', type=int)
def overriding(x, y=4):
    '''Yoo an override test.

    :param x <str>: This is so wroong!!! Let's hope it gets overridden by @arg.
    :param -y <metavar>: This too!!'''

    return x - y


@program.command
def dashes(a, b=5):
    '''Usual command help.

    :param a <int>: A help obviously.
    :param b <int>: Yooo.'''
    return a ** b


@program.command
def append(acc=[]):
    return acc


GENERIC_COMMANDS_CASES = [
    ('goo 2', [['2', False, None]]),
    ('goo 2 --verbose', [['2', True, None]]),
    ('goo 2 --bar 9', [['2', False, '9']]),
    ('goo 2 --verbose --bar 8', [['2', True, '8']]),
    ('vara 2 3', [['2', '3', 24]]),
    ('vara 2 3 --spam 8', [['2', '3', 8]]),
    # Unfortunately this is an argparse "bug". See:
    # http://bugs.python.org/issue15112
    # You cannot intermix positional and optional arguments for now.
    #('vara 1 2 --spam 8 9 8', ['1', '2', 8, '9', '8']),
    ('vara 1 2 4 5 --spam 8', [['1', '2', 8, '4', '5']]),
    ('vara --spam 8 1 2 4 5', [['1', '2', 8, '4', '5']]),
    ('vara 9 8 1 2 3 4', [['9', '8', 24, '1', '2', '3', '4']]),
    ('another 2', [['2', 42, False, None]]),
    ('another 2 -j', [['2', 42, True, None]]),
    ('another 2 -t 1 -o 3', [['2', 3, False, '1']]),
    ('another 2 --owl 89 --tomawk 98', [['2', 89, False, '98']]),
    ('another 2 --json -o 1', [['2', 1, True, None]]),
    ('another 3 --owl 8 --json --tomawk 8', [['3', 8, True, '8']]),
    ('alias 5 -b 9', [['5', 9], 'analiased']),
    ('more-power 9 -y 2', [[9, 2], 'more_power']),
    ('more-powerful 9 -y 3', [[9, 3], 'more_power_2']),
    ('more-powerful 9 --epsilon 3', [[9, 3], 'more_power_2']),
    ('overriding 2', [[2, 4]]),
    ('overriding 2 -y 7', [[2, 7]]),
    ('dashes 2', [[2, 5]]),
    ('dashes 8 -b 7', [[8, 7]]),
    ('append', [[[]]]),
    ('append --acc 2', [[['2']]]),
    ('append --acc 2 --acc 3', [[['2', '3']]]),
]


@pytest.mark.parametrize('args,rest', GENERIC_COMMANDS_CASES)
def test_generic_commands(args, rest):
    args = args.split()
    if len(rest) == 1:
        to_args = rest[0]
        real_name = args[0]
    else:
        to_args = rest[0]
        real_name = rest[1]
    parsed = program.parse(args)
    assert real_name == parsed[0].__name__
    assert to_args == parsed[1]


PROGRAM_EXECUTE_CASES = [
    ('power 2', 4),
    ('power 2 -y 4', 16),
    ('more-power 3', 9),
    ('more-power 3 -y 4', 81),
    ('more-powerful 4 -y 2', 16),
    ('more-powerful 4 --epsilon 2', 16),
    ('overriding 2', -2),
    ('overriding 2 -y 7', -5),
    ('dashes 2', 32),
    ('dashes 7 -b 3', 343),
]


@pytest.mark.parametrize('args,result', PROGRAM_EXECUTE_CASES)
def test_program_execute(args, result):
    args = args.split()
    assert result == program.execute(args)
    assert program.parse(args)[0].__name__ == program._current_command


PROGRAM_OPTIONS_CASES = [
    ('          getopt foo', 'bar'),
    ('   -f xyz getopt foo', 'xyz'),
    ('--foo xyz getopt foo', 'xyz'),
    ('          sub         powOfSub  2 3', 8),
    ('   -f xyz sub    -i 1 powOfSub  2 3', 9),
    ('--foo xyz sub --inc 2 powOfSub  2 3', 10),
    ('          sub         powOfSub2 2 3', 8),
    ('   -f xyz sub    -i 1 powOfSub2 2 3', 7),
    ('--foo xyz sub --inc 2 powOfSub2 2 3', 6),
]


@pytest.mark.parametrize('args,result', PROGRAM_OPTIONS_CASES)
def test_program_options(args, result):
    args = args.split()
    assert "example.py" == program.name
    assert result == program.execute(args)

