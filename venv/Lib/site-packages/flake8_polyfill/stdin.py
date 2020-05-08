"""Monkey-patching for pep8 and pycodestyle."""
try:
    import pep8
except ImportError:
    pep8 = None

try:
    import pycodestyle
except ImportError:
    pycodestyle = None

from flake8_polyfill import version

__all__ = ('monkey_patch',)

modules = {
    'pep8': [pep8],
    'pycodestyle': [pycodestyle],
    'all': [pep8, pycodestyle],
}


def monkey_patch(which):
    """Monkey-patch the specified module with the appropriate stdin.

    On Flake8 2.5 and lower, Flake8 would would monkey-patch
    ``pep8.stdin_get_value`` for everyone. This avoided problems where
    stdin might be exhausted.

    On Flake8 2.6, Flake8 stopped patching ``pep8`` and started
    monkey-patching ``pycodestyle.stdin_get_value``.

    On Flake8 3.x, Flake8 has no need to monkey patch either ``pep8`` or
    ``pycodestyle``.

    This function accepts three parameters:

    - pep8
    - pycodestyle
    - all

    "all" is a special value that will monkey-patch both "pep8" and
    "pycodestyle".

    :param str which:
        The name of the module to patch.
    :returns:
        Nothing.
    :rtype:
        NoneType
    """
    if (2, 0) <= version.version_info < (3, 0):
        from flake8.engine import pep8 as _pep8
        stdin_get_value = _pep8.stdin_get_value
    elif (3, 0) <= version.version_info < (4, 0):
        from flake8 import utils
        stdin_get_value = utils.stdin_get_value

    for module in modules[which]:
        if module is None:
            continue
        module.stdin_get_value = stdin_get_value
