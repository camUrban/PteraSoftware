'''Main module containing the class Program(), which allows the conversion from
ordinary Python functions into commands for the command line. It uses
:py:module:``argparse`` behind the scenes.'''

import sys
import inspect
import argparse
try:
    getfullargspec = inspect.getfullargspec
except AttributeError:
    getfullargspec = inspect.getargspec
try:
    from itertools import izip_longest
except ImportError:  # pragma: no cover
    from itertools import zip_longest as izip_longest

from mando.napoleon import Config, GoogleDocstring, NumpyDocstring

from mando.utils import (purify_doc, action_by_type, find_param_docs,
                         split_doc, ensure_dashes, purify_kwargs)


_POSITIONAL = type('_positional', (object,), {})
_DISPATCH_TO = '_dispatch_to'


class SubProgram(object):

    def __init__(self, parser, argspecs):
        self.parser = parser
        self._subparsers = self.parser.add_subparsers()
        self._argspecs = argspecs

    @property
    def name(self):
        return self.parser.prog

    # Add global script options.
    def option(self, *args, **kwd):
        assert args and all(arg.startswith('-') for arg in args), \
            "Positional arguments not supported here"
        completer = kwd.pop('completer', None)
        arg = self.parser.add_argument(*args, **kwd)
        if completer is not None:
            arg.completer = completer
        # do not attempt to shadow existing attributes
        assert not hasattr(self, arg.dest), "Invalid option name: " + arg.dest
        return arg

    def add_subprog(self, name, **kwd):
        # also always provide help= to fix missing entry in command list
        help = kwd.pop('help', "{} subcommand".format(name))
        prog = SubProgram(self._subparsers.add_parser(name, help=help, **kwd),
                          self._argspecs)
        # do not attempt to overwrite existing attributes
        assert not hasattr(self, name), "Invalid sub-prog name: " + name
        setattr(self, name, prog)
        return prog

    def command(self, *args, **kwargs):
        '''A decorator to convert a function into a command. It can be applied
        as ``@command`` or as ``@command(new_name)``, specifying an alternative
        name for the command (default one is ``func.__name__``).'''
        if len(args) == 1 and hasattr(args[0], '__call__'):
            return self._generate_command(args[0])
        else:
            def _command(func):
                return self._generate_command(func, *args, **kwargs)
            return _command

    def arg(self, param, *args, **kwargs):
        '''A decorator to override the parameters extracted from the docstring
        or to add new ones.

        :param param: The parameter's name. It must be among the function's
            arguments names.'''
        def wrapper(func):
            if not hasattr(func, '_argopts'):
                func._argopts = {}
            func._argopts[param] = (args, kwargs)
            return func
        return wrapper

    def _generate_command(self, func, name=None, doctype='rest',
                          *args, **kwargs):
        '''Generate argparse's subparser.

        :param func: The function to analyze.
        :param name: If given, a different name for the command. The default
            one is ``func.__name__``.'''
        func_name = func.__name__
        name = func_name if name is None else name
        argspec = getfullargspec(func)
        self._argspecs[func_name] = argspec
        argz = izip_longest(reversed(argspec.args),
                            reversed(argspec.defaults or []),
                            fillvalue=_POSITIONAL())
        argz = reversed(list(argz))
        doc = (inspect.getdoc(func) or '').strip() + '\n'
        if doctype == 'numpy':
            config = Config(napoleon_google_docstring=False,
                            napoleon_use_rtype=False)
            doc = str(NumpyDocstring(doc, config))
        elif doctype == 'google':
            config = Config(napoleon_numpy_docstring=False,
                            napoleon_use_rtype=False)
            doc = str(GoogleDocstring(doc, config))
        elif doctype == 'rest':
            pass
        else:
            raise ValueError('doctype must be one of "numpy", "google", '
                             'or "rest"')
        cmd_help, cmd_desc = split_doc(purify_doc(doc))
        subparser = self._subparsers.add_parser(name,
                                                help=cmd_help or None,
                                                description=cmd_desc or None,
                                                **kwargs)
        params = find_param_docs(doc)
        for a, kw in self._analyze_func(func, params, argz, argspec.varargs):
            completer = kw.pop('completer', None)
            arg = subparser.add_argument(*a, **purify_kwargs(kw))
            if completer is not None:
                arg.completer = completer

        subparser.set_defaults(**{_DISPATCH_TO: func})
        return func

    def _analyze_func(self, func, params, argz, varargs_name):
        '''Analyze the given function, merging default arguments, overridden
        arguments (with @arg) and parameters extracted from the docstring.

        :param func: The function to analyze.
        :param params: Parameters extracted from docstring.
        :param argz: A list of the form (arg, default), containing arguments
            and their default value.
        :param varargs_name: The name of the variable arguments, if present,
            otherwise ``None``.'''
        for arg, default in argz:
            override = getattr(func, '_argopts', {}).get(arg, ((), {}))
            yield merge(arg, default, override, *params.get(arg, ([], {})))
        if varargs_name is not None:
            kwargs = {'nargs': '*'}
            kwargs.update(params.get(varargs_name, (None, {}))[1])
            yield ([varargs_name], kwargs)


class Program(SubProgram):

    def __init__(self, prog=None, version=None, **kwargs):
        parser = argparse.ArgumentParser(prog, **kwargs)
        if version is not None:
            parser.add_argument('-v', '--version', action='version',
                                version=version)

        super(Program, self).__init__(parser, dict())
        self._options = None
        self._current_command = None

    # Attribute lookup fallback redirecting to (internal) options instance.
    def __getattr__(self, attr):
        return getattr(self._options, attr)

    def parse(self, args):
        '''Parse the given arguments and return a tuple ``(command, args)``,
        where ``args`` is a list consisting of all arguments. The command can
        then be called as ``command(*args)``.

        :param args: The arguments to parse.'''
        try:
            # run completion handler before parsing
            import argcomplete
            argcomplete.autocomplete(self.parser)
        except ImportError:  # pragma: no cover
            # ignore error if not installed
            pass

        self._options = self.parser.parse_args(args)
        arg_map = self._options.__dict__
        if _DISPATCH_TO not in arg_map:  # pragma: no cover
            self.parser.error("too few arguments")

        command = arg_map.pop(_DISPATCH_TO)
        argspec = self._argspecs[command.__name__]
        real_args = []
        for arg in argspec.args:
            real_args.append(arg_map.pop(arg))
        if arg_map and arg_map.get(argspec.varargs):
            real_args.extend(arg_map.pop(argspec.varargs))
        return command, real_args

    def execute(self, args):
        '''Parse the arguments and execute the resulting command.

        :param args: The arguments to parse.'''
        command, a = self.parse(args)
        self._current_command = command.__name__
        return command(*a)

    def __call__(self):  # pragma: no cover
        '''Parse ``sys.argv`` and execute the resulting command.'''
        return self.execute(sys.argv[1:])


def merge(arg, default, override, args, kwargs):
    '''Merge all the possible arguments into a tuple and a dictionary.

    :param arg: The argument's name.
    :param default: The argument's default value or an instance of _POSITIONAL.
    :param override: A tuple containing (args, kwargs) given to @arg.
    :param args: The arguments extracted from the docstring.
    :param kwargs: The keyword arguments extracted from the docstring.'''
    opts = [arg]
    if not isinstance(default, _POSITIONAL):
        opts = list(ensure_dashes(args or opts))
        kwargs.update({'default': default, 'dest': arg})
        kwargs.update(action_by_type(default))
    else:
        # positionals can't have a metavar, otherwise the help is screwed
        # if one really wants the metavar, it can be added with @arg
        kwargs['metavar'] = None
    kwargs.update(override[1])
    return override[0] or opts, kwargs
