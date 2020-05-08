import re
import textwrap

SPHINX_RE = re.compile(
    r'^([\t ]*):'
    r'(?P<field>param|type|returns|rtype|parameter|arg|argument|key|keyword)'
    r' ?(?P<var1>[-\w_]+,?)?'
    r' ?(?P<var2>[-<>\w_]+)?'
    r' ?(?P<var3>[<>\w_]+)?:'
    r'(?P<help>[^\n]*\n+((\1[ \t]+[^\n]*\n)|\n)*)',
    re.MULTILINE)
ARG_RE = re.compile(
    r'-(?P<long>-)?'
    r'(?P<key>(?(long)[^ =,]+|.))[ =]?'
    r'(?P<meta>[^ ,]+)?')
POS_RE = re.compile(
    r'(?P<meta>[^ ,]+)?')
ARG_TYPE_MAP = {
    'n': int, 'num': int, 'number': int,
    'i': int, 'int': int, 'integer': int,
    's': str, 'str': str, 'string': str,
    'f': float, 'float': float,
    None: None, '': None,
}


def purify_doc(string):
    '''Remove Sphinx's :param: and :type: lines from the docstring.'''
    return SPHINX_RE.sub('', string).rstrip()


def split_doc(string):
    '''Split the documentation into help and description.

    A two-value list is returned, of the form ``[help, desc]``. If no
    description is provided, the help is duplicated.'''
    parts = [part.strip() for part in string.split('\n\n', 1)]
    if len(parts) == 1:
        return parts * 2
    return parts


def purify_kwargs(kwargs):
    '''If type or metavar are set to None, they are removed from kwargs.'''
    for key, value in kwargs.copy().items():
        if key in set(['type', 'metavar']) and value is None:
            del kwargs[key]
    return kwargs


def find_param_docs(docstring):
    '''Find Sphinx's :param:, :type:, :returns:, and :rtype: lines and return
       a dictionary of the form:
       ``param: (opts, {metavar: meta, type: type, help: help})``.'''
    paramdocs = {}
    typedocs = {}
    for m in SPHINX_RE.finditer(docstring + '\n'):
        if m.group('field') in ['param',
                                'parameter',
                                'arg',
                                'argument',
                                'key',
                                'keyword']:
            # mando
            #     :param name: Help text.               name   None   None    0
            #     :param name <type>: Help text.        name   <type> None    1
            #     :param -n: Help text.                 -n     None   None    2
            #     :param -n <type>: Help text.          -n     <type> None    3
            #     :param --name: Help text.             --name None   None    4
            #     :param --name <type>: Help text.      --name <type> None    5
            #     :param -n, --name: Help text.         -n,    --name None    6
            #     :param -n, --name <type>: Help text.  -n,    --name <type>  7
            # sphinx
            #     :param name: Help text.               name   None   None    8
            #     :param type name: Help text.          type   name   None    9
            #     :type name: str

            # The following is ugly, but it allows for backward compatibility

            if m.group('var2') is None:  # 0, 2, 4, 8
                vname = m.group('var1')
                vtype = None
            # 1, 3, 5
            elif m.group('var2') is not None and '<' in m.group('var2'):
                vname = m.group('var1')
                vtype = m.group('var2')
            elif '-' in m.group('var1') and '-' in m.group('var2'):  # 6, 7
                vname = '{0} {1}'.format(m.group('var1'), m.group('var2'))
                vtype = m.group('var3')
            else:                        # 9
                vname = m.group('var2')
                vtype = m.group('var1')

            name, opts, meta = get_opts('{0} {1}'.format(vname.strip(),
                                                         vtype or ''))
            name = name.replace('-', '_')

            helpdoc = m.group('help').strip()
            helpdoc = helpdoc.splitlines(True)
            if len(helpdoc) > 1:
                helpdoc = helpdoc[0] + textwrap.dedent(''.join(helpdoc[1:]))
            else:
                helpdoc = helpdoc[0]
            paramdocs[name] = (opts, {
                'metavar': meta or None,
                'type': ARG_TYPE_MAP.get(meta.strip('<>')),
                'help': helpdoc,
            })
        elif m.group('field') == 'type':
            typedocs[m.group('var1').strip()] = m.group('help').strip()
    for key in typedocs:
        paramdocs[key][1]['type'] = ARG_TYPE_MAP.get(typedocs[key])
    return paramdocs


def get_opts(param):
    '''Extract options from a parameter name.'''
    if param.startswith('-'):
        opts = []
        names = []
        meta = None
        for long, name, meta in ARG_RE.findall(param):
            prefix = ['-', '--'][len(long)]
            opts.append('{0}{1}'.format(prefix, name))
            names.append(name)
        return max(names, key=len), opts, meta
    opt, meta = (list(filter(None, POS_RE.findall(param))) + [''])[:2]
    return opt, [opt], meta


def action_by_type(obj):
    '''Determine an action and a type for the given object if possible.'''
    kw = {}
    if isinstance(obj, bool):
        return {'action': ['store_true', 'store_false'][obj]}
    elif isinstance(obj, list):
        kw = {'action': 'append'}
    kw.update(get_type(obj))
    return kw


def get_type(obj):
    '''Determine the type of the object if among some of the built-in ones.'''
    otype = type(obj)
    if any(otype is t for t in set([int, float, str, bool])):
        return {'type': otype}
    return {}


def ensure_dashes(opts):
    '''Ensure that the options have the right number of dashes.'''
    for opt in opts:
        if opt.startswith('-'):
            yield opt
        else:
            yield '-' * (1 + 1 * (len(opt) > 1)) + opt
