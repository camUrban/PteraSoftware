"""Option handling polyfill for Flake8 2.x and 3.x."""
import optparse
import os


def register(parser, *args, **kwargs):
    r"""Register an option for the Option Parser provided by Flake8.

    :param parser:
        The option parser being used by Flake8 to handle command-line options.
    :param \*args:
        Positional arguments that you might otherwise pass to ``add_option``.
    :param \*\*kwargs:
        Keyword arguments you might otherwise pass to ``add_option``.
    """
    try:
        # Flake8 3.x registration
        parser.add_option(*args, **kwargs)
    except (optparse.OptionError, TypeError):
        # Flake8 2.x registration
        # Pop Flake8 3 parameters out of the kwargs so they don't cause a
        # conflict.
        parse_from_config = kwargs.pop('parse_from_config', False)
        comma_separated_list = kwargs.pop('comma_separated_list', False)
        normalize_paths = kwargs.pop('normalize_paths', False)
        # In the unlikely event that the developer has specified their own
        # callback, let's pop that and deal with that as well.
        base_callback = kwargs.pop('callback', store_callback)
        callback = generate_callback_from(comma_separated_list,
                                          normalize_paths,
                                          base_callback)
        kwargs['callback'] = callback
        kwargs['action'] = 'callback'

        # We've updated our args and kwargs and can now rather confidently
        # call add_option.
        option = parser.add_option(*args, **kwargs)
        if parse_from_config:
            parser.config_options.append(option.get_opt_string().lstrip('-'))


def parse_comma_separated_list(value):
    """Parse a comma-separated list.

    :param value:
        String or list of strings to be parsed and normalized.
    :returns:
        List of values with whitespace stripped.
    :rtype:
        list
    """
    if not value:
        return []

    if not isinstance(value, (list, tuple)):
        value = value.split(',')

    return [item.strip() for item in value]


def normalize_path(path, parent=os.curdir):
    """Normalize a single-path.

    :returns:
        The normalized path.
    :rtype:
        str
    """
    # NOTE(sigmavirus24): Using os.path.sep allows for Windows paths to
    # be specified and work appropriately.
    separator = os.path.sep
    if separator in path:
        path = os.path.abspath(os.path.join(parent, path))
    return path.rstrip(separator)


def parse_normalized_paths(value):
    """Normalize the path(s) value."""
    if isinstance(value, list):
        normalized = [normalize_path(s) for s in value]
    else:
        normalized = normalize_path(value)
    return normalized


def store_callback(option, opt_str, value, parser, *args, **kwargs):
    """Implement optparse's "store" action as a callback."""
    setattr(parser.values, option.dest, value)


def generate_callback_from(comma_separated_list, normalize_paths,
                           base_callback):
    """Generate a callback from parameters provided for the option."""
    def _callback(option, opt_str, value, parser, *args, **kwargs):
        """Wrap `base_callback` by transforming `value` for option params."""
        if comma_separated_list:
            value = parse_comma_separated_list(value)
        if normalize_paths:
            value = parse_normalized_paths(value)
        base_callback(option, opt_str, value, parser, *args, **kwargs)

    return _callback
