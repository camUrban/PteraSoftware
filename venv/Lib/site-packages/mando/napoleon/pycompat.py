from six import PY3


# UnicodeMixin
if PY3:
    class UnicodeMixin(object):
        """Mixin class to handle defining the proper __str__/__unicode__
        methods in Python 2 or 3."""

        def __str__(self):
            return self.__unicode__()
else:
    class UnicodeMixin(object):
        """Mixin class to handle defining the proper __str__/__unicode__
        methods in Python 2 or 3."""

        def __str__(self):
            return self.__unicode__().encode('utf8')
