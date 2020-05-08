# This is important: it will make all literals unicode under 2.x
from __future__ import unicode_literals

import unittest

from mando import Program


program = Program('example.py', '1.0.10')


class Test_unicode_docstring_on_py2(unittest.TestCase):

    def test_py2_unicode_literals(self):
        @program.command
        def some_command():
            'this is a unicode doc-string!'

        assert not isinstance(some_command.__doc__, bytes)
        # TODO: check that the generated help is correct
