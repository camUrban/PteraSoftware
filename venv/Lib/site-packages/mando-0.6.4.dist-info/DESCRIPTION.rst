mando: CLI interfaces for Humans!
=================================

.. image:: https://img.shields.io/travis/rubik/mando/master.svg
    :alt: Travis-CI badge
    :target: https://travis-ci.org/rubik/mando

.. image:: https://img.shields.io/coveralls/rubik/mando/master.svg
    :alt: Coveralls badge
    :target: https://coveralls.io/r/rubik/mando

.. image:: https://img.shields.io/pypi/v/mando.svg
    :alt: Latest release
    :target: https://pypi.python.org/pypi/mando

.. image:: https://img.shields.io/pypi/dm/mando.svg
    :alt: PyPI downloads count
    :target: https://pypi.python.org/pypi/mando

.. image:: https://img.shields.io/pypi/format/mando.svg
    :alt: Download format
    :target: http://pythonwheels.com/

.. image:: https://img.shields.io/pypi/l/mando.svg
    :alt: Mando license
    :target: https://pypi.python.org/pypi/mando/

mando is a wrapper around ``argparse``, and allows you to write complete CLI
applications in seconds while maintaining all the flexibility.

Installation
------------

Mando is tested across all Python versions from **Python 2.6** to **Python
3.6** and also on **Pypy**. You can install it with Pip::

    $ pip install mando

The problem
-----------

While ``argparse`` is great for simple command line applications with only
one, default command, when you have to add multiple commands and manage them
things get really messy and long. But don't worry, mando comes to help!

Quickstart
----------

.. code-block:: python

    from mando import command, main

    @command
    def echo(text, capitalize=False):
        if capitalize:
            text = text.upper()
        print(text)

    if __name__ == '__main__':
        main()

Generated help:

.. code-block:: console

    $ python example.py -h
    usage: example.py [-h] {echo} ...

    positional arguments:
      {echo}
        echo      Echo the given text.

    optional arguments:
      -h, --help  show this help message and exit

    $ python example.py echo -h
    usage: example.py echo [-h] [--capitalize] text

    Echo the given text.

    positional arguments:
      text

    optional arguments:
      -h, --help    show this help message and exit
      --capitalize

Actual usage:

.. code-block:: console

    $ python example.py echo spam
    spam
    $ python example.py echo --capitalize spam
    SPAM


A *real* example
----------------

Something more complex and real-world-*ish*. The code:

.. code-block:: python

    from mando import command, main


    @command
    def push(repository, all=False, dry_run=False, force=False, thin=False):
        '''Update remote refs along with associated objects.

        :param repository: Repository to push to.
        :param --all: Push all refs.
        :param -n, --dry-run: Dry run.
        :param -f, --force: Force updates.
        :param --thin: Use thin pack.'''

        print ('Pushing to {0}. All: {1}, dry run: {2}, force: {3}, thin: {4}'
               .format(repository, all, dry_run, force, thin))


    if __name__ == '__main__':
        main()

mando understands Sphinx-style ``:param:``'s in the docstring, so it creates
short options and their help for you.

.. code-block:: console

    $ python git.py push -h
    usage: git.py push [-h] [--all] [-n] [-f] [--thin] repository

    Update remote refs along with associated objects.

    positional arguments:
      repository     Repository to push to.

    optional arguments:
      -h, --help     show this help message and exit
      --all          Push all refs.
      -n, --dry-run  Dry run.
      -f, --force    Force updates.
      --thin         Use thin pack.

Let's try it!

.. code-block:: console

    $ python git.py push --all myrepo
    Pushing to myrepo. All: True, dry run: False, force: False, thin: False
    $ python git.py push --all -f myrepo
    Pushing to myrepo. All: True, dry run: False, force: True, thin: False
    $ python git.py push --all -fn myrepo
    Pushing to myrepo. All: True, dry run: True, force: True, thin: False
    $ python git.py push --thin -fn myrepo
    Pushing to myrepo. All: False, dry run: True, force: True, thin: True
    $ python git.py push --thin
    usage: git.py push [-h] [--all] [-n] [-f] [--thin] repository
    git.py push: error: too few arguments

Amazed uh? Yes, mando got the short options and the help from the docstring!
You can put much more in the docstring, and if that isn't enough, there's an
``@arg`` decorator to customize the arguments that get passed to argparse.

Mando has lots of other options. For example, it supports different docstring
styes (Sphinx, Google and NumPy), supports shell autocompletion via the
``argcomplete`` package and supports custom format classes. For a complete
documentation, visit https://mando.readthedocs.org/.


