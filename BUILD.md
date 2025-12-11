# Building the GUI Installer

This document describes how to build the **Windows installer for the GUI**. This
bundles the entire project (pterasoftware package, examples, and GUI) into a
standalone executable.

**Note:** For building the PyPI package (wheel), use `python -m build` from the project
root. That process is separate and uses `setup.cfg` and `pyproject.toml`.

When packaging the installer, don't forget to check versions in both `setup.cfg` and
`make_installer.iss` and ensure they match!

## Creating the Windows Installer

**Currently, this is only supported on Windows.**

The installer uses `PyInstaller` to extract all the dependencies and `InnoSetup` to
create a Windows installer package
out of them.

1. Ensure you have installed `InnoSetup` [here](https://jrsoftware.org/isdl.php).

2. Inside the root repository directory with the virtualenv active:

   ```commandline
   python -O -m PyInstaller --noconfirm "pterasoftware.spec"
   ```
   
   We run `python` with the first level optimize flag `-O` to slim down some now
   unnecessary debug code. *Do not use second
   level optimisation `-OO`, as this removes some docstrings that break dependencies.*

3. Open `make_installer.iss` in InnoSetup and run the packaging process. This can take a
   while, but should output an
   installer in the `Output` folder.
