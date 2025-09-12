# Packaging

There are 2 available packaging formats: Installer or wheel. The former is what is
distributed, while the latter is
required for use of the CLI for instance.

When packaging don't forget to check versions in both `setup.py` and
`make_installer.iss` and ensure they match!

### Installer

**Currently, this is only supported on Windows.**

The installer uses `PyInstaller` to extract all the dependencies and `InnoSetup` to
create a Windows installer package
out of them.

1. Ensure you have installed `InnoSetup` [here](https://jrsoftware.org/isdl.php).

2. Inside the root repository directory with the virtualenv active:

   ```commandline
   python -O -m PyInstaller --noconfirm "pterasoftware.spec"
   ```
   
   We run `python` with the first level optimise flag `-O` to slim down some now
   unnecessary debug code. *Do not use second
   level optimisation `-OO`, as this removes some docstrings that break dependencies.*

3. Open `make_installer.iss` in InnoSetup and run the packaging process. This can take a
   while, but should output an
   installer in the `Output` folder.
