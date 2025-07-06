# -*- mode: python ; coding: utf-8 -*-


import os
import sys

import PySide6
import shiboken2
import cmocean
from PyInstaller.building.api import PYZ, EXE, COLLECT
from PyInstaller.building.build_main import Analysis

block_cipher = None

one_dir_mode = True

binaries = []
if sys.platform.startswith('win'):
    qt_plugins_path = os.path.join(PySide6.__path__[0], "Qt", "plugins")
    binaries = [
        (os.path.join(PySide6.__path__[0], "Qt", "plugins"), "PySide6")
    ]
elif sys.platform.startswith('linux'):
    qt_plugins_path = os.path.join(PySide6.__path__[0], "Qt", "plugins", "platforms")
    binaries = [
        (os.path.join(sys.base_prefix, "lib", "libspatialindex_c.so"), '.'),
        # (os.path.join(PySide6.__path__[0], "Qt", "plugins", "platforms"), '.')
    ]

upx = False  # UPX does not play with anything Qt
upx_exclude = [
    'PySide6',
    'shiboken2',
    'qwindows.dll'
]
a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=[],
    datas=[('docs/Logo.ico', 'docs'),
                 ('docs/Black_Text_Logo.ico', 'docs'),
                 ('docs/Logo.png', 'docs'),
                 ('docs/Black_Text_Logo.png', 'docs'),
                 ("README.md", '.'),
                 ('examples/analyze_steady_trim_example.py', 'examples'),
                 ('examples/analyze_unsteady_trim_example.py', 'examples'),
                 ('examples/steady_convergence_example.py', 'examples'),
                 ('examples/steady_horseshoe_vortex_lattice_method_solver.py', 'examples'),
                 ('examples/steady_ring_vortex_lattice_method_solver.py', 'examples'),
                 ('examples/unsteady_ring_vortex_lattice_method_solver_static.py', 'examples'),
                 ('examples/unsteady_ring_vortex_lattice_method_solver_variable.py', 'examples'),
                 ('examples/unsteady_ring_vortex_lattice_method_solver_variable_formation.py', 'examples'),
                 ('examples/unsteady_static_convergence_example.py', 'examples'),
                 ('examples/unsteady_variable_convergence_example.py', 'examples'),
                 (shiboken2.__path__[0], "shiboken2"),
                 (cmocean.__path__[0], "cmocean")
                 ],
    hiddenimports=['vtkmodules','vtkmodules.all','vtkmodules.qt.QVTKRenderWindowInteractor','vtkmodules.util','vtkmodules.util.numpy_support', 'examples'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='PteraSoftware',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='docs/Logo.ico'
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='main',
)
