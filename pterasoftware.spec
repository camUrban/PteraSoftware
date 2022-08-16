# -*- mode: python ; coding: utf-8 -*-


import os
import sys

import PySide2
from PyInstaller.building.api import PYZ, EXE, COLLECT
from PyInstaller.building.build_main import Analysis

block_cipher = None

one_dir_mode = True

binaries = []
if sys.platform.startswith('win'):
    qt_plugins_path = os.path.join(PySide2.__path__[0], "plugins")
    binaries = [
        (os.path.join(PySide2.__path__[0], "plugins"), 'PySide2')
    ]
elif sys.platform.startswith('linux'):
    qt_plugins_path = os.path.join(PySide2.__path__[0], "Qt", "plugins", "platforms")
    binaries = [
        (os.path.join(sys.base_prefix, "lib", "libspatialindex_c.so"), '.'),
        # (os.path.join(PySide2.__path__[0], "Qt", "plugins", "platforms"), '.')
    ]

upx = False  # UPX does not play with anything Qt
upx_exclude = [
    'PySide2',
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
                 ("README.md", '.')
                 ],
    hiddenimports=['vtkmodules','vtkmodules.all','vtkmodules.qt.QVTKRenderWindowInteractor','vtkmodules.util','vtkmodules.util.numpy_support'],
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
    name='main',
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
