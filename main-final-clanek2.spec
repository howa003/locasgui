# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(
    ['main-final-clanek2.py'],
    pathex=[],
    binaries=[],
    datas=[('C:\\Users\\Admin\\PycharmProjects\\LOCASGUI\\venv\\lib\\site-packages\\eel\\eel.js', 'eel'), ('static_web_folder', 'static_web_folder')],
    hiddenimports=['bottle_websocket'],
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
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='main-final-clanek2',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
