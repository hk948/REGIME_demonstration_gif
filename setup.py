from cx_Freeze import setup, Executable

build_exe_options = {
    "packages": ["numpy", "matplotlib"],
    "includes": ["tkinter"],
    "include_files": [],
}

# base = "Win32GUI"

setup(
    name="REGIME potential",
    version="1.0.0",
    description="Simulates impacts of disturbance regimes on carbon storage dynamics",
    options={"build_exe": build_exe_options},
    executables=[Executable("main.py", icon=None)]
)
