from pathlib import Path
import os
import shutil
import platform
import subprocess
from setuptools import setup, Distribution
from setuptools.command.build_py import build_py
from wheel.bdist_wheel import bdist_wheel


class BuildRust(build_py):
    def run(self):
        root = Path(__file__).resolve().parent.parent
        subprocess.run(["cargo", "build", "--release", "--locked", "--bin", "rustdock-vina-server"], cwd=root, check=True)
        super().run()
        target = Path(os.environ.get("CARGO_TARGET_DIR", root / "target"))
        name = "rustdock-vina-server.exe" if os.name == "nt" else "rustdock-vina-server"
        destination = Path(self.build_lib) / "rustdock_vina" / "bin"
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target / "release" / name, destination / name)


class PlatformDistribution(Distribution):
    def has_ext_modules(self):
        return True


class NativeWheel(bdist_wheel):
    def get_tag(self):
        _, _, tag = super().get_tag()
        if platform.system() == "Darwin":
            tag = "macosx_11_0_" + platform.machine()
        return "py3", "none", tag


setup(cmdclass={"build_py": BuildRust, "bdist_wheel": NativeWheel}, distclass=PlatformDistribution)
