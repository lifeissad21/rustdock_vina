"""The upstream Python API backed by a persistent native Rust engine process."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading


class Vina:
    def __init__(self, sf_name, cpu, seed, verbosity, no_refine):
        root = Path(__file__).resolve().parents[2]
        candidates = [
            os.environ.get("RUSTDOCK_VINA_SERVER"),
            Path(__file__).with_name("bin") / ("rustdock-vina-server.exe" if os.name == "nt" else "rustdock-vina-server"),
            root / "target" / "release" / ("rustdock-vina-server.exe" if os.name == "nt" else "rustdock-vina-server"),
            shutil.which("rustdock-vina-server"),
        ]
        binary = next((str(p) for p in candidates if p and Path(p).is_file()), None)
        if binary is None:
            raise RuntimeError("Build the Rust engine with cargo build --release --workspace, or set RUSTDOCK_VINA_SERVER")
        self._lock = threading.Lock()
        self._process = subprocess.Popen([binary], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, encoding="utf-8")
        try:
            self._call("new", sf_name, cpu, seed, verbosity, no_refine)
        except Exception:
            self.close()
            raise

    def _call(self, method, *args):
        payload = json.dumps({"method": method, "args": args}, allow_nan=False)
        with self._lock:
            process = self._process
            if process is None or process.poll() is not None:
                raise RuntimeError("Rust docking engine is closed or exited")
            try:
                process.stdin.write(payload + "\n")
                process.stdin.flush()
                line = process.stdout.readline()
            except (BrokenPipeError, OSError) as error:
                raise RuntimeError("Rust docking engine communication failed") from error
            if not line:
                raise RuntimeError("Rust docking engine exited without a response")
            response = json.loads(line)
            if "error" in response:
                raise RuntimeError(response["error"])
            return response["result"]

    def cite(self):
        print("AutoDock Vina: Trott & Olson (2010), DOI 10.1002/jcc.21334; Eberhardt et al. (2021), DOI 10.1021/acs.jcim.1c00203")

    def __getattr__(self, method):
        if method.startswith("_"):
            raise AttributeError(method)
        return lambda *args: self._call(method, *args)

    def close(self):
        process = getattr(self, "_process", None)
        if process is None:
            return
        self._process = None
        process.stdin.close()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        process.stdout.close()

    def __del__(self):
        self.close()
