"""Isolated Lumerical v252 launch support; importing this module never runs a solver.

The user's existing CAD process is never attached to or terminated.  Each context
creates a new hidden CAD process, restores resource settings, and closes it.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import importlib
import json
import os
import subprocess
import sys
import threading
import time
import warnings
import ctypes


LUMERICAL_ROOT = Path(os.environ.get("LUMERICAL_ROOT", r"D:\Ansys Inc\v252\Lumerical"))
LUMAPI_DIRECTORY = LUMERICAL_ROOT / "api" / "python"
LUMERICAL_PYTHON = LUMERICAL_ROOT / "python" / "python.exe"
CPU_ENGINE = LUMERICAL_ROOT / "bin" / "fdtd-engine.exe"
LICENSE_SERVER = "1055@localhost"
DEFAULT_THREADS = 8
_LAUNCH_LOCK = threading.RLock()

def load_lumapi():
    """Import the installed API without launching CAD or modifying license files."""
    api_path = str(LUMAPI_DIRECTORY)
    if api_path not in sys.path:
        sys.path.insert(0, api_path)
    return importlib.import_module("lumapi")


def child_environment():
    """Return a child-only environment, preserving existing valid license setup.

    License.ini may have precedence over this environment variable.  No ini file
    or persistent Windows user/system environment is changed here.
    """
    environment = os.environ.copy()
    environment.setdefault("ANSYSLMD_LICENSE_FILE", LICENSE_SERVER)
    environment["PATH"] = str(LUMERICAL_ROOT / "bin") + os.pathsep + environment.get("PATH", "")
    return environment


def _resource_snapshot(fdtd, row):
    properties = ("active", "device type", "processes", "threads", "capacity", "job launching preset")
    snapshot = {}
    for prop in properties:
        try:
            snapshot[prop] = fdtd.getresource("FDTD", row, prop)
        except Exception:
            pass  # e.g. capacity is unavailable for GPU rows; never invent it.
    return snapshot


def _local_cpu_row(fdtd):
    for row in range(1, int(fdtd.getresource("FDTD")) + 1):
        if str(fdtd.getresource("FDTD", row, "device type")).upper() == "CPU":
            host = str(fdtd.getresource("FDTD", row, "hostname")).lower()
            if host in ("localhost", "127.0.0.1", "", os.environ.get("COMPUTERNAME", "").lower()):
                return row
    raise RuntimeError("No existing local CPU resource; refusing to alter remote resources")


def configure_cpu(fdtd, threads=DEFAULT_THREADS, resource_row=None):
    """Configure only this CAD session's local CPU row; return observed settings.

    Local Computer uses the non-MPI engine and therefore one solver process.
    The hidden_fdtd context restores these properties before closing CAD.
    """
    if isinstance(threads, bool) or int(threads) != threads or not 1 <= int(threads) <= 16:
        raise ValueError("threads must be an integer in [1, 16]")
    row = _local_cpu_row(fdtd) if resource_row is None else int(resource_row)
    for prop, value in (("device type", "CPU"), ("processes", "1"),
                        ("threads", str(int(threads))), ("capacity", "1"),
                        ("job launching preset", "Local Computer"), ("active", "1")):
        fdtd.setresource("FDTD", row, prop, value)
    observed = {"resource_row": row, "name": fdtd.getresource("FDTD", row, "name"),
                **_resource_snapshot(fdtd, row)}
    if int(float(observed["threads"])) != int(threads):
        raise RuntimeError(f"CPU thread configuration did not stick: {observed}")
    if str(observed["job launching preset"]) != "Local Computer":
        raise RuntimeError(f"Expected non-MPI Local Computer resource: {observed}")
    if "processes" in observed and int(float(observed["processes"])) != 1:
        raise RuntimeError(f"Expected one solver process: {observed}")
    return observed


@contextmanager
def hidden_fdtd(threads=DEFAULT_THREADS):
    """Yield a fresh hidden blank CAD session, close it even on caller failure.

    No model, source, or solver object is created by this helper.  Callers should
    use fdtd.run('FDTD', 'CPU', resource_name) if running through CAD, or close CAD
    after saving layout and call run_cpu_engine for separate stdout/stderr logs.
    """
    fdtd = None
    row = None
    snapshot = None
    wake_request = False
    try:
        if os.name == 'nt':
            # Process/thread-scoped idle-sleep inhibition for an authorized
            # overnight solve. No persistent power/display setting is changed.
            result=ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
            if not result:raise RuntimeError('Could not request system-awake state during FDTD solve')
            wake_request=True
        # lumapi launches through the parent environment; only missing license
        # information is temporarily supplied, then immediately restored.
        with _LAUNCH_LOCK:
            original = os.environ.get("ANSYSLMD_LICENSE_FILE")
            if original is None:
                os.environ["ANSYSLMD_LICENSE_FILE"] = LICENSE_SERVER
            try:
                fdtd = load_lumapi().FDTD(hide=True, serverArgs={"threads": str(int(threads))})
            finally:
                if original is None:
                    os.environ.pop("ANSYSLMD_LICENSE_FILE", None)
        row = _local_cpu_row(fdtd)
        snapshot = _resource_snapshot(fdtd, row)
        configure_cpu(fdtd, threads, row)
        yield fdtd
    finally:
        try:
            if fdtd is not None:
                try:
                    if snapshot is not None:
                        # Restore launch preset first, since it controls which other
                        # properties are editable in the resource manager.
                        ordered = ["job launching preset"] + [p for p in snapshot if p != "job launching preset"]
                        for prop in ordered:
                            if prop in snapshot:
                                try:
                                    fdtd.setresource("FDTD", row, prop, snapshot[prop])
                                except Exception as exc:
                                    warnings.warn(f"Could not restore own CAD resource {row}/{prop}: {exc}")
                finally:
                    fdtd.close()
        finally:
            if wake_request:ctypes.windll.kernel32.SetThreadExecutionState(0x80000000)


def cpu_engine_command(project_path, threads=DEFAULT_THREADS):
    """One non-MPI CPU process; shell=False preserves spaces and Unicode paths."""
    project = Path(project_path).resolve()
    if project.suffix.lower() != ".fsp":
        raise ValueError("Expected an .fsp project")
    if not 1 <= int(threads) <= 16:
        raise ValueError("threads must be in [1, 16]")
    return [str(CPU_ENGINE), "-t", str(int(threads)), "-fullinfo", "-log-stdout", str(project)]


def run_cpu_engine(project_path, log_dir, threads=DEFAULT_THREADS, timeout=None):
    """Run an existing authorized project, capturing distinct stdout/stderr logs.

    This is opt-in: it is never called at import or by the environment probe.
    A timeout terminates only the child process started by subprocess.run.
    """
    project = Path(project_path).resolve(strict=True)
    logs = Path(log_dir).resolve()
    logs.mkdir(parents=True, exist_ok=True)
    command = cpu_engine_command(project, threads)
    stdout_path, stderr_path = logs / "engine_stdout.log", logs / "engine_stderr.log"
    environment = child_environment()
    environment["OMP_NUM_THREADS"] = str(int(threads))
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        result = subprocess.run(command, cwd=project.parent, env=environment,
                                stdout=stdout, stderr=stderr, shell=False,
                                creationflags=flags, timeout=timeout, check=False)
    return {"command": command, "returncode": result.returncode,
            "started_utc": started_utc, "elapsed_seconds": time.perf_counter() - started,
            "stdout": str(stdout_path), "stderr": str(stderr_path),
            "threads": int(threads), "processes": 1, "device": "CPU"}
