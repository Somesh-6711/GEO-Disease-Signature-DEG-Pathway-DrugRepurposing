import subprocess
from pathlib import Path

def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("\n>>", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)
