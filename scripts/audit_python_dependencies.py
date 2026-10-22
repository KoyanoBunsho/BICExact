"""Audit every pinned Poetry dependency, including platform-specific packages."""
from pathlib import Path
import subprocess
import tempfile
import tomllib


def main():
    root = Path(__file__).resolve().parents[1]
    with (root / "src/poetry.lock").open("rb") as stream:
        lock = tomllib.load(stream)
    with tempfile.TemporaryDirectory(prefix="bicexact-audit-") as directory:
        requirements = Path(directory) / "requirements.txt"
        requirements.write_text("".join(
            f"{package['name']}=={package['version']}\n"
            for package in lock["package"]
        ))
        return subprocess.run([
            "pip-audit", "--disable-pip", "--no-deps", "-r", str(requirements),
        ], check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
