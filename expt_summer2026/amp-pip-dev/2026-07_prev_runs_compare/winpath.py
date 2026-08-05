#!/usr/bin/env python3
"""Convert a Windows path to its WSL/Linux mount-point equivalent.

Usage:
  python winpath.py "Z:\\2026-05 Cindy\\P#5 Amplicon NGS\\file.fa"

Prints:
  /mnt/z/2026-05 Cindy/P#5 Amplicon NGS/file.fa

Reconstructed Aug 2026 from shell history after accidental deletion.
"""
import sys


def win_to_wsl(winpath):
    """Convert 'Z:\\foo\\bar' -> '/mnt/z/foo/bar'."""
    p = winpath.replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        drive = p[0].lower()
        p = f"/mnt/{drive}{p[2:]}"
    return p


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <windows_path>", file=sys.stderr)
        return 1
    print(win_to_wsl(sys.argv[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
