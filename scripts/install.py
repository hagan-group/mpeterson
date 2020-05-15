#! /usr/bin/env python3

import argparse
import shutil
import sys

from pathlib import Path


SUCCESS = 0
FAILURE = 1


def parse_args():
    parser = argparse.ArgumentParser(
        description="Installs the Python scripts in this directory."
    )

    parser.add_argument("--path",
        type=Path,
        default=Path('~/.local/bin'),
        help="The installation directory. Defaults to '$HOME/.local/bin'"
    )

    parser.add_argument("scripts",
        nargs="*",
        help="A list of scripts to install. Defaults to installing all of them."
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    cwd = Path.cwd()
    install_path = args.path.expanduser().resolve()

    if not args.scripts:
        scripts = [
            cwd/"make-movie-grid.py",
            cwd/"rename.py"
        ]
    else:
        for script in args.scripts:
            path = (cwd / script).with_suffix('.py')
            scripts.append(path)

    for script in scripts:
        shutil.copy(script, install_path)


if __name__ == "__main__":
    try:
        main()
        sys.exit(SUCCESS)
    except Exception as e:
        print("ERROR:", e)
        sys.exit(FAILURE)
