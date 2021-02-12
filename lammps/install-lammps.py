#! /usr/bin/env python3

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


HOME = Path.home()
EXTENSIONS = (Path(__file__)/"../extensions").resolve()
LAMMPS_REPO = "https://github.com/lammps/lammps"


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--build-type",
        choices=["Debug", "Release", "RelWithDebInfo", "MinSizeRel"],
        default="Release",
        help="Set the build type for LAMMPS. Default is 'Release'"
    )
    parser.add_argument("-j", "--jobs",
        default=os.cpu_count(),
        help="""\
            Sets the number of jobs (processes) to use when building LAMMPS.\
            Defaults to the number of available CPU cores.\
        """
    )
    parser.add_argument("--packages",
        nargs='+',
        default=[],
        help="A space-separated list of packages to enable when building LAMMPS."
    )
    parser.add_argument("--path",
        type=Path,
        default=HOME/".lammps",
        help="""\
            The path to the LAMMPS repository. LAMMPS will be cloned to this\
            location if it doesn't already exist.\
        """
    )
    parser.add_argument("--prefix",
        type=Path,
        default=HOME/".local",
        help="""\
            Where LAMMPS will be installed. Defaults to $HOME/.local, rather\
            than /usr/local, in order to avoid needing escalated privileges.
        """
    )
    parser.add_argument("--tune-flags",
        default="",
        help="""\
            Additional compiler flags to use when compiling LAMMPS. These can\
            be used to tune performance on certain architectures. If you plan\
            to run LAMMPS on the host machine, consider passing '-march=native'\
        """
    )

    return parser.parse_args()


def clone_or_pull(path):
    if path.is_dir():
        if (path/".git").is_dir():
            cmd = shlex.split(f'git -C "{path}" pull')
            subprocess.run(cmd)
        else:
            msg = f"{path} exists but is not a git repository!"
            raise RuntimeError(msg)
    else:
        cmd = shlex.split(f'git clone "{LAMMPS_REPO}" "{path}"')
        subprocess.run(cmd)

    os.chdir(path)
    subprocess.run(["git", "checkout", "stable"])


def enable_extensions():
    dst=Path.cwd()/"src"
    for extension in EXTENSIONS.glob("ext-*"):
        name = extension.name[4:].replace('-', ' ').title()
        print(f"... enabling extension: '{name}'")
        for src in extension.iterdir():
            if src.match("*.h") or src.match("*.cpp"):
                shutil.copy(src, dst)


def configure(args):
    build_path = Path("build")
    try:
        build_path.mkdir()
    except FileExistsError:
        if (build_path/"CMakeCache.txt").is_file():
            (build_path/"CMakeCache.txt").unlink()
    finally:
        os.chdir(build_path)

    cmd = f"""\
        cmake ../cmake
            -DCMAKE_BUILD_TYPE={args.build_type}
            -DCMAKE_INSTALL_PREFIX={args.prefix}
            -DCMAKE_TUNE_FLAGS="{args.tune_flags}"
    """
    
    for pkg in args.packages:
        cmd += f" -DPKG_{pkg.upper()}=ON"

    cmd = shlex.split(cmd)
    subprocess.run(cmd, env=os.environ)


def build(args):
    subprocess.run(shlex.split(f"make -j {args.jobs}"))


def install():
    subprocess.run(["make", "install"])


def main():
    args = parse_args()

    print("Cloning LAMMPS repository")
    clone_or_pull(args.path)
    
    print("Enabling extensions")
    enable_extensions()

    print("Configuring")
    configure(args)
    
    print("Building")
    build(args)
    
    print("Installing")
    install()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", e)
        sys.exit(1)