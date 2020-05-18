#! /usr/bin/env python3

import argparse
from pathlib import Path
from subprocess import call

parser = argparse.ArgumentParser(
    description="Merge a collection of movies into a single grid of movies."
)

parser.add_argument("movies",
    nargs='+',
    type=Path,
    help="A list of paths to movie files to use."
)

parser.add_argument("-s", "--shape",
    metavar=('ROWS', 'COLUMNS'),
    nargs=2,
    type=int,
    required=True,
    help="The number of rows and columns of the output movie grid."
)

parser.add_argument("-r", "--resolution",
    metavar=('WIDTH', 'HEIGHT'),
    nargs=2,
    type=int,
    default=[128, 128],
    help="The desired resolution for each individual movie in the grid. "\
         "Defaults to 128x128."
)

parser.add_argument("-c", "--clip-length",
    action="store_true",
    help="If set, will clip the length of the output movie to be the same as "\
         "the shortest input movie."
)

parser.add_argument("output",
    type=Path,
    help="The destination of the output movie."
)

args = parser.parse_args()


# constants

MOVIES = [movie.expanduser().resolve() for movie in args.movies]
NROWS = args.shape[0]
NCOLS = args.shape[1]
CLIP = int(args.clip_length)
OUTPUT = args.output.expanduser().resolve()

assert(len(MOVIES) == NROWS * NCOLS)


def pos2ind(pos):
    m, n = pos
    return m * NCOLS + n


def ind2pos(index):
    row, col = index // NROWS, index % NCOLS
    return (row, col)


def build_cmd():
    filt = ""
    
    for row in range(NROWS):
        for col in range(NCOLS):
            i = pos2ind((row, col))
            filt += f"[{i}:v]"
            filt += f"scale=h=128:w=128"
            filt += f"[a{i}];"

    for row in range(NROWS):
        for col in range(NCOLS):
            i = pos2ind((row, col))
            filt += f"[a{i}]"

        filt += f"hstack=inputs={NCOLS}:shortest={CLIP}"
        filt += f"[row{row}];"

    for row in range(NROWS):
        filt += f"[row{row}]"

    filt += f"vstack=inputs={NROWS}:shortest={CLIP}"
    filt += f"[out]"

    cmd = ["ffmpeg"]
    for movie in MOVIES:
        cmd += ["-i", f"{movie}"]

    cmd += ["-filter_complex", f"{filt}"]
    cmd += ["-map", "[out]"]
    cmd += ["-c:v", "libx264"]
    cmd += ["-pix_fmt", "yuv420p"]
    cmd += [f"{OUTPUT}"]

    return cmd


if __name__ == "__main__":
    cmd = build_cmd()

    print(cmd)
    call(cmd)
