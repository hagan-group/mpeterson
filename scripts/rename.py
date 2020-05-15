#! /usr/bin/env python3

import argparse
import os
import re
import shlex
import subprocess
import sys

from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Bulk rename files and directories, with support for '\
                    'specific naming conventions and replacement rules.'
    )
    
    parser.add_argument('-f', '--files',
        action='store_true',
        help='Also rename files'
    )
    parser.add_argument('-n', '--dry-run',
        action='store_true',
        help='Performs a dry run. Prints the changes that would be made '\
             'without making them.'
    )
    parser.add_argument('-p', '--positive',
        action='store_true',
        help='Treat all numbers as positive when deciding how to rename '\
             'possible kebab-case files and directories'
    )
    parser.add_argument('-r', '--recursive',
        action='store_true',
        help='Whether to recursively enter directories'
    )

    parser.add_argument('-R', '--replace',
        default=[],
        action='append',
        type=make_rule,
        help='Replacement rules in the form "from:to[:count]". The first '\
             '`count` substrings `from` are replaced with `to` (defaults to '\
             'replacing all matched substrings). Many such rules can given, '\
             'and replacements are made in the order they appear.'
    )
    parser.add_argument('--style',
        choices=[
            'MACRO_CASE',
            'kebab-case',
            'snake_case',
            'PascalCase',
            'camelCase',
            'Title Case',
            'Sentence case',
        ],
        help='Converts the name to the desired style.'
    )
    parser.add_argument('src',
        nargs='+',
        help='A directory or file to rename with the given styles and rules.'
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-u', '--upper',
        action='store_true',
        help='Maintain the current format, but force all letters to uppercase.'
    )
    group.add_argument('-l', '--lower',
        action='store_true',
        help='Maintain the current format, but force all letters to lowercase.'
    )
    
    return parser.parse_args()


def make_rule(rule):
    rule = rule.split(':')
    if len(rule) not in (2, 3):
        raise TypeError(f"Invalid rule {rule}")
    else:
        old = rule[0]
        new = rule[1]
        try:
            cnt = int(rule[2])
        except:
            cnt = -1

    return (old, new, cnt)


def rename_to_canonical(name, opts):
    """
    Converts a string to snake_case.
    """
    
    number_regex = '' if opts.positive else r'(?!\d|\.)'
    tokens = re.split(rf"""(?x)  # enable verbose mode
        \s                      # split at any whitespace
        | _                     # or at an underscore
        | \B(?=[A-Z])(?<![A-Z]) # or before a mid-word capitalized letter
        | -{number_regex}       # or at a hyphen not followed by a number
        | (?=-(?:\d|\.))        # or before any negative number
        | (?<=[a-zA-Z])(?=\d)   # or between a letter and number
        """,
        name,
    )
    tokens = filter(lambda tok: len(tok) > 0, tokens)
    name = '_'.join(tok.lower() for tok in tokens)

    # macro case
    if name == name.upper():
        return name.lower();

    return name


def rename_from_canonical(name: str, style: str):
    """
    Assumes that `name` is in canonical form; that is, snake_case. This function
    converts `name` to the desired style.
    """
    
    if style == 'snake_case':
        return name
    elif style == 'MACRO_CASE':
        return name.upper()
    elif style == 'camelCase':
        words = name.split('_')
        return words[0] + ''.join(word.capitalize() for word in words[1:])
    elif style == 'PascalCase':
        return ''.join(word.capitalize() for word in name.split('_'))
    elif style == 'kebab-case':
        return name.replace('_', '-')
    elif style == 'Title Case':
        return ' '.join(word.capitalize() for word in name.split('_'))
    elif style == 'Sentence case':
        return ' '.join(name.split('_')).capitalize()


def renamed(name: str, opts):
    for rule in opts.replace:
        name = name.replace(*rule)

    if opts.style is not None:
        name = rename_from_canonical(rename_to_canonical(name, opts), opts.style)

    if opts.upper:
        name = name.upper()
    elif opts.lower:
        name = name.lower()

    return name


def rename_once(src, opts):
    dst = src.parent/renamed(src.name, opts)
    if src != dst and not dst.exists():
        print(src, '-->', dst)
        if not opts.dry_run:
            subprocess.run(shlex.split(f'mv --no-clobber "{src}" "{dst}"'),
                           check=True)
    elif src != dst and dst.exists():
        print(f'WARNING: not moving {src} to {dst} because {dst} already exists')
    return src if opts.dry_run else dst


def rename_recursive(start, opts):
    root = Path(start).expanduser()
    
    items = root.glob('*')
    if opts.files:
        items = filter(lambda item: item.is_dir() or item.is_file(), items)
    else:
        items = filter(lambda item: item.is_dir(), items)
    
    for item in items:
        item = rename_once(item, opts)
        if item.is_dir():
            rename_recursive(item, opts)


def main():
    opts = parse_args()

    if opts.dry_run:
        print("----- Begin Dry Run -----")

    for i, src in enumerate(opts.src):
        if opts.recursive:
            rename_recursive(src, opts)
        else:
            rename_once(Path(src), opts)

    if opts.dry_run:
        print("----- End Dry Run -----")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except subprocess.CalledProcessError as e:
        print("ERROR:", e)
        sys.exit(e.returncode)
    except Exception as e:
        print("ERROR:", e)
        sys.exit(128)
