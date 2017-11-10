#!/usr/bin/python
 
import sys
import argparse
import fileinput

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', action='store_true')
    parser.add_argument('-r', action='store_true')
    args = parser.parse_args()

    if args.a:
        for line in sys.stdin:
            if line[0:1] == "#":
                sys.stdout.write(line)
            else:
                sys.stdout.write('chr' + line)
    elif args.r:
        for line in sys.stdin:
            if line[0:1] == "#":
                sys.stdout.write(line)
            else:
                sys.stdout.write(line[3:])