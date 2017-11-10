#!/usr/bin/python -u
 
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', action='store_true')
    parser.add_argument('-r', action='store_true')
    args = parser.parse_args()

    sys.stderr.write('Begin streaming VCF.')

    if args.a:
        counter = 0
        counter_write = 0
        while 1:
            line = sys.stdin.readline()
            if not line:
                break
            counter += 1
            if line[0:1] == "#":
                sys.stdout.write(line)
                counter_write += 1
            else:
                sys.stdout.write('chr' + line)
                sys.stdout.flush()
                counter_write += 1
        sys.stderr.write(str(counter) + '\n')
        sys.stderr.write(str(counter_write) + '\n')
    elif args.r:
        for line in sys.stdin:
            if line[0:1] == "#":
                sys.stdout.write(line)
            else:
                sys.stdout.write(line[3:])
    sys.stderr.write('Streaming of VCF file complete.')
