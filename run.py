#!/usr/bin/env python3
from datetime import datetime
from os import chdir, getcwd, makedirs
from os.path import isdir, isfile

# useful constants
VERSION = '0.0.1'
LOGFILE = None

# check if user is just printing version
if '--version' in argv:
    print("FAVITES-COVID-Lite version %s" % VERSION); exit()

# print to the log (None implies stdout only)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    print(tmp, end=end); stdout.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('--gzip_output', action='store_true', help="Gzip Compress Output Files")
    parser.add_argument('--version', action='store_true', help="Display Version")
    return parser.parse_args()

# main execution
if __name__ == "__main__":
    # parse and check user args
    args = parse_args()
    if isdir(args.output) or isfile(args.output):
        raise ValueError("Output directory exists: %s" % args.output)

    # set up output directory
    makedirs(args.output, exist_ok=True)
    LOGFILE_fn = "%s/log.txt" % args.output
    LOGFILE = open(LOGFILE_fn, 'w')

    # print run info to log
    print_log("===== RUN INFORMATION =====")
    print_log("FAVITES-ABM-HIV-SD-Lite Version: %s" % VERSION)
    print_log("FAVITES-ABM-HIV-SD-Lite Command: %s" % ' '.join(argv))
    print_log("Output Directory: %s" % args.output)
    print_log()
