import sys
import os
import argparse
import shilofue
import shilofue.TwoDSubduction as _TwoDSubduction

def main():
    """
    todo
    """
    # parse options and arguments
    parser = argparse.ArgumentParser(description='todo')
    parser.add_argument('-p', '--project', type=str,
                        default ='todo',
                        help='todo')
    parser.add_argument('-t', '--type', type=str,
                        default = 'todo',
                        help='todo')
    parser.add_argument('-i', '--inputfile', type=str,
                        default = None,
                        help='todo')
    parser.add_argument('-o', '--outputfile', type=str,
                        default = None,
                        help='todo')
    parser.add_argument('-c', '--configfile', type=str,
                        default = None,
                        help='todo')
    # only one argument
    arg = parser.parse_args()
    # todo
    if arg.project == 'twoDSubduction':
        if arg.type == 'parse':
            _TwoDSubduction.Parse(arg.inputfile, arg.outputfile)


if __name__ == "__main__":
    main()