#! /usr/bin/env python

import os
import sys

BIN_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.dirname(BIN_DIR)
SIM_DIR = os.path.join(PROJECT_DIR, 'simulations')
VAL_DIR = os.path.join(SIM_DIR, 'validation')
CONFIG_DIR = os.path.join(PROJECT_DIR, 'configs')

def main():
    sys.stdout.write("{0}\n".format(PROJECT_DIR))

if __name__ == '__main__':
    main()

