#!/bin/bash
./promolden filename.inp
gprof ./promolden | ./gprof2dot.py | dot -Tpng -o filename.png
