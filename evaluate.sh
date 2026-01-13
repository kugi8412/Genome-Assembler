#!/bin/bash

bowtie2 -a --local --mp 2,2 --rdg 10,2 --rfg 10,2 -f -x ./training/reference/reference -U $1 | ./evaluate.py
