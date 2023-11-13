#!/bin/bash

PWM=$1

# Determine highest LOD score, echo consensus

awk -F"\t" '{if ($1>$2 && $1>$3 && $1>$4) print "A"; else if ($2>$1 && $2>$3 && $2>$4) print "C"; else if ($3>$1 && $3>$2 && $3>$4) print "G"; else if ($4>$1 && $4>$2 && $4>$3) print "T"; else print "N"}' ${PWM} | awk '{printf $1}'
