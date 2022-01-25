#!/bin/bash
octave-cli --eval 'mpver'
octave-cli --eval "warning('off', 'Octave:singular-matrix'); test_matpower"
