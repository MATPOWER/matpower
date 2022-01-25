#!/bin/bash
octave-cli --eval 'mpver'
octave-cli --eval 'test_matpower'
octave-cli --eval 'test_most'
