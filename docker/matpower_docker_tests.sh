#!/bin/bash
octave-cli --eval 'mpver'
octave-cli --eval 'test_matpower && test_most && test_mll_main && test_syngrid'
