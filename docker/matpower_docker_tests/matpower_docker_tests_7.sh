#!/bin/bash
octave-cli --eval 'mpver'
octave-cli --eval "if strcmp(mpver, '7.0') && have_fcn('octave', 'vnum') >= 6, warning('off', 'Octave:empty-index'); end, test_matpower && test_most && test_mll_main && test_syngrid"
