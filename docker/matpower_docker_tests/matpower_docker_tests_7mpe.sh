#!/bin/bash
octave-cli --eval 'mpver'
octave-cli --eval "f = @mp_table; test_matpower && test_most && ~have_feature('mp_element', 0) && test_mll_main && test_syngrid && have_feature('mp_element', 1) && test_mp_element"
