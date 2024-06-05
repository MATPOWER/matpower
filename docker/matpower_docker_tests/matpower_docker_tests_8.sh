#!/bin/bash
octave-cli --eval 'mpver'
octave-cli --eval "test_matpower && test_most && ~have_feature('mp_core', 0) &&  ~have_feature('sdp_pf', 0) && test_matpower && test_mll_main && test_syngrid"
