MP-Element
==========

MP-Element is a new, generalized network and element modeling layer for
[MATPOWER][1]. It is currently **under active development** in its own
separate repository, with the intention of eventually being incorporated
into the core MATPOWER code.


System Requirements
-------------------
*   [MATLAB][2] version 9.1 (R2016b) or later, or
*   [GNU Octave][3] version 6.2 or later
*   A version of [MATPOWER][1] that explicitly includes support for
    MP-Element  
    _The `master` branch of the [MATPOWER GitHub repository][1] now includes
    the required support. See the [CHANGES][8] file._


Installation
------------

Installation and use of MP-Element requires familiarity with the basic operation
of MATLAB or Octave, including setting up your MATLAB/Octave path.

1.  Clone the repository or download and extract the zip file of the MP-Element
    distribution from the [MP-Element project page][4] to the location of your
    choice.

2.  Either ...
  - Move the resulting `mp-element` directory to the directory
    containing MATPOWER (must be a version that supports MP-Element),
    and re-run `install_matpower` (the directory must be named
    `mp-element` for the installer to recognize its presence),  
  ... _or_ ...
  - Add the following directories to your MATLAB/Octave path:
    * `<MPELEMENT>/lib`
    * `<MPELEMENT>/lib/t`

    where `<MPELEMENT>` is used to denote the path to the `mp-element`
    directory you cloned or unzipped in step 1.

3.  At the MATLAB/Octave prompt, type `test_mp_element` to run the test
    suite and verify that MP-Element is properly installed and functioning.
    The result should resemble the following:
```
  >> test_mp_element
  t_mp_mapped_array.........ok
  t_mp_table................ok
  t_mp_data_model...........ok
  t_mp_dm_converter_mpc2....ok
  t_nm_element..............ok
  t_port_inj_current_acc....ok
  t_port_inj_current_acp....ok
  t_port_inj_power_acc......ok
  t_port_inj_power_acp......ok
  t_node_test...............ok
  t_run_mp..................ok
  t_run_mp_3p...............ok
  All tests successful (2075 of 2075)
  Elapsed time 29.94 seconds.
```


Getting Started
---------------

#### Default Behavior

With MP-Element installed, its modeling is used by default by MATPOWER's
`runpf`, `runcpf` and `runopf` for the following:
  - DC power flow
  - DC optimal power flow
  - AC power flow for all except radial and hybrid Newton-Raphson
    formulations/solvers
    - a new `'FSOLVE'` option based on `fsolve()` function is available for
      AC power flow
  - AC continuation power flow
  - AC OPF for solvers MIPS, `fmincon`, IPOPT, and Artelys Knitro, for
    all formulations

MP-Element modeling can be turned off in favor of the legacy MATPOWER modeling
with `have_feature('mp_element', 0)`.

Note: The MP-Opt-Model object is used for power flow and continuation power
flow as well as OPF and is added as `om` field to power flow and CPF `results`
struct.

#### New Functions

In addition to the old `runpf`, `runcpf` and `runopf` functions which are
backward compatible when using the MP-Element modeling, there are 3 new
corresponding functions, namely `run_pf`, `run_cpf` and `run_opf`. They are
simple wrappers around a new _(unfinished)_ function, `run_mp`, which does
not convert the MATPOWER case struct to internal numbering before use, and
will eventually provide a much more flexible environment for customization.


Documentation
-------------

Given that this code is **under active development**, the documentation is
very incomplete and sometimes missing entirely (i.e. the code _is_ the
documentation).

#### Technical Note

The following MATPOWER Technical Note describes some of the design of
MP-Element.

- R. D. Zimmerman, ["MP-Element: A Unified MATPOWER Element Model, with
  Corresponding Functions and Derivatives,"][5] _MATPOWER Technical Note 5_,
  October 2020.  
  Available: [https://matpower.org/docs/TN5-MP-Element.pdf][5]  
  doi: [10.5281/zenodo.3237866][6].


License
-------

MP-Element is distributed under the [3-clause BSD license][7].


Acknowledgments
---------------

This material is based upon work supported in part by the National Science
Foundation under Grant Nos. 1642341 and 1931421. Any opinions, findings, and
conclusions or recommendations expressed in this material are those of the
author(s) and do not necessarily reflect the views of the funding agencies.

----
[1]: https://github.com/MATPOWER/matpower
[2]: https://www.mathworks.com/
[3]: https://www.gnu.org/software/octave/
[4]: https://github.com/MATPOWER/mp-element
[5]: https://matpower.org/docs/TN5-MP-Element.pdf
[6]: https://doi.org/10.5281/zenodo.4110676
[7]: LICENSE
[8]: https://github.com/MATPOWER/matpower/blob/master/CHANGES.md
