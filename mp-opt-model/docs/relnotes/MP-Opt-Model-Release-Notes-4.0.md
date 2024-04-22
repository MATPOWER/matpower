What's New in MP-Opt-Model 4.0
------------------------------

#### Released October 18, 2021

Below is a summary of the changes since version 3.0 of MP-Opt-Model. See
the [`CHANGES.md`][1] file for all the gory details. For release notes
for previous versions, see Appendix C of the [MP-Opt-Model User's
Manual][2].


#### New Features
  - Support for new class of problems -- parameterized nonlinear equations
    (PNE). Either create a model with only equality constraints (no
    inequalities or costs) and with number of variables equal to 1 more than
    number of constraints, _or_ call `pnes_master()` directly. See Section 4.5
    of User's Manual for details.
    *Thanks to Shrirang Abhyankar and Alexander Flueck for contributions to this
    feature, which is based on the continuation power flow code in MATPOWER 7.1.*
    - Predictor/corrector numerical continuation method for tracing solution
      curves for PNE problems.
    - Plotting of solution curves.
    - User-defined event functions and callback functions.
    - Warm-start capabilities.
  - Optional threshold for detecting failure of LEQ solve, by setting the
    `leq_opt.thresh` option. If the absolute value of any element of the
    solution vector exceeds the threshold, `exitflag` is set to 0, indicating
    failure.
  - New functions:
      - `pnes_master()` provides unified interface for parameterized nonlinear
         equation (PNE) solvers.
      - `pne_callback_default()` collects PNE results and optionally plots
         solution curve.
      - `pne_callback_nose()` handles event signaling a nose point or limit
         has been reached.
      - `pne_callback_target_lam()` handles event signaling a target value
         of parameter &#955; has been reached.
      - `pne_detect_events()` detects events from event function values.
      - `pne_detected_event()` returns detected event details for events
        with a particular name.
      - `pne_event_nose()` detects the limit or nose point.
      - `pne_event_target_lam()` detects a target &#955; value.
      - `pne_pfcn_arc_length()` implements arc length parameterization.
      - `pne_pfcn_natural()` implements natural parameterization.
      - `pne_pfcn_pseudo_arc_length()` implements pseudo arc length
        parameterization.
      - `pne_register_callbacks()` registers callback functions.
      - `pne_register_events()` registers event functions.
      - `mp_idx_manager/set_type_idx_map()` method returns information about
        mapping of indices for a given set type back to the corresponding
        named (and possibly indexed) sets.
      - `mpopt2pneopt()` creates or modifies an options struct for
        `pnes_master()` from a MATPOWER options struct.

#### Bugs Fixed:
  - Calling the `problem_type()` or `is_mixed_integer()` method on an empty
    model no longer causes a fatal error.

#### Other Changes
  - Labels from the `set_types` property are now used as headers for
    `opt_model.display()` to simplify things facilitate use by subclasses.
  - Refactored `describe_idx` into a new method, `set_type_idx_map`, that
    returns in information in a programmatically usable form, and an updated
    `describe_idx` that calls the new method, then formats the results in
    the expected char array(s).


[1]: ../../CHANGES.md
[2]: ../MP-Opt-Model-manual.pdf
