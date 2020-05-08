Change history for MP-Opt-Model
===============================


Version 1.0 - *May 8, 2020*
---------------------------

#### 5/8/20
  - Release 1.0.

#### 5/7/20
  - Add MP-Opt-Model User's Manual in `docs`, with LaTeX source in
    `docs/src`.
  - Add usage examples to `README.md`.

#### 4/30/20
  - Add `README.md`, `CHANGES.md`, `AUTHORS`, `CONTRIBUTING.md`, `LICENSE`.
  - Refactor `opt_model` class to inherit from new abstract base class
    `mp_idx_manager`, which can be used to manage the indexing of other sets of
    parameters, etc. in other contexts.


Version 0.8 - *Apr 29, 2020*
----------------------------

#### 4/29/20
  - Version 0.8.
  - Add `mpomver()` to define MP-Opt-Model version information.
  - **INCOMPATIBLE CHANGE:** Renamed the following functions and
    modified the order of their input args so that the MP-Opt-Model
    object appears first. Ideally, these would be defined as methods
    of the `opt_model` class, but Octave 4.2 and earlier is not
    able to find them via a function handle (as used in the `solve()`
    method) if they are inherited by a sub-class.
    - `opf_consfcn()` --> `nlp_consfcn()`
    - `opf_costfcn()` --> `nlp_costfcn()`
    - `opf_hessfcn()` --> `nlp_hessfcn()`
  - Add GitHub Actions CI workflow and [Travis-CI][3] configuration.
  - Add `test_mp_opt_model()` to run all tests.
  - Remove MATPOWER dependencies.
  - Move code related to solver interfaces, `opt_model` and a
    few other functions like `have_fcn()` and `nested_struct_copy()`
    from main [MATPOWER][1] repository to new [MP-Opt-Model][2]
    repository.


⬆ _In [MP-Opt-Model][2] repository_ ⬆

-----------------------------------

⬇ _In [MATPOWER][1] repository_ ⬇


#### 4/28/20
  - Move deprecated `opt_model` methods and code related to legacy
    user-defined OPF costs from `@opt_model` to `@opf_model`.
  - **INCOMPATIBLE CHANGE:** Modify order of default output arguments of
    `opt_model/get_idx()` (again), removing the one related to legacy
    costs.

#### 3/18/20
  - Add `nlpopf_solver()` based on the new `solver()` method of
    `opt_model`. This single function replaces `mipsopf_solver()`,
    `fmincopf_solver()`, `ipoptopf_solver()`, and `ktropf_solver()`.
  - Convert `dcopf_solver()` to use the new `solver()` method of
    `opt_model` instead of calling `qps_matpower()` manually.
  - Add new top-level wrapper function `nlps_matpower()` to provide
    a standard interface for MATPOWER's nonlinear program (NLP)
    solvers (`fmincon`, IPOPT, Artelys Knitro, and MIPS), with
    calling syntax similar to `mips()`. It includes the ability to
    pass in solver-specific input options.
  - Add `nlps_fmincon()`, `nlps_ipopt()` and `nlps_knitro()`, with
    interface that matches `nlps_matpower()` to handle implementation
    for `fmincon`, IPOPT, and Artelys Knitro solvers, respectively.
  - Add `mpopt2nlpopt()` to set up an options struct for
    `nlps_matpower()` based on a MATPOWER options struct.
  - Add three new methods to `opt_model` class:
    - `is_mixed_integer()` - returns true if the model includes any binary
      or integer variables
    - `problem_type()` - returns one of the following strings, based on
      the characteristics of the variables, costs and constraints in the
      model:
      - `'NLP'` - nonlinear program
      - `'LP'` - linear program
      - `'QP'` - quadratic program
      - `'MILP'` - mixed-integer linear program
      - `'MIQP'` - mixed-integer quadratic program
    - `solve()` - solves the model using `qps_matpower()`,
      `miqps_matpower()`, or `nlps_matpower()`, depending on the problem
      type (`'MINLP'` problems are not yet implemented)

#### 3/12/20
  - Fix bug in `ktropf_solver()` where Artelys Knitro was still using
    `fmincon` options.

#### 3/6/20
  - Fix issue with missing objective function value from `miqps_mosek()`
    and `qps_mosek()` when return status is "Stalled at or near optimal
    solution."

#### 3/4/20
  - Remove unused input arguments from `opf_consfcn()` and `opf_hessfcn()`.

#### 2/27/20
  - Add `copy()` method to `opt_model` class to get around issues
    with inheritance in constructors that was preventing copy constructor
    from working in Octave 5.2 and earlier (see also [Octave bug
    52614](https://savannah.gnu.org/bugs/?52614).

#### 2/26/20
  - Significant performance improvement for CPLEX on small problems by
    eliminating call to `cplexoptimset()`, which was a huge bottleneck.
  - Fix CPLEX 12.10 compatibility [issue #90][6].

#### 2/18/20
  - Artelys Knitro 12.1 compatibility fix.

#### 8/15/19
  - Improve performance of `opt_model/add_named_set()`.
    (See [issue #79][5].)
    *Thanks to Baraa Mohandes.*
  - Refactor code in `opt_model/params_lin_constraint()` and
    `opt_model/params_quad_cost()` to speed up sparse matrix construction
    when there are lots of constraint or cost sets. Results in significant
    speedups for some problems during problem setup in MOST.
    (See [pull request #70][4].)
    *Thanks to Daniel Muldrew.*


Version 0.7.0 - *Jun 20, 2019*
------------------------------

#### 6/20/19
  - This change history begins with the code that was part of the
    MATPOWER 7.0 release, which is tagged as version 0.7.0 in the
    MP-Opt-Model repository.

----
[1]: https://github.com/MATPOWER/matpower
[2]: https://github.com/MATPOWER/mp-opt-model
[3]: https://travis-ci.org
[4]: https://github.com/MATPOWER/matpower/pull/70
[5]: https://github.com/MATPOWER/matpower/issues/79
[6]: https://github.com/MATPOWER/matpower/issues/90
