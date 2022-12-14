What's New in MOST 1.2
----------------------

#### Released Dec 13, 2022

Below is a summary of the changes since version 1.1 of MOST. See the
[`CHANGES.md`][1] file for all the gory details. For release notes for
previous versions, see Appendix B of the [MOST User's Manual][2].

#### Changes:
  - Ramping reserves and constraints are now included for the
    transition from the initial state into period 1, except for
    single-period problems.
  - Added calculation of expected TLMP (temporal locational marginal
    price) based on work by Guo, Chen, Tong in [[1](#references),
    [2](#references), [3](#references)]. For generators, these are
    returned in `mdo.results.GenTLMP` and `mdo.results.CondGenTLMP`.
    For storage units they are returned in `mdo.results.StorageTLMPc`,
    `mdo.results.StorageTLMPd`, `mdo.results.CondStorageTLMPc`, and
    `mdo.results.CondStorageTLMPd`. See Table 5-13 in the
    [MOST User's Manual][2].
  - For deterministic cases with storage where `ForceCyclicStorage` is 0,
    ensure that initial storage bounds are equal to initial storage and
    output a warning if they are modified. Fix deterministic UC tests
    where this was causing results to change depending on value of `rho`.

#### Bugs Fixed:
  - Plotting of commitment schedule using `plot_uc()` did not work
    properly in a subplot, such as in `t_most_uc()`. *Thanks to Lim Han.*
  - Fix tests that were failing under Octave 7.x.
  - Fix issue [#29][3] where a typo caused a check on `md.UC.MinDown` > 1
    to be skipped. *Thanks to Talha Iqbal.*

#### Incompatible Changes:
  - Modified definition of ramping reserves for period *t* (and all
    corresponding input and output parameters) to refer to the transition
    from *t*-1 to *t*, not *t* to *t*+1. This means that the ramping
    reserves for the transition into the first period are now optimization
    variables and the corresponding constraints are explicit. This is for
    multiperiod problems only. Ramping reserves and contraints are explicitly
    excluded for single-period problems.  
    *Note:* This change also corrects an error in (4.11), where \gamma^t
    is now correct. Previously it should have been \gamma^{t+1}, as it was
    in the code.

#### References

1. Y. Guo, C. Chen and L. Tong, [“Pricing Multi-Interval Dispatch Under
   Uncertainty Part I: Dispatch-Following Incentives,”][4] *IEEE Transactions
   on Power Systems*, vol. 36, no. 5, pp. 3865–3877, Sept. 2021,
   doi: [10.1109/TPWRS.2021.3055730][4].

2. C. Chen, Y. Guo and L. Tong, [“Pricing Multi-Interval Dispatch Under
   Uncertainty Part II: Generalization and Performance,”][5] *IEEE Transactions
   on Power Systems*, vol. 36, no. 5, pp. 3878–3886, Sept. 2021,
   doi: [10.1109/TPWRS.2020.3045162][5].

3. C. Chen and L. Tong, [“Pricing Real-time Stochastic Storage Operations,”][6]
   *2022 Power Systems Computation Conference (PSCC)*, June 27–July 1, 2022,
   doi: [10.48550/arXiv.2204.08140][6].


[1]: ../../CHANGES.md
[2]: ../MOST-manual.pdf
[3]: https://github.com/MATPOWER/most/issues/29
[4]: https://doi.org/10.1109/TPWRS.2021.3055730
[5]: https://doi.org/10.1109/TPWRS.2020.3045162
[6]: https://doi.org/10.48550/arXiv.2204.08140

