MATPOWER Developer's Guide
--------------------------

Consider this to be a bare-bones, rough draft of an evolving
document. I'm sure there are many helpful things that are missing.
Please feel free to discuss additions, modifications, etc. on the
[MATPOWER-DEV-L email list][1].

One of the primary goals for the MATPOWER code is to **keep things
as simple as possible**. That is, we want to keep the code as simple
to understand and maintain as possible without sacrificing
performance.

The following are some guidelines, in no particular order.

- All functions should include a help section that fully describes
  the inputs, outputs and what the function does. You can use the
  help section from [`runpf`](../lib/runpf) as a template.

- Following the help section, should be a blank line, followed by the
  copyright and license lines, like:
     ```matlab
     %   MATPOWER
     %   Copyright (c) 2016-2017, Power Systems Engineering Research Center (PSERC)
     %   by <Your Name>, <Your Affiliation>
     %
     %   This file is part of MATPOWER.
     %   Covered by the 3-clause BSD License (see LICENSE file for details).
     %   See http://www.pserc.cornell.edu/matpower/ for more info.
     ```

- Write tests to cover all functionality.

- Use the named constants defined by `idx_bus`, `idx_gen`, etc. rather
  than the numerical values.

- As much as possible, for each function the input arguments should be
  limited to what is needed to do the job. The output arguments should
  only include things the function creates or modifies.

- All input and output arguments **must** be well defined and documented.
  **No** using structs (or any other data structure) to allow passing of
  arbitrary, undocumented data.
  
- Every piece of data **must** be well-defined somewhere.

- Use consistent variable names. That is, use the same name to refer to
  the same piece of data in all contexts, unless there is a good reason
  to use a different name.

- Use ...
    - 4 spaces to indent (no tabs)
    - a single space but no parens following `if`, `for`, `while`, etc,
      except where necessary
    - no `end` or `return` at the end of functions
    - spaces after commas in argument lists and index lists for non-trivial
      array indexing
        - e.g `matrix(r,c)` is ok
        - but use a space for `matrix(row+3, col-7)`
    - lower case function names, with `_` to separate words
    - lower case variable names
    - upper case constant names
    - no `return` in the middle of a function
    - comments:
        - `%% descriptive text` to document the code
        - `% ` to temporarily comment out unused code
        - `%%-----  heading text  -----%%` for section headings in the code
        - *Why? For consistency ... once upon a time I began doing it
           that way for some reason which no longer matters.*

- Fatal errors should generally be reserved for invalid inputs or other
  misuse of a function.


[1]: http://www.pserc.cornell.edu/matpower/mailinglists.html#devlist
