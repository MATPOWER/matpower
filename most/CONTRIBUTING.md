MOST Contributors Guide
=======================


Getting Involved
----------------

[MOST][1] and [MATPOWER][2] are community efforts and your involvement
is greatly appreciated.
We are always looking for help in identifying and fixing bugs, writing test
cases, improving documentation, answering queries on the mailing lists,
enhancing existing functionality and implementing new features that fit
within the general scope of the project.

Please take a moment to review this document before contributing, as it
helps to communicate your respect for the time and efforts of the MATPOWER
developer community managing and developing this open source project.


Repository Organization
-----------------------

This repository, for the MOST distribution, is organized in the
following directories:
- `docs` - User's Manual
- `lib` - software that implements MOST's core functions

This repository is also included in the core MATPOWER repository as a
subrepository using [`git subrepo`][4].


Mailing Lists
-------------

Your participation is welcomed on the two mailing lists for MATPOWER and MOST
users and developers, [MATPOWER-L][5] and [MATPOWER-DEV-L][6], respectively.
We can always use help in answering questions from newer users.


Reporting a Bug
---------------

The [MOST issue tracker][7] is the preferred channel for reporting bugs or
submitting code changes (pull requests). A good bug report is extremely
helpful and benefits the entire community, so if you find a bug, including
a mistake in the documentation, please report it.

1. **Confirm it is a bug.** You should be able to demonstrate that it
is an error caused by the code in this repository. If it is simply that
you do not understand a result you are getting, ask a question on the
[MATPOWER discussion mailing list][5] instead of submitting an issue
*(after searching [the archives][8] to see if your question has already
been answered, of course)*.

2. **Check if it has already been reported or fixed.** Make sure the bug
still exists by attempting to reproduce it using the latest `master`
branch in the repository. Search the issues on GitHub to make sure it has
not already been reported.

3. **Isolate the problem.** Create a reduced test case with no external
dependencies, if possible, that demonstrates the problem, so others can
easily reproduce it. The simpler the case the better.

4. **Submit an issue that includes a detailed report.** A good bug report
will avoid the need for the developers to track you down for more information.
    - Use an accurate, descriptive title for your issue.
    - Describe your environment - better yet, include the output of `mpver`
      on your system.
    - Include a few lines of code or attach a script file that reproduces
      the bug.
    - Describe the result you got and what you expected.
    - Select "bug" under "Labels".

_**Note:** Bugs or issues related to core MATPOWER functions should be
submitted to the main [MATPOWER issue tracker][9]._


Reviewing Issues
----------------

If you see an issue or bug report submitted by someone else, please consider
helping out by attempting to reproduce the issue on your system. Report your
experience in a comment on the issue. Even if you are not comfortable
submitting a patch to fix it, any information you can add to help the
developers locate a solution is greatly appreciated. Simply leave your
comments on the issue. The same goes for review of pull requests, which
are discussed below.


Submitting Additions or Modifications to the Code
-------------------------------------------------

Code contributions are a great help and are always welcome. This includes
bug fixes, enhancements to existing functionality, new features or tests,
and even edits to the included documentation. It is always a good idea to
discuss your ideas first on the [MATPOWER developer mailing list][6], especially
for larger or more complex contributions.

Contributions should be submitted as pull requests, as described below.
Before submitting your pull request, please make sure you have tested
your changes and that they follow the [MATPOWER developer guidelines][10].

### Getting started

#### Step 1 : [Set Up Git][11].

Make sure git knows your name and email address:

  ```bash
  git config --global user.name "Random Citizen"
  git config --global user.email "random.citizen@example.com"
  ```

#### Step 2 : [Fork][12] the repository.

Click "Fork" on the [repository page][1] on GitHub to create your own fork
of the project.

#### Step 3 : Clone

Check out your copy locally and configure the remotes:

  ```bash
  # clone your fork of the repo into the current directory
  git clone https://github.com/<your-username>/most.git
  # go to the newly cloned directory
  cd most
  # assign the original repo to a remote called "upstream"
  git remote add upstream https://github.com/MATPOWER/most.git
  ```

### Making your changes

#### Step 4 : Update

If it has been a while since you first made your clone, get all of the
latest changes from the upstream repository:

  ```bash
  git checkout master
  git pull upstream master
  ```

#### Step 5 : Branch

Create a new topic branch where you can work on your new feature, change
or fix. Always create it from an up-to-date `master` branch:

  ```bash
  git checkout -b <topic-branch-name>
  ```

#### Step 6 : Commit

Edit your local copy of the files to implement your feature, change or fix.
Then commit your changes in logical chunks. Do not combine multiple logical
changes in a single commit. And please adhere to the guidlines below for
commit messages.

Add and commit:

  ```bash
  git add my/changed/files
  git commit
  ```

Writing good commit messages is important. A commit message should describe
what changed and why. Follow these guidelines when writing one:

1. The first line should be 50 characters or less and contain a short
   description of the change. Begin with a capitalized imperative verb,
   for example, "Fix issue #4", not "Fixed issue #4" or "Fixes issue #4".
   All other words in this description should be in lowercase with the
   exception of proper nouns, acronyms, and references to code, such as
   function or variable names.
2. Keep the second line blank.
3. Wrap all other lines at 72 columns.

If your patch fixes an open issue, the issue should be referenced in the
first line of the commit message with the issue number preceded by a hash
symbol (`#`), e.g. `Fix #4, Q limit violations` and at the end of the
message with the full URL. Use the `Fixes:` prefix for bug fixes. For
other references use `Refs:`. For example, a good commit message might
look something like:

  ```text
  Fix issue #4, Q limit violations in CPF.
  
  More detailed explanatory text, if necessary.  Wrap it to about 72
  characters or so.  In some contexts, the first line is treated as the
  subject of an email and the rest of the text as the body.  The blank
  line separating the summary from the body is critical (unless you omit
  the body entirely); tools like rebase can get confused if you run the
  two together.

  Further paragraphs come after blank lines.

  - Bullet points are okay, too.

  - Typically a hyphen or asterisk is used for the bullet, followed by a
    single space, with blank lines in between, but conventions vary here.

  - Use a hanging indent

  Fixes: https://github.com/MATPOWER/most/issues/4
  Refs: http://www.mail-archive.com/matpower-l@cornell.edu/msg05557.html
  Refs: https://github.com/MATPOWER/most/pull/5
  ```

#### Step 7 : Test

Bug fixes and features **should come with tests**, either added to the
appropriate existing test function in `lib/t`, or in a new test function
with a descriptive name beginning with `t_`, in which case it should also
be added to `test_most.m`. See the documentation for [MP-Test][3] and
the existing MOST test files (e.g. [`t_most_sp`](lib/t/t_most_sp.m)) for
examples of how to write tests.

You can run your tests by typing the name of your test function at the
MATLAB or Octave prompt, or `test_most` to run the entire test suite.

### Sharing your changes

#### Step 8 : Rebase

Use `git rebase` (not `git merge`) to sync your work with the upstream
development branch from time to time, and especially before pushing.

  ```bash
  git fetch upstream
  git rebase upstream/master
  ```

And use Git's [interactive rebase][13] feature to tidy up your commits, if
necessary, *before* making them public. See [this article][14] for some
helpful background on `git rebase` vs. `git merge`.

#### Step 9 : Push

Push your topic branch up to your fork on GitHub.

  ```bash
  git push origin <topic-branch-name>
  ```

#### Step 10 : [Open a Pull Request][15]

Open a pull request against the `master` branch, by going to the GitHub page
for your fork (`https://github.com/<your-username>/most`), and selecting
your topic branch. Click the "Pull Request" button and fill out the form
using a clear, accurate title and description.

_**IMPORTANT:** By submitting a pull request, you represent that you have
the right to license your contribution to PSERC and the community, and
agree by submitting the patch that your contributions are licensed under
the [3-clause BSD license][16]._

#### Step 11 : Discuss and update

You will probably get feedback or requests for changes to your pull request.
This is a normal part of the submission process, so don't be surprised or
discouraged.

To make changes to an existing pull request, make the changes to your branch.
When you push that branch to your fork, GitHub will automatically update the
pull request. Each time the pull request is updated it triggers a [CI
(continuous integration) test run][17] which results in a check mark or an
X next to the commit on the pull request page and next to the pull request
name in the issue tracker, indicating whether or not all tests passed.

After your pull request has been reviewed and approved by a MATPOWER
Collaborator, it can be merged into the upstream MOST repository. Once
this has happened the pull request is closed and your contributions are
part of the project and available to the world! Congratulations and many
thanks!


### Git Workflow and Branching Model

We use [GitHub Flow][18], meaning that the `master` branch should always
be ready for release and all new work is done in descriptively named
branches off of `master`, which are then reviewed via pull requests. The
only addition is that we do still have a `release` branch that always
points to the latest versioned release. We also use tags like `6.0` to
tag each release.


Licensing
---------

By submitting a pull request, you represent that you have the right to
license your contribution to PSERC and the community, and agree by
submitting the patch that your contributions are licensed under the
[3-clause BSD license][16].


Thanks
------

We would like to express our appreciation to everyone who has contributed to
the MOST project and community for helping to make it a useful tool for
the power systems community. **Thank you!**

--

[1]: https://github.com/MATPOWER/most
[2]: https://github.com/MATPOWER/matpower
[3]: https://github.com/MATPOWER/mptest
[4]: https://github.com/ingydotnet/git-subrepo
[5]: http://www.pserc.cornell.edu/matpower/mailinglists.html#discusslist
[6]: http://www.pserc.cornell.edu/matpower/mailinglists.html#devlist
[7]: https://github.com/MATPOWER/most/issues
[8]: http://www.mail-archive.com/matpower-l@cornell.edu/
[9]: https://github.com/MATPOWER/matpower/issues
[10]: https://github.com/MATPOWER/matpower/blob/master/docs/MATPOWER-dev-guide.md
[11]: https://help.github.com/articles/set-up-git/
[12]: https://help.github.com/articles/fork-a-repo/
[13]: https://help.github.com/articles/interactive-rebase
[14]: https://medium.com/@porteneuve/getting-solid-at-git-rebase-vs-merge-4fa1a48c53aa#.7gmldhj6m
[15]: https://help.github.com/articles/about-pull-requests/
[16]: LICENSE
[17]: https://travis-ci.org/MATPOWER/most
[18]: http://scottchacon.com/2011/08/31/github-flow.html
