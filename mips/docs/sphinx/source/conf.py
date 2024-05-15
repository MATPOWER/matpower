# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sphinx_rtd_theme
import re


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MIPS'
copyright = '2009-2024, Power Systems Engineering Research Center (PSERC)'
author = 'Ray D. Zimmerman, Hongye Wang'

# The full version, including alpha/beta/rc tags
release = '1.5.1'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinxcontrib.matlab',
    'sphinx.ext.autodoc',
    'sphinx.ext.extlinks',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
    'sphinx_tabs.tabs',
]
matlab_src_dir = 'matlab-source'
primary_domain = 'mat'

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['css/matpower.css']


# -- Options for Sphinx RTD Theme ------------------------------------------------
# See https://pypi.org/project/sphinx-rtd-theme/
#     https://github.com/readthedocs/sphinx_rtd_theme
#     https://sphinx-rtd-theme.readthedocs.io/
#     https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'both',
    'style_nav_header_background': '#568085',   # medium dark teal
}


# -- Other Options -----------------------------------------------------------

matlab_short_links = True
# matlab_auto_link = "basic"
matlab_auto_link = "all"
matlab_show_property_default_value = True
# autoclass_content = 'both'         # 'class', 'init', 'both'
autodoc_member_order = 'bysource'   # 'alphabetical', 'groupwise', 'bysource'
napoleon_use_param = False
napoleon_use_rtype = False
napoleon_custom_sections = [
    ('Input', 'params_style'),
    ('Inputs', 'params_style'),
    ('Output', 'params_style'),
    ('Outputs', 'params_style'),
]

rst_prolog = """
.. include:: /mp-docs-shared/prolog.rst.txt
"""
extlinks = {
#     'mostman': ('https://matpower.org/docs/MOST-manual-%s.pdf', 'MOST User''s Manual'),
    'doi': ('https://doi.org/%s', 'doi: %s'),
}

numfig = True               # enable numbered references
numfig_format = {
    'figure': 'Figure %s',
    'code-block': 'Listing %s',
    'section': 'Section %s' }
highlight_language = 'octave'   # 'matlab' has issues with '...'
                                # e.g. y = fcn(x, ...)
math_number_all = True
math_eqref_format = '({number})'


# -- MathJax Macro Setup -----------------------------------------------------
# -- based on https://stackoverflow.com/questions/9728292/creating-latex-math-macros-within-sphinx#comment124804327_60497853
mathjax3_config = {
    'loader': {'load': ['[tex]/upgreek']},
    'tex': {'packages': {'[+]': ['upgreek']},
            'macros': {}}   # create empty
}

with open('mp-docs-shared/mathCmds.tex.txt', 'r') as f:
    for line in f:
        macros = re.findall(r'\\(re)?newcommand{\\(.*?)}(\[(\d)\])?{(.+)}', line)
        for macro in macros:
            if len(macro[2]) == 0:
                mathjax3_config['tex']['macros'][macro[1]] = macro[4]
            else:
                mathjax3_config['tex']['macros'][macro[1]] = [macro[4], int(macro[3])]
