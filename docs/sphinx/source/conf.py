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
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sphinx_rtd_theme
import re


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MATPOWER Documentation'
copyright = '1996-2024, Power Systems Engineering Research Center (PSERC)'
author = 'Ray D. Zimmerman, Carlos E. Murillo-Sánchez, Hongye Wang, et. al.'

# The full version, including alpha/beta/rc tags
release = '8.0b1'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinxcontrib.matlab',
    'sphinx.ext.autodoc',
#    'sphinx.ext.autosectionlabel', # (links break if you edit a section title)
    'sphinx.ext.extlinks',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
    'sphinx_tabs.tabs',
#    'sphinx.ext.autosummary',  # seems work only with Python since it attempts
                                # to load the module it's going to summarize
]
matlab_src_dir = 'matlab-source'
primary_domain = 'mat'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'basic'
#html_theme = 'alabaster'
#html_theme = 'classic'      # ugly
#html_theme = 'traditional'  # ugly
#html_theme = 'haiku'
#html_theme = 'nature'       # nice
#html_theme = 'agogo'
#html_theme = 'scrolls'
#html_theme = 'sphinxdoc'
#html_theme = 'pyramid'      # ok
#html_theme = 'bizstyle'     # nice
html_theme = 'sphinx_rtd_theme' # favorie, installed via pip
#html_theme = 'epub'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'css/matpower.css',
]

# -- Options for Sphinx RTD Theme ------------------------------------------------
# See https://pypi.org/project/sphinx-rtd-theme/
#     https://github.com/readthedocs/sphinx_rtd_theme
#     https://sphinx-rtd-theme.readthedocs.io/
#     https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html

#html_logo = 'MATPOWER-md.png'
html_theme_options = {
#     'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
#     'analytics_anonymize_ip': False,
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'both',
#    'style_nav_header_background': 'white',   # white
#    'style_nav_header_background': '#faf1d9',   # light yellow
#    'style_nav_header_background': '#013f53',   # dark teal
    'style_nav_header_background': '#568085',   # medium dark teal
#   'style_nav_header_background': '#8eaca6',   # medium teal
#    'style_nav_header_background': '#d2e9e1',   # light teal
#     'style_external_links': False,
#     'vcs_pageview_mode': '',
#     'style_nav_header_background': 'white',
#     # Toc options
#     'collapse_navigation': True,
#     'sticky_navigation': True,
#     'navigation_depth': 4,
#     'includehidden': True,
#     'titles_only': False
}


# -- Options for LaTeX output ------------------------------------------------

latex_logo = 'MATPOWER-md.png'
latex_documents = [
#     (   'users-manual/index',
#         'matpower_users_manual.tex',
#         '\\matpower{} \\textrm{User\'s Manual}',
#         'Ray D. Zimmerman, Carlos E. Murillo-Sánchez, Hongye Wang, et. al.',
#         'manual'),
    (   'dev-manual/index',
        'matpower_dev_manual.tex',
        '\\matpower{} \\textrm{Developers\'s Manual}',
        'Ray D. Zimmerman, et. al.',\
        'manual'),
    (   'ref-manual/index',
        'matpower_ref_manual.tex',
        '\\matpower{} \\textrm{Reference Manual}',
        'Ray D. Zimmerman, et. al.',\
        'manual'),
#     (   'users-manual-legacy/index',
#         'matpower_users_manual_legacy.tex',
#         '\\matpower{} \\textrm{User\'s Manual (Legacy)}',
#         'Ray D. Zimmerman, Carlos E. Murillo-Sánchez',
#         'manual'),
]
latex_additional_files = ["mp-docs-shared/preamble.tex.txt", "mp-docs-shared/mathCmds.tex.txt"]
latex_show_pagerefs = True          # False is default
latex_elements = {
    'preamble': r'\input{preamble.tex.txt}',
    'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    'printindex': r'\footnotesize\raggedright\printindex',
    'maxlistdepth': '8',
    'extraclassoptions': 'openany',     # chapters can start on any page
    'maketitle': r'\date{May 7, 2024}\newcommand\sphinxbackoftitlepage{\begin{center}{~\vfill\copyright~1996--2024~\PSERC\\All Rights Reserved}\end{center}}\sphinxmaketitle',
    'extrapackages': r"""\usepackage{bm}
\usepackage{upgreek}
\usepackage{amsfonts}
\usepackage{amstext}
\usepackage{soul}
\usepackage{relsize}
""", 
}
# latex_show_urls = 'footnote'


# -- Other Options -----------------------------------------------------------

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = False
# matlab_keep_package_prefix = False
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

# Make sure autosection targets are unique
#autosectionlabel_prefix_document = True

numfig = True               # enable numbered references
numfig_format = {
    'figure': 'Figure %s',
    'code-block': 'Listing %s',
    'section': 'Section %s' }
# numfig_secnum_depth = 2
# highlight_language = 'matlab'
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
#                mathjax3_config['tex']['macros'][macro[1]] = "{"+macro[4]+"}"
            else:
                mathjax3_config['tex']['macros'][macro[1]] = [macro[4], int(macro[3])]
#                mathjax3_config['tex']['macros'][macro[1]] = ["{"+macro[4]+"}", int(macro[3])]

# print("MathJax Macros: ", mathjax3_config['tex']['macros'])
