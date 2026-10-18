import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'SPECTRA'
copyright = '2026, BioSystems Engineering and control (BiSECt) lab, IIT Madras'
author = 'Pavan Kumar S'

release = '0.1'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
    'sphinx.ext.mathjax',
    'sphinxcontrib.matlab',
    'sphinx.ext.linkcode',
]


def linkcode_resolve(domain, info):
    """Point each documented function at its source on GitHub."""
    filename = info['module'].replace('.', '/')
    return ('https://github.com/bisect-group/SPECTRA/blob/main/'
            + filename + '/' + info['fullname'] + '.m')


napoleon_google_docstring = True

# The MATLAB sources are documented in the COBRA Toolbox house style, whose
# section headers are not the ones napoleon knows. Mapping them to
# params_style is what turns "INPUTS:" and the rest into parameter tables
# rather than loose prose, so the headings in src/*.m have to stay exactly
# as they are written.
napoleon_custom_sections = [
    ('INPUTS', 'params_style'),
    ('INPUT', 'params_style'),
    ('OUTPUTS', 'params_style'),
    ('OUTPUT', 'params_style'),
    ('REQUIRED INPUTS', 'params_style'),
    ('REQUIRED INPUT', 'params_style'),
    ('OPTIONAL INPUTS', 'params_style'),
    ('OPTIONAL INPUT', 'params_style'),
    ('USAGE', 'params_style'),
    'Authors',
]

# -- MATLAB domain -----------------------------------------------------------

matlab_src_dir = os.path.abspath(os.path.join('..'))
primary_domain = 'mat'

add_module_names = False

source_suffix = ['.rst']
source_encoding = 'utf-8'
root_doc = 'index'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'src/depricated']

# src/depricated holds superseded code that is kept for reference only, and
# some of it does not satisfy the matlab domain's file/function name rule.
matlab_src_encoding = 'utf-8'

# -- HTML output -------------------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

epub_show_urls = 'footnote'
