# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

import algait2de
import pygait2d

DOCS_CONF_PATH = os.path.realpath(__file__)
DOCS_DIR = os.path.dirname(DOCS_CONF_PATH)
REPO_DIR = os.path.realpath(os.path.join(DOCS_DIR, '..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Gait2D'
copyright = '2014-2026, Jason K. Moore and Ton van den Bogert'
author = 'Jason K. Moore and Ton van den Bogert'
version = pygait2d.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_gallery.gen_gallery',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

intersphinx_mapping = {
    'sympy': ('https://docs.sympy.org/latest/', None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    'github_repo': 'gait2d',
    'github_type': 'star',
    'github_user': 'csu-hmc',
    'page_width': '1080px',  # 960 doesn't show 79 linewidth examples
}

# Display long function signatures better.
maximum_signature_line_length = 50

# -- sphinx-gallery settings --------------------------------------------------
sphinx_gallery_conf = {
    'examples_dirs': os.path.join(REPO_DIR, 'examples'),
    'gallery_dirs': 'examples',
    'matplotlib_animations': True,
    'remove_config_comments': True,
}
