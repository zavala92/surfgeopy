# Sphinx Documentation Configuration File

# This configuration file includes essential options for the Sphinx documentation builder.
# For a comprehensive list of options, refer to the official documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path Setup --------------------------------------------------------------

import os
import sys

# If extensions (or modules for autodoc) are in another directory,
# add them to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute.
sys.path.insert(0, os.path.abspath('..'))


# -- Project Information -----------------------------------------------------

project = 'surfgeopy'
copyright = '2024, Gentian Zavalani'
author = 'Gentian Zavalani'

# Version information, including alpha/beta/rc tags.
release = '1.0.0-alpha'


# -- General Configuration ---------------------------------------------------

# List of Sphinx extension modules.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
]

# Add paths containing templates, relative to this directory.
templates_path = ['_templates']

# Specify the suffix(es) of source filenames.
# You can specify multiple suffixes as a list of strings:
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# Language for content autogenerated by Sphinx.
# This is also used for content translation via gettext catalogs.
language = 'en'

# Ignore patterns for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML Output -------------------------------------------------

# Choose the theme for HTML and HTML Help pages. See the documentation for
# a list of built-in themes.
html_theme = 'sphinx_rtd_theme'

# Theme-specific options to customize the look and feel.
html_theme_options = {
    "collapse_navigation": False,
}

# Paths containing custom static files (e.g., style sheets).
# These files are copied after the built-in static files.
html_static_path = ['_static']

# Additional CSS files.
html_css_files = ['css/custom.css']

# Project logo.
html_logo = '../images/surfgeopy_logo.png'