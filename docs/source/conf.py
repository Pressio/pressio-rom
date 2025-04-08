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
import datetime
import os
import sys
sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------
# The master toctree document.
master_doc = "index"

project = 'Pressio-ROM'
# copyright = f'{datetime.datetime.now().year}, 2021, National Technology & Engineering Solutions of Sandia, LLC (NTESS)'
copyright = u"2021, National Technology & Engineering Solutions of Sandia, LLC (NTESS)"
author = 'Francesco Rizzi'

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
def get_version():
    version_splits = []
    with open("../../include/pressio/pressio_macros.hpp") as version_file:
        for line in version_file:
            splits = line.strip().split()
            if len(splits) == 3 and "VERSION" in splits[1]:
                version_splits.append(splits[2])
    return ".".join(version_splits)

# The full version, including alpha/beta/rc tags
release = get_version()

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
  "myst_parser", \
  "sphinx.ext.autodoc", \
  "sphinx.ext.viewcode", \
  "sphinx.ext.intersphinx", \
  "sphinx_copybutton"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = "%B %d, %Y"

# -- Options for HTML output -------------------------------------------------
# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "tango"
pygments_dark_style = "monokai"
#from pygments.styles import get_all_styles
#styles = list(get_all_styles())
#print(styles)

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_title = project + "-"+release

html_theme = 'furo'

html_theme_options = {
  "sidebar_hide_name": False,
  "light_css_variables": {
    "color-brand-primary": "#336790",  # blue
    "color-brand-content": "#336790"
  },
  "dark_css_variables": {
    "color-brand-primary": "#E5B62F",  # yellow1
    "color-brand-content": "#E5B62F"
    #"color-brand-primary": "#F4EB4D",  # yellow2
    #"color-brand-content": "#F4EB4D"
    #"color-brand-primary": "#F39C12",  # orange
    #"color-brand-content": "#F39C12",
    #"color-brand-primary": "#00FF40",  # lime
    #"color-brand-content": "#00FF40"
    #"color-brand-primary": "#6ed0f5",  # cyan
    #"color-brand-content": "#6ed0f5"
  },
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = ['hacks.css']

source_suffix = {
  '.rst': 'restructuredtext',
  '.md': 'markdown',
}

# If not '', a "Last updated on:" timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = "%b %d, %Y"

# need to figure out why this does not work
# rst_prolog = """
# .. include:: special.rst
# """

html_sidebars = {
    "**": [
        "sidebar/scroll-start.html",
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
    ]
}

html_logo = "_static/logo.png"
