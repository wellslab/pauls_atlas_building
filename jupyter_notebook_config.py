"""
This file runs automatically after saving a jupyter notebook, if placed in 
~username/.jupyter/jupyter_notebook_config.py. The code here can be used to run conversion of the notebook
to html file for example. Modify the parameters below to suit your needs.

Adapted from: https://www.svds.com/jupyter-notebook-best-practices-for-data-science/
"""

import os
from subprocess import check_call

# Customise the following parameters to suit your needs. --------------------------------

# This is the subfolder (under the current notebook) where converted files will be placed. 
subfolder = "output"

# These are the converted file formats. Set to True for format you want to save.
formats = {"script":True, "html":True, "markdown":False}

# If you have plotly plots in the notebook, saving as html requires plotly javascript library, which
# isn't included with the html file by default. So copy this file into the subfolder also.
# (Adding a reference to the cdn in the html file didn't seem to work - requires local physical copy of .js file.)
addPlotlyFile = True

# ---------------------------------------------------------------------------------------

def post_save(model, os_path, contents_manager):
    """post-save hook for converting notebooks.
    """
    if model['type'] != 'notebook':
        return # only do this for notebooks

    d, fname = os.path.split(os_path)  # get current directory and filename

    # save all output_files in subfolder
    subfolder_path = os.path.join(d, subfolder)
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)    
    
    for key,value in formats.items():
        if value:
            check_call(['jupyter', 'nbconvert', '--output-dir', subfolder_path, '--to', key, fname], cwd=d)

    if addPlotlyFile:  # add plotly.js file in the subfolder
        plotlyFile = os.path.join(subfolder_path, "plotly.js")
        if not os.path.exists(plotlyFile):
            import requests
            r = requests.get("https://cdn.plot.ly/plotly-latest.min.js", allow_redirects=True)
            f = open(plotlyFile, 'wb')
            f.write(r.content)
            f.close()

c.FileContentsManager.post_save_hook = post_save
