#!/usr/bin/env python
# coding: utf-8

# An example jupyter notebook, which includes how to read gene expression files and make a plot using plotly.

# In[2]:


# Set up the environment.
import pandas, numpy, os
import plotly
import plotly.graph_objs as go


# In[9]:


# Read input files and check some key metrics. These input files are downloaded from stemformatics.org:
#
# http://data.stemformatics.org/files/notta_expression.tsv
# http://data.stemformatics.org/files/notta_samples.tsv
#
# Then saved under received/Stemformatics/
#
def readDataFiles():
    # Read expression and samples.
    df = pandas.read_csv("../received/Stemformatics/notta_expression.tsv", sep="\t", index_col=0)
    samples = pandas.read_csv("../received/Stemformatics/notta_samples.tsv", sep="\t", index_col=0)
        
    # Show first 5 lines of these matrices
    display(df.head())
    display(samples.head())
    
    # Min/max values (shows that this matrix must have been logged already).
    print("min", df.min().min(), "max", df.max().max())
    
    # Check that columns of df match row ids of samples
    print(set(df.columns)==set(samples.index))
    
    return df, samples

exp,samples = readDataFiles()


# In[10]:


# Make a plot of library size
def libPlot():
    traces = [go.Bar(x=exp.columns, y=exp.sum(), name="library size")]
    go.Figure(data=traces, layout=dict(width=900, height=400, template="plotly_white")).show()
    
libPlot()


# In[3]:


# Show conda environment which was used to run this notebook, as well as versions of key packages
get_ipython().system('conda info')
print("pandas version", pandas.__version__)
print("plotly version", plotly.__version__)

