import pandas as pd
import numpy as np

import requests, json

username = 'pwangel@unimelb.edu.au'
password = 'YcY64N'

base_URL = 'https://api.stemformatics.org/'
response = requests.get(base_URL, auth=(username, password))

print(response.json())

meta_URL = 'https://api.stemformatics.org/samples / metadata/7284/metadata'
response = requests.get(meta_URL, auth=(username, password))

if response.status_code != 200:
    pass

meta_URL = 'https://api.stemformatics.org/samples / metadata/7284/metadata'
response = requests.get(meta_URL, auth=(username, password))

#example 
json.loads(response.json())[0]['dataset_id']  