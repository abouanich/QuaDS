  
import yaml
with open("config_file.yml", "r") as yamlfile:
    config = yaml.load(yamlfile, Loader=yaml.FullLoader)

import pandas as pd
import plotly.express as px
import os

# System library to manipulate the file system
from utils import write_excel
import shutil
from quads import *
data_path = config["directory_management"]["data_path"]
result_path = config["directory_management"]["result_path"]
#name the files
if config["file_management"]["table"] == "excel" :
  df=pd.read_excel(data_path+config["file_management"]["original_data_file"],\
                     index_col=0)
elif config["file_management"]["table"] == "csv" :
  df=pd.read_csv(data_path+config["file_management"]["original_data_file"],\
                   index_col=0)

file_name_qualitative = result_path+\
	                    config["file_management"]["qualitative_enrichment"]

factor= config["variable_management"]["factor_variable"]

#make the visualisations
col = {'over-represented' : config["figure_management"]["over_represented"], \
	   'under-represented' : config["figure_management"]["under_represented"], \
	   'Not significant': config["figure_management"]["not_significant"]}
general_color = config["figure_management"]["general_color"]

if config["file_management"]["table"] == "excel" :
  df = pd.read_excel(file_name_qualitative,index_col=0)
elif config["file_management"]["table"] == "csv" :
  df = pd.read_csv(file_name_qualitative,index_col=0)

sunburst = px.sunburst(df, path=[factor, 'variables', 'modalities'],\
                       values=config["figure_management"]["statistic"], \
                       color = 'interpretation',\
                       color_discrete_map=col,\
                       color_discrete_sequence = [general_color])

sunburst.add_annotation(x=0.2,y=1,text= 'Over-represented',\
font = dict(color=config["figure_management"]["over_represented"],size=14),\
showarrow=False)
sunburst.add_annotation(x=0.2,y=0.95,text= 'Under-represented',\
font = dict(color=config["figure_management"]["under_represented"],size=14),\
showarrow=False)
sunburst.add_annotation(x=0.2,y=0.9,text= 'Not significant',\
font = dict(color=config["figure_management"]["not_significant"],size=14),\
showarrow=False)
sunburst.show()

