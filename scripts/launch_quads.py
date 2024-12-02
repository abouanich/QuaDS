#!/usr/bin/python3

import yaml
with open("config_file.yml", "r") as yamlfile:
    config = yaml.load(yamlfile, Loader=yaml.FullLoader)

import pandas as pd
import plotly.express as px
import os
import logging
logger = logging.getLogger(__name__) 
logging.basicConfig(filename=config["logging"]["log_file"], encoding='utf-8',\
                    level=logging.DEBUG)
# System library to manipulate the file system
from utils import write_excel
import shutil
from quads import *
import sys

###############################################################################
#name variables
###############################################################################

data_path = config["directory_management"]["data_path"]
result_path = config["directory_management"]["result_path"]
file_name_x2 = config["directory_management"]["result_path"]+\
	           config["file_management"]["qualitative_test"]
file_name_qualitative = result_path+\
	                    config["file_management"]["qualitative_enrichment"]
file_name_weight = result_path+\
	               config["file_management"]["variable_weight"]
file_name_anova = result_path+\
	              config["file_management"]["quantitative_test"]
file_name_quantitative = result_path+\
	                     config["file_management"]["quantitative_enrichment"]
separator = config["file_management"]["sep"]
index = config["file_management"]["index"]
tab_type = config["file_management"]["table"]
cluster= config["variable_management"]["cluster_variable"]
if type(cluster) != str :
  if config["logging"]["log_level"]=="twice":
    print("Your cluster variable have to be repertoried as a string")
    logger.warning("Your cluster variable have to be repertoried as a string")
  elif config["logging"]["log_level"]== "console" :
    print("Your cluster variable have to be repertoried as a string")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("Your cluster variable have to be repertoried as a string")

  sys.exit()

############################################################################### 
#open the file
###############################################################################

if tab_type == "xslx" and index == True:
  try :
    df=pd.read_excel(data_path+config["file_management"]["original_data_file"],\
                   sep=separator, index_col=0)
  except FileNotFoundError: 
    if config["logging"]["log_level"]=="twice":
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "console" :
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    sys.exit()

elif tab_type == "xslx" and index == False:
  try :
    df=pd.read_excel(data_path+config["file_management"]["original_data_file"]\
                     ,sep=separator)
  except FileNotFoundError: 
    if config["logging"]["log_level"]=="twice":
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "console" :
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    sys.exit()

elif tab_type == "csv" and index == True:
  try :
    df=pd.read_csv(data_path+config["file_management"]["original_data_file"],\
                 sep = separator, index_col=0)
  except FileNotFoundError: 
    if config["logging"]["log_level"]=="twice":
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "console" :
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    sys.exit()


elif tab_type == "csv" and index == False:
  try :
    df=pd.read_csv(data_path+config["file_management"]["original_data_file"],\
                 sep = separator)
  except FileNotFoundError: 
    if config["logging"]["log_level"]=="twice":
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "console" :
      print("The file is not present in the repository",data_path)
      print("Or the name of your file is not the good one.")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("The file is not present in the repository",data_path)
      logger.warning("Or the name of your file is not the good one.")
    sys.exit()

###############################################################################
#quantitatives variables
###############################################################################

quantitative =config["variable_management"]["quantitative_variables"]
if type(quantitative) != list :
  if config["logging"]["log_level"]=="twice":
    print("Your quantative variables have to be repertoried in a list")
    logger.warning("Your quantative variables have to be repertoried in a list")
  elif config["logging"]["log_level"]== "console" :
    print("Your quantative variables have to be repertoried in a list")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("Your quantative variables have to be repertoried in a list")
  sys.exit()
try :
  df_quantitative = df[quantitative+\
				    [config["variable_management"]["cluster_variable"]]]
except KeyError:
  if config["logging"]["log_level"]=="twice":
    print("One or more quantitative variable(s) is/are not in the header table.")
    print("or cluster variable is not in the header table.")
    logger.warning("One or more quantitative variable(s) is/are not in the header table.")
    logger.warning("or cluster variable is not in the header table.")
  elif config["logging"]["log_level"]== "console" :
    print("One or more quantitative variable(s) is/are not in the header table.")
    print("or cluster variable is not in the header table.")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("One or more quantitative variable(s) is/are not in the header table.")
    logger.warning("or cluster variable is not in the header table.")
  sys.exit()

try :
  df_quantitative = df_quantitative.infer_objects()
  df_quantitative = df_quantitative.fillna(0)
except ValueError :
  if config["logging"]["log_level"]=="twice":
    print("One/or more of your quantitative variable(s) is/are not quantitative")
    logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
  elif config["logging"]["log_level"]== "console" :
    print("One/or more of your quantitative variable(s) is/are not quantitative")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
  sys.exit()

###############################################################################
#make the quantitative analysis for each quantitative variable
###############################################################################

sd = sdquanti(df_quantitative, \
			  config["variable_management"]["quantitative_variables"], \
			  config["variable_management"]["cluster_variable"],\
				config["thresholds_management"]["anova_threshold"])
quanti_a = quanti_analysis(sd, \
					df_quantitative, \
					config["variable_management"]["quantitative_variables"],\
					config["variable_management"]["cluster_variable"], \
					config["thresholds_management"]["anova_threshold"], \
          config["thresholds_management"]["gaussian_threshold"])

if config["logging"]["log_level"]=="twice":
  print("quantitative analysis done.")
  logger.info("quantitative analysis done.")
elif config["logging"]["log_level"]== "console" :
  print("quantitative analysis done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("quantitative analysis done.")

###############################################################################
#qualitatives variables
###############################################################################
qualitative = config["variable_management"]["qualitative_variables"]
if qualitative != [] :
  if config["logging"]["log_level"]=="twice":
    print("!! WARNING !!")
    print("Take care, be sure your qualitative variables are indeed qualitative")
    logger.warning("Take care, be sure your qualitative variables are indeed qualitative")
  elif config["logging"]["log_level"]== "console" :
    print("!! WARNING !!")
    print("Take care, be sure your qualitative variables are indeed qualitative")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("Take care, be sure your qualitative variables are indeed qualitative")

for variable in qualitative :
  variable_modalities = df[variable].to_list()
  if variable in variable_modalities :
    if config["logging"]["log_level"]=="twice":
      print("Take care, a modality have the same name as its variable!")
      print("Change the name of the variable or its modalities in your table")
      logger.warning("Take care, a modality have the same name as its variable!")
      logger.warning("Change the name of the variable or its modalities in your table")
    elif config["logging"]["log_level"]== "console" :
      print("Take care, a modality have the same name as its variable!")
      print("Change the name of the variable or its modalities in your table")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("Take care, a modality have the same name as its variable!")
      logger.warning("Change the name of the variable or its modalities in your table")
    sys.exit()

if type(qualitative) != list :
  if config["logging"]["log_level"]=="twice":
    print("Your qualitative variables have to be repertoried in a list")
    logger.warning("Your qualitative variables have to be repertoried in a list")
  elif config["logging"]["log_level"]== "console" :
    print("Your qualitative variables have to be repertoried in a list")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("Your qualitative variables have to be repertoried in a list")
  sys.exit()
try :
  df_qualitative = df[qualitative+\
					[config["variable_management"]["cluster_variable"]]]
except KeyError:
  if config["logging"]["log_level"]=="twice":
    print("One or more qualitative variable(s) is/are not in the header table.")  
    print("or cluster variable is not in the header table.")
    logger.warning("One or more qualitative variable(s) is/are not in the header table.")
    logger.warning("or cluster variable is not in the header table.")
  elif config["logging"]["log_level"]== "console" :
    print("One or more qualitative variable(s) is/are not in the header table.")  
    print("or cluster variable is not in the header table.")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("One or more qualitative variable(s) is/are not in the header table.")
    logger.warning("or cluster variable is not in the header table.")
  sys.exit()

df_qualitative = df_qualitative.astype(str)
df_qualitative = df_qualitative.fillna('missing values')

###############################################################################
#make the qualitative analysis
###############################################################################

sdqualitative = sdquali(df_qualitative, \
						qualitative, \
						cluster, \
						config["thresholds_management"]["x2_threshold"])
quali_a = quali_analysis(cluster)
if config["logging"]["log_level"]=="twice":
  print("qualitative analysis done.")
  logger.info("qualitative analysis done.")
elif config["logging"]["log_level"]== "console" :
  print("qualitative analysis done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("qualitative analysis done.")

#cla/mod
cla_mod = clamod(quali_a,cluster)
if config["logging"]["log_level"]=="twice":
  print("cla/mod calcul done.")
  logger.info("cla/mod calcul done.")
elif config["logging"]["log_level"]== "console" :
  print("cla/mod calcul done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("cla/mod calcul done.")

#mod/cla
mod_cla = modcla(quali_a,cluster)
if config["logging"]["log_level"]=="twice":
  print("mod/cla calcul done")
  logger.info("mod/cla calcul done")
elif config["logging"]["log_level"]== "console" :
  print("mod/cla calcul done")
elif config["logging"]["log_level"]== "logger": 
  logger.info("mod/cla calcul done")

#global
global_stat = globa(quali_a)
if config["logging"]["log_level"]=="twice":
  print("global calcul done.")
  logger.info("global calcul done.")
elif config["logging"]["log_level"]== "console" :
  print("global calcul done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("global calcul done.")

#pvalue
p_value = pvalue(quali_a,cluster)
if config["logging"]["log_level"]=="twice":
  print('pvalue done.')
  logger.info('pvalue done.')
elif config["logging"]["log_level"]== "console" :
  print('pvalue done.')
elif config["logging"]["log_level"]== "logger": 
  logger.info('pvalue done.')

#test value
test_value = vtest(quali_a,\
				config["thresholds_management"]["hypergeometric_threshold"])
if config["logging"]["log_level"]=="twice":
  print("hypergeometric distribution done.")
  logger.info("hypergeometric distribution done.")
elif config["logging"]["log_level"]== "console" :
  print("hypergeometric distribution done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("hypergeometric distribution done.")

#variable weight  
weight = variable_weight(quali_a)
if config["logging"]["log_level"]=="twice":
  print("variable weight done.")
  logger.info("variable weight done.")
elif config["logging"]["log_level"]== "console" :
  print("variable weight done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("variable weight done.")

###############################################################################
#out :
#create the new path for the result
###############################################################################

if not os.path.exists(config["directory_management"]["result_path"]) :
  os.makedirs(config["directory_management"]["result_path"])

if config["file_management"]["table"] == "excel" :
  write_excel(file_name_x2, sdqualitative, idx=True)
  write_excel(file_name_qualitative, test_value,idx=True)
  write_excel(file_name_weight, weight,idx=True)
  write_excel(file_name_anova, sd,idx=True)
  write_excel(file_name_quantitative, quanti_a, idx=True)
elif config["file_management"]["table"] == "csv" :
  sdqualitative.to_csv(file_name_x2)
  test_value.to_csv(file_name_qualitative)
  weight.to_csv(file_name_weight)
  sd.to_csv(file_name_anova)
  quanti_a.to_csv(file_name_quantitative)

###############################################################################
#make the visualisations
###############################################################################

col = {'over-represented' : config["figure_management"]["over_represented"], \
	   'under-represented' : config["figure_management"]["under_represented"],\
	   'Not significant': config["figure_management"]["not_significant"]}
general_color = config["figure_management"]["general_color"]

if config["file_management"]["table"] == "excel" :
  df = pd.read_excel(file_name_qualitative,index_col=0)
elif config["file_management"]["table"] == "csv" :
  df = pd.read_csv(file_name_qualitative,index_col=0)

sunburst = px.sunburst(df, path=[cluster, 'variables', 'modalities'],\
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

