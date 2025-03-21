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
import time 

###############################################################################
#start timer
###############################################################################
start_time = time.time()

###############################################################################
#name variables
###############################################################################

data_path = config["directory_management"]["data_path"]
result_path = config["directory_management"]["result_path"]
file_name_x2 = result_path+\
	           config["file_management"]["chi2_test"]
file_name_fisher = result_path+\
	           config["file_management"]["fisher_test"]
file_name_qualitative = result_path+\
	                    config["file_management"]["qualitative_enrichment"]
file_name_weight = result_path+\
	               config["file_management"]["variable_weight"]
file_name_normality = result_path+\
	              config["file_management"]["normality_test"]
file_name_homoscedasticity = result_path+\
	              config["file_management"]["homoscedasticity_test"]
file_name_anova = result_path+\
	              config["file_management"]["anova_test"]
file_name_kruskal_wallis = result_path+\
	              config["file_management"]["kruskal_wallis_test"]
file_name_quantitative = result_path+\
	                     config["file_management"]["quantitative_enrichment"]
separator = config["file_management"]["sep"]
index = config["file_management"]["index"]
tab_type = config["file_management"]["table"]
factor= config["variable_management"]["factor_variable"]
if type(factor) != str :
  if config["logging"]["log_level"]=="twice":
    print("Your factor variable have to be repertoried as a string")
    logger.warning("Your factor variable have to be repertoried as a string")
  elif config["logging"]["log_level"]== "console" :
    print("Your factor variable have to be repertoried as a string")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("Your factor variable have to be repertoried as a string")

  sys.exit()

############################################################################### 
#open the file
###############################################################################

if tab_type == "xlsx" and index == True:
  try :
    df=pd.read_excel(data_path+config["file_management"]["original_data_file"],\
                    index_col=0)
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

elif tab_type == "xlsx" and index == False:
  try :
    df=pd.read_excel(data_path+config["file_management"]["original_data_file"])
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
ms_quanti=config["missing_data_management"]["quanti_missing_data"]
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
				    [config["variable_management"]["factor_variable"]]]
except KeyError:
  if config["logging"]["log_level"]=="twice":
    print("One or more quantitative variable(s) is/are not in the header table.")
    print("or factor variable is not in the header table.")
    logger.warning("One or more quantitative variable(s) is/are not in the header table.")
    logger.warning("or factor variable is not in the header table.")
  elif config["logging"]["log_level"]== "console" :
    print("One or more quantitative variable(s) is/are not in the header table.")
    print("or factor variable is not in the header table.")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("One or more quantitative variable(s) is/are not in the header table.")
    logger.warning("or factor variable is not in the header table.")
  sys.exit()

if ms_quanti=="drop":
  try :
    df_quantitative = df_quantitative.infer_objects()
    for col in quantitative:
      na_count = df_quantitative[col].isnull().values.sum()
      if na_count != 0:
        if config["logging"]["log_level"]=="twice":
          print(na_count, "missing values are in the column",col,\
                "and the line containing these missing values are delete")
          logger.info(str(na_count)+ " missing values are in the column "+\
            col+" and the line containing these missing values are delete")
        elif config["logging"]["log_level"]== "console" :
          print(na_count, "missing values are in the column",col,\
                "and the line containing these missing values are delete")
        elif config["logging"]["log_level"]== "logger": 
          logger.info(str(na_count)+ " missing values are in the column "+\
            col+" and the line containing these missing values are delete")
    df_quantitative = df_quantitative.dropna()
  except ValueError :
    if config["logging"]["log_level"]=="twice":
      print("One/or more of your quantitative variable(s) is/are not quantitative")
      logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
    elif config["logging"]["log_level"]== "console" :
      print("One/or more of your quantitative variable(s) is/are not quantitative")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
    sys.exit()

elif ms_quanti=="zero":
  try :
    df_quantitative = df_quantitative.infer_objects()
    for col in quantitative:
      na_count = df_quantitative[col].isnull().values.sum()
      if na_count != 0:
        df_quantitative[col] = df_quantitative[col].fillna(0)
        if config["logging"]["log_level"]=="twice":
          print(na_count, "missing values are in the column",col,\
                "and the missing values are replaced by 0")
          logger.info(str(na_count)+" missing values are in the column "+\
            col+" and the missing values are replaced by 0")
        elif config["logging"]["log_level"]== "console" :
          print(na_count, "missing values are in the column",col,\
                "and the missing values are replaced by 0")
        elif config["logging"]["log_level"]== "logger": 
          logger.info(str(na_count)+" missing values are in the column "+\
            col+" and the missing values are replaced by 0")
  except ValueError :
    if config["logging"]["log_level"]=="twice":
      print("One/or more of your quantitative variable(s) is/are not quantitative")
      logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
    elif config["logging"]["log_level"]== "console" :
      print("One/or more of your quantitative variable(s) is/are not quantitative")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
    sys.exit()

elif ms_quanti=="mean":
  try :
    df_quantitative = df_quantitative.infer_objects()
    for col in quantitative:
      na_count = df_quantitative[col].isnull().values.sum()
      if na_count != 0 : 
        df_quantitative[col] = df_quantitative[col].fillna(df_quantitative[col].mean())
        if config["logging"]["log_level"]=="twice":
          print(na_count, "missing values are in the column",col,\
                "and the missing values are replaced by the mean of the column",col)
          logger.info(str(na_count)+" missing values are in the column "+\
            col+" and the missing values are replaced by the mean of the column "+col)
        elif config["logging"]["log_level"]== "console" :
          print(na_count, "missing values are in the column",col,\
                "and the missing values are replaced by the mean of the column",col)
        elif config["logging"]["log_level"]== "logger": 
          logger.info(str(na_count)+" missing values are in the column "+\
            col+" and the missing values are replaced by the mean of the column "+col)
  except ValueError :
    if config["logging"]["log_level"]=="twice":
      print("One/or more of your quantitative variable(s) is/are not quantitative")
      logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
    elif config["logging"]["log_level"]== "console" :
      print("One/or more of your quantitative variable(s) is/are not quantitative")
    elif config["logging"]["log_level"]== "logger": 
      logger.warning("One/or more of your quantitative variable(s) is/are not quantitative")
    sys.exit()

else:
  if config["logging"]["log_level"]=="twice":
    print("Your parameter in config file is not known for missing data quantitative management")
    logger.warning("Your parameter in config file is not known for missing data quantitative management")
  elif config["logging"]["log_level"]== "console" :
    print("Your parameter in config file is not known for missing data quantitative management")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("Your parameter in config file is not known for missing data quantitative management")
  sys.exit()


###############################################################################
#make the quantitative analysis for each quantitative variable
###############################################################################
homosc_calcul, homosc_var, non_homos_var = quanti_homoscedasticity(df_quantitative,\
			  config["variable_management"]["quantitative_variables"], \
			  config["variable_management"]["factor_variable"],\
        config["thresholds_management"]["bartlett_threshold"])
normality_calcul, normal_var, non_normal_var =quanti_normality(df_quantitative,\
			    config["variable_management"]["quantitative_variables"],\
			    config["variable_management"]["factor_variable"],\
          config["thresholds_management"]["shapiro_threshold"]) 
var_anova = []
for i in normal_var :
  if i in homosc_var :
    var_anova.append(i)
sd_anova, anova_var = anova(df_quantitative, \
			  var_anova, \
			  config["variable_management"]["factor_variable"],\
				config["thresholds_management"]["anova_threshold"])
sd_kruskal_wallis, kw_var = kruskal_wallis(df_quantitative, \
			  non_homos_var, \
        non_normal_var, \
			  config["variable_management"]["factor_variable"],\
				config["thresholds_management"]["kruskal_wallis_threshold"])

var_quanti_desc = anova_var+kw_var

quanti_a = quanti_analysis(df_quantitative, \
			    config["variable_management"]["quantitative_variables"],\
          var_quanti_desc,\
					config["variable_management"]["factor_variable"], \
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
ms_quali = config["missing_data_management"]["quali_missing_data"]
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
					[config["variable_management"]["factor_variable"]]]
except KeyError:
  if config["logging"]["log_level"]=="twice":
    print("One or more qualitative variable(s) is/are not in the header table.")  
    print("or factor variable is not in the header table.")
    logger.warning("One or more qualitative variable(s) is/are not in the header table.")
    logger.warning("or factor variable is not in the header table.")
  elif config["logging"]["log_level"]== "console" :
    print("One or more qualitative variable(s) is/are not in the header table.")  
    print("or factor variable is not in the header table.")
  elif config["logging"]["log_level"]== "logger": 
    logger.warning("One or more qualitative variable(s) is/are not in the header table.")
    logger.warning("or factor variable is not in the header table.")
  sys.exit()

#df_qualitative = df_qualitative.astype(str)

for col in qualitative:
  na_count = df_qualitative[col].isnull().values.sum()
  if na_count != 0 :
    if ms_quali=="drop":
      if config["logging"]["log_level"]=="twice":
        print(na_count, "missing values are in the column",col,\
              "and the line containing these missing values are delete")
        logger.info(str(na_count)+ " missing values are in the column "+\
          col+" and the line containing these missing values are delete")
      elif config["logging"]["log_level"]== "console" :
        print(na_count, "missing values are in the column",col,\
              "and the line containing these missing values are delete")
      elif config["logging"]["log_level"]== "logger": 
        logger.info(str(na_count)+ " missing values are in the column "+\
          col+" and the line containing these missing values are delete")
      df_qualitative  = df_qualitative.dropna()
    else :
      if config["logging"]["log_level"]=="twice":
        print(na_count, "missing values are in the column",col,\
              "and the missing values are replaced by the modality you choose:"\
                ,ms_quali,"for the column",col)
        logger.info(str(na_count)+ " missing values are in the column "+\
          col+" and the missing values are replaced by the modality you choose: "\
          +ms_quali+" for the column "+col)
      elif config["logging"]["log_level"]== "console" :
        print(na_count, "missing values are in the column",col,\
              "and the missing values are replaced by the modality you choose:"\
                ,ms_quali,"for the column",col)
      elif config["logging"]["log_level"]== "logger": 
        logger.info(str(na_count)+ " missing values are in the column "+\
          col+" and the missing values are replaced by the modality you choose: "\
          +ms_quali+" for the column "+col)
    df_qualitative = df_qualitative.fillna(ms_quali)

###############################################################################
#make the qualitative analysis
###############################################################################

chi2, fisher = sdquali(df_qualitative, \
						qualitative, \
						factor, \
						config["thresholds_management"]["x2_threshold"],\
            config["thresholds_management"]["fisher_threshold"])
quali_a = quali_analysis(factor)
if config["logging"]["log_level"]=="twice":
  print("qualitative analysis done.")
  logger.info("qualitative analysis done.")
elif config["logging"]["log_level"]== "console" :
  print("qualitative analysis done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("qualitative analysis done.")

#cla/mod
cla_mod = clamod(quali_a,factor)
if config["logging"]["log_level"]=="twice":
  print("cla/mod calcul done.")
  logger.info("cla/mod calcul done.")
elif config["logging"]["log_level"]== "console" :
  print("cla/mod calcul done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("cla/mod calcul done.")

#mod/cla
mod_cla = modcla(quali_a,factor)
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
p_value = pvalue(quali_a,factor)
if config["logging"]["log_level"]=="twice":
  print('pvalue done.')
  logger.info('pvalue done.')
elif config["logging"]["log_level"]== "console" :
  print('pvalue done.')
elif config["logging"]["log_level"]== "logger": 
  logger.info('pvalue done.')

#test value
test_value = vtest(quali_a,\
				config["thresholds_management"]["hypergeometric_threshold"],\
        config["variable_management"]["factor_variable"])
if config["logging"]["log_level"]=="twice":
  print("hypergeometric distribution done.")
  logger.info("hypergeometric distribution done.")
elif config["logging"]["log_level"]== "console" :
  print("hypergeometric distribution done.")
elif config["logging"]["log_level"]== "logger": 
  logger.info("hypergeometric distribution done.")

#variable weight  
weight = variable_weight(test_value)
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

if tab_type == "xlsx" :
  if len(chi2) != 0 :
    write_excel(file_name_x2, chi2, idx=True)
  if len(fisher) != 0 :
    write_excel(file_name_fisher, fisher, idx=True)
  if len(test_value) != 0 :  
    write_excel(file_name_qualitative, test_value,idx=True)
  if len(weight) != 0 :
    write_excel(file_name_weight, weight,idx=True)
  if len(normality_calcul) != 0 :
    write_excel(file_name_normality, normality_calcul,idx=True)
  if len(homosc_calcul) != 0 :
    write_excel(file_name_homoscedasticity, homosc_calcul,idx=True)
  if len(sd_anova) != 0 :
    write_excel(file_name_anova, sd_anova,idx=True)
  if len(sd_kruskal_wallis) != 0 :
    write_excel(file_name_kruskal_wallis, sd_kruskal_wallis,idx=True)
  if len(quanti_a) != 0 :
    write_excel(file_name_quantitative, quanti_a, idx=True)
if tab_type == "csv" :
  if len(chi2) != 0:
    chi2.to_csv(file_name_x2)
  if len(fisher) != 0:
    fisher.to_csv(file_name_fisher)
  if len(test_value) != 0:
    test_value.to_csv(file_name_qualitative)
  if len(weight) != 0:
    weight.to_csv(file_name_weight)
  if len(normality_calcul) != 0:
    normality_calcul.to_csv(file_name_normality)
  if len(homosc_calcul) != 0:
    homosc_calcul.to_csv(file_name_homoscedasticity)
  if len(sd_anova) != 0:
    sd_anova.to_csv(file_name_anova)
  if len(sd_kruskal_wallis) != 0:
    sd_kruskal_wallis.to_csv(file_name_kruskal_wallis)
  if len(quanti_a) != 0 :
    quanti_a.to_csv(file_name_quantitative)

###############################################################################
#make the visualisations
###############################################################################

col = {'over-represented' : config["figure_management"]["over_represented"], \
	   'under-represented' : config["figure_management"]["under_represented"],\
	   'Not significant': config["figure_management"]["not_significant"]}
general_color = config["figure_management"]["general_color"]

if tab_type == "xlsx" :
  df = pd.read_excel(file_name_qualitative,index_col=0)
if tab_type == "csv" :
  df = pd.read_csv(file_name_qualitative,index_col=0)

sunburst = px.sunburst(df, path=[factor, 'variables', 'modalities'],\
                       values=config["figure_management"]["statistic"], \
                       color = 'interpretation',\
                       color_discrete_map=col,\
                       color_discrete_sequence = [general_color])

sunburst.add_annotation(x=0.2,y=1,text= 'Over-represented',\
font = dict(color=config["figure_management"]["over_represented"],size=20),\
showarrow=False)
sunburst.add_annotation(x=0.2,y=0.95,text= 'Under-represented',\
font = dict(color=config["figure_management"]["under_represented"],size=20),\
showarrow=False)
sunburst.add_annotation(x=0.2,y=0.9,text= 'Not significant',\
font = dict(color=config["figure_management"]["not_significant"],size=20),\
showarrow=False)
sunburst.show()

###############################################################################
#end timer
###############################################################################
end_time = time.time()
if config["logging"]["log_level"]=="twice":
  print("QuaDS time: "+str(round((end_time-start_time),2))+" seconds")
  logger.info("QuaDS time: "+str(round((end_time-start_time),2))+" seconds")
elif config["logging"]["log_level"]== "console" :
  print("QuaDS time: "+str(round((end_time-start_time),2))+" seconds")
elif config["logging"]["log_level"]== "logger": 
  logger.info("QuaDS time: "+str(round((end_time-start_time),2))+" seconds")
