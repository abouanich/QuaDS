
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
cluster= config["variable_management"]["cluster_variable"]

#quantitatives variables
quantitative =config["variable_management"]["quantitative_variables"]
df_quantitative = df[quantitative+\
				    [config["variable_management"]["cluster_variable"]]]
df_quantitative = df_quantitative.fillna(0)

#make the quantitative analysis for each quantitative variable
sd = sdquanti(df_quantitative, \
			  config["variable_management"]["quantitative_variables"], \
			  config["variable_management"]["cluster_variable"])
quanti_a = quanti_analysis(sd, \
					df_quantitative, \
					config["variable_management"]["quantitative_variables"],\
					config["variable_management"]["cluster_variable"], \
					config["thresholds_management"]["anova_threshold"])
print("quantitative analysis done.")

#qualitatives variables
qualitative = config["variable_management"]["qualitative_variables"]
df_qualitative = df[qualitative+\
					[config["variable_management"]["cluster_variable"]]]
df_qualitative = df_qualitative.fillna('missing values')


#make the qualitative analysis
sdqualitative = sdquali(df_qualitative, \
						qualitative, \
						cluster, \
						config["thresholds_management"]["x2_threshold"])
quali_a = quali_analysis(cluster)
print("qualitative analysis done.")
cla_mod = clamod(quali_a,cluster)
print("cla/mod calcul done.")
mod_cla = modcla(quali_a,cluster)
print("mod/cla calcul done")
global_stat = globa(quali_a)
print("global calcul done.")
p_value = pvalue(quali_a,cluster)
print('pvalue done.')
test_value = vtest(quali_a,\
				config["thresholds_management"]["hypergeometric_threshold"])
print("hypergeometric distribution done.")
weight = variable_weight(quali_a)
print("variable weight done.")


#out :
#create the new path for the result
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



#make the visualisations
col = {'overrepresented' : config["figure_management"]["over_represented"], \
	   'underrepresented' : config["figure_management"]["under_represented"],\
	   'Not significant': config["figure_management"]["not_significant"]}


if config["file_management"]["table"] == "excel" :
  df = pd.read_excel(file_name_qualitative,index_col=0)
elif config["file_management"]["table"] == "csv" :
  df = pd.read_csv(file_name_qualitative,index_col=0)

sunburst = px.sunburst(df, path=[cluster, 'variables', 'modalities'],\
                       values=config["figure_management"]["statistic"], \
                       color = 'signification',\
                       color_discrete_map=col)
sunburst.add_annotation(x=0,y=1.1,\
                        font=dict(color='black',size=14),\
                        showarrow=False)
sunburst.add_annotation(x=0.2,y=1,text= 'Overrepresented',\
font = dict(color=config["figure_management"]["over_represented"],size=14),\
showarrow=False)
sunburst.add_annotation(x=0.2,y=0.95,text= 'Underrepresented',\
font = dict(color=config["figure_management"]["under_represented"],size=14),\
showarrow=False)
sunburst.add_annotation(x=0.2,y=0.9,text= 'Not significant',\
font = dict(color=config["figure_management"]["not_significant"],size=14),\
showarrow=False)
sunburst.show()

