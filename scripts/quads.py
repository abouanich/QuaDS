#!/usr/bin/python3
"""
:author: Andrea Bouanich
:version: 3.0
:email: andrea.bouanich@inrae.fr
:email: julie.bourbeillon@institut-agro.fr
"""

################################################################################
# Import Libraries
################################################################################

# Data analysis and manipulation library
import pandas as pd
import numpy as np

#Data manipulation library
from collections import OrderedDict

#Data statistics analysis library
from scipy.stats import chi2_contingency, fisher_exact, norm
from scipy.stats import shapiro, bartlett, kruskal
import math
from math import *
from itertools import combinations
import statsmodels.formula.api as smf
from statsmodels.formula.api import ols
import statsmodels.api as sm
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
#from fichier import files_enter

################################################################################
# catdes' traduction code
################################################################################

def bi(n,k) :
  """
    Function used to calculate the binomial coefficient

	  Args:
        n and k
        
      Returns:
        the binomial coeffcient of (n,k)
  """
  return math.comb(n,k)

def sdquali(df, columns, variable_cat, threshold_chi2, threshold_fisher_exact) :
  """
    Function used to select which modalities are over and under represented
    in the different groups of the variable_cat

    Actions performed:

      * Chi2 test for the variables with the variable_cat
      * Separate the pandas a dictionary of pandas for each group of variable_cat
      * Make the v-test and final statistics
	
	Args:
      df: a pandas DataFrame 
      columns : the selected columns
      variable_cat : the variable to test the chi2
      chi2_p_value : p-value for the chi2 test
        
    Returns:
      Dataframe : a panda DataFrame containing the chi2 statistics analysis 
          	     for the variables
  """

  #each columns make a chi2 test with the variable variable_cat
  global column
  column = []
  chi_p_value = []
  chi = []
  chi_column = []
  fisher_variables = []
  fisher_pvalue = []
  chi_significative = []
  fisher_interpretation = []
  for col in df[columns]:
    cont = pd.crosstab(df[col],df[variable_cat])
    cat_modalities = cont.columns.tolist()
	# Chi-square test of independence
    chi2, p_chi2, dof, expected = chi2_contingency(cont)
    if np.all(expected >= 5) :
      chi_column.append(col)
      if p_chi2 < 0.000001 :
        chi_p_value.append("<10-6")
      else:
        chi_p_value.append(round(p_chi2,7))
      chi.append(round(chi2,2))

      #if the p-value < threshold_chi2, the variable is variable chi2's
      #  dependent we can move the variable that are independant from 
      # the variable chi2 for the next steps
      if p_chi2 < threshold_chi2 :
        column.append(col)
        chi_significative.append('Significant')
      else :
        chi_significative.append('Not significant')
    else :
      fisher_variables.append(col)
      var_modalities = list(cont.index)
      cont = np.array(cont)
      stats = importr('stats')
      res = stats.fisher_test(cont)
      p_fisher = res[0][0]
      if p_fisher < 0.000001 :
        fisher_pvalue.append("<10-6")
      else:
        fisher_pvalue.append(round(p_fisher,6))
      if p_fisher < threshold_fisher_exact :
        fisher_interpretation.append('Significant')
        column.append(col)
      else : 
        fisher_interpretation.append("Not significant")
  global new_df
  new_df = df[column]
  new_df.insert(len(column),variable_cat,df[variable_cat].to_list())
  #generate the table from the chi2 test with the variables and their p_value

  X2 = pd.DataFrame({'Variables' : chi_column,
                     'Chi2 Statistic' : chi, 
                     'p-value' : chi_p_value,
			               'interpretation' : chi_significative})
  FISHER = pd.DataFrame({'Variables' : fisher_variables,
                         'p-value' : fisher_pvalue,
                         'interpretation':fisher_interpretation})
  return X2, FISHER

def quali_analysis(variable_cat):
  """
    Function used to select which modalities are over and under represented
    in the different groups of the variable_cat

    Actions performed:

      * Chi2 test for the variables with the variable_cat
      * Separate the pandas a dictionary of pandas for each group of variable_cat
      * Make the v-test and final statistics
	
	Args:
      variable_cat : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing statistics analysis for each
          	     variable_cat, variables and modalities
  """
  #separation of the dataframe for each variable chi2 with 
  #a dictionary of variable chi2
  global dictio
  dictio = {}
  for i in new_df[variable_cat].unique():
    dictio[i] = new_df[  new_df[variable_cat] == i ]
	#order the dictionary by the number of variable_cat
    dictio = OrderedDict(sorted(dictio.items(), key=lambda x: x[0]))
	
  #column : modalities
  global modality
  index = []
  nb_mod_by_var = []
  for col in new_df:
    unique_elements = new_df[col].unique()
    nb_mod_by_var.append(len(unique_elements))  
    index.extend(unique_elements)  
  index_chi2 = np.sort(new_df[variable_cat].unique())  
  index = [elem for elem in index if elem not in index_chi2]  
  modality = index * len(index_chi2) + list(index_chi2)

  #column : variable_cat	
  global chi2_var
  el = np.sort(new_df[variable_cat].unique())
  variable = [str(element) for element in el for _ in range(len(index))]
  chi2_var = variable + list(index_chi2)

  #column : variables
  global variables
  global variables_2
  var = [([i]*j) for i,j in zip(column,nb_mod_by_var)]
  nb_variable_cat = [variable_cat] * len(new_df[variable_cat].unique())
  variables = []
  for i in var:
    variables.extend(i)
  variables = variables * len(new_df[variable_cat].unique()) + nb_variable_cat
  variables_2 = []
  for i in var:
    variables_2.extend(i)
  variables_2.extend(nb_variable_cat)

  global NA
  NA = []
  for element in chi2_var :
    NA.append('Not present')
	
  result = pd.DataFrame({
		variable_cat : chi2_var,
		'variables' : variables ,
		'modalities':modality, 
		'cla/mod' : NA, 
		'mod/cla' : NA, 
		'global' : NA, 
		'p-value' : NA, 
		'v-test' : NA,
		'interpretation': NA})
  return result

#column : cla/mod
def clamod(result, variable_cat):
  """
    Function used to make the cla/mod column for the result table

    Actions performed:

      * create the cla/mod statistics for each variable_cat, variables and modalities
	
    Args:
      result : the pandas realised in the sdquali's function
      variable_cat : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing cla/mod statistics analysis
          			 for each variable_cat, variables and modalities
  """
  #add the cla/mod column from clamod to cla/mod column from result
  clamod = pd.DataFrame(columns = [variable_cat,'variables','modalities','cla/mod'])
  i_clamod = []
  l_clamod = []
  variable = []
  var_chi2 = []
  for c in new_df :
    for v in dictio :
      #chi2var/modality named : cla_mod 
      # % of the modality in each chi2var 
      # exemple : x% of the modality pink color are in the chi2var 'y'
      #global clamod_value
      clamod_value=pd.value_counts(dictio[v][c])/pd.value_counts(new_df[c])*100
      for i in clamod_value :
        l_clamod.append(i)
        variable.append(clamod_value.name)
        var_chi2.append(v)
      li_clamod = clamod_value.index.tolist()
      i_clamod.append(li_clamod)
      modalities = []
      for i in i_clamod :
        modalities.extend(i)
  for g,h,i,j in zip(var_chi2,variable,modalities,l_clamod) :
    clamod = clamod.append({variable_cat : g, 
								'variables' : h,
								'modalities' : i,
								'cla/mod' : j },
								ignore_index=True)
    clamod = clamod.fillna(0)	
  
  clamod_dict = {
    (str(row[variable_cat]), \
     str(row['variables']), \
     str(row['modalities'])): round(row['cla/mod'], 7)
    for index, row in clamod.iterrows()
  }
  for i, row in result.iterrows():
    key = (str(row[variable_cat]), \
           str(row['variables']), \
           str(row['modalities']))
    if key in clamod_dict:
      result.at[i, 'cla/mod'] = clamod_dict[key]
  return result

	
#column : mod/cla
def modcla(result,variable_cat):
  """
    Function used to make the mod/cla column for the result table

    Actions performed:

    * create the mod/cla statistics for each variable_cat, variables and modalities
	
    Args:
      result : the pandas realised in the sdquali's function
      variable_cat : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing mod/cla statistics analysis 
          			 for each variable_cat, variables and modalities
  """
  modcla = pd.DataFrame(columns = [variable_cat,'variables','modalities','mod/cla'])
  i_modcla = []
  l_modcla = []
  variable = []
  var_chi2 = []
	
  for c in new_df :
    for v in dictio :
      #modality/chi2var named : mean_category
      # % of the chi2var in the modality
      #exemple : x% of the colors in the chi2var 'y' are pink color
      #global modcla_value
      modcla_value = pd.value_counts(dictio[v][c])/dictio[v].shape[0]*100
      for i in modcla_value :
        l_modcla.append(i)
        variable.append(modcla_value.name)
        var_chi2.append(v)
      li_modcla = modcla_value.index.tolist()
      i_modcla.append(li_modcla)
      modalities = []
      for i in i_modcla :
        modalities.extend(i)
  for g,h,i,j in zip(var_chi2,variable,modalities,l_modcla) :
    modcla = modcla.append({variable_cat : g, 
								'variables' : h,
								'modalities' : i,
								'mod/cla' : j },
								ignore_index=True)
  modcla_dict = {
    (str(row[variable_cat]), \
     str(row['variables']), \
     str(row['modalities'])): round(row['mod/cla'], 7)
    for index, row in modcla.iterrows()
  }
  for i, row in result.iterrows():
    key = (str(row[variable_cat]), \
           str(row['variables']), \
           str(row['modalities']))
    if key in modcla_dict:
      result.at[i, 'mod/cla'] = modcla_dict[key]
  result["mod/cla"] = result["mod/cla"].replace('Not present',0)
  return result

#column : global	
def globa(result):
  """
    Function used to make global column for the result table

    Actions performed:

    * create the global statistics for each variable_cat, variables and modalities
	
    Args:
      result : the pandas realised in the sdquali's function
        
    Returns:
      DataFrame: a pandas DataFrame containing global statistics 
          			 analysis for each variable_cat, variables and modalities
  """

  glo = pd.DataFrame(columns = ['variables','modalities','global'])
  i_global = []
  l_global = []
  for col in new_df :
  #global named : glob
  # % of the global modality
  #exemple : in total, x% of the samples have pink color
    global glob
    glob = pd.value_counts(new_df[col])/(new_df.shape[0])*100
    for i in glob :
      l_global.append(i)
    list_modalities = glob.index.tolist()	
    i_global.append(list_modalities)
    modalities = []
    for i in i_global :
      modalities.extend(i)
  for h,i,j in zip(variables_2,modalities,l_global) :		
    glo = glo.append({'variables' : h,
				   		  'modalities' : i,
				   		  'global' : j},
				   		  ignore_index=True)
  global_dict = { 
     (str(row['variables']), \
     str(row['modalities'])): round(row['global'], 7)
    for index, row in glo.iterrows()
  }
  for i, row in result.iterrows():
    key = (str(row['variables']), \
           str(row['modalities']))
    if key in global_dict:
      result.at[i, 'global'] = global_dict[key]	
  
  return result

#column : p_value
def pvalue(result,variable_cat):
  """
    Function used to make the p-value column for the result table

    Actions performed:

    * create a table of different effectives
    * calculate the p-value
	
    Args:
      result : the pandas realised in the sdquali's function
      variable_cat : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing the p-value for each variable_cat, 
          		     variables and modalities
  """
  global list_pval
  list_pval = []  

  global table
  table = pd.DataFrame({
    variable_cat : chi2_var,
    'variables' : variables,
    'modalities': modality, 
    'nj' : NA, 
    'nkj' : NA,
    'nk' : NA
  })

  n = new_df.shape[0]

  nj_table = []
  for c in new_df.columns:
    nj = new_df[c].value_counts().reset_index()
    nj.columns = ['modalities', 'nj']
    nj['variables'] = c
    nj_table.append(nj)

  table_nj = pd.concat(nj_table, ignore_index=True)
  table_nj = table_nj.fillna(0)

  table_nk = pd.DataFrame(columns=[variable_cat,'nk'])
  for v, df in dictio.items():
    nk = df.shape[0]
    table_nk = table_nk.append({variable_cat: v,'nk':nk}, ignore_index=True)


  nkj_data = []
  for c in new_df.columns:
    for v, df in dictio.items():
      nkj = df[c].value_counts().reset_index()
      nkj.columns = ['modalities','nkj']
      nkj['variables'] = c
      nkj[variable_cat] = v
      nkj_data.append(nkj)

  table_nkj = pd.concat(nkj_data, ignore_index=True)
  table_nkj = table_nkj.fillna(0)
  table.merge(table_nj[['variables', 'modalities', 'nj']], \
                      on=['variables', 'modalities'], how='left')
  table.merge(table_nkj[[variable_cat, 'variables', 'modalities', 'nkj']], \
                      on=[variable_cat, 'variables', 'modalities'], how='left')
  table.merge(table_nk[[variable_cat, 'nk']],\
                      on=[variable_cat], how='left')

  table_nkj_dict = {
    (str(row[variable_cat]), str(row['variables']), str(row['modalities'])): row['nkj']
    for _, row in table_nkj.iterrows()
  }
  table_nj_dict = {
    (str(row['variables']), str(row['modalities'])): row['nj']
    for _, row in table_nj.iterrows()
  }
  table_nk_dict = {
    (str(row[variable_cat])): row['nk']
    for _, row in table_nk.iterrows()
  }
  for i, row in table.iterrows():
    key_nkj = (str(row[variable_cat]), str(row['variables']), str(row['modalities']))
    if key_nkj in table_nkj_dict:
        table.at[i, 'nkj'] = table_nkj_dict[key_nkj]
    key_nj = (str(row['variables']), str(row['modalities']))
    if key_nj in table_nj_dict:
        table.at[i, 'nj'] = table_nj_dict[key_nj]
    key_nk = str(row[variable_cat])
    if key_nk in table_nk_dict:
        table.at[i, 'nk'] = table_nk_dict[key_nk]
    if table.at[i, 'nkj'] == 'Not present':
        table.at[i, 'nkj'] = 0
  

  table_dict = {
    (str(row[variable_cat]), str(row['variables']), str(row['modalities'])): {
        'nj': row['nj'],
        'nk': row['nk'],
        'nkj': row['nkj']
    }
    for _, row in table.iterrows()
  }
  for i, row in result.iterrows():
    key=(str(row[variable_cat]),str(row['variables']),str(row['modalities']))
    if key in table_dict:
      nj=table_dict[key]['nj']
      nk=table_dict[key]['nk']
      nkj=table_dict[key]['nkj']
      n_nk=n-nk
      s1=0
      s2=0

      for y in range(nkj+1,nj):
        nm=nj-y
        p=(bi(nk,y)*bi(n_nk,nm))/bi(n,nj)
        s1+=p
        if s1==s1+p:
          break

      for y in range(nkj,nj):
        nm=nj-y
        p=(bi(nk,y)*bi(n_nk,nm))/bi(n,nj)
        s2+=p
        if s2==s2+p:
          break
        
      p_value=2*min(s1,(1-round(s2,7)))+(bi(nk,nkj)*bi(n_nk,nj-nkj))/bi(n,nj)
        
      list_pval.append(p_value)
        
      if p_value < 0.000001:
        result.at[i,'p-value']="<10-6"
      else:
        result.at[i,'p-value']=round(p_value,7)
  return result

def vtest(result, v_p_value,cluster) :
  """
    Function used to make the v-test column for the result table

    Actions performed:

    * calcule the v-test on each variable_cat, variables and modalities
	
    Args:
          
    Returns:
      DataFrame: a pandas DataFrame containing v-test for each variable_cat, 
          			 variables and modalities
  """

  for i in range(len(result)) :
    if result['mod/cla'].iloc[i] > result['global'].iloc[i] :
      val = 1
    else :
      val = 0 
    if list_pval[i] < v_p_value :
      result['v-test'].iloc[i] = round((1-2*val)*\
                                  norm.ppf(list_pval[i]/2),7)
    else :
      result['v-test'].iloc[i] = 'Not significant'
		
    if result['v-test'].iloc[i] != 'Not significant':
      if result['v-test'].iloc[i] < 0 :
        result['interpretation'].iloc[i] ='under-represented'
      elif result['v-test'].iloc[i] > 0 :
        result['interpretation'].iloc[i] ='over-represented'
    elif result['mod/cla'].iloc[i] == 0 :
      result['interpretation'].iloc[i] = 'Not present'
    elif result['v-test'].iloc[i] == 'Not significant' :
      result['interpretation'].iloc[i] ='Not significant'
  result = result[(result["variables"] != cluster)].reset_index(drop=True)
  return result

def variable_weight(result):
  """
    Function used to know the weight of each variable in the clustering

    Actions performed:

    * calcule the weight of a variable counting the number of modality
      over-represented in this variable divided by the number of modality
      of the variable
	
    Args: the qualitative_analysis table
          
    Returns:
      Dataframe : a panda DataFrame containing the chi2 statistics analysis 
          	         for the variables and the weight of these variables.	
	
  """
  over = []
  under = []
  over_under = []
  mod = []
  for i in range(len(result)) :
    if result['interpretation'][i] == 'over-represented' :
      over.append(result['variables'][i])
    if result['interpretation'][i] == 'under-represented' :
      under.append(result['variables'][i])
    if result['interpretation'][i] == 'over-represented' or\
        result['interpretation'][i] == 'under-represented' :
      over_under.append(result['variables'][i])
    mod.append(result['variables'][i])
  var = set(over)
  var = list(var)
  var.sort()
  over.sort()
  over_under.sort()

  nb_var = []
  nb_mod_over = []
  nb_mod_under = []
  nb_mod_overunder = []
  for element in var :
    nb_var.append(mod.count(element))
    nb_mod_over.append(over.count(element))
    nb_mod_under.append(under.count(element))	
    nb_mod_overunder.append(over_under.count(element))
  ratio1 = []
  for i,j in zip (range(len(nb_mod_over)),range(len(nb_var))) :
      ratio1.append(round(nb_mod_over[i]/nb_var[j]*100,2))	
  m = ratio1.copy()
  ranking1 = ratio1.copy()
  m.sort()
  r1 = np.array(ratio1)
  i=1
  for num in m :
    maxi = max(m)
    ind = np.where(r1 == maxi)
    for element in ind[0] :
      ranking1[element] = i
    m=m[0:-1]
    i=i+1	
			
  ratio2 = []
  for i,j in zip(range(len(nb_mod_under)),range(len(nb_var))):
    ratio2.append(round(nb_mod_under[i]/nb_var[j]*100,2))
  n = ratio2.copy()
  ranking2 = ratio2.copy()
  n.sort()
  r2 = np.array(ratio2)
  i=1
  for num in n :
    maximum = max(n)
    ind = np.where(r2 == maximum)
    for element in ind[0] :
      ranking2[element] = i
    n=n[0:-1]
    i=i+1	
	
  ratio3 = []
  for i,j in zip(range(len(nb_mod_overunder)),range(len(nb_var))):
    ratio3.append(round(nb_mod_overunder[i]/nb_var[j]*100,2))	
  o = ratio3.copy()
  ranking3 = ratio3.copy()
  o.sort()
  r3 = np.array(ratio3)
  i=1
  for num in o :
    maxim = max(o)
    ind = np.where(r3 == maxim)
    for element in ind[0] :
      ranking3[element] = i
    o=o[0:-1]
    i=i+1
				
  weight = pd.DataFrame({'variables' : var,
						            'sum of all mod of all groups' : nb_var,
		    			          'number mod over' : nb_mod_over, 
						            'ratio over mod' : ratio1,
						            'contribution over mod' : ranking1,
		    			          'number mod under' : nb_mod_under,
						            'ratio under mod' : ratio2,
						            'contribution under mod' : ranking2,
						            'number mod over&under' : nb_mod_overunder,
						            'ratio over&under mod' : ratio3,
						            'contribution over&under mod' : ranking3})
  return weight

def quanti_normality(df,quanti_var,variable_cat,shapiro_pvalue):
  """
    Actions performed:
    * Make the normality test on each quantitative variable
	
    Args:
      df: a pandas DataFrame containing only the quantitatives variables
      quanti_var: name of the quantitative variable
      threshold_normality: threshold choose by the user
        
    Returns:
      DataFrame: a pandas DataFrame containing the normality results
      List: a list containing the normal variables
      List: a list containing the non normal variables
  """
  list_stat=[]
  list_variable = []
  list_factor=[]
  list_pvalue = []
  normal_variables=[]
  non_normal_variables = []
  factor_modalities = list(set(df[variable_cat].to_list()))
  for variable in quanti_var :
    count=0
    for factor_modality in factor_modalities :
      list_variable.append(variable)
      list_factor.append(factor_modality)
      df_cat = df[df[variable_cat]==factor_modality]
      stat, p_value = shapiro(df_cat[variable])
      list_stat.append(round(stat,6))
      if p_value < 0.000001 :
        list_pvalue.append("<10-6")
      else: 
        list_pvalue.append(round(p_value,6))
      if p_value < shapiro_pvalue :
        count = count
      else :
        count +=1
    if count != len(factor_modalities) :
        non_normal_variables.append(variable)
    else:
      normal_variables.append(variable)
  output_shapiro = pd.DataFrame({"variable":list_variable,\
                                 "Factor modality":list_factor,\
                                 "statistic":list_stat,\
                                 "p-value":list_pvalue})
  return output_shapiro, normal_variables, non_normal_variables

def quanti_homoscedasticity(df,quanti_var, variable_cat,homoscedasticity_pvalue):
  """
    Actions performed:
    * Make the homoscedasticity test on each quantitative variable
	
    Args:
      df: a pandas DataFrame containing only the quantitatives variables
      quanti_var: name of the quantitative variable
      variable_cat: the categorial variable
      homoscedasticity_normality: threshold choose by the user
        
    Returns:
      DataFrame: a pandas DataFrame containing the homoscedasticity results
      List: a list containing the homoscedastic variables
      List: a list containing the non homoscedastic variables
  """
  list_stat = []
  list_pvalue = []
  homoscedasticity_variables = []
  non_homoscedasticity_variables = []
  for var in quanti_var :
    df_cat = [df[df[variable_cat] == cat][var] for cat in df[variable_cat].unique()]
    stat, p_value = bartlett(*df_cat)
    list_stat.append(round(stat,6))
    if p_value < 0.000001 :
      list_pvalue.append("<10-6")
    else: 
      list_pvalue.append(round(p_value,6))
    if p_value > homoscedasticity_pvalue : 
      homoscedasticity_variables.append(var)
    else :
      non_homoscedasticity_variables.append(var)
  output_bartlett = pd.DataFrame({"variable":quanti_var,\
                                  "statistic":list_stat,\
                                  "p-value":list_pvalue})
  return output_bartlett, homoscedasticity_variables, non_homoscedasticity_variables


def anova(df, var, variable_cat, threshold_anova):
  """
    Actions performed:
    * Make the anova on each quantitative variable with variable_cat
	
    Args:
      df: a pandas DataFrame containing only the quantitatives variables
      var : the quantitative variable
      variable_cat : the variable to test
      threshold_anova: threshold choose by the user
        
    Returns:
      DataFrame: a pandas DataFrame containing the anova results
      List: a list containing the significative variables to the anova
  """
  #separation of the dataframe for each variable_cat with 
  #a dictionary of variable_cat
  global dictionary
  dictionary = {}
  for i in df[variable_cat].unique():
    dictionary[i] = df[  df[variable_cat] == i ]
    #order the dictionary by the number of variable_cat
    dictionary = OrderedDict(sorted(dictionary.items(),key=lambda x: x[0]))
	
  #a dictionary that don't contain missing values
  global dictionary_1
  dictionary_1 = {}
  df_na = df.dropna()#remove the missing values
  for i in df_na[variable_cat].unique():
    dictionary_1[i] = df_na[  df_na[variable_cat] == i ]
    #order the dictionary_1 by the number of variable_cat
    dictionary_1 = OrderedDict(sorted(dictionary_1.items(),key=lambda x: x[0]))
	
  #anova
  list_var = []
  list_eta2 = []
  list_pvalue = []
  info_interpretation = []
  signi_anova_var = []
  for varia in var :
    lm = ols('df_na[varia] ~ C(df_na[variable_cat])', data = df_na).fit()
    table = sm.stats.anova_lm(lm)
  #eta2 : Etat squared ranges from 0 to 1, where values close to 1 indicate a
  #higher proportion of variance that can be explained by a given variable in 
  #the model
    eta2 = table['sum_sq'][0]/(table['sum_sq'][0]+table['sum_sq'][1])
    pval = table['PR(>F)'][0]
    list_var.append(varia)
    list_eta2.append(round(eta2,6))
    if pval < 0.000001 :
      list_pvalue.append("<10-6")
    else: 
      list_pvalue.append(round(pval,6))
    if pval < threshold_anova :
      signi_anova_var.append(varia)
      info_interpretation.append("Significant")
    else :
      info_interpretation.append("Not significant")
  anova = pd.DataFrame({'variable' : list_var,
						  'eta-squared' : list_eta2,
						  'p-value' : list_pvalue,
              'interpretation' : info_interpretation})
  return anova, signi_anova_var

def kruskal_wallis(df,var_non_homos,var_non_normal,variable_cat,threshold_kw):
  """
    Actions performed:
    * Make the kruskal wallis on each quantitative variable non homoscedatic
     and/or non normal with variable_cat
	
    Args:
      df: a pandas DataFrame containing only the quantitatives variables
      var_non_homos: the quantitative variable non homoscedastic
      var_non_normal: the quantitative variable non normal
      variable_cat: the variable to test
      threshold_kw: threshold choose by the user
        
    Returns:
      DataFrame: a pandas DataFrame containing the kruskal wallis results
      List: a list containing the significative variables to the kw and not
            contained in the var_non_homos list
  """
  quanti_var= []
  for i in var_non_homos:
    if i not in quanti_var :
      quanti_var.append(i)
  for i in var_non_normal :
    if i not in quanti_var:
      quanti_var.append(i)
  list_stat = []
  list_pvalue = []
  list_interpretation = []
  signi_kw_var = []
  for var in quanti_var :
    df_cat=[df[df[variable_cat]==cat][var] for cat in df[variable_cat].unique()]
    stat, p_value = kruskal(*df_cat)
    list_stat.append(stat)
    if p_value < 0.000001 :
      list_pvalue.append("<10-6")
    else: 
      list_pvalue.append(round(p_value,6))
    if p_value < threshold_kw and var not in var_non_homos:
      signi_kw_var.append(var)
      list_interpretation.append("Significant")
    else : 
      list_interpretation.append("Not significant")

  output_kruskal_wallis = pd.DataFrame({"variable":quanti_var,\
                                        "statistic":list_stat,\
                                        "p-value":list_pvalue,\
                                        "interpretation":list_interpretation})
  return output_kruskal_wallis, signi_kw_var

def quanti_analysis(df, var, signi_variable, variable_cat, thres_gaussian):
  """
    Actions performed:
    * Make the v-test and final statistics
	
    Args:
      df: a pandas DataFrame containing only quantitative variable
      var : the quantitative variables
      signi_variable: a list containing the significative variables to the anova
                      and kruskal wallis and homoscedastic 
      variable_cat : the variable to test
      thres_gaussian : the gaussian threshold for the distribution
        
    Returns:
      DataFrame: a pandas DataFrame containing statistics analysis for each
          	         variable_cat and variables
  """
  df_na = df.dropna()#dataframe without missing values
  variable_cluster = []
  variable = []
  v_test = []
  mean_category = []
  overall_mean = []
  category_sd = [] 
  overall_sd = [] 
  pv = []
  info_interpretation = []
  for varia in var : 
    I = len(df) #number of individuals
    Im = len(df_na)#number of individuals that don't contain missing values
    x = round(df[varia].mean(),6)#mean of the modalities
    #check if all the data in the variable_cat are missing values
    ms = []
    for element in dictionary :
      variable_cluster.append(element)
      variable.append(varia)
      ms.append(dictionary[element][varia].isna().all())
    if True in ms :
      mean_category.append('Not significant')
      overall_mean.append('Not significant')
      category_sd.append('Not significant')
      overall_sd.append('Not significant')
      pv.append('Not significant')
      v_test.append('Not significant')
      info_interpretation.append('Not significant')

    else : 
      sous_ov = np.empty((0,2), float)
      df_na = df_na.reset_index(drop=True)
      for el in range(len(df_na)) :
        df_na=df_na.astype(type('float'))
        sous_ov=np.append(sous_ov,np.array([[float(df_na[varia][el]),\
                                             float(x)]]),axis=0)
      som_ov = 0
      for i in sous_ov :
        som_ov = som_ov+(i[0]-i[1])*(i[0]-i[1])
      et_overall = round(sqrt(som_ov/Im),6)
      if varia not in signi_variable : 
        for i,j in zip(dictionary_1,dictionary) :
          v_test.append('Not significant')
          pv.append("Not significant")
          info_interpretation.append("Not significant")
          Iq_df_na = len(dictionary_1[i])
          Iq = len(dictionary[j][varia])
          xq = round(dictionary[j][varia].mean(),6)
          #mean of the modalities of the cluster 
          data = dictionary_1[i]
          mean_category.append(float(xq))
          overall_mean.append(float(x))
          sous_cat = np.empty((0,2), float)
          for el in range(len(data)) :
            data = data.reset_index(drop=True)
            data=data.astype(type('float'))
            sous_cat=np.append(sous_cat,np.array([[float(data[varia][el]), \
                                                   float(xq)]]),axis=0)
          som_cat = 0
          for k in sous_cat :
            som_cat = som_cat+(k[0]-k[1])*(k[0]-k[1])
          et_category = round(sqrt(som_cat/Iq_df_na),6)
          category_sd.append(et_category)
          overall_sd.append(et_overall)
      else : 
        for i,j in zip(dictionary_1,dictionary) :
          Iq_df_na = len(dictionary_1[i])
          Iq = len(dictionary[j][varia])
          #mean of the modalities of the cluster
          xq = round(dictionary[j][varia].mean(),6) 
          data = dictionary_1[i]
          mean_category.append(float(xq))
          overall_mean.append(float(x))
          sous_cat = np.empty((0,2), float)
          for el in range(len(data)) :
            data = data.reset_index(drop=True)
            data=data.astype(type('float'))
            sous_cat=np.append(sous_cat,np.array([[float(data[varia][el]),\
                                                 float(xq)]]),axis=0)
          som_cat = 0
          for k in sous_cat :
            som_cat = som_cat+(k[0]-k[1])*(k[0]-k[1])
          et_category = round(sqrt(som_cat/Iq_df_na),6)
          category_sd.append(et_category)
          overall_sd.append(et_overall)
          vtest= round((float(xq)-float(x))/sqrt(((et_overall**2)/Iq)*\
                                                 ((I-Iq)/(I-1))),6)
          v_test.append(vtest)
          pvalue = (1-norm.cdf(abs(vtest)))*2
          if pvalue < 0.000001 :
            pv.append("<10-6")
          else: 
            pv.append(round(pvalue,6))
          if pvalue < thres_gaussian and vtest < 0 :
            info_interpretation.append("below average")
          elif pvalue < thres_gaussian and vtest > 0 :
            info_interpretation.append("above average")
          elif pvalue > thres_gaussian :
            info_interpretation.append("Not significant")

  quantitative = pd.DataFrame({variable_cat : variable_cluster,
							    'variable' : variable,
								'Mean in category':mean_category,
								'Overall mean': overall_mean,
								'Sd in category':category_sd,
								'Overall sd':overall_sd,
								'p-value':pv,
								'v-test':v_test,
                'interpretation':info_interpretation})
  return quantitative
