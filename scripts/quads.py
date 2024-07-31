#!/usr/bin/python3
"""
:author: Andrea Bouanich
:version: 1.0
:email: andrea.bouanich@inrae.fr
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
from scipy.stats import chi2_contingency
from scipy.stats import norm
import math
from math import *
import statsmodels.formula.api as smf
from statsmodels.formula.api import ols
import statsmodels.api as sm
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

def sdquali (df, columns, vchi2, chi2_p_value) :
  """
    Function used to select which modalities are over and under represented
    in the different groups of the vchi2

    Actions performed:

      * Chi2 test for the variables with the vchi2
      * Separate the pandas a dictionary of pandas for each group of vchi2
      * Make the v-test and final statistics
	
	Args:
      df: a pandas DataFrame 
      columns : the selected columns
      vchi2 : the variable to test the chi2
      chi2_p_value : p-value for the chi2 test
        
    Returns:
      Dataframe : a panda DataFrame containing the chi2 statistics analysis 
          	     for the variables
  """

  #each columns make a chi2 test with the variable vchi2
  global column
  column = []
  p_value = []
  chi = []
  significative = []
  for col in df[columns]:
    cont = pd.crosstab(df[col],df[vchi2])
	# Chi-square test of independence
    chi2, p, dof, expected = chi2_contingency(cont)
    if p < 0.000001 :
      p_value.append("<10-6")
    else:
      p_value.append(round(p,7))

    chi.append(round(chi2,2))

    #if the p-value < chi2_p_value, the variable is variable chi2's dependent
	#we can move the variable that are independant from the variable chi2 
	# for the next steps
    if p < chi2_p_value :
      column.append(col)
      significative.append('Significant')
    else :
      significative.append('Not significant')
  global new_df
  new_df = df[column]
  new_df.insert(len(columns)-1,vchi2,df[vchi2].to_list())
	
  #generate the table from the chi2 test with the variables and their p_value
  v = []
  for element in columns :
    v.append('Not significant')
  global X2

  X2 = pd.DataFrame({'Variables' : columns,
                     'Chi2 Statistic' : chi, 
                     'p-value' : p_value,
			               'interpretation' : significative})
  return X2

def quali_analysis(vchi2):
  """
    Function used to select which modalities are over and under represented
    in the different groups of the vchi2

    Actions performed:

      * Chi2 test for the variables with the vchi2
      * Separate the pandas a dictionary of pandas for each group of vchi2
      * Make the v-test and final statistics
	
	Args:
      vchi2 : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing statistics analysis for each
          	     vchi2, variables and modalities
  """
  #separation of the dataframe for each variable chi2 with 
  #a dictionary of variable chi2
  global dictio
  dictio = {}
  for i in new_df[vchi2].unique():
    dictio[i] = new_df[  new_df[vchi2] == i ]
	#order the dictionary by the number of vchi2
    dictio = OrderedDict(sorted(dictio.items(), key=lambda x: x[0]))
	

  #column : modalities
  index = []
  nb_mod_by_var=[]
  for col in new_df :
    nb = 0
    for element in new_df[col].unique() :
      index.append(element)
      nb +=1
    nb_mod_by_var.append(nb)
  nb_mod_by_var = nb_mod_by_var[0:-1]
  index_chi2 = index[-len(new_df[vchi2].unique()):]
  index_chi2.sort()
  index = index[0:-len(new_df[vchi2].unique())]#remove the vchi2 index
  global modality
  modality = index*len(new_df[vchi2].unique())
  for element in index_chi2 :
    modality.append(element)


  #column : vchi2	
  variable = []
  el = (new_df[vchi2].unique())
  el.sort()
  for element in el:
    t = (' '+str(element))*len(index)
    t =t.split()
    variable.append(t)
  global chi2_var
  chi2_var =[]
  for i in variable :
    chi2_var = chi2_var+i
  for j in index_chi2 :
    chi2_var.append(j)


  #column : variables
  var = []
  for i,j in zip(column,nb_mod_by_var) :
    u = (i+';')*j
    u = u.split(';')
    var.append(u)
  nb_vchi2 = []
  for element in new_df[vchi2].unique():
    nb_vchi2.append(vchi2)
  global variables
  variables = []
  for i in var :
    variables = variables+i
  while '' in variables :
    del variables[variables.index('')]
  variables = variables*len(new_df[vchi2].unique())
  for element in nb_vchi2 :
    variables.append(element)		
  global variables_2
  variables_2 =[]
  for i in var : 
    variables_2 = variables_2+i
  while '' in variables_2 :
    del variables_2[variables_2.index('')]
  variables_2 = variables_2 + nb_vchi2

  global NA
  NA = []
  for element in chi2_var :
    NA.append('Not present')
	
  result = pd.DataFrame({
		vchi2 : chi2_var,
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
def clamod(result, vchi2):
  """
    Function used to make the cla/mod column for the result table

    Actions performed:

      * create the cla/mod statistics for each vchi2, variables and modalities
	
    Args:
      result : the pandas realised in the sdquali's function
      vchi2 : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing cla/mod statistics analysis
          			 for each vchi2, variables and modalities
  """
  #add the cla/mod column from clamod to cla/mod column from result
  clamod = pd.DataFrame(columns = [vchi2,'variables','modalities','cla/mod'])
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
      l_modalities = []
      for i in i_clamod :
        for j in i :
          v = str(j)+(';')
          v = v.split(';')
          l_modalities.append(v)
	
        modalities = []
        for i in l_modalities :
          modalities = modalities+i
        while '' in modalities :
          del modalities[modalities.index('')]
			
  for g,h,i,j in zip(var_chi2,variable,modalities,l_clamod) :
    clamod = clamod.append({vchi2 : g, 
								'variables' : h,
								'modalities' : i,
								'cla/mod' : j },
								ignore_index=True)
    clamod = clamod.fillna(0)	
  for i in range(len(result)) :
    for j in range(len(clamod)):
      if str(result[vchi2][i]) == str(clamod[vchi2][j]) \
		and str(result['variables'][i]) == str(clamod['variables'][j]) \
			and str(result['modalities'][i]) == str(clamod['modalities'][j]):
        result['cla/mod'][i] = round(clamod['cla/mod'][j],7)
  return result

	
#column : mod/cla
def modcla(result,vchi2):
  """
    Function used to make the mod/cla column for the result table

    Actions performed:

    * create the mod/cla statistics for each vchi2, variables and modalities
	
    Args:
      result : the pandas realised in the sdquali's function
      vchi2 : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing mod/cla statistics analysis 
          			 for each vchi2, variables and modalities
  """
  modcla = pd.DataFrame(columns = [vchi2,'variables','modalities','mod/cla'])
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
      l_modalities = []
      for i in i_modcla :
        for j in i :
          v = str(j)+(';')
          v = v.split(';')
          l_modalities.append(v)
        modalities = []
        for i in l_modalities :
          modalities = modalities+i
        while '' in modalities :
          del modalities[modalities.index('')]
				
  for g,h,i,j in zip(var_chi2,variable,modalities,l_modcla) :
    modcla = modcla.append({vchi2 : g, 
								'variables' : h,
								'modalities' : i,
								'mod/cla' : j },
								ignore_index=True)
	
  for i in range(len(result)) :
  #add the mod/cla column from modcla to mod/cla column from result		
    for j in range(len(modcla)) :	
      if str(result[vchi2][i]) == str(modcla[vchi2][j]) and \
		str(result['variables'][i]) == str(modcla['variables'][j]) and \
		str(result['modalities'][i]) == str(modcla['modalities'][j]):
        result['mod/cla'][i] = round(modcla['mod/cla'][j],7)
      if result['mod/cla'][i] == 'Not present' :
        result['mod/cla'][i] = 0
  return result

#column : global	
def globa(result):
  """
    Function used to make global column for the result table

    Actions performed:

    * create the global statistics for each vchi2, variables and modalities
	
    Args:
      result : the pandas realised in the sdquali's function
        
    Returns:
      DataFrame: a pandas DataFrame containing global statistics 
          			 analysis for each vchi2, variables and modalities
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
    l_modalities = []
    for i in i_global :
      for j in i :
        v = str(j)+(';')
        v = v.split(';')
        l_modalities.append(v)
    modalities = []
    for i in l_modalities :
      modalities = modalities+i
    while '' in modalities :
      del modalities[modalities.index('')]
	
  for h,i,j in zip(variables_2,modalities,l_global) :		
    glo = glo.append({'variables' : h,
				   		  'modalities' : i,
				   		  'global' : j},
				   		  ignore_index=True)
	
  for i in range(len(result)) :	
  #add the global column from glo to global column from result		
    for j in range(len(glo)):
      if str(result['variables'][i]) == str(glo['variables'][j]) and \
		str(result['modalities'][i]) == str(glo['modalities'][j]):
        result['global'][i] = round(glo['global'][j],7)
  return result

#column : p_value
def pvalue(result,vchi2):
  """
    Function used to make the p-value column for the result table

    Actions performed:

    * create a table of different effectives
    * calculate the p-value
	
    Args:
      result : the pandas realised in the sdquali's function
      vchi2 : the variable to test the chi2
        
    Returns:
      DataFrame: a pandas DataFrame containing the p-value for each vchi2, 
          		     variables and modalities
  """
  global list_pval
  list_pval = []  
  #table is a df that contain all the numbers we need for the p-value calcul
  global table
  table = pd.DataFrame({vchi2 : chi2_var,
						  'variables' : variables ,
						  'modalities':modality, 
						  'nj' : NA, 
						  'nkj' : NA,
						  'nk' : NA})
  n = new_df.shape[0]#number of individuals
	
  #make the nj table
  #nj number is the size of the modality for all individuals	
  nj_values = [] 
  nj_mod = []
  table_nj = pd.DataFrame(columns = ['variables','modalities','nj'])
  for c in new_df :
    nj = pd.value_counts(new_df[c])
    for i in nj :
      nj_values.append(i)
    list_modalities_nj = nj.index.tolist()
    nj_mod.append(list_modalities_nj)
    l_modalities_nj = []
    for i in nj_mod :
      for j in i :
        v = str(j)+(';')
        v = v.split(';')
        l_modalities_nj.append(v)
    modalities_nj = []
    for i in l_modalities_nj :
      modalities_nj = modalities_nj+i
    while '' in modalities_nj :
      del modalities_nj[modalities_nj.index('')]
  for h,i,j in zip(variables_2,modalities_nj, nj_values) :	
    table_nj = table_nj.append({'variables' : h,
				   			  		'modalities' : i,
				   			  		'nj' : j},
				   			  		ignore_index=True)
    table_nj = table_nj.fillna(0)

  #make the nk table
  #nk number is the size of the cluster
  nk_values = []
  name_chi2var = []
  table_nk = pd.DataFrame(columns = [vchi2,'nk'])	
  for v in dictio :
    name_chi2var.append(v)
    nk_values.append(dictio[v].shape[0])
  for i,j in zip(name_chi2var, nk_values) :
    table_nk = table_nk.append({vchi2 : i, 
									'nk' : j},
									ignore_index=True) 
  #make the nkj table
  # #nkj number is the size of modality per clusters
  nkj_modalities = []
  values_nkj = []
  variables_nkj = []
  chi2var_nkj = []
  table_nkj = pd.DataFrame(columns = [vchi2,'variables','modalities','nkj'])	
  for c in new_df :
    for v in dictio :
      nkj = pd.value_counts(dictio[v][c]) 
      nk = dictio[v].shape[0]
      for i in nkj :
        values_nkj.append(i)
        variables_nkj.append(nkj.name)
        chi2var_nkj.append(v)
      mod_nkj = nkj.index.tolist()
      nkj_modalities.append(mod_nkj)
      l_modalities_nkj = []
      for i in nkj_modalities :
        for j in i :
          v = str(j)+(';')
          v = v.split(';')
          l_modalities_nkj.append(v)
        modalities_nkj = []
        for i in l_modalities_nkj :
          modalities_nkj = modalities_nkj+i
        while '' in modalities_nkj :
          del modalities_nkj[modalities_nkj.index('')]
			
  for g,h,i,j in zip(chi2var_nkj,variables_nkj,modalities_nkj,values_nkj) :
    table_nkj = table_nkj.append({vchi2 : g, 
						              'variables' : h,
						              'modalities' : i,
						              'nkj' : j },
						               ignore_index=True)
    table_nkj = table_nkj.fillna(0)
	
  #attribute the nkj number from table_nkj to table
  for i in range(len(table)) :	
    for j in range(len(table_nkj)) :	
      if str(table[vchi2][i]) == str(table_nkj[vchi2][j]) and \
          str(table['variables'][i]) == str(table_nkj['variables'][j]) and \
          str(table['modalities'][i]) == str(table_nkj['modalities'][j]):
        table['nkj'][i] = table_nkj['nkj'][j]
	
    #attribute the nj number from table_nj to table		
    for j in range(len(table_nj)):
      if str(table['variables'][i]) == str(table_nj['variables'][j]) and \
          str(table['modalities'][i]) == str(table_nj['modalities'][j]):
        table['nj'][i] = table_nj['nj'][j]
	
    #attribute the nk number from table_nk to table	
    for j in range(len(table_nk)):	
      if str(table[vchi2][i]) == str(table_nk[vchi2][j]) :
        table['nk'][i] = table_nk['nk'][j]
    if table['nkj'][i] == 'Not present' :
      table['nkj'][i] = 0

  #p-value calcul
  for i in range (len(result)) : 
    for j in range(len(table)) :
      if str(result[vchi2][i]) == str(table[vchi2][j]) and \
        str(result['variables'][i]) == str(table['variables'][j]) and \
        str(result['modalities'][i]) == str(table['modalities'][j]):
        nj = table['nj'][j]
        nk = table['nk'][j]
        nkj = table['nkj'][j]
        n_nk = n-nk
        s1=0
        s2=0
        #if aux2 > aux3
        for y in range(nkj+1,nj):
          nm = table['nj'][j]-y
          p = (bi(nk,y)*bi(n_nk,nm))/bi(n,nj)
          if s1 == s1+p :
            break
          else :
            s1 = s1+p
					
        #if aux2 < aux3
        for y in range(nkj,nj):
          nm = table['nj'][j]-y
          p = (bi(nk,y)*bi(n_nk,nm))/bi(n,nj)
          if s2 == s2+p :
            break
          else :
            s2 = s2+p
        p_value = 2*min(s1,(1-round(s2,7)))+\
            (bi(nk,nkj)*bi(n_nk,nj-nkj))/bi(n,nj)
        list_pval.append(p_value)
        if p_value < 0.000001 :
          result['p-value'][i] = "<10-6"
        else : 
          result['p-value'][i] = round(p_value,7)
  return result

def vtest(result, v_p_value) :
  """
    Function used to make the v-test column for the result table

    Actions performed:

    * calcule the v-test on each vchi2, variables and modalities
	
    Args:
          
    Returns:
      DataFrame: a pandas DataFrame containing v-test for each vchi2, 
          			 variables and modalities
  """
	
  for i in range(len(result)) :
    for j in range(len(table)) :
      if result['mod/cla'][i] > result['global'][i] :
        val = 1
      else :
        val = 0 
      if list_pval[i] < v_p_value :
        result['v-test'][i] = round((1-2*val)*\
                                    norm.ppf(list_pval[i]/2),7)
      else :
        result['v-test'][i] = 'Not significant'
		
    if result['v-test'][i] != 'Not significant':
      if result['v-test'][i] < 0 :
        result['interpretation'][i] ='under-represented'
      elif result['v-test'][i] > 0 :
        result['interpretation'][i] ='over-represented'
    elif result['mod/cla'][i] == 0 :
      result['interpretation'][i] = 'Not present'
    elif result['v-test'][i] == 'Not significant' :
      result['interpretation'][i] ='Not significant'
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
  m = ratio1[0:-1].copy()
  ranking1 = ratio1[0:-1].copy()
  m.sort()
  r1 = np.array(ratio1[0:-1])
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
  n = ratio2[0:-1].copy()
  ranking2 = ratio2[0:-1].copy()
  n.sort()
  r2 = np.array(ratio2[0:-1])
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
  o = ratio3[0:-1].copy()
  ranking3 = ratio3[0:-1].copy()
  o.sort()
  r3 = np.array(ratio3[0:-1])
  i=1
  for num in o :
    maxim = max(o)
    ind = np.where(r3 == maxim)
    for element in ind[0] :
      ranking3[element] = i
    o=o[0:-1]
    i=i+1
				
  weight = pd.DataFrame({'variables' : var[0:-1],
		    			   'number mod over' : nb_mod_over[0:-1], 
		    			   'number mod under' : nb_mod_under[0:-1],
						   'number mod over&under' : nb_mod_overunder[0:-1],
						   'sum of all mod of all groups' : nb_var[0:-1],
						   'ratio over/mod' : ratio1[0:-1],
						   'contribution over/mod' : ranking1,
						   'ratio under/mod' : ratio2[0:-1],
						   'contribution under/mod' : ranking2,
						   'ratio over&under/mod' : ratio3[0:-1],
						   'contribution over&under/mod' : ranking3})
  return weight

def sdquanti(df, var, variable_cat, threshold_anova):
  """
    Actions performed:
    * Make the anova on each quantitative variable with variable_cat
	
    Args:
      df: a pandas DataFrame containing only the quantitatives variables
      var : the quantitative variable
      variable_cat : the variable to test
        
    Returns:
      DataFrame: a pandas DataFrame containing the anova results
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
      info_interpretation.append("Significant")
    else :
      info_interpretation.append("Not significant")
  anova = pd.DataFrame({'variable' : list_var,
						  'eta-squared' : list_eta2,
						  'p-value' : list_pvalue,
              'interpretation' : info_interpretation})
  return anova

def quanti_analysis(anova, df, var, variable_cat, thres_anova,thres_gaussian):
  """
    Actions performed:
    * Make the v-test and final statistics
	
    Args:
      anova : the dataframe from the anova analysis
      df: a pandas DataFrame containing only quantitative variable 
      var : the quantitative variables
      variable_cat : the variable to test
      thres_anova : the anova threshold to continue the statistics
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
  pv = []#
  info_interpretation = []
  for varia in var : 
    I = len(df) #number of individuals
    Im = len(df_na)#number of individuals that don't contain missing values
    x = round(df[varia].mean(),6)#mean of the modalities
    #check if all the data in the vchi2 are missing values
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
      index_pvalue = anova.index[anova["variable"] == varia].to_list()
      for ind in index_pvalue :
        index_pv = ind
      if str(anova["p-value"][index_pv]) != "<10-6" :
        if anova["p-value"][index_pv] > thres_anova :
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
