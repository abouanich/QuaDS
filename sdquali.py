"""
Import input rose data and preprocess the dataset

:project: DIVIS: Biological ontology Integration and Visualisation on RoseData
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

################################################################################
# catdes' traduction code
################################################################################


def sdquali (df, columns, vchi2, chi2_p_value) :
	"""
        Function used to select which modalities are over and under represented
        in the different groups of the vchi2

        Actions performed:

        * Chi2 test for the variables with the vchi2
        * Separate the pandas a dictionary of pandas for each group of 
          vchi2
        * Make the v-test and final statistics
	
		Args:
        	df: a pandas DataFrame 
        	columns : the selected columns
        	vchi2 : the variable to test the chi2
        	chi2_p_value : p-value for the chi2 test
        
        Returns:
          Dataframe : a panda DataFrame containing the chi2 statistics analysis 
          	          for the variables
          DataFrame: a pandas DataFrame containing statistics analysis for each
          	         vchi2, variables and modalities
        """

	#each columns make a chi2 test with the variable vchi2
	global column
	column = []
	p_value = []
	chi = []
	significative = []
	for col in df:
		cont = pd.crosstab(df[col],df[vchi2])
		# Chi-square test of independence
		chi2, p, dof, expected = chi2_contingency(cont)
		p_value.append(p)
		chi.append(round(chi2,2))

		#if the p-value < chi2_p_value, the variable is variable 
		#chi2's dependent
		#we can move the variable that are independant from the 
		#variable chi2 for the next steps
		if p < chi2_p_value :
			column.append(col)
			significative.append('Significant')
		else :
			significative.append('Not significant')
	global new_df
	new_df = df[column]
	
	#generate the table from the chi2 test with the variables 
	#and their p_value
	v = []
	for element in columns :
		v.append('Not significant')
	global X2
	X2 = pd.DataFrame({'Variables' : columns,
			   'Chi2 Statistic' : chi[0:-1], 
			   'p_value' : p_value[0:-1],
			   'significance' : significative[0:-1]})
	return X2

def quali_analysis(df, columns, vchi2):
	"""
        Function used to select which modalities are over and under represented
        in the different groups of the vchi2

        Actions performed:

        * Chi2 test for the variables with the vchi2
        * Separate the pandas a dictionary of pandas for each group of 
          vchi2
        * Make the v-test and final statistics
	
		Args:
        	df: a pandas DataFrame 
        	columns : the selected columns
        	vchi2 : the variable to test the chi2
        	chi2_p_value : p-value for the chi2 test
        
        Returns:
          Dataframe : a panda DataFrame containing the chi2 statistics analysis 
          	          for the variables
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
		'signification': NA})
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
	i_cm = []
	l_cm = []
	ln_cm = []
	lv_cm = []
	for c in new_df :
		for v in dictio :
			#chi2var/modality named : cla_mod 
			# % of the modality in each chi2var 
			# exemple : x% of the modality pink color are in the chi2var 'y'
			global cm
			cm = pd.value_counts(dictio[v][c])/pd.value_counts(new_df[c])*100
			for i in cm :
				l_cm.append(i)
				ln_cm.append(cm.name)
				lv_cm.append(v)
			li_clamod = cm.index.tolist()
			i_cm.append(li_clamod)
			l_index_cm = []
			for i in i_cm :
				for j in i :
					v = str(j)+(';')
					v = v.split(';')
					l_index_cm.append(v)
	
				index_cm = []
				for i in l_index_cm :
					index_cm = index_cm+i
				while '' in index_cm :
					del index_cm[index_cm.index('')]
			
	for g,h,i,j in zip(lv_cm,ln_cm,index_cm,l_cm) :
		clamod = clamod.append({vchi2 : g, 
								'variables' : h,
								'modalities' : i,
								'cla/mod' : j },
								ignore_index=True)
		clamod = clamod.fillna(0)	
	for i in range(len(result)) :
		for j in range(len(clamod)):
			if str(result[vchi2][i]) == str(clamod[vchi2][j]) and str(result['variables'][i]) == str(clamod['variables'][j]) and str(result['modalities'][i]) == str(clamod['modalities'][j]):
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
	i_mc = []
	l_mc = []
	ln_mc = []
	lv_mc = []
	
	for c in new_df :
		for v in dictio :
			#modality/chi2var named : mc
			# % of the chi2var in the modality
			#exemple : x% of the colors in the chi2var 'y' are pink color
			global mc
			mc = pd.value_counts(dictio[v][c])/dictio[v].shape[0]*100
			for i in mc :
				l_mc.append(i)
				ln_mc.append(mc.name)
				lv_mc.append(v)
			li_modcla = mc.index.tolist()
			i_mc.append(li_modcla)
			li_mc = []
			for i in i_mc :
				for j in i :
					v = str(j)+(';')
					v = v.split(';')
					li_mc.append(v)
				i_modcla = []
				for i in li_mc :
					i_modcla = i_modcla+i
				while '' in i_modcla :
					del i_modcla[i_modcla.index('')]
				
	for g,h,i,j in zip(lv_mc,ln_mc,i_modcla,l_mc) :
		modcla = modcla.append({vchi2 : g, 
								'variables' : h,
								'modalities' : i,
								'mod/cla' : j },
								ignore_index=True)
	
	for i in range(len(result)) :
	#add the mod/cla column from modcla to mod/cla column from result		
		for j in range(len(modcla)) :	
			if str(result[vchi2][i]) == str(modcla[vchi2][j]) and str(result['variables'][i]) == str(modcla['variables'][j]) and str(result['modalities'][i]) == str(modcla['modalities'][j]):
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
	i_g = []
	l_g = []
	for col in new_df :
		#global named : glob
		# % of the global modality
		#exemple : in total, x% of the samples have pink color
		global glob
		glob = pd.value_counts(new_df[col])/(new_df.shape[0])*100
		for i in glob :
			l_g.append(i)
		list_index_global = glob.index.tolist()
		
		i_g.append(list_index_global)
		l_index_global = []
		for i in i_g :
			for j in i :
				v = str(j)+(';')
				v = v.split(';')
				l_index_global.append(v)
		index_g = []
		for i in l_index_global :
			index_g = index_g+i
		while '' in index_g :
			del index_g[index_g.index('')]
	
	for h,i,j in zip(variables_2,index_g,l_g) :	
		
		glo = glo.append({'variables' : h,
				   		  'modalities' : i,
				   		  'global' : j},
				   		  ignore_index=True)
	
	for i in range(len(result)) :	
	#add the global column from glo to global column from result		
		for j in range(len(glo)):
			if str(result['variables'][i]) == str(glo['variables'][j]) and str(result['modalities'][i]) == str(glo['modalities'][j]):
				result['global'][i] = round(glo['global'][j],7)
	return(result)


def bi(n,k) :
	"""
        Function used to calculate the binomial coefficient

		Args:
        	n and k
        
        Returns:
          the binomial coeffcient of (n,k)
	"""
	return math.comb(n,k)

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
	l_nj = []
	i_nj = []
	i_nkj = []
	i_nk = []
	l_nkj = []
	l_nk = []
	ln_nkj = []
	lv_nkj = []
	table_nj = pd.DataFrame(columns = ['variables','modalities','nj'])
	table_nkj = pd.DataFrame(columns = [vchi2,'variables','modalities','nkj'])
	table_nk = pd.DataFrame(columns = [vchi2,'nk'])
	#table is a dataframe that contain all the numbers we need for the p-value 
	#calcul
	global table
	table = pd.DataFrame({
		vchi2 : chi2_var,
		'variables' : variables ,
		'modalities':modality, 
		'nj' : NA, 
		'nkj' : NA,
		'nk' : NA})
	n = new_df.shape[0]#number of individuals
	for c in new_df :
		nj = pd.value_counts(new_df[c]) #size modality for all individuals
		for i in nj :
			l_nj.append(i)
		list_index_nj = nj.index.tolist()
		i_nj.append(list_index_nj)
		l_index_nj = []
		for i in i_nj :
			for j in i :
				v = str(j)+(';')
				v = v.split(';')
				l_index_nj.append(v)
		index_nj = []
		for i in l_index_nj :
			index_nj = index_nj+i
		while '' in index_nj :
			del index_nj[index_nj.index('')]
	for h,i,j in zip(variables_2,index_nj,l_nj) :	
		
		table_nj = table_nj.append({'variables' : h,
				   			  'modalities' : i,
				   			  'nj' : j},
				   			  ignore_index=True)
		table_nj = table_nj.fillna(0)
	for v in dictio :
		i_nk.append(v)
		l_nk.append(dictio[v].shape[0])
	for i,j in zip(i_nk, l_nk) :
		table_nk = table_nk.append({vchi2 : i, 
									'nk' : j},
									ignore_index=True) 
	for c in new_df :
		for v in dictio :
			nkj = pd.value_counts(dictio[v][c]) #size modality for the clusters
			nk = dictio[v].shape[0]#size of the clusters
			for i in nkj :
				l_nkj.append(i)
				ln_nkj.append(nkj.name)
				lv_nkj.append(v)
			li_nkj = nkj.index.tolist()
			i_nkj.append(li_nkj)
			l_index_nkj = []
			for i in i_nkj :
				for j in i :
					v = str(j)+(';')
					v = v.split(';')
					l_index_nkj.append(v)
	
				index_nkj = []
				for i in l_index_nkj :
					index_nkj = index_nkj+i
				while '' in index_nkj :
					del index_nkj[index_nkj.index('')]
			
	for g,h,i,j in zip(lv_nkj,ln_nkj,index_nkj,l_nkj) :
		table_nkj = table_nkj.append({vchi2 : g, 
						              'variables' : h,
						              'modalities' : i,
						              'nkj' : j },
						               ignore_index=True)
		table_nkj = table_nkj.fillna(0)
	
	#attribute the nkj number from table_nkj to table
	for i in range(len(table)) :	
		for j in range(len(table_nkj)) :	
			if str(table[vchi2][i]) == str(table_nkj[vchi2][j]) and str(table['variables'][i]) == str(table_nkj['variables'][j]) and str(table['modalities'][i]) == str(table_nkj['modalities'][j]):
				table['nkj'][i] = table_nkj['nkj'][j]
	
		#attribute the nj number from table_nj to table		
		for j in range(len(table_nj)):
			if str(table['variables'][i]) == str(table_nj['variables'][j]) and str(table['modalities'][i]) == str(table_nj['modalities'][j]):
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
			if str(result[vchi2][i]) == str(table[vchi2][j]) and str(result['variables'][i]) == str(table['variables'][j]) and str(result['modalities'][i]) == str(table['modalities'][j]):
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
					s1 = s1+p
					
				#if aux2 < aux3
				for y in range(nkj,nj):
					nm = table['nj'][j]-y
					p = (bi(nk,y)*bi(n_nk,nm))/bi(n,nj)
					s2 = s2+p
				result['p-value'][i] = 2*min(s1,(1-round(s2,7))) + (bi(nk,nkj)*bi(n_nk,nj-nkj))/bi(n,nj)
	return result

def vtest(result, vchi2, v_p_value) :
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
			#if result['p-value'][i] != 'Not present':
			if result['p-value'][i] < v_p_value :
				result['v-test'][i] = round((1-2*val)*norm.ppf(result['p-value'][i]/2),7)
			else :
				result['v-test'][i] = 'Not significant'
		
		if result['v-test'][i] != 'Not significant':
			if result['v-test'][i] < 0 :
				result['signification'][i] ='underrepresented'
			elif result['v-test'][i] > 0 :
				result['signification'][i] ='overrepresented'
		elif result['mod/cla'][i] == 0 :
			result['signification'][i] = 'Not present'
		elif result['v-test'][i] == 'Not significant' :
			result['signification'][i] ='Not significant'
	return result
	

def variable_weight(result):
	"""
    	Function used to know the weight of each variable in the clustering

    	Actions performed:

        * calcule the weight of a variable counting the number of modality
          overrepresented in this variable divided by the number of modality
          of the variable
	
		Args: the qualitative_analysis table
          
        Returns:
  		 Dataframe : a panda DataFrame containing the chi2 statistics analysis 
          	         for the variables and the weight of these variables.	
	
	"""
	over = []
	under = []
	ov_un = []
	mod = []
	for i in range(len(result)) :
		if result['signification'][i] == 'overrepresented' :
			over.append(result['variables'][i])
		if result['signification'][i] == 'underrepresented' :
			under.append(result['variables'][i])
		if result['signification'][i] == 'overrepresented' or result['signification'][i] == 'underrepresented' :
			ov_un.append(result['variables'][i])
		mod.append(result['variables'][i])
	var = set(over)
	var = list(var)
	var.sort()
	over.sort()
	ov_un.sort()
	nb_mod_over = []
	nb_mod_under = []
	nb_mod_overunder = []
	nb_var = []
	
	for element in var :
		nb = over.count(element)
		nb_mod_over.append(nb)
	
	for element in var :
		nb = under.count(element)
		nb_mod_under.append(nb)

	for element in var :
		nb = ov_un.count(element)
		nb_mod_overunder.append(nb)
	
	for element in var :
		nb = mod.count(element)
		nb_var.append(nb)
	
	ratio1 = []
	for i,j in zip (range(len(nb_mod_over)),range(len(nb_var))) :
			ratio = nb_mod_over[i]/nb_var[j]*100
			ratio1.append(round(ratio,2))
	
	ratio2 = []
	for i,j in zip(range(len(nb_mod_under)),range(len(nb_var))):
		ratio = nb_mod_under[i]/nb_var[j]*100
		ratio2.append(round(ratio,2))
		
	ratio3 = []
	for i,j in zip(range(len(nb_mod_overunder)),range(len(nb_var))):
		ratio = nb_mod_overunder[i]/nb_var[j]*100
		ratio3.append(round(ratio,2))
	
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
						   'sum of all mod of all clusters' : nb_var[0:-1],
						   'ratio over/mod' : ratio1[0:-1],
						   'contribution over/mod' : ranking1,
						   'ratio under/mod' : ratio2[0:-1],
						   'contribution under/mod' : ranking2,
						   'ratio over&under/mod' : ratio3[0:-1],
						   'contribution over&under/mod' : ranking3})
	return weight

def sdquanti(df, var, variable_cat, pvalue_var_cat):
	"""
        Actions performed:
        * Make the v-test and final statistics
	
		Args:
        	df: a pandas DataFrame 
        	var : the quantitative variables
        	variable_cat : the variable to test
        
        Returns:
          DataFrame: a pandas DataFrame containing statistics analysis for each
          	         variable_cat and variables
        """
	#separation of the dataframe for each variable chi2 with 
	#a dictionary of variable_cat
	df1 = df.dropna()
	global dictionary
	dictionary = {}
	for i in df[variable_cat].unique():
		dictionary[i] = df[  df[variable_cat] == i ]
		#order the dictionary by the number of variable_cat
		dictionary = OrderedDict(sorted(dictionary.items(), key=lambda x: x[0]))
	
	#a dictionary that don't contain missing values
	global dictionary_1
	dictionary_1 = {}
	for i in df1[variable_cat].unique():
		dictionary_1[i] = df1[  df1[variable_cat] == i ]
		#order the dictionary_1 by the number of variable_cat
		dictionary_1 = OrderedDict(sorted(dictionary_1.items(), key=lambda x: x[0]))
	
	#anova
	lm = ols('df1[var] ~ C(df1[variable_cat])', data = df1).fit()
	table = sm.stats.anova_lm(lm)
	table.to_csv('table.csv')
	eta = table['sum_sq'][0]/(table['sum_sq'][0]+table['sum_sq'][1])
	anova = pd.DataFrame({'variable' : var,
						  'Eta2' : eta,
						  'p-value' : table['PR(>F)'][0]},index=[0])
	#print(anova)
	
	
	
	
	
	
	#Overall mean : x
	I = len(df)
	Im = len(df1)
	x= round(df[var].mean(),6)
	somme=0
	for i,j in zip(dictionary_1,dictionary) :
		Iq_df1 = len(dictionary_1[i])
		xq = round(dictionary[j].mean(),2)
		xq=xq[0].astype(type('float'))
		data = dictionary_1[i]
		
		sous = np.empty((0,2), float)
		for el in range(len(data)) :
			data = data.reset_index(drop=True)
			data=data.astype(type('float'))
			sous=np.append(sous,np.array([[float(data[var][el]), float(xq)]]),axis=0)
		som = 0
		for i in sous :
			som = som+(i[0]-i[1])*(i[0]-i[1])
		somme = somme+som
		et_category = som/Iq_df1
		print(et_category)
		
	et_overall = somme/Im
	for i in dictionary :
		Iq = len(dictionary[i])
		xq = round(dictionary[i].mean(),6)
		xq=xq[0].astype(type('float'))
		data = dictionary[i].dropna()
		sous = np.empty((0,2), float)
		for el in range(len(data)) :
			data = data.reset_index(drop=True)
			data=data.astype(type('float'))
			sous=np.append(sous,np.array([[float(data[var][el]), float(x)]]),axis=0)
		for i in sous:
			vtest = (i[0]-i[1])/sqrt(((et_overall**2)/Iq)*((I-Iq)/(I-1)))
			pvalue = (1-norm.cdf(abs(vtest)))*2
		#print(vtest)
		#print(pvalue)
		
	
if __name__ =="__main__" : 
	#make a pandas for the input file
	#form of the pandas
	#variable1    variable2    variable 3   variable x    vchi2

	#take the different variable from the first table
	df1 = pd.read_csv(r'input_data_file.csv', sep =';', encoding='latin-1')
	
	#quantitatives variables
	quantitative =['Name (original)','Number of flowers by volume']
	#qualitatives variables
	qualitative = ['Name (original)',
		      	   'Breeding period',
		      	   'Geographic origin',
		       	   'Horticultural group',
		      	   'Ploidy',
		       	   'Bush height',
		       	   'Type',
		           'Quantity of prickles',
		       	   'Perfume intensity',
		       	   'Continious flowering',
		           'Duplicature',
   		       	   'Petal color']
   	
	df_quali = df1[qualitative]
	df_quanti = df1[quantitative]
	
	#take the variable cluster from the second table 
	#with a merge from df1 to df2
	#merge : df1 = Name(original) ; df2 = Unnamed: 0
	df2 = pd.read_csv(r'semantic_cluster_coordinates7.csv', 
	sep =';', encoding='latin-1')
	#rename the df2 columns from Unnamed: 0 to Name(original) to make the merge
	df2.rename(columns={'Unnamed: 0' : 'Name (original)'}, inplace = True)
	columns_df2 = ['Name (original)','cluster']
	df2 = df2[columns_df2]
	
	#make the dataframe that contain only the qualitatives variables
	dataframe_quali = df_quali.merge(df2)
	dataframe_quali = dataframe_quali.fillna('missing values')
	dataframe_quali.to_csv('df_qualitative.csv')
	
	#make the dataframe that contain only the quantitatives variables
	dataframe_quanti = df_quanti.merge(df2)
	dataframe_quanti = dataframe_quanti.rename(columns = {'Number of flowers by volume' : 'Number_of_flowers_by_volume'})
	dataframe_quanti.to_csv('df_quantitative.csv')
	
	#make the qualitative analysis
	#sdquali = sdquali(dataframe_quali, qualitative, 'cluster', 0.05)
	#qa = quali_analysis(dataframe_quali, qualitative, 'cluster')
	#cm = clamod(qa,'cluster')
	#mc = modcla(qa,'cluster')
	#g = globa(qa)
	#pv = pvalue(qa,'cluster')
	#vtest = vtest(qa,'cluster',0.05)
	#w = variable_weight(qa)

	#out :
	#sdquali=sdquali.rename_axis('file : florhige_synthese_english_7.05.22, code : sdquali_10.05.22')
	#sdquali.to_csv('x2_semantic_cluster7.csv')	
	#vtest=vtest.rename_axis('file : florhige_synthese_english_7.05.22, code : sdquali_10.05.22')
	#vtest.to_csv("semantic_cluster7_sdquali.csv")
	#w=w.rename_axis('file : florhige_synthese_english_7.05.22, code : sdquali_10.05.22')
	#w.to_csv('semantic_cluster7_weight.csv')	
		
	#make the quantitative analysis for each quantitative variable
	quanti = sdquanti(dataframe_quanti,'Number_of_flowers_by_volume', 'cluster',0.05)
