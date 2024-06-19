# Data analysis and manipulation library
import pandas as pd
import plotly.express as px
import os

# System library to manipulate the file system
from os import path
from scripts.utils import write_excel
import shutil
from tqdm import tqdm
from scripts.quads import *

df = pd.ExcelFile('data/gower_cluster_coordinates5.xlsx')
sheets = df.sheet_names
df1 = pd.read_excel('data/input_data_file.xlsx')
for sheet in tqdm(sheets) :
    df2 = pd.read_excel(df, sheet)
    #quantitatives variables
    quantitative =['Name (original)','Number of flowers per inflorescence']
    #qualitatives variables
    qualitative = ['Name (original)',
                   'Breeding period',
                   'Geographic origin',
                   'Horticultural group',
                   'Ploidy',
                   'Bush height',
                   'Shape',
                   'Quantity of prickles',
                   'Perfume intensity',
                   'Repeat flowering',
                   'Quantity of bristles by branch',
                   'Shine of upper face',
                   'Corolla form',
                   'Corolla size',
                   'Color repartition',
                   'Duplicature',
                   'Petal color']
   	
    df_quali = df1[qualitative]
    df_quanti = df1[quantitative]

    #take the variable cluster from the second table 
    #with a merge from df1 to df2
    #merge : df1 = Name(original) ; df2 = Unnamed: 0
    #rename the df2 columns from Unnamed: 0 to Name(original) to make the merge
    df2.rename(columns={'Unnamed: 0' : 'Name (original)'}, inplace = True)
    columns_df2 = ['Name (original)','cluster']
    df2 = df2[columns_df2]
	
    #make the dataframe that contain only the qualitatives variables
    dataframe_quali = df_quali.merge(df2)
    dataframe_quali = dataframe_quali.fillna('missing values')
	
    #make the dataframe that contain only the quantitatives variables
    dataframe_quanti = df_quanti.merge(df2)
    dataframe_quanti = dataframe_quanti.rename(columns = {'Number of flowers per inflorescence' : 'Number_of_flowers_per_inflorescence'})

    #make the qualitative analysis
    sdqualitative = sdquali(dataframe_quali, qualitative, 'cluster', 0.05)
    sdqualitative=sdqualitative.rename_axis('file : 20220615_florhige_synthese_english, code : 20220615_quads')
    quali_a = quali_analysis(dataframe_quali, qualitative, 'cluster')
    cm = clamod(quali_a,'cluster')
    mc = modcla(quali_a,'cluster')
    g = globa(quali_a)
    pv = pvalue(quali_a,'cluster')
    test_value = vtest(quali_a,'cluster',0.05)
    test_value=test_value.rename_axis('file : 20220615_florhige_synthese_english, code : 20220615_quads')
    w = variable_weight(quali_a)
    w=w.rename_axis('file : 20220615_florhige_synthese_english, code : 20220615_quads')	
	
    #make the quantitative analysis for each quantitative variable
    sd = sdquanti(dataframe_quanti,'Number_of_flowers_per_inflorescence', 'cluster')
    sd = sd.rename_axis('file : 20220615_florhige_synthese_english, code : 20220615_quads')
    quanti_a = quanti_analysis(sd, dataframe_quanti,'Number_of_flowers_per_inflorescence', 'cluster',0.05,0.05)
    quanti_a = quanti_a.rename_axis('file : 20220615_florhige_synthese_english, code : 20220615_quads')
			
    #out :
    #create the new path for the result
    if not os.path.exists('results/gower/cluster5') :
        os.makedirs('results/gower/cluster5')
    path = 'results/gower/cluster5/'
	
    #name the files
    file_name_x2 = 'x2_gower_cluster5.xlsx'
    file_name_qualitative = 'qualitative_analysis_gower_cluster5.xlsx'
    file_name_weight = 'weight_gower_cluster5.xlsx'
    file_name_anova = 'anova_gower_cluster5.xlsx'
    file_name_quantitative = 'quantitative_analysis_gower_cluster5.xlsx'
	
    #create the excel files
    write_excel(file_name_x2, sheet, sdqualitative, idx=True)
    write_excel(file_name_qualitative, sheet, test_value,idx=True)
    write_excel(file_name_weight, sheet, w,idx=True)
    write_excel(file_name_anova, sheet, sd,idx=True)
    write_excel(file_name_quantitative, sheet, quanti_a, idx=True)

#make the visualisations
data = pd.ExcelFile(file_name_qualitative)
sheets = data.sheet_names
col = {'overrepresented' : 'red', 'underrepresented' : 'blue', 'Not significant': 'grey'}
for sheet in sheets :
	title = 'Proportions of modalities in each clusters with Gower distance and '+sheet+' method'
	df = pd.read_excel(data, sheet)
	legend=''
	for i in range (len(df)):
		if legend == '' :
			pass
		else : 
			legend = legend+' ; '
		if df['variables'][i] =='cluster' :
			legend= legend+ str(df['cluster'][i])+' : '+str(round(df['global'][i],2))+'%'
	sunburst = px.sunburst(df, path=['cluster', 'variables', 'modalities'],values='mod/cla',title=title, color = 'signification',color_discrete_map=col)
	sunburst.add_annotation(x=0,y=1.1,text=legend,font = dict(color='black',size=14),showarrow=False)
	sunburst.add_annotation(x=0.2,y=1,text= 'Overrepresented',font = dict(color='red',size=14),showarrow=False)
	sunburst.add_annotation(x=0.2,y=0.95,text= 'Underrepresented',font = dict(color='blue',size=14),showarrow=False)
	sunburst.add_annotation(x=0.2,y=0.9,text= 'Not significant',font = dict(color='grey',size=14),showarrow=False)
	sunburst.show()

#move the files in the good directory	
shutil.move(file_name_x2,path+file_name_x2)
shutil.move(file_name_qualitative,path+file_name_qualitative)
shutil.move(file_name_weight,path+file_name_weight)
shutil.move(file_name_anova,path+file_name_anova)
shutil.move(file_name_quantitative,path+file_name_quantitative)
