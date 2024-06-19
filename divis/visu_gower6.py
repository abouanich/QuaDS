import plotly.express as px
import pandas as pd

#ind_ov = df[df['signification'] == 'overrepresented'].index
#ind_und = df[df['signification'] == 'underrepresented'].index	
#ind_np = df[df['signification'] == 'Not present'].index
#ind_ns = df[df['signification'] == 'Not significant'].index

#Overrepresented modalities data
#df_ov = df.copy()
#df_ov.drop(ind_und,inplace=True)
#df_ov.drop(ind_np,inplace=True)
#df_ov.drop(ind_ns,inplace=True)
#sunburst_over = px.sunburst(df_ov, path=['cluster', 'variables', 'modalities'], values='cla/mod')
#sunburst_over.show()

#Underrepresented modalities data
#df_und = df.copy()
#df_und.drop(ind_ov,inplace=True)
#df_und.drop(ind_np,inplace=True)
#df_und.drop(ind_ns,inplace=True)
#sunburst_under = px.sunburst(df_und, path=['cluster', 'variables', 'modalities'], values='cla/mod')
#sunburst_under.show()


file_name_qualitative = 'results/gower/cluster6/qualitative_analysis_gower_cluster6.xlsx'	
data = pd.ExcelFile(file_name_qualitative)
col = {'overrepresented' : 'red', 'underrepresented' : 'blue', 'Not significant': 'grey'}
sheets = data.sheet_names
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
	
	#treemap = px.treemap(df, path=['cluster', 'variables', 'modalities'],values='mod/cla',title=legend, color = 'signification')
	#treemap.update_traces(root_color="lightgrey")
	#treemap.update_layout(margin = dict(t=50, l=25, r=25, b=25))
	#treemap.show()


