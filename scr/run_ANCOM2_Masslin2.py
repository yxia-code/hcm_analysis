import os
import pandas as pd 
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import Formula
import os

def pd2ri(df):
	with localconverter(ro.default_converter + pandas2ri.converter):
		r_df = ro.conversion.py2rpy(df)
	return r_df
def ri2pd(df):
	with localconverter(ro.default_converter + pandas2ri.converter):
		r_df = ro.conversion.rpy2py(df)
	return r_df
def lib_rlib():
	global base, maaslin2
	maaslin2 = importr('Maaslin2')
	base = importr('base')
	return 

lib_rlib()


class Dataset():
	def __init__(self, asv_path: str, genus_path: str, meta_path: str) -> None:
		try:
			self.asv_df = pd.read_csv(asv_path, index_col=0, header=1, sep='\t')
			self.genus_df = pd.read_csv(genus_path, index_col=0, header=1, sep='\t')
			self.meta = pd.read_csv(meta_path, index_col=0)
		except:
			print('No tables')

		self.metircs = ['shannon', 'chao1']
		self.bins = [('M01', 'M02', 'M03'), ('M04', 'M05', 'M06'), ('M07', 'M08', 'M09', \
			'M10', 'M11', 'M12')]

		self.cov_0 = ['delivery', 'edu', 'excl_bf_maxtimepoint',\
							'gestation2', 'HAZ12_2', 'HIV', 'income2',\
							'mom_smoke', 'pets', 'sex', 'siblings2', 'WAZ0_2',\
								'WAZ12_2', 'collection_season2']
		self.adj_0 = ['age_days_numeric', 'collection_season2']
		
		self.cov_1 = ['delivery', 'edu', 'excl_bf_maxtimepoint',\
							'gestation2', 'HAZ12_2', 'HIV', 'income2',\
							'mom_smoke', 'pets', 'sex', 'siblings2', 'WAZ0_2',\
								'WAZ12_2', 'collection_season2', 'daycare_prior_sampling', 'first_antibiotics_upto_100days_post']
		self.adj_1 = ['age_days_numeric', 'collection_season2']
		
		self.cov_2 = ['delivery', 'edu', 'excl_bf_maxtimepoint',\
							'gestation2', 'HAZ12_2', 'HIV', 'income2',\
							'mom_smoke', 'pets', 'sex', 'siblings2', 'WAZ0_2',\
								'WAZ12_2', 'collection_season2', 'daycare_prior_sampling', 'first_antibiotics_upto_100days_post',\
									'first_hosp_upto_100days_post', 'first_INH_upto_100days_post']
		self.adj_2 = ['age_days_numeric', 'collection_season2']
		
		self.add_cols = ['shannon', 'chao1', 'short_PID', 'timepoint_month']

	def gen(self, bin: int, table_type = 'asv') -> tuple:
		covs = eval('self.cov_{}'.format(bin)) 
		adjs = eval('self.adj_{}'.format(bin))
		
		sel_cols = covs + adjs + self.add_cols
		sel_cols = list(set(sel_cols))
		meta = self.meta[sel_cols] 

		sel_tps = self.bins[bin]
		meta = meta[meta.timepoint_month.isin(sel_tps)]

		sel_samples = meta.index.tolist()

		if table_type == 'asv':
			table = self.asv_df
			table = table[sel_samples]
			tax = pd.read_csv('./data/taxonomy.tsv', sep='\t')
			dic = {}
			for t in tax.index:
				ser = tax.loc[t]
				asv_id = ser['Feature ID'] 
				tax_name = ser['Taxon']
				tax_name = tax_name.split('g__')[-1]
				tax_name = tax_name.split('; s__')
				tax_name = '_'.join(tax_name)
				dic[asv_id] = tax_name
			idx = table.index.tolist()
			idx = [x + '_' + dic[x] for x in idx]
			table.index = idx
			table.index.name = '#OTU ID'

		elif table_type == 'genus':
			table = self.genus_df
			table = table[sel_samples]
			idx = table.index.tolist()
			idx = [x.split('g__')[-1] for x in idx]
			table.index = idx
			table.index.name = '#OTU ID'

		else:
			print('No table!')
			table = None

		return (table, meta, covs, adjs)

class DEanalysis():
	def __init__(self) -> None:
		self.ds = Dataset('./data/asv/feature-table.tsv', './data/genus/feature-table.tsv', './data/labdata_metadata_Bv3.csv')

	def run(self, bin: int, table_type: str, model: str) -> None:
		table_df, meta_df, covs, adjs = self.ds.gen(bin, table_type)
		p_cutoff = 0.05
		meta_df.index.name = 'SampleID'
		s_pid = meta_df['short_PID'].tolist()
		s_pid = ['N' + str(x) for x in s_pid]
		meta_df['short_PID'] = s_pid
		
		## Filter pids
		pids = meta_df.short_PID.value_counts()
		if bin == 2: 
			pids = pids[pids==6].index
		else:
			pids = pids[pids==3].index
		meta_df = meta_df[meta_df.short_PID.isin(pids)]
		
		meta_df.drop(['shannon', 'chao1'], axis=1, inplace=True)

		table_df = table_df[meta_df.index]

		## Filter table
		for idx in table_df.index:
			tmp = table_df.loc[idx]
			if (tmp==0.0).sum() > 0.9*len(tmp):
				table_df.drop(idx, axis=0, inplace=True)
		if model == 'maaslin2':
			output_root_path = os.path.join('Results', 'Maaslin2_revised')
		if model == 'maaslin2_zicp':
			output_root_path = os.path.join('Results', 'Maaslin2_ZICP')
		elif model == 'ancom2':
			output_root_path = os.path.join('Results', 'ANCOM2')
		
		os.system('mkdir {}'.format(output_root_path))
		output_root_path = os.path.join(output_root_path, str(table_type))
		os.system('mkdir {}'.format(output_root_path))
		output_root_path = os.path.join(output_root_path, 'Bin_{}'.format(bin+1))
		os.system('mkdir {}'.format(output_root_path))
		
		table_out_path = os.path.join(output_root_path, 'table.tsv')
		meta_out_path = os.path.join(output_root_path, 'meta.tsv')
		table_df.to_csv(table_out_path, sep='\t')
		meta_df.to_csv(meta_out_path, sep='\t')
		
		table = table_out_path
		meta = meta_out_path

		for cov in ['timepoint_month','collection_season2', 'excl_bf_maxtimepoint']+covs:
			if model == 'maaslin2':
				meta_df = pd.read_csv(meta,sep='\t', index_col=0)
				feature_df = pd.read_csv(table, sep='\t',index_col=0)
				## Transpose asv table
				feature_df = feature_df.T
				meta_df_py = meta_df.copy()

				meta_df = pd2ri(meta_df)
				feature_df = pd2ri(feature_df)
				
				if cov in adjs:
					all_vars = adjs 
				else:
					all_vars = adjs + [cov]

				ref_lst = []
				for var in all_vars:
					if var == 'age_days_numeric':
						continue
					labels = meta_df_py[var].unique().tolist()
					labels.sort()
					ref_lab = labels[0]
					ref_lst.append(var)
					ref_lst.append(ref_lab)

				## Run maaslin2
				output_path = os.path.join(output_root_path, cov)
				maaslin2.Maaslin2(input_data = feature_df,
								input_metadata = meta_df,
								output = output_path,
								normalization='CLR',
								transform='NONE',
								fixed_effects = ro.StrVector(all_vars),
    							random_effects = ro.StrVector(["short_PID"]),
								reference = ro.StrVector(ref_lst))
				res_df = pd.read_csv(os.path.join(output_path,'significant_results.tsv'),sep='\t',header=0,index_col=0)
				res_df = res_df[res_df.metadata==cov]
				res_df = res_df[res_df.qval<0.05]
				res_df.to_csv(os.path.join(output_path,'final_results.csv'),index=True)
			if model == 'maaslin2_zicp':
				meta_df = pd.read_csv(meta,sep='\t', index_col=0)
				feature_df = pd.read_csv(table, sep='\t',index_col=0)
				## Transpose asv table
				feature_df = feature_df.T
				meta_df_py = meta_df.copy()

				meta_df = pd2ri(meta_df)
				feature_df = pd2ri(feature_df)
				
				if cov in adjs:
					all_vars = adjs 
				else:
					all_vars = adjs + [cov]

				ref_lst = []
				for var in all_vars:
					labels = meta_df_py[var].unique().tolist()
					labels.sort()
					ref_lab = labels[0]
					ref_lst.append(var)
					ref_lst.append(ref_lab)
				## Run maaslin2
				output_path = os.path.join(output_root_path, cov)
				maaslin2.Maaslin2(analysis_method = 'ZICP',
								input_data = feature_df,
								input_metadata = meta_df,
								output = output_path,
								normalization='TSS',
								transform='LOG',
								fixed_effects = ro.StrVector(all_vars),
    							random_effects = ro.StrVector(["short_PID"]),
								reference = ro.StrVector(ref_lst))
			elif model == 'ancom2':
				output_path = os.path.join(output_root_path, cov)
				os.system('mkdir {}'.format(output_path))

				rand_adj_scr = 'run_random_inter_adj.R'
				adj_formula = '+'.join(adjs)
				cmd = 'Rscript ./R_scr/{} {} {} {} {} {} {}'.format(rand_adj_scr, table, \
				meta, cov, p_cutoff, output_path, adj_formula)
				print(cmd)
				os.system(cmd)
			else:
				print('Wrong model name')
		
		return



os.system('mkdir Results')
func = DEanalysis()
for bin in [0,1,2]:
	for type in ['asv']:
		func.run(bin, type, 'maaslin2')
		func.run(bin, type, 'ancom2')
