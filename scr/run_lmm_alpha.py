# %%
import pandas as pd
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf

class Dataset():
	def __init__(self, asv_path: str, genus_path: str, meta_path: str) -> None:
		try:
			self.asv_df = pd.read_csv(asv_path, index_col=0, header=1, sep='\t')
			self.genus_df = pd.read_csv(genus_path, index_col=0, header=1, sep='\t')
			self.meta = pd.read_csv(meta_path, index_col=0)
			self.meta.replace('none','T00000_none',inplace=True)
		except:
			print('No tables')
		
		meta_cols = self.meta.columns.tolist()
		meta_cols = ['mom_smoke' if x == 'mom.smoke' else x for x in meta_cols]
		self.meta.columns = meta_cols

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
		
		meta_collect_season2 = meta.collection_season2.tolist()
		meta_collect_season2 = ['aasummer' if x == 'summer' else x for x in meta_collect_season2]
		meta.collection_season2 = meta_collect_season2
		
		sel_tps = self.bins[bin]
		meta = meta[meta.timepoint_month.isin(sel_tps)]

		sel_samples = meta.index.tolist()

		if table_type == 'asv':
			table = self.asv_df
			table = table[sel_samples]
			tax = pd.read_csv('/data/ANCOM2_Shantelle/data/taxonomy.tsv', sep='\t')
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


ds = Dataset('./data/asv/feature-table.tsv', './data/genus/feature-table.tsv', './data/labdata_metadata_Bv4_intervalPIDs_complete.csv')

OUTPUT_PATH = 'Results_alpha_lmm'
os.system(f'mkdir {OUTPUT_PATH}')

for i in [0,1,2]:
	bin_name = str(i+1)
	_, meta_df, covs, adjs = ds.gen(i)
	
	for div in ['shannon','chao1']:
		output_df = pd.DataFrame()
		output_path = os.path.join(OUTPUT_PATH,f'{div}_bin_{bin_name}.csv')
		for co in covs:
			ad_covs = '+'.join(adjs)

			form = f"{div} ~ {co}+timepoint_month+{ad_covs}"
			groups = "short_PID"
			
			md = smf.mixedlm(form, meta_df, groups=meta_df[groups])
			mdf = md.fit()
			out = mdf.summary().tables[1]

			out = out.loc[out.index.str.startswith(co)]
			out['formula'] = [form for x in range(out.shape[0])]
			out['groups'] = [groups for x in range(out.shape[0])] 
			
			output_df = pd.concat([output_df, out])
		output_df.to_csv(output_path)	
	

