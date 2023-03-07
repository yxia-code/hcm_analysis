# %%
import pandas as pd
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf

class Dataset():
	def __init__(self, asv_path: str, genus_path: str, meta_path: str, bi_path) -> None:
		try:
			self.asv_df = pd.read_csv(asv_path, index_col=0, header=1, sep='\t')
			self.genus_df = pd.read_csv(genus_path, index_col=0, header=1, sep='\t')
			self.meta = pd.read_csv(meta_path, index_col=0,sep='\t')
			self.meta['short_PID'] = [x.split('.')[1] for x in self.meta.index]
			self.meta['timepoint_month'] = [x.split('.')[5].split('n')[0] for x in self.meta.index]
			self.meta.replace('none','T00000_none',inplace=True)
			self.bi = pd.read_csv(bi_path,index_col=0,sep='\t',header=None)
			self.bi.index = [str(int(x)) for x in self.bi.index]
		except:
			print('No tables')
		meta_cols = self.meta.columns.tolist()
		meta_cols = ['mom_smoke' if x == 'mom.smoke' else x for x in meta_cols]
		self.meta.columns = meta_cols

		self.metircs = ['']
		self.bins = [('1m')]

		self.cov_0 = ['delivery', 'edu', 'excl_bf_maxtimepoint',\
							'gestation2', 'HAZ12_2', 'HIV', 'income2',\
							'mom_smoke', 'pets', 'sex', 'siblings2', 'WAZ0_2',\
								'WAZ12_2']
		self.adj_0 = ['']
		
		self.cov_1 = ['delivery', 'edu', 'excl_bf_maxtimepoint',\
							'gestation2', 'HAZ12_2', 'HIV', 'income2',\
							'mom_smoke', 'pets', 'sex', 'siblings2', 'WAZ0_2',\
								'WAZ12_2']
		self.adj_1 = ['']
		
		self.cov_2 = ['delivery', 'edu', 'excl_bf_maxtimepoint',\
							'gestation2', 'HAZ12_2', 'HIV', 'income2',\
							'mom_smoke', 'pets', 'sex', 'siblings2', 'WAZ0_2',\
								'WAZ12_2']
		self.adj_2 = ['']
		
		self.add_cols = ['short_PID']

	def gen(self, bin: int, table_type = 'genus') -> tuple:
		covs = eval('self.cov_{}'.format(bin)) 
		adjs = eval('self.adj_{}'.format(bin))
		
		sel_cols = covs + self.add_cols 
		sel_cols = list(set(sel_cols))
		meta = self.meta[sel_cols] 
		
		#meta_collect_season2 = meta.collection_season2.tolist()
		#meta_collect_season2 = ['aasummer' if x == 'summer' else x for x in meta_collect_season2]
		#meta.collection_season2 = meta_collect_season2
		
		sel_tps = self.bins[bin]
		#meta = meta[meta.timepoint_month.isin(sel_tps)]

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

		return (table, meta, covs, '',self.bi)


ds = Dataset('../data/asv/feature-table.tsv', '../data/genus/feature-table.tsv','./trajectory.tsv' ,'./biplots_1to12m.tsv')

OUTPUT_PATH = 'Results_alpha_lm_PCs_1to12m'
os.system(f'mkdir {OUTPUT_PATH}')

for i in [0]:
	bin_name = str(i+1)
	_, meta_df, covs, adjs, bi= ds.gen(i)
	meta_df.index = meta_df.short_PID
	meta_df = meta_df.drop_duplicates()
	meta_df = meta_df[~meta_df.index.duplicated()]
	bi.columns = ['PC1','PC2','PC3']
	meta_df = pd.concat([meta_df,bi],axis=1)
	for div in ['PC1','PC2','PC3']:
		output_df = pd.DataFrame()
		output_path = os.path.join(OUTPUT_PATH,f'{div}.csv')
		for co in covs:
			form = f"{div} ~ {co}"
			groups = "short_PID"
			try:	
				md = smf.mixedlm(form, meta_df, groups=meta_df[groups])
			except:
				continue
			mdf = md.fit()
			out = mdf.summary().tables[1]

			out = out.loc[out.index.str.startswith(co)]
			out['formula'] = [form for x in range(out.shape[0])]
			#out['groups'] = [groups for x in range(out.shape[0])] 
			
			output_df = pd.concat([output_df, out])
		output_df.to_csv(output_path)	
	

