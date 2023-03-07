import pandas as pd
import os

class Dataset():
	def __init__(self, asv_path: str, genus_path: str, meta_path: str) -> None:
		try:
			self.asv_df = pd.read_csv(asv_path, index_col=0, header=1, sep='\t')
			self.genus_df = pd.read_csv(genus_path, index_col=0, header=1, sep='\t')
			self.meta = pd.read_csv(meta_path, index_col=0)
		except:
			print('No tables')

		self.metircs = ['shannon', 'chao1']
		self.bins = [('M01', 'M02', 'M03', 'M04', 'M05', 'M06','M07','M08','M09','M10','M11','M12')]

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

ds = Dataset('./data/asv/feature-table.tsv', './data/genus/feature-table.tsv', './data/labdata_metadata_Bv3.csv')

df,meta,cov_lst,d = ds.gen(0,'asv')
meta.index.name = '#SampleID'
os.system('mkdir TMP')

dff = df.copy()
for i in dff.columns:
	ser = dff[i]
	ser_sum = ser.sum()
	ser = [x/ser_sum for x in ser]
	dff[i] = ser 
dff = dff.sum(axis=1)
dff.sort_values(inplace=True,ascending=False)
dff = dff[:20]
df = df.loc[dff.index]
df.index = [x.split('_')[1] for x in df.index]

df.to_csv('TMP/asv.tsv',sep='\t')
meta=meta.replace('M01',1)
meta=meta.replace('M02',2)
meta=meta.replace('M03',3)
meta=meta.replace('M04',4)
meta=meta.replace('M05',5)
meta=meta.replace('M06',6)
meta=meta.replace('M07',7)
meta=meta.replace('M08',8)
meta=meta.replace('M09',9)
meta=meta.replace('M10',10)
meta=meta.replace('M11',11)
meta=meta.replace('M12',12)
meta.to_csv('TMP/meta.tsv',sep='\t')

os.system('biom convert -i TMP/asv.tsv -o TMP/asv.biom --table-type=\'OTU table\' --to-hdf5')
os.system('qiime tools import --type FeatureTable[Frequency] --input-path TMP/asv.biom --output-path TMP/asv.qza')

cmd = 'qiime gemelli ctf --i-table TMP/asv.qza --m-sample-metadata-file TMP/meta.tsv \
--p-state-column timepoint_month \
--p-individual-id-column short_PID \
--output-dir Results_CTF'
os.system(cmd)

cmd = 'qiime longitudinal volatility \
--m-metadata-file \'Results_CTF/state_subject_ordination.qza\' \
--p-state-column timepoint_month --p-individual-id-column subject_id \
--p-default-group-column \'delivery\' --p-default-metric PC1 \
--o-visualization \'Results_CTF/PCs_plots.qzv\''
os.system(cmd)

COV = 'group'
PATH_COV = f'TMP/{COV}.tsv'
mf = meta.copy()
mf = mf.groupby('short_PID').agg('first')
mf.index = [float(x) for x in mf.index]
mf.index.name = '#SampleID'
mf.to_csv(PATH_COV,sep='\t')
OUT_COV = f'TMP/group_biplot.qzv'
cmd='qiime emperor biplot --i-biplot Results_CTF/subject_biplot.qza \
--m-sample-metadata-file TMP/group.tsv \
--p-number-of-features 20 --o-visualization Results_CTF/group_biplot.qzv'
os.system(cmd)

cmd = 'qiime qurro loading-plot --i-table TMP/asv.qza \
--i-ranks Results_CTF/subject_biplot.qza \
--m-sample-metadata-file TMP/meta.tsv --o-visualization Results_CTF/qurro.qzv'
os.system(cmd)

os.system('rm -rf TMP')
