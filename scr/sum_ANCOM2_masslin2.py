import pandas as pd
import os
import glob 

def gen(lev,bin):
    p11 = f'./Results/ANCOM2/{lev}/{bin}/*/ancom_results.csv'
    p21 = f'./Results/Maaslin2_revised/{lev}/{bin}/*/all_results.tsv'
    OUT = './Results/sum_results'
    os.system(f'mkdir {OUT}')
    OUT1 = f'./Results/sum_results/{lev}'
    os.system(f'mkdir {OUT1}')
    OUT1 = f'./Results/sum_results/{lev}/{bin}'
    os.system(f'mkdir {OUT1}')
    return p11, p21, OUT1


def run(lev, bin):
    p11, p21, OUT1 = gen(lev,bin)
    
    lst = glob.glob(p11)
    lst2 = glob.glob(p21)
    for i in lst:
        cov = i.split('/')[-2]
        df1 = pd.read_csv(i,index_col=0)
        df1 = df1[['detected_0.6','W']]
        df1.columns=['ANCOM_W_cutoff_0.6','W']
        df11 = df1.copy()
        
        i2 = [x for x in lst2 if cov in x][0]
        df2 = pd.read_csv(i2,index_col=0,sep='\t')
        df2 = df2[df2.metadata==cov]
        df22 = df2.copy()
        df2 = df2[df2.qval<0.05]

        df1 = df1.loc[df2.index]

        df = pd.concat([df2,df1],axis=1) 
        df11 = df11[df11['ANCOM_W_cutoff_0.6']==True]
        df22 = df22.loc[df11.index]
        if df11.shape[0]==0 or df22.shape[0]==0:
            dfff = df
        else:
            
            dff = pd.DataFrame(columns=['feature']+df22.columns.tolist()+df11.columns.tolist())
            for j in df22.index:
                s2 = df22.loc[j]
                s1 = df11.loc[j]
                sum = 0
                if len(s2.shape)==2 and s2.shape[0] > 1:
                    for n,k in enumerate(s2.index):
                        l1=s2.iloc[n].tolist()
                        l2=s1.tolist()
                        l = [k]+l1+l2
                        dff.loc[sum] = l
                        sum += 1
                        
                elif len(s2.shape)==1: 
                    s1s2 = pd.concat([s2,s1])
                    s1s2 = s1s2.tolist()
                    l = [s1.name] + s1s2
                    dff.loc[sum] = l
                    sum += 1

                else:
                    pass
            dff.index=dff.feature
            dff.drop('feature',axis=1,inplace=True)
            dfff = pd.concat([df,dff],axis=0)
        dfff.to_csv(os.path.join(OUT1, f'{cov}.csv'))
    return 

run('asv','Bin_1')
run('asv','Bin_2')
run('asv','Bin_3')
run('genus','Bin_1')
run('genus','Bin_2')
run('genus','Bin_3')

