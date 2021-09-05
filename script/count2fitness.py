import os
import glob
import pandas as pd


def cal_enri(count_DF,column_name,input):
    count_DF[column_name+'_enri'] = ((count_DF[column_name]) / count_DF.groupby('Amplicon')[column_name].transform('sum'))/(((count_DF[input])+1)/count_DF.groupby('Amplicon')[input].transform('sum'))
    return count_DF

def cal_fitness(count_DF,enri_ls,input_DNA):
    for col in enri_ls:
        count_DF=cal_enri(count_DF, col, input_DNA)
    #calculate the fitness and rename
    count_DF.loc[:, enri_ls[0]+'_enri':] = count_DF.loc[:, enri_ls[0]+'_enri':].div(count_DF.iloc[-1][enri_ls[0]+'_enri':])
    count_DF=count_DF.rename(columns={col+'_enri': col+'_fit' for col in enri_ls})

    return count_DF


def cal_mean(count_DF,mean,column1,column2):
    count_DF[mean] = (count_DF[column1]+count_DF[column2])/2
    return count_DF

def grab_files_as_Dict(path,suffix):
    #path = 'result/*/*'
    #pattern = '_count.tsv'
    files = glob.glob(path+suffix)
    #filenames = [os.path.basename(i).split(suffix, 1)[0] for i in files]
    file_dict={}
    for file in files:
        filename=os.path.basename(file).split(suffix, 1)[0]
        file_dict[filename]=file
    return file_dict

def main():
    path = 'result/*/*'
    suffix = '_count.tsv'
    file_dict=grab_files_as_Dict(path,suffix)
    DF=pd.DataFrame(columns=['Mutation'])
    for filename,file in file_dict.items():
        df= pd.read_csv(file, header=0, sep= '\t')
        df.columns=['Mutation',filename]
        DF=pd.merge(DF,df,how='outer',on='Mutation')
    DF.to_csv('result/Mos99_all_count.csv')
    # split amplicon and mutation
    DF[['Amplicon','Mutation']] = DF.Mutation.str.split("|", n=1, expand=True)
    DF=DF.fillna(0)
    DF['WT_plasmid' + '_norm'] = ((DF['WT_plasmid']) / DF.groupby('Amplicon')['WT_plasmid'].transform('sum'))
    DF['Rep1' + '_norm'] = ((DF['Rep1']) / DF.groupby('Amplicon')['Rep1'].transform('sum'))
    DF['Rep2' + '_norm'] = ((DF['Rep2']) / DF.groupby('Amplicon')['Rep2'].transform('sum'))
    DF['Input' + '_norm'] = ((DF['Input']) / DF.groupby('Amplicon')['Input'].transform('sum'))
    DF['Rep1' + '_enri'] = DF['Rep1' + '_norm']/DF['Input' + '_norm']
    DF['Rep2' + '_enri'] = DF['Rep2' + '_norm'] / DF['Input' + '_norm']
    # enri_ls=['R_1_DMS_RNA','R_2_DMS_RNA','R_3_DMS_RNA','Y_1_DMS_RNA','Y_2_DMS_RNA','Y_3_DMS_RNA','WT_1_DMS_RNA','WT_2_DMS_RNA','WT_3_DMS_RNA']
    # DF=cal_fitness(DF, enri_ls, 'HA_DMS_Plasmid')
    DF.to_csv('result/Moss99_fit_raw.csv')
    # onemut_df = DF[~DF["Mutation"].str.contains('-')]
    # pos_df = onemut_df.Mutation.str.extract('(\d+)')
    # one_mut_sortdf = onemut_df.join(pos_df, lsuffix='_caller', rsuffix='_other') \
    #     .set_index(0) \
    #     .fillna(0) \
    #     .sort_index()
    #
    # correlation = one_mut_sortdf.corr(method='pearson')
    # print(one_mut_sortdf.head(5))
    # correlation.to_csv('result/HA1_mut_correlation.csv')
    # one_mut_sortdf.to_csv('result/HA1_fit.csv')


if __name__ == "__main__":
    main()