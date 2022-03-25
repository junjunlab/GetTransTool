def main():
    """
    Extract longest transcript from gtf format annotation file based on gencode/ensembl/ucsc database.
    """

    # 导入模块
    import pandas as pd
    import gzip
    import time
    import argparse
    from pyfaidx import Fasta

    parser = argparse.ArgumentParser(usage="GetLongestTransFromGTF --database ensembl --gtffile Homo_sapiens.GRCh38.101.gtf.gz --genome Homo_sapiens.GRCh38.dna.primary_assembly.fa --outfile longest_trans.fa",
                                    description="Extract longest transcript from gtf format annotation file based on gencode/ensembl/ucsc database.",
                                    epilog="Thank your for your support, if you have any questions or suggestions please contact me: 3219030654@stu.cpu.edu.cn.")
    # version
    parser.add_argument('-v','--version', action='version', version='%(prog)s 0.0.7')
    # 读取注释类型文件
    parser.add_argument('-d','--database',type=str,action="store",dest="database",metavar="databse",choices=['ucsc','ensembl','gencode'],default="ensembl",
                        help='which annotation database you choose. (default="ensembl", ucsc/ensembl/gencode)')
    # 读取gtf文件
    parser.add_argument('-g','--gtffile', type=str,action="store",dest="gtffile",metavar="gtffile",
                        help='input your GTF file with ".gz" format.')
    # 读取基因组fasta文件
    parser.add_argument('-fa','--genome',type=str,action="store",dest="genome",metavar="genome",
                        help='your genome fasta file matched with your GTF file with ".fa/.fasta" format. (Homo_sapiens.GRCh38.dna.primary_assembly.fa)')
    # 导出文件名称
    parser.add_argument('-o','--outfile', type=str,action="store",dest="longestfile",metavar="longestfile",
                        help='output your longest transcript file. (longest_trans.fa)')
    # 解析参数
    args = parser.parse_args()

    # 获取参数
    db = args.database
    gtffile = args.gtffile
    genomefile = args.genome
    outfile = args.longestfile

    # main fuction
    print("Your job is running, please wait...")
    ######################################################################
    job_start = time.time()
    #######################
    if db == 'ensembl':
    ####################################################################################
        # 信息保存在字典里
        info = {}
        # 打开测试 gtf 文件
        with gzip.open(gtffile,'rt') as gtf:
            for line in gtf:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                # 分割
                fields = line.split()
                # 类型
                type = fields[2]
                if type == 'exon':
                    # 名称
                    gene_name = fields[19].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[13].replace('"','').replace(';','')
                    biotype = fields[23].replace('"','').replace(';','')
                    # 连接名称
                    key = '|'.join([gene_name,gene_id,trans_id,biotype])
                    # 计算多个外显子长度
                    start = int(fields[3])
                    end = int(fields[4])
                    length = end - start + 1
                    # 累计求和
                    info.setdefault(key,0)
                    info[key] += length
        ######################################
        # 转为数据框
        res = pd.DataFrame(pd.Series(info), columns = ['transcript_length'])

        # 添加基因名列
        res['gene_name'] = [line.split(sep='|')[0] for line in list(res.index[:])]

        # 排序
        res_sorted = res.sort_values(by = ['gene_name','transcript_length'],ascending=False)

        # 筛选最长转录本id
        longest_id = res_sorted.drop_duplicates(subset=['gene_name'],keep='first').index.values.tolist()

        # 筛选最长转录本表格
        longest_data = res_sorted.loc[res_sorted.index.isin(longest_id)]

        longest_data['ID'] = longest_data.index
        # order columns
        longest_data = longest_data[['ID','gene_name','transcript_length']]
        
        # 保存
        longest_data.to_csv(r'longest_transcripts_info.csv', index=False)

        ###########################################################
        # 储存最长转录本id  
        transid = {line.split(sep='_')[2]:line for line in list(longest_data.ID)}

        infolist = []
        with gzip.open(gtffile) as gtf:
            for line in gtf:
                # skip
                if line.startswith('#'):
                    continue
                # split
                fields = line.split()
                # feature type
                type = fields[2]
                if type == 'exon':
                    # pos
                    chr = fields[0]
                    start = fields[3]
                    end = fields[4]
                    strand = fields[6]
                    # name
                    gene_name = fields[19].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[13].replace('"','').replace(';','')
                    if trans_id in transid:
                        infolist.append([chr,start,end,strand,type,gene_name,gene_id,trans_id,transid[trans_id]])
                    else:
                        pass
                else:
                    pass

        # to dataframe
        dfinfo = pd.DataFrame(infolist,columns=['chr','start','end','strand','type','gene_name','gene_id','trans_id','id'])
        dfinfo_1_strand = dfinfo[dfinfo['strand'] == '+']

        # descrese coord by - strand gene
        dfinfo_2_strand = dfinfo[dfinfo['strand'] == '-']
        dfinfo_2_strand = dfinfo_2_strand.sort_values(by = ['trans_id','start','end'],ascending = False)

        # merge
        df_fianl = pd.concat([dfinfo_1_strand,dfinfo_2_strand],axis=0)

        ###########################################################
        # extact  sequnece from genome

        # load genome
        genome = Fasta(genomefile)

        # chrmosome info
        chrmosome_list = genome.keys()

        # save in dict
        res = {}
        for line in range(0,df_fianl.shape[0]):
            
            # chromosome strand
            fileds = df_fianl.iloc[line]
            chrom = fileds['chr']
            strand = fileds['strand']
            start = int(fileds['start'])
            end = int(fileds['end'])
            # key
            key = fileds['id']
            # filter chromoseome
            if chrom in chrmosome_list:
                # extarct sequence
                if strand == '+':
                    seq = genome[chrom][(start-1):end].seq
                elif strand == '-':
                    seq = genome[chrom][(start-1):end].complement.reverse.seq
                else:
                    pass
                # save in dict
                res.setdefault(key,'')
                res[key] += seq
            else:
                pass
        
        ###########################################################
        # 输出序列
        outputfile = open(outfile,'w')

        # fasta序列分割长度
        my_length = 60
        ###########################################################
        # 输出
        for key,val in res.items():
            outputfile.write('>' + key + '\n')
            while len(val) > my_length:
                outputfile.write(val[0:my_length] + '\n')
                val = val[my_length:len(val)]
            outputfile.write(val + '\n')

        # 关闭文件
        outputfile.close()
    ####################################################################################
    elif db == 'gencode':
        # 信息保存在字典里
        info = {}
        # 打开测试 gtf 文件
        with gzip.open(gtffile,'rt') as gtf:
            for line in gtf:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                # 分割
                fields = line.split()
                # 类型
                type = fields[2]
                if type == 'exon':
                    # 名称
                    gene_name = fields[15].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[11].replace('"','').replace(';','')
                    biotype = fields[13].replace('"','').replace(';','')
                    # 连接名称
                    key = '|'.join([gene_name,gene_id,trans_id,biotype])
                    # 计算多个外显子长度
                    start = int(fields[3])
                    end = int(fields[4])
                    length = end - start + 1
                    # 累计求和
                    info.setdefault(key,0)
                    info[key] += length

        ######################################
        # 转为数据框
        res = pd.DataFrame(pd.Series(info), columns = ['transcript_length'])

        # 添加基因名列
        res['gene_name'] = [line.split(sep='|')[0] for line in list(res.index[:])]

        # 排序
        res_sorted = res.sort_values(by = ['gene_name','transcript_length'],ascending=False)

        # 筛选最长转录本id
        longest_id = res_sorted.drop_duplicates(subset=['gene_name'],keep='first').index.values.tolist()

        # 筛选最长转录本表格
        longest_data = res_sorted.loc[res_sorted.index.isin(longest_id)]

        longest_data['ID'] = longest_data.index
        # order columns
        longest_data = longest_data[['ID','gene_name','transcript_length']]
        
        # 保存
        longest_data.to_csv(r'longest_transcripts_info.csv', index=True)

        ###########################################################
        # 储存最长转录本id  
        transid = {line.split(sep='|')[2]:line for line in list(longest_data.ID)}

        infolist = []
        with gzip.open(gtffile,'rt') as gtf:
            for line in gtf:
                # skip
                if line.startswith('#'):
                    continue
                # split
                fields = line.split()
                # feature type
                type = fields[2]
                if type == 'exon':
                    # pos
                    chr = fields[0]
                    start = fields[3]
                    end = fields[4]
                    strand = fields[6]
                    # name
                    gene_name = fields[15].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[11].replace('"','').replace(';','')
                    if trans_id in transid:
                        infolist.append([chr,start,end,strand,type,gene_name,gene_id,trans_id,transid[trans_id]])
                    else:
                        pass
                else:
                    pass

        # to dataframe
        dfinfo = pd.DataFrame(infolist,columns=['chr','start','end','strand','type','gene_name','gene_id','trans_id','id'])
        dfinfo_1_strand = dfinfo[dfinfo['strand'] == '+']

        # descrese coord by - strand gene
        dfinfo_2_strand = dfinfo[dfinfo['strand'] == '-']
        dfinfo_2_strand = dfinfo_2_strand.sort_values(by = ['trans_id','start','end'],ascending = False)

        # merge
        df_fianl = pd.concat([dfinfo_1_strand,dfinfo_2_strand],axis=0)
        ###########################################################

        # extact  sequnece from genome

        # load genome
        genome = Fasta(genomefile)

        # chrmosome info
        chrmosome_list = genome.keys()

        # save in dict
        res = {}
        for line in range(0,df_fianl.shape[0]):
            
            # chromosome strand
            fileds = df_fianl.iloc[line]
            chrom = fileds['chr']
            strand = fileds['strand']
            start = int(fileds['start'])
            end = int(fileds['end'])
            # key
            key = fileds['id']
            # filter chromoseome
            if chrom in chrmosome_list:
                # extarct sequence
                if strand == '+':
                    seq = genome[chrom][(start-1):end].seq
                elif strand == '-':
                    seq = genome[chrom][(start-1):end].complement.reverse.seq
                else:
                    pass
                # save in dict
                res.setdefault(key,'')
                res[key] += seq
            else:
                pass

        ###########################################################
        # 输出序列
        outputfile = open(outfile,'w')

        # fasta序列分割长度
        my_length = 60

        # 输出
        for key,val in res.items():
            outputfile.write('>' + key + '\n')
            while len(val) > my_length:
                outputfile.write(val[0:my_length] + '\n')
                val = val[my_length:len(val)]
            outputfile.write(val + '\n')

        # 关闭文件
        outputfile.close()
    ####################################################################################
    elif db == 'ucsc':
        # 信息保存在字典里
        info = {}
        # 打开测试 gtf 文件
        with gzip.open(gtffile,'rt') as gtf:
            for line in gtf:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                # 分割
                fields = line.split()
                # 类型
                type = fields[2]
                if type == 'exon':
                    # 名称
                    gene_name = fields[17].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[11].replace('"','').replace(';','')
                    # 连接名称
                    key = '|'.join([gene_name,gene_id,trans_id])
                    # 计算多个外显子长度
                    start = int(fields[3])
                    end = int(fields[4])
                    length = end - start + 1
                    # 累计求和
                    info.setdefault(key,0)
                    info[key] += length
        ######################################
        # 转为数据框
        res = pd.DataFrame(pd.Series(info), columns = ['transcript_length'])

        # 添加基因名列
        res['gene_name'] = [line.split(sep='|')[0] for line in list(res.index[:])]

        # 排序
        res_sorted = res.sort_values(by = ['gene_name','transcript_length'],ascending=False)

        # 筛选最长转录本id
        longest_id = res_sorted.drop_duplicates(subset=['gene_name'],keep='first').index.values.tolist()

        # 筛选最长转录本表格
        longest_data = res_sorted.loc[res_sorted.index.isin(longest_id)]

        longest_data['ID'] = longest_data.index
        # order columns
        longest_data = longest_data[['ID','gene_name','transcript_length']]
        
        # 保存
        longest_data.to_csv(r'longest_transcripts_info.csv', index=True)

        ###########################################################
        # 储存最长转录本id  
        transid = {line.split(sep='|')[2]:line for line in list(longest_data.ID)}

        infolist = []
        with gzip.open(gtffile) as gtf:
            for line in gtf:
                # skip
                if line.startswith('#'):
                    continue
                # split
                fields = line.split()
                # feature type
                type = fields[2]
                if type == 'exon':
                    # pos
                    chr = fields[0]
                    start = fields[3]
                    end = fields[4]
                    strand = fields[6]
                    # name
                    gene_name = fields[17].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[11].replace('"','').replace(';','')
                    if trans_id in transid:
                        infolist.append([chr,start,end,strand,type,gene_name,gene_id,trans_id,transid[trans_id]])
                    else:
                        pass
                else:
                    pass

        # to dataframe
        dfinfo = pd.DataFrame(infolist,columns=['chr','start','end','strand','type','gene_name','gene_id','trans_id','id'])
        dfinfo_1_strand = dfinfo[dfinfo['strand'] == '+']

        # descrese coord by - strand gene
        dfinfo_2_strand = dfinfo[dfinfo['strand'] == '-']
        dfinfo_2_strand = dfinfo_2_strand.sort_values(by = ['trans_id','start','end'],ascending = False)

        # merge
        df_fianl = pd.concat([dfinfo_1_strand,dfinfo_2_strand],axis=0)

        ###########################################################
        # extact  sequnece from genome

        # load genome
        genome = Fasta(genomefile)

        # chrmosome info
        chrmosome_list = genome.keys()

        # save in dict
        res = {}
        for line in range(0,df_fianl.shape[0]):
            
            # chromosome strand
            fileds = df_fianl.iloc[line]
            chrom = fileds['chr']
            strand = fileds['strand']
            start = int(fileds['start'])
            end = int(fileds['end'])
            # key
            key = fileds['id']
            # filter chromoseome
            if chrom in chrmosome_list:
                # extarct sequence
                if strand == '+':
                    seq = genome[chrom][(start-1):end].seq
                elif strand == '-':
                    seq = genome[chrom][(start-1):end].complement.reverse.seq
                else:
                    pass
                # save in dict
                res.setdefault(key,'')
                res[key] += seq
            else:
                pass
        
        ###########################################################

        # 输出序列
        outputfile = open(outfile,'w')

        # fasta序列分割长度
        my_length = 60

        # 输出
        for key,val in res.items():
            outputfile.write('>' + key + '\n')
            while len(val) > my_length:
                outputfile.write(val[0:my_length] + '\n')
                val = val[my_length:len(val)]
            outputfile.write(val + '\n')

        # 关闭文件
        outputfile.close()

    ####################################################################################
    job_stop = time.time()
    print("Your job is done! ")
    print("Running with " + str(round(job_stop - job_start,2)) + " seconds!")

    if __name__=="__main__":
	    main()