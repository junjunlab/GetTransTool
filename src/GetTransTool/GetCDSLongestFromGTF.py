def main():
    """
    Extract longest CDS regeion with longest transcript from gtf format annotation file based on ensembl/ucsc database.
    """

    # 导入模块
    import pandas as pd
    import pyfastx
    import gzip
    import time
    import warnings
    import argparse

    warnings.filterwarnings('ignore')

    parser = argparse.ArgumentParser(usage="GetCDSLongestFromGTF --database ensembl --gtffile Homo_sapiens.GRCh38.101.gtf.gz --genome Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --outfile longest_cds_trans.fa",
                                    description="Extract longest CDS regeion with longest transcript from gtf format annotation file based on ensembl/ucsc database.",
                                    epilog="Thank your for your support, if you have any questions or suggestions please contact me: 3219030654@stu.cpu.edu.cn.")
    # version
    parser.add_argument('-v','--version', action='version', version='%(prog)s 0.0.3')
    # 读取注释类型文件
    parser.add_argument('-d','--database',type=str,action="store",dest="database",metavar="databse",choices=['ucsc','ensembl'],default="ensembl",
                        help='which annotation database you choose. (default="ensembl", ucsc/ensembl)')
    # 读取gtf文件
    parser.add_argument('-g','--gtffile', type=str,action="store",dest="gtffile",metavar="gtffile",
                        help='input your GTF file with ".gz" format.')
    # 读取基因组fasta文件
    parser.add_argument('-fa','--genome',type=str,action="store",dest="genome",metavar="genome",
                        help='your genome fasta file matched with your GTF file with ".gz" format. (Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)')
    # 导出文件名称
    parser.add_argument('-o','--outfile', type=str,action="store",dest="cdslongestfile",metavar="cdslongestfile",
                        help='output your longest transcript file. (longest_cds_trans.fa)')
    # 解析参数
    args = parser.parse_args()

    # 获取参数
    db = args.database
    gtffile = args.gtffile
    genomefile = args.genome
    outfile = args.cdslongestfile

    # main fuction
    print("Your job is running, please wait...")
    ######################################################################
    job_start = time.time()
    ################################################################################################
    if db == 'ensembl':
        # 打开 gtf 文件
        with gzip.open(gtffile,'rt') as gtf:
            # 信息保存在字典里
            trans_len = {}
            utr5_len = {}
            cds_len = {}
            for line in gtf:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                # 分割
                fields = line.split()
                # 类型
                type = fields[2]
                if len(fields) > 24:
                    biotype = fields[23].replace('"','').replace(';','')
                    utr_biotype = fields[21].replace('"','').replace(';','')
                    if biotype == 'protein_coding' and type == 'exon':
                        # 名称
                        gene_name = fields[19].replace('"','').replace(';','')
                        gene_id = fields[9].replace('"','').replace(';','')
                        trans_id = fields[13].replace('"','').replace(';','')
                        # 连接名称
                        key = '_'.join([gene_name,gene_id,trans_id])
                        # 计算多个外显子长度
                        length = int(fields[4]) - int(fields[3]) + 1
                        # 累计求和
                        trans_len.setdefault(key,0)
                        trans_len[key] += length
                    elif biotype == 'protein_coding' and type == 'CDS':
                        # 名称
                        gene_name = fields[19].replace('"','').replace(';','')
                        gene_id = fields[9].replace('"','').replace(';','')
                        trans_id = fields[13].replace('"','').replace(';','')
                        # 连接名称
                        key = '_'.join([gene_name,gene_id,trans_id])
                        # 计算多个CDS长度
                        length = int(fields[4]) - int(fields[3]) + 1
                        # 累计求和
                        cds_len.setdefault(key,0)
                        cds_len[key] += length
                    elif utr_biotype == 'protein_coding':
                        # 名称
                        gene_name = fields[17].replace('"','').replace(';','')
                        gene_id = fields[9].replace('"','').replace(';','')
                        trans_id = fields[13].replace('"','').replace(';','')
                        # 连接名称
                        key = '_'.join([gene_name,gene_id,trans_id])
                        if type == 'five_prime_utr':
                            # 计算多个5'UTR长度
                            length = int(fields[4]) - int(fields[3]) + 1
                            # 累计求和
                            utr5_len.setdefault(key,0)
                            utr5_len[key] += length
                        else:
                            # 若无则为 0
                            utr5_len.setdefault(key,0) 
                else:
                    pass
        
        # transorm into dataframe and merge by id
        df_tran = pd.DataFrame.from_dict(trans_len,orient='index',columns=['translength'])
        df_tran['ID'] = df_tran.index

        df_cds = pd.DataFrame.from_dict(cds_len,orient='index',columns=['cdslength'])
        df_cds['ID'] = df_cds.index

        df_5utr = pd.DataFrame.from_dict(utr5_len,orient='index',columns=['utr5length'])
        df_5utr['ID'] = df_5utr.index

        # 按id合并表格
        data_info = df_cds.merge(df_tran.merge(df_5utr,on='ID'),on='ID')

        # 添加基因名列
        data_info['gene_name'] = [i.split(sep='_')[0] for i in data_info['ID']]

        # 按gen_name cdslength translength 降序排序
        data_infonew = data_info.sort_values(by = ['gene_name','cdslength','translength'],ascending = False,inplace=False)

        # 保存
        data_infonew.to_csv(r'All_transcripts_cds_info.csv', index=False)

        ############################
        # 筛选最长转录本id
        longest_id = list(data_infonew.drop_duplicates(subset=['gene_name'],keep='first')['ID'])

        # 筛选最长转录本表格
        longest_data = data_infonew.loc[data_infonew.ID.isin(longest_id)]

        # 保存
        longest_data.to_csv(r'longest_cds_transcripts_info.csv', index=False)

        # 给 ID 添加 CDS 位置信息和转录本长度信息
        longest_data['ID'] = longest_data.ID + '_' + (longest_data.utr5length + 1).map(str) + '_' + \
                                (longest_data.utr5length + longest_data.cdslength).map(str) + '_' + \
                                (longest_data.translength).map(str)

        # 储存最长转录本id  
        trans_id = {line.split(sep='_')[1]:line.split(sep='_')[2] for line in list(longest_data.ID)}

        gid_dict = {key: 0 for key,val in trans_id.items()}
        transid_dict = {val: 0 for key,val in trans_id.items()}

        ###########################################
        # 新建导出最长转录本 gtf 文件
        out_gtf = open('CDS_longest_trans.gtf','w')

        # loop export
        with gzip.open(gtffile,'rt') as gtf:
            for line in gtf:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                # 分割
                fields = line.split()
                # 类型
                id_type = fields[2].replace('"','').replace(';','')
                gid = fields[9].replace('"','').replace(';','')
                tranid = fields[13].replace('"','').replace(';','')
                if id_type == 'gene':
                    if gid in gid_dict:
                        out_gtf.write(line)
                    else:
                        pass
                elif tranid in transid_dict:
                    out_gtf.write(line)
                else:
                    pass
                
        # 关闭文件
        out_gtf.close()
        ###########################################
        # prepare filter id  
        filter_id = {line.split(sep='_')[2]:line for line in list(longest_data.ID)}

        # 读取基因组
        genome = pyfastx.Fasta(genomefile)
        # chrmosome info
        chrmosome_list = [name for name, seq in pyfastx.Fasta(genomefile, build_index=False)]
        
        # 保存结果字典
        res = {}

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
                    trans_id = fields[13].replace('"','').replace(';','')
                    if trans_id in filter_id:
                        # 名称
                        fa_id = filter_id[trans_id]
                        # chromosome strand
                        chrom = str(fields[0])
                        strand = str(fields[6])
                        # 储存外显子区间
                        interval = (int(fields[3]), int(fields[4]))
                        if chrom in chrmosome_list:
                            # 提取序列
                            seq = genome.fetch(chrom, interval, strand=strand)
                            # 储存
                            res.setdefault(fa_id,'')
                            res[fa_id] += seq
                        else:
                            pass
                    else:
                        pass
                else:
                    pass

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

    ################################################################################################
    elif db == 'ucsc':
        # 打开 gtf 文件
        with gzip.open(gtffile,'rt') as gtf:
            # 信息保存在字典里
            trans_len = {}
            utr5_len = {}
            cds_len = {}
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
                    length = int(fields[4]) - int(fields[3]) + 1
                    # 累计求和
                    trans_len.setdefault(key,0)
                    trans_len[key] += length
                elif type == 'CDS':
                    # 名称
                    gene_name = fields[17].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[11].replace('"','').replace(';','')
                    # 连接名称
                    key = '|'.join([gene_name,gene_id,trans_id])
                    # 计算多个CDS长度
                    length = int(fields[4]) - int(fields[3]) + 1
                    # 累计求和
                    cds_len.setdefault(key,0)
                    cds_len[key] += length
                elif type == '5UTR':
                    # 名称
                    gene_name = fields[17].replace('"','').replace(';','')
                    gene_id = fields[9].replace('"','').replace(';','')
                    trans_id = fields[11].replace('"','').replace(';','')
                    # 连接名称
                    key = '|'.join([gene_name,gene_id,trans_id])
                    # 计算多个5'UTR长度
                    length = int(fields[4]) - int(fields[3]) + 1
                    # 累计求和
                    utr5_len.setdefault(key,0)
                    utr5_len[key] += length
                else:
                    pass
            else:
                pass
        
        # fillwith no 5UTR genes
        new_utr5_len = {key:utr5_len.get(key,0) for key,val in cds_len.items()}

        # transorm into dataframe and merge by id
        df_tran = pd.DataFrame.from_dict(trans_len,orient='index',columns=['translength'])
        df_tran['ID'] = df_tran.index

        df_cds = pd.DataFrame.from_dict(cds_len,orient='index',columns=['cdslength'])
        df_cds['ID'] = df_cds.index

        df_5utr = pd.DataFrame.from_dict(new_utr5_len,orient='index',columns=['utr5length'])
        df_5utr['ID'] = df_5utr.index

        # 按id合并表格
        data_info = df_cds.merge(df_tran.merge(df_5utr,on='ID'),on='ID')

        # 添加基因名列
        data_info['gene_name'] = [i.split(sep='|')[0] for i in data_info['ID']]

        # 按gen_name cdslength translength 降序排序
        data_infonew = data_info.sort_values(by = ['gene_name','cdslength','translength'],ascending = False,inplace=False)

        # 保存
        data_infonew.to_csv(r'All_transcripts_cds_info.csv', index=False)

        ############################
        # 筛选最长转录本id
        longest_id = list(data_infonew.drop_duplicates(subset=['gene_name'],keep='first')['ID'])

        # 筛选最长转录本表格
        longest_data = data_infonew.loc[data_infonew.ID.isin(longest_id)]

        # 保存
        longest_data.to_csv(r'longest_cds_transcripts_info.csv', index=False)

        # 给 ID 添加 CDS 位置信息和转录本长度信息
        longest_data['ID'] = longest_data.ID + '|' + (longest_data.utr5length + 1).map(str) + '|' + \
                                (longest_data.utr5length + longest_data.cdslength).map(str) + '|' + \
                                (longest_data.translength).map(str)

        # 储存最长转录本id  
        trans_id = {line.split(sep='|')[1]:line.split(sep='|')[2] for line in list(longest_data.ID)}

        transid_dict = {val: 0 for key,val in trans_id.items()}

        # 新建导出最长转录本 gtf 文件
        out_gtf = open('CDS_longest_trans.gtf','w')

        # loop export
        with gzip.open(gtffile,'rt') as gtf:
            for line in gtf:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                # 分割
                fields = line.split()
                # 类型
                tranid = fields[11].replace('"','').replace(';','')
                if tranid in transid_dict:
                    out_gtf.write(line)
                else:
                    pass
                
        # 关闭文件
        out_gtf.close()
        ###########################################
        # prepare filter id  
        filter_id = {line.split(sep='|')[2]:line for line in list(longest_data.ID)}

        # 读取基因组
        genome = pyfastx.Fasta(genomefile)

        # chrmosome info
        chrmosome_list = [name for name, seq in pyfastx.Fasta(genomefile, build_index=False)]

        # 保存结果字典
        res = {}

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
                    trans_id = fields[11].replace('"','').replace(';','')
                    if trans_id in filter_id:
                        fa_id = filter_id[trans_id]
                        # chromosome strand
                        chrom = str(fields[0])
                        strand = str(fields[6])
                        # 储存外显子区间
                        interval = (int(fields[3]), int(fields[4]))
                        if chrom in chrmosome_list:
                            # 提取序列
                            seq = genome.fetch(chrom, interval, strand=strand)
                            # 储存
                            res.setdefault(fa_id,'')
                            res[fa_id] += seq
                        else:
                            pass
                    else:
                        pass
                else:
                    pass

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