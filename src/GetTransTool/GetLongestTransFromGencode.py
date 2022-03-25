def main():
    """
    Extract longest transcript from gencode database transcripts fasta file.
    Only for human and mouse transcripts fasta file.
    """

    # 导入模块
    import pandas as pd
    import gzip
    import time
    import argparse

    parser = argparse.ArgumentParser(usage="GetLongestTransFromGencode --file gencode.vM28.transcripts.fa.gz --outfile longest_trans.fa",
                                    description="Get longest transcripts from gencode transcripts fasta file.",
                                    epilog="Thank your for your support, if you have any questions or suggestions please contact me: 3219030654@stu.cpu.edu.cn.")
    parser.add_argument('-v','--version', action='version', version='%(prog)s 0.0.3')

    # 读取转录本文件
    parser.add_argument('-f','--file', type=str,action="store",dest="transfile",metavar="transfile",help='input your transcripts file with ".gz" format. (gencode.vM28.transcripts.fa.gz)')
    # 导出文件名称
    parser.add_argument('-o','--outfile', type=str,action="store",dest="longestfile",metavar="longestfile",help='output your longest transcript file. (longest_trans.fa)')
    # 解析参数
    args = parser.parse_args()

    # 获取参数
    inputfile = args.transfile
    outfile = args.longestfile

    # main fuction
    print("Your job is running, please wait...")
    ######################################################################
    job_start = time.time()
    #######################
    # 保存ID列表
    info = []
    with gzip.open(inputfile,'rt') as tx:
        for line in tx:
            if line.startswith('>'):
                info_name = line.replace('\n','')
                info.append(info_name)

    # 更改名称            
    comb = []
    for line in info:
        # 分割取出需要信息
        gene_name = line.split('|')[5]
        gene_id = line.split('|')[1]
        transcript_id = line.split('|')[0].replace('>','')
        length = line.split('|')[6]
        
        # 连接完整名称
        # finame = gene_name + '_' + gene_id + '_' + transcript_id + '_' + length
        finame = '|'.join([gene_name,gene_id,transcript_id,length])
        # 保存到列表
        comb.append([finame,gene_name,int(length)])

    # 转为数据框
    data_info = pd.DataFrame(comb,columns=['fullname','gene_name','translength'])

    ######################################################################
    # 新建输出改名结果文件
    output_fa = open('name_changed.fa','w')

    # loop change ID
    with gzip.open(inputfile,'rt') as tx:
        for line in tx:
            if line.startswith('>'):
                info_name = line.replace('\n','')
                # 分割取出需要信息
                gene_name = info_name.split('|')[5]
                gene_id = info_name.split('|')[1]
                transcript_id = info_name.split('|')[0].replace('>','')
                length = info_name.split('|')[6]
                # 连接完整名称
                finame = '>' + gene_name + '|' + gene_id + '|' + transcript_id + '|' + length + '\n'
                # 命名
                output_fa.write(finame)
            else:
                seq = line
                output_fa.write(seq)

    # 关闭文件
    output_fa.close()

    ######################################################################
    # 按转录本序列降序排序
    data_infonew = data_info.sort_values(by = ['gene_name','translength'],ascending = False,inplace=False)

    # 筛选最长转录本id
    longest_id = list('>' + data_infonew.drop_duplicates(subset=['gene_name'],keep='first')['fullname'])

    ####################### save csv
    # remove '>'
    longest_idClean = [id.replace('>','') for id in longest_id]
        
    # 筛选最长转录本表格
    longest_data = data_infonew.loc[data_infonew.fullname.isin(longest_idClean)]

    # 保存
    longest_data.to_csv(r'longest_transcripts_info.csv', index=False)

    ######################################################################
    # 保存筛选的id为字典
    filter_id = {id:0 for id in longest_id}

    # 读取 fasta 文件保存为字典
    with open('name_changed.fa') as fa:
        fa_dict = {}
        for line in fa:
            if line.startswith('>'):
                seq_name = line.strip()
                fa_dict[seq_name] = ''
            else:
                # 序列
                fa_dict[seq_name] += line.replace('\n','')
                
    # 新建输出结果文件
    output_fa = open(outfile,'w')

    # fasta序列分割长度
    my_length = 60

    # 输出
    for key,val in fa_dict.items():
        if key in filter_id:
            output_fa.write(key + '\n')
            while len(val) > my_length:
                output_fa.write(val[0:my_length] + '\n')
                val = val[my_length:len(val)]
            output_fa.write(val + '\n')
        
    # 关闭文件
    output_fa.close()

    ####################################################################################
    job_stop = time.time()
    print("Your job is done! ")
    print("Running with " + str(round(job_stop - job_start,2)) + " seconds!")

    if __name__=="__main__":
	    main()