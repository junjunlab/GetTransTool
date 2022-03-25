def main():
    """
    Extract longest CDS regeion with longest transcript from gencode database transcripts fasta file.
    Only for human and mouse transcripts fasta file.
    """

    # 导入模块
    import pandas as pd
    import gzip
    import time
    import argparse

    parser = argparse.ArgumentParser(usage="GetCDSLongestFromGencode --file gencode.vM28.pc_transcripts.fa.gz --outfile longest_cds_trans.fa",
                                    description="Extract longest CDS regeion with longest transcript from gencode database transcripts fasta file.",
                                    epilog="Thank your for your support, if you have any questions or suggestions please contact me: 3219030654@stu.cpu.edu.cn.")
    parser.add_argument('-v','--version', action='version', version='%(prog)s 0.0.3')

    # 读取转录本文件
    parser.add_argument('-f','--file', type=str,action="store",dest="transfile",metavar="transfile",help='input your protein-coding transcripts file with ".gz" format. (gencode.vM28.pc_transcripts.fa.gz)')
    # 导出文件名称
    parser.add_argument('-o','--outfile', type=str,action="store",dest="longestfile",metavar="longestfile",help='output your longest transcript file. (longest_cds_trans.fa)')
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
    # 储存id
    cdsinfo = {}

    # 储存改名文件
    tmpfile = open('name_changed.fa','w')

    # main code
    with gzip.open(inputfile,'rt') as pc:
        for line in pc:
            if line.startswith('>'):
                # split id
                fileds = line.split(sep='|')
                gene_name = '>' + fileds[5]
                gene_id = fileds[1]
                trans_id = fileds[0].replace('>','')
                trans_length = fileds[6]
                # include 5UTR+CDS+3UTR or 5UTR+CDS/
                if fileds[8].startswith('CDS'):
                    cds_range = fileds[8].split(sep=':')[1].split(sep='-')
                    cds_start = cds_range[0]
                    cds_end = cds_range[1]
                    # cds length
                    cds_len = int(cds_end) - int(cds_start) + 1
                    fullname = '|'.join([gene_name,gene_id,trans_id,cds_start,cds_end,trans_length])
                    # save
                    tmpfile.write(fullname + '\n')
                    # save in dict
                    cdsinfo[fullname] = str(cds_len)
                # include CDS+3UTR or only CDS
                else:
                    cds_range = fileds[7].split(sep=':')[1].split(sep='-')
                    cds_start = cds_range[0]
                    cds_end = cds_range[1]
                    # cds length
                    cds_len = int(cds_end) - int(cds_start) + 1
                    fullname = '|'.join([gene_name,gene_id,trans_id,cds_start,cds_end,trans_length])
                    # save
                    tmpfile.write(fullname + '\n')
                    # save in dict
                    cdsinfo[fullname] = str(cds_len)
            else:
                # write seq
                tmpfile.write(line)

    # close file
    tmpfile.close()

    #################################################
    # transform into datafarme
    tmp = [[key,key.split(sep='|')[0],int(key.split(sep='|')[5]),int(val)] for key,val in cdsinfo.items()]

    # 转为数据框
    data_info = pd.DataFrame(tmp,columns=['fullname','gene_name','translength','cdslength'])

    # 按gen_name cdslength translength 降序排序
    data_infonew = data_info.sort_values(by = ['gene_name','cdslength','translength'],ascending = False,inplace=False)

    # 保存
    data_infonew.to_csv(r'All_transcripts_cds_info.csv', index=False)

    # 筛选最长转录本id
    longest_id = list(data_infonew.drop_duplicates(subset=['gene_name'],keep='first')['fullname'])

    # 筛选最长转录本表格
    longest_data = data_infonew.loc[data_infonew.fullname.isin(longest_id)]

    # 保存
    longest_data.to_csv(r'longest_cds_transcripts_info.csv', index=False)

    #################################################
    # prepare filter id  
    filter_id = {id:0 for id in list(longest_data.fullname)}

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