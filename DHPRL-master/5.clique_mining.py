##############Protein Complex identification algorithm based on DPPN########################
import sys
from math import sqrt

from optparse import OptionParser
from scipy import *
import scipy as sp

from scipy.linalg import *

from numpy import *
import numpy as np

from numpy.linalg import *
import string
import os




protein_out_file="protein.temp"
cliques_file="cliques"
ppi_pair_file="ppi.pair"
ppi_matrix_file="ppi.matrix"
new_PPI_file="dataset/DIP_new.txt"
new_TIMEPPI_file="f_timeout_temp.txt"

def f_key(a):
    return(a[-1])


def getOptions():
    optparser = OptionParser(usage="%prog [options]\n-h for help")
    optparser.add_option("-p", "--PPI_file", dest="PPI_file", help="Input file of PPI data",default="dataset/DIP.txt")
    
    optparser.add_option("-e", "--Gene_expression_file", dest="EXPRESSION_file", help="Input file of gene expression data",default="dataset/series_matrix.txt")
    
    optparser.add_option("-o", "--output", dest="output", help="Output file for writing the identified complexes",default="result/result_DIP.txt",)
   
    optparser.add_option("-c", "--core_thes_parameter", dest="Core_thes_parameter", help="expand_thes_parameter: range from 0 to 1", default="0.8")# 0.8

    optparser.add_option("-a", "--attach_thes_parameter", dest="Attach_thes_parameter", help="expand_thes_parameter: range from 0 to 1", default="0.2")# 0.2

    
    
    (options, args) = optparser.parse_args()
    if not options.PPI_file:
        optparser.error("No input PPI data file")
    if not options.EXPRESSION_file:
        optparser.error("No input Gene expression data file")
    if not options.output:
        optparser.error("No output file defined")
    return options, args



if __name__ == "__main__":
    
    options, args = getOptions()

    Dic_map={}
    index=0
    All_node=set([])
    All_node_index=set([])
    go_set_list=[]
    expression_list=[]
    All_node_go=set([])
    neighbor_list=[]
    protein_time_list=[]
    time_protein_list = []
    Seed_edge_list=[]
    Cliques_list=[]
    Core_list=[]
    Complex_list=[]
    Core_protein_set=set([])
    
    Core_edge_set=set([])
    Protein_SD_weight={}
    
    Node1=[]
    Node2=[]
    
    Node_PPI_1=[]
    Node_PPI_2=[]
    Weight = []
    neighbor_PPI_list=[]

    
    Time_num=12  # 时段分为12时段
    Times_for_SD=3  # 每个时段有三个活性水平


    
    Seed_weight_thresh=float(options.Core_thes_parameter)  # 0.09 种子边得分阈值
    Attach_thresh=float(options.Attach_thes_parameter)  # 0.05
    
   
    ########### input high-throughput PPI data（高通量PPI数据） ############
    f = open(options.PPI_file,"r")
    f_protein_out=open(protein_out_file,"w")  # protein.temp
    

    
    for line in f:
        line = line.strip().split()
        if len(line)==2:
            
            
            if line[0] not in Dic_map:
                
                Dic_map[line[0]]=index
                f_protein_out.write(line[0]+"\n")
                
                index+=1
                go_set_list.append(set([]))
                expression_list.append([])
                neighbor_list.append(set([]))
                neighbor_PPI_list.append(set([]))
                protein_time_list.append(set([])) 
                
                
            if line[1] not in Dic_map:
                
                Dic_map[line[1]]=index
                f_protein_out.write(line[1]+"\n")
                
                index+=1
                go_set_list.append(set([]))
                expression_list.append([])
                neighbor_list.append(set([]))
                neighbor_PPI_list.append(set([]))
                protein_time_list.append(set([]))

           
                
    Protein_num=index
    f.close()
    f_protein_out.close()  # 蛋白质名称文件

    ###########input Gene expression data（基因表达数据）##################
   
    f=open(options.EXPRESSION_file,"r")
    
    
    for line in f:
        line = line.strip().split()
        if len(line)==38:
            if line [1] in Dic_map:  # 如果基因表达数据中的蛋白质 出现在 数据集蛋白质集中
                
                for i in range(2,Time_num+2):  # 2~13
                    expression_list[Dic_map[line[1]]].append((float(line[i])+float(line[i+12])+float(line[i+24]))/3)
                    # 计算蛋白质节点在12个时间节点上的活性表达值（三个周期的平均值）
           
    f.close()


    ###############compute active time attribute for proteins 计算蛋白质活性时间属性###################（包括活跃时间点及对应的活性概率）
    for t in range(0,Time_num):
        time_protein_list.append(set([]))  # 建立12个时间集合


    for instance in Dic_map:
        
        if len(expression_list[Dic_map[instance]])>=12:
            Temp_mean_value= 0.
            Temp_sd_value=0.
            Temp_thresh_3SD=0.
            Temp_thresh_2SD=0.
            Temp_thresh_1SD=0.
            
            for j in range(0,Time_num):  # 12时刻表达值平均数
                Temp_mean_value+=expression_list[Dic_map[instance]][j]
            Temp_mean_value/=Time_num

            for j in range(0,Time_num):# 12时刻表达值标准差
                Temp_sd_value+=(expression_list[Dic_map[instance]][j]-Temp_mean_value)**2
            Temp_sd_value/=(Time_num-1)

            k = 1
            Temp_thresh_kSD = Temp_mean_value + k * (Temp_sd_value ** 0.5) * (Temp_sd_value / (1 + Temp_sd_value))

            for j in range(0, Time_num):
                if expression_list[Dic_map[instance]][j] >= Temp_thresh_kSD:
                    protein_time_list[Dic_map[instance]].add(j)  # 记录蛋白质节点的活跃时间节点
                    time_protein_list[j].add(Dic_map[instance])  # 记录各时间活跃的蛋白质
    
    ###########input PPI data 去除掉一些在任何时间段都不共同活跃的PPI ############（有共活跃时间的PPI构建的新PPI网络）


    all_cliques_file = open("all_cliques.txt", "w")
    for t in range(0, Time_num):
        f = open(options.PPI_file, "r")
        f_timeout_temp=open("f_timeout_temp.txt","w")
        for line in f:
            line = line.strip().split()
            if len(line) == 2:
                if line[0] != line[1]:
                    if Dic_map[line[0]] in time_protein_list[t] and Dic_map[line[1]] in time_protein_list[t]:
                        f_timeout_temp.write(line[0]+" "+line[1]+"\n")
        f_timeout_temp.close()
        f.close()
        os.system("ConvertPPI.exe " + new_TIMEPPI_file + " " + protein_out_file + " " + ppi_pair_file + " " + ppi_matrix_file)
        os.system("Mining_Cliques.exe " + ppi_matrix_file + " " + "1" + " " + "3" + " " + str(Protein_num) + " " + cliques_file)

        with open(cliques_file, 'r') as f:
            cliques_lines = f.readlines()
            for i in cliques_lines:
                all_cliques_file.write(i)

    all_cliques_file.close()



