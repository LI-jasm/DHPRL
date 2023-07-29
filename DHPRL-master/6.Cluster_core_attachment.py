from optparse import OptionParser
from numpy import *

from numpy.linalg import *
import string
import os


def f_key(a):
    return(a[-1])


def density_score(temp_set, matrix):
    temp_density_score = 0.
    for m in temp_set:
        for n in temp_set:
            if n != m and matrix[m, n] != 0:
                temp_density_score += matrix[m, n]

    temp_density_score = temp_density_score / (len(temp_set) * (len(temp_set) - 1))
    return temp_density_score


def merge_cliques(new_cliques_set, matrix):
    seed_clique = []

    while (True):
        temp_cliques_set = []
        if len(new_cliques_set) >= 2:
            seed_clique.append(new_cliques_set[0])  # 把候选核集合中的第一个（分数较大的）添加到种子核中

            for i in range(0, len(new_cliques_set)):
                if len(new_cliques_set[i].intersection(new_cliques_set[0])) == 0:  # 如果和候选核集合中第一个的重复节点数为0
                    temp_cliques_set.append(new_cliques_set[i])
                elif len(new_cliques_set[i].difference(new_cliques_set[0])) >= 3:  # 如果和候选核集合中第一个的不同节点数大于等于3
                    temp_cliques_set.append(new_cliques_set[i].difference(new_cliques_set[0]))

            cliques_set = []

            for i in temp_cliques_set:

                clique_score = density_score(i, matrix)
                temp_list = []
                for j in i:
                    temp_list.append(j)

                temp_list.append(clique_score)
                cliques_set.append(temp_list)

            cliques_set.sort(key=f_key, reverse=True)

            new_cliques_set = []
            for i in range(len(cliques_set)):
                temp_set = set([])
                for j in range(len(cliques_set[i]) - 1):
                    temp_set.add(cliques_set[i][j])
                new_cliques_set.append(temp_set)

        elif len(new_cliques_set) == 1:
            seed_clique.append(new_cliques_set[0])
            break
        else:
            break

    return seed_clique


def expand_cluster(seed_clique, all_protein_set, matrix, expand_thres):
    expand_set = []
    complex_set = []

    for instance in seed_clique:
        avg_node_score = density_score(instance, matrix)

        temp_set = set([])
        for j in all_protein_set.difference(instance):
            temp_score = 0.
            for n in instance:
                temp_score += matrix[n, j]
            temp_score /= len(instance)

            if (temp_score ) >= expand_thres:
                temp_set.add(j)
        expand_set.append(temp_set)
    for i in range(len(seed_clique)):
        complex_set.append(seed_clique[i].union(expand_set[i]))

    return (complex_set)


protein_out_file="protein.temp"  #for each protein a index
cliques_file="cliques"
ppi_pair_file="ppi.pair"
ppi_matrix_file= "ppi.matrix"

if __name__ == "__main__":

    #options, args = getOptions()

    ###########input PPI data############
    #f = open(options.PPI_file, "r")
    Time_num = 12

    f = open("dataset/krogan2006core_attr_sim.txt", "r")  # 在新PPI上的节点相似度（实际包含的节点数量少）
    f1 = open("dataset/krogan2006core.txt", "r")
    f_protein_out = open(protein_out_file, "w")
    Dic_map = {}
    index = 0
    Node1 = []
    Node2 = []
    Weight = []
    All_node = set([])
    All_node_index = set([])


    for line in f:
        line = line.strip().split()
        if len(line) == 3:
            Node1.append(line[0])
            Node2.append(line[1])
            Weight.append(float(line[2]))

    f.close()

    #print Dic_map
    #----------------------------------------------------------------------------------

    for line in f1:
        line = line.strip().split()
        if len(line) == 2:
            All_node.add(line[0])  # all_node重新获取
            All_node.add(line[1])
            if line[0] not in Dic_map:  # 以下重新获取
                Dic_map[line[0]] = index
                f_protein_out.write(line[0] + "\n")
                All_node_index.add(index)
                index += 1
            if line[1] not in Dic_map:
                Dic_map[line[1]] = index
                f_protein_out.write(line[1] + "\n")
                All_node_index.add(index)
                index += 1
    Node_count = index  # 节点是PPI中的节点
    f1.close()
    f_protein_out.close()
    #print Dic_map


    ######dic_map to map_dic###########
    Map_dic = {}
    for key in Dic_map.keys():
        Map_dic[Dic_map[key]] = key

    #print Map_dic

    ######bulid Adj_matrix###########

    Adj_Matrix = mat(zeros((Node_count, Node_count), dtype=float))

    if len(Node1) == len(Node2):

        for i in range(len(Node1)):
            if Node1[i] in Dic_map and Node2[i] in Dic_map:
                Adj_Matrix[Dic_map[Node1[i]], Dic_map[Node2[i]]] = Weight[i]
                Adj_Matrix[Dic_map[Node2[i]], Dic_map[Node1[i]]] = Weight[i]
    #print Adj_Matrix.shape[0]

    # 在source_code中生成all_cliques.txt
    cliques_set = []
    # f = open(cliques_file, "r")
    f = open("all_cliques.txt", "r")
    for line in f:
        temp_set = []
        line = line.strip().split()
        for i in range(1, len(line)):
            temp_set.append(int(line[i]))
        cliques_set.append(temp_set)  # 检测到的候选核集合

    f.close()
    avg_clique_score = 0.

    for instance in cliques_set:
        clique_score = density_score(instance, Adj_Matrix)  # 候选核的密集度得分
        avg_clique_score += clique_score
        instance.append(clique_score)
    avg_clique_score /= len(cliques_set)  # 所有候选核的平均得分
    cliques_set.sort(key=f_key, reverse=True)  # 根据核的密集度排序
    #print cliques_set

    new_cliques_set = []
    for i in range(len(cliques_set)):
        temp_set = set([])
        for j in range(len(cliques_set[i]) - 1):
            temp_set.add(cliques_set[i][j])
        new_cliques_set.append(temp_set)  # 根据密度得分排序后的候选核
    #print new_cliques_set
    #print len(new_cliques_set)

    seed_clique = merge_cliques(new_cliques_set, Adj_Matrix)  # 得到重复度小的种子核
    #print seed_clique
    #print len(seed_clique)

    expand_thres = 0.3  # 0.3
    complex_set = expand_cluster(seed_clique, All_node_index, Adj_Matrix, expand_thres)  # 预测到的复合物集合
    print("##########output predicted complexes##########\n")
    final_file = open("result/result_krogan2006core", "w")

    for i in range(len(complex_set)):

        line = ""
        for m in complex_set[i]:
            line += Map_dic[m] + " "
        line += "\n"

        final_file.write(line)
    final_file.close()

    print(len(complex_set))

    print("##########COAN completes############")