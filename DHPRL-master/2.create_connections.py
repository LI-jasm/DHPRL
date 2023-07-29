import numpy as np
import re
str1="dataset/krogan2006core.txt"

file1=open(str1)
filetemp=open("krogan2006core_node.txt",'w')
print("create topological network!")
node=[]
for j in file1:
    temp1=j.split('	')[0]
    temp2=j.split('	')[1].rstrip('\n')
    if temp1 not in node:
        node.append(temp1)
        filetemp.write(temp1)
        filetemp.write('\n')
    if temp2 not in node:
        node.append(temp2)
        filetemp.write(temp2)
        filetemp.write('\n')
file1.close()
filetemp.close()
print(node)

# -----------建立蛋白质节点索引，提取protein-protein关系-----------
f = open("krogan2006core_node.txt", "r")
str4="dataset/krogan_new.txt"
nfile=open(str4)

Dic_map = {}
index = 0
All_node_index = set([])

for line in f:
    line = line.strip()
    if line not in Dic_map:
        Dic_map[line] = index
        All_node_index.add(index)
        index += 1
Node_count = index  # 至此为蛋白质节点建立索引Dic_map
print("no_protein:", Node_count)

file5=open("krogan2006core_protein-protein_relation.txt",'w')
for line in nfile:
    temp1 = line.split(' ')[0]
    temp2 = line.split(' ')[1].rstrip('\n')
    file5.write(str(Dic_map[temp1]))
    file5.write(' ')
    file5.write(str(Dic_map[temp2]))
    file5.write('\n')  # 至此protein-protein关系建立完成
file5.close()
f.close()
#-----------建立P节点索引，提取protein-P关系-----------
f6=open("krogan2006core_protein-P_temp.txt",'w')
f = open("krogan2006core_node.txt", "r")
Dic_map_P = {}
noP=0
for i in f:
    i=i.rstrip('\n')
    print(i)
    f6.write(i)
    f6.write(' ')
    file = open("dataset/go_slim_mapping.tab.txt")
    for j in file:
        node_name = j.split('	')[0]
        go_tag=j.split('	')[3]
        go = j.split('	')[5]
        if node_name==i:
            if go_tag=="P":
                if go!='' :
                    go = go.strip('GO:')
                    f6.write(go)
                    f6.write(' ')
                    if go not in Dic_map_P:
                        Dic_map_P[go] = noP
                        noP += 1
    file.close()
    f6.write('\n')
f.close()
f6.close()
print("noP:",noP)

f6=open("krogan2006core_protein-P_temp.txt")
file6=open("krogan2006core_protein-P_relation.txt",'w')
gon = []
for line in f6:
    n_name = line.split(' ', 1)[0]
    n_go = line.split(' ', 1)[1].rstrip('\n').rstrip(' ')
    n_go = n_go.split(' ')
    one = {}
    one['node_name'] = n_name
    one['node_go'] = n_go
    gon.append(one)
for i in node:
    for j in gon:
        if i==j['node_name']:
            if j['node_go']!=['']:
                for z in j['node_go']:
                    file6.write(str(Dic_map[i]))
                    file6.write(' ')
                    file6.write(str(Dic_map_P[z]))
                    file6.write('\n')
file6.close()
f6.close()
#-----------建立F节点索引，提取protein-F关系-----------
f7=open("krogan2006core_protein-F_temp.txt",'w')
f = open("krogan2006core_node.txt", "r")
Dic_map_F = {}
noF=0
for i in f:
    i=i.rstrip('\n')
    print(i)
    f7.write(i)
    f7.write(' ')
    file = open("dataset/go_slim_mapping.tab.txt")
    for j in file:
        node_name = j.split('	')[0]
        go_tag=j.split('	')[3]
        go = j.split('	')[5]
        if node_name==i:
            if go_tag=="F":
                if go!='' :
                    go = go.strip('GO:')
                    f7.write(go)
                    f7.write(' ')
                    if go not in Dic_map_F:
                        Dic_map_F[go] = noF
                        noF += 1
    file.close()
    f7.write('\n')
f.close()
f7.close()
print("noF:",noF)

f7=open("krogan2006core_protein-F_temp.txt")
file7=open("krogan2006core_protein-F_relation.txt",'w')
gom = []
for line in f7:
    node_name = line.split(' ', 1)[0]
    node_go = line.split(' ', 1)[1].rstrip('\n').rstrip(' ')
    node_go = node_go.split(' ')
    one = {}
    one['node_name'] = node_name
    one['node_go'] = node_go
    gom.append(one)
for i in node:
    for j in gom:
        if i==j['node_name']:
            if j['node_go'] != ['']:
                for z in j['node_go']:
                    file7.write(str(Dic_map[i]))
                    file7.write(' ')
                    file7.write(str(Dic_map_F[z]))
                    file7.write('\n')
file7.close()
f7.close()