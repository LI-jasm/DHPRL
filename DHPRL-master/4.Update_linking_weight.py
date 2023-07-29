import numpy as np
import pickle

# str1="dataset/krogan2006core.txt"
str2="dataset/krogan2006core_attr_vector.txt"  # 记录节点向量表示
str3="dataset/krogan2006core_attr_sim.txt"  # 加权新的PPI网络
str4="dataset/krogan_new.txt"  # 新的PPI网络

path="0.pkl"
data=open(path, 'rb')
file= pickle.load(data)  # 学习到的节点嵌入  # krogan 2708个节点表示√
file1=open(str4)  # 新的PPI网络
file2=open(str2,'w')  # 记录节点向量表示
print("get the vector representation: ")
node=[]
for j in file1:
    temp1=j.split(' ')[0]
    temp2=j.split(' ')[1].rstrip('\n')
    if temp1 not in node:
        node.append(temp1)
    if temp2 not in node:
        node.append(temp2)
file1.close()  # node list（新的PPI网络中涉及到的所有节点）


with open("krogan2006core_node.txt", 'r') as f:
    node = f.readlines()
    for i in range(len(node)):
        s1= node[i].rstrip('\n')
        file2.write(s1)
        file2.write(' ')
        s2 = str(file[i, :]).replace('[', '')
        s2 = s2.replace(']', '')
        s2 = s2.replace('\n', ' ')
        file2.write(s2)
        file2.write('\n')
    file2.close()  # file2 首列是蛋白质名称name；次列及之后为蛋白质节点嵌入表示


#-----------------------------计算节点之间相似度-----------------------------------
print("calculate the similarity between two nodes:")
file=open(str2)  # 向量
file1=open(str4)  # 新的PPI邻接边
file2=open(str3,'w')  # 记录加权PPI
file3=open("krogan2006core_node.txt", "r")

def cos_sim(vector1,vector2):
    dot_product = 0.0
    normA = 0.0
    normB = 0.0
    for a, b in zip(vector1, vector2):
        dot_product += a * b
        normA += a ** 2
        normB += b ** 2
    result=dot_product / ((normA * normB) ** 0.5)
    return result

edge_name_name=[]
for i in file1:
    node_name1=i.split(' ')[0]
    node_name2=i.split(' ')[1]
    node_name2 = node_name2.split('\n')[0]
    d={}
    d['node_name1']=node_name1
    d['node_name2']=node_name2
    edge_name_name.append(d)  # node_name: ,node_name:

vector=[]
for i in file:  # 节点向量文件
    if not i.strip(): continue
    node_name=i.split(' ',1)[0]
    node_vector=i.split(' ',1)[1].rstrip('\n')
    node_vector = node_vector.split()
    # node_vector = ' '.join(node_vector.split())
    node_vector = list(map(float, node_vector))
    d = {}
    d['node_name'] = node_name
    d['node_vector'] = node_vector
    vector.append(d)  # node_name: ,node_vector:

v1=[]
v2=[]
for i in edge_name_name:
    temp1=0
    temp2=0
    for j in vector:
        if(i['node_name1']==j['node_name']):
            v1=np.array(j['node_vector'])
            temp1=1
    for z in vector:
        if(i['node_name2']==z['node_name']):
            v2=np.array(z['node_vector'])
            temp2=1
    if(temp1==1)and(temp2==1):
        result=cos_sim(v1,v2)
        file2.write(i['node_name1'])
        file2.write(' ')
        file2.write(i['node_name2'])
        file2.write(' ')
        file2.write(str(result))
        file2.write('\n')

file.close()
file1.close()
file2.close()