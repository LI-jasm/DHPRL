# DHPRL


**This is the data and code for our paper** ` Temporal Protein Complex Identification Based on Dynamic Heterogeneous Protein Information Network Representation Learning`.

## Requirements

- python==3.8.5
- scipy==1.5.4
- torch==1.7.0
- numpy==1.19.2
- scikit_learn==0.24.2


Run
------------
    Step 1: 1.dynamisation.py: Introduce time-course gene expression data of proteins from series_matrix.txt to obtain the temporal sub-network.
    Step 2: 2.create_connections.py: Establish three types of node connections in Dynamic Heterogeneous Protein Information Network (DHPIN). In addition to protein-protein interactions with common active time, protein-molecule functional connections and protein-biological process connections are obtained according to go_slim_mapping.tab.txt.
    Step 3: 3.protein_embedding: On top of the DHPIN network constructed in Step 2, learn protein node representations from the 1-hop relation view and meta-path view, and optimise the protein vector representations with a dual-view co-comparison mechanism.
    Step 4: 4.Update_linking_weight.py: Calculate the cosine similarity based on the protein vector representation and update the weight of each edge in the temporal sub-network.
    Step 5: 5.clique_mining.py: Generate a set of seed cores of proteins that are structurally closely related in each temporal sub-network separately with cluster mining algorithm, and group them together.
    Step 6: 6.Cluster_core_attachment.py: Generate complexes based on the core-attachment structure.
    Step 7: 7.Compare_performance.py: Compare the performance of DHPRL and other classical methods based on different evaluation metrics.

## Usage of protein_embedding
Fisrt, You need to make protein-protein interactions(protein-protein), protein-molecule functional connectivity (protein-F), and protein-biological process connectivity (protein-P) as inputs to generate 1-hop relations and meta-paths with neibor.py and mp_gen.py, respectively. And then you can use the following commend to run our model: 

> python main.py DIP --gpu=0

Here, "DIP" can be replaced with the name of your dataset.



## Contact

If you have any questions, please feel free to contact us!

**Zeqian Li** @github.com/LI-jasm \
**Email:** lizeqian@dlmu.edu.cn \
**Site:** [GitHub](https://github.com/LI-jasm)

