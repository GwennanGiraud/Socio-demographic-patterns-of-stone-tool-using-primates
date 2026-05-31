################Python script for homophily correction#######################

#To run into R

#Based on paper of Karimi & Oliveira 2023

#First, run R_Script_to_prepare_for_Python !!

#Make sure your environnement is setup for Python and not for R, 
#but it should change by itself when you run the code 
#(try to run it twice if there is an error the first time).

import os # to handle file paths
import networkx as nx # networkx library for network analysis
import pandas as pd # pandas library for data manipulation
import numpy as np # numpy library for numerical operations
from assortativity import* #function manuscript Karimi & Oliveira 2023 (src folder)
import network_generator.er_mixing # function to generate a network with assortativity
from scipy.optimize import minimize # to perform optimization
import scipy.optimize as opt # optimization functions
import matplotlib.pyplot as plt # matplotlib for plotting

###Import data
m=pd.read_csv("mat.csv", index_col=0) # read the matrix from a CSV file #m=mGroom here
m=np.array(m)# Vector transformation
df=pd.read_csv("df.csv") # read the dataframe from a CSV file
char=df['Tool.user'].values #Extract individual attributes
char = (char == 'Yes').astype(int) # Convert 'Yes'/'No' to 1/0 for binary attributes
G = nx.from_numpy_array(m) #Transform matrix in networkx object
attributes_dict = {node_id: char[node_id] for node_id in G.nodes()} # create a dictionary mapping node IDs to their attributes


###Compute adjusted assortativity (correction 1)---------------------
nx.set_node_attributes(G, attributes_dict, "tool") # set the node attributes in the graph
#Calculate the adjusted assortativity coefficient with groups of unequal sizes
nx_adjusted_assortativity(G, "tool") 
      #=> close to 0; so no homophily



###Compute asymmetric mixing (correction 2)-----------------------
def get_f0(G, attribute):
    counts, unique = nx_group_fraction(G, attribute)
    f_0 = min(counts / counts.sum())
    return f_0
f_0 = get_f0(G, 'tool') #Here 'tool' is the binary variable

def estimate_h_er_numerical(nx_graph, attribute="color"):
    import networkx as nx
    counts, unique = nx_group_fraction(nx_graph, attribute)
    f_0 = min(counts / counts.sum())
    group_map = dict([(unique[i], i) for i in range(len(unique))])

    absolute_mixing_matrix = nx.attribute_mixing_matrix(nx_graph, attribute, normalized=False, mapping=group_map)
    absolute_mixing_matrix /= 2
    E = absolute_mixing_matrix.sum()
    E_00 = absolute_mixing_matrix[0, 0]
    E_11 = absolute_mixing_matrix[1, 1]
    e_00 = E_00 / E
    e_11 = E_11 / E
    f_1 = 1 - f_0
    f = lambda h_11: np.abs(
        h_11 - (e_11 / e_00) * (f_0 ** 2 / f_1 ** 2) * (f_0 * f_1 + f_1 ** 2 * h_11 + f_0 * f_1 * (1 - h_11)) / (
                    f_0 ** 2 / e_00 - f_0 ** 2 + f_0 * f_1))
    optimization = minimize(f, 0.5, method="L-BFGS-B", bounds=[(0, 1)])
    h_11 = optimization['x'][0]
    h_00 = (e_00 / e_11) * (f_1 ** 2 / f_0 ** 2) * h_11
    return h_00, h_11  #h=probability of same-group nodes being connected

def estimate_h_er_analytical(nx_graph, attribute="color"):
    import networkx as nx
    counts, unique = nx_group_fraction(nx_graph, attribute)
    f_0 = min(counts / counts.sum())
    group_map = dict([(unique[i], i) for i in range(len(unique))])

    absolute_mixing_matrix = nx.attribute_mixing_matrix(nx_graph, attribute, normalized=False, mapping=group_map)
    absolute_mixing_matrix /= 2
    E = absolute_mixing_matrix.sum()
    E_00 = absolute_mixing_matrix[0, 0]
    E_11 = absolute_mixing_matrix[1, 1]
    e_00 = E_00 / E
    e_11 = E_11 / E
    f_1 = 1 - f_0

    sum_p_ij = 2 * f_0 * f_1 / (1 - e_00 * (1 - f_1 / f_0) - e_11 * (1 - f_0 / f_1))
    h_11 = e_11 * sum_p_ij / f_1 ** 2
    h_00 = (e_00 / e_11) * (f_1 ** 2 / f_0 ** 2) * h_11
    return h_00, h_11
  
#estimate_h_er_numerical(G, "tool")
    
#h_00, h_11 = estimate_h_ba(G, 'tool')
h_00, h_11 = estimate_h_er_numerical(G, 'tool')
print(h_00, h_11)

#Plot the distribution of the mixing matrix
plt.hist(m)
plt.show()

#Based on links distribution we selected power low estimation
h_00  #h=probability of same-group nodes being connected
h_11

###Assessing asymmetric mixing patterns in networks (correction 2)
#Analytical assortativity calculation: asymmetric mixing patterns
analytical_assortativity(f_0, h_00, h_11, model="er")
        #=> close to 0; so no asymmetric mixing patterns



###Mixing matrix------------------------------------------
# Create the mapping so that the indices are correct
counts, unique = nx_group_fraction(G, 'tool')
group_map = {unique[i]: i for i in range(len(unique))}

# Calculate the raw mixing matrix (absolute, unnormalized)
mixing_matrix_absolute = nx.attribute_mixing_matrix(G, 'tool', normalized=False, mapping=group_map)
print(mixing_matrix_absolute) # Real connexion in the network (raw factions)

# Calculate the normalized mixing matrix
mixing_matrix_normalized = nx.attribute_mixing_matrix(G, 'tool', normalized=True, mapping=group_map)
print(mixing_matrix_normalized) # Proportion of connexion in the network (to visualize patterns)



###Adjusted mixing matrix---------------------------------
#Fractions of nodes in each group
f0 = counts.min() / counts.sum()
f1 = 1 - f0
#Randomly expected fractions of edges
E = np.sum(mixing_matrix_absolute)  # Total number of edges in the graph

# Adjusted matrix (Normalized)
adjusted_matrix_normalized = np.array([[f0**2 * h_00, f0*f1*(1 - h_00)], # NTU -> (NTU, TU)
                            [f0*f1*(1 - h_11), f1**2 * h_11]]) # TU -> (NTU, TU)
print(adjusted_matrix_normalized) # Adjusted mixing matrix (Proportion of connexion in the network corrected for asymmetric mixing)

# Adjusted matrix (Brut)
adjusted_matrix_absolute = adjusted_matrix_normalized * E
print(adjusted_matrix_absolute) # Adjusted mixing matrix in absolute fractions (Real connexion in the network corrected for asymmetric mixing)
