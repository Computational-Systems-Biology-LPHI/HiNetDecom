import os, libsbml, sys, numpy, networkx, tkinter
from tkinter import filedialog, Tk, simpledialog
import matplotlib.pyplot as plt
#=======================================================================================================================
#=======================================================================================================================
#b2AR model
# file_path = 'b2AR_PKA_v5_removed_reaction.xml'
#EGFR model
file_path = 'D4_model_EGFR_v13b_removed_reaction.xml'
#=======================================================================================================================
#=======================================================================================================================
file_name = os.path.basename(file_path).split('.')[0] # Get the file names
reader = libsbml.SBMLReader()
document = reader.readSBML(file_path)
model = document.getModel()
#=======================================================================================================================
#=======================================================================================================================
def logical_test(reaction , species_id: str) -> bool:
        for reactant in reaction.getListOfReactants():
            if reactant.getSpecies() == species_id:
                return True
        for product in reaction.getListOfProducts():
            if product.getSpecies() == species_id:
                return True
        return False
#-------------------
#Species and rates
species = [species.getId() for species in model.getListOfSpecies()]
print(f"\tSpecies: {species}.\n")
print()
rates = [reaction.getKineticLaw().getFormula() if reaction.getKineticLaw() else None for reaction in model.getListOfReactions()]
#-------------------
#Interaction graph
G = networkx.DiGraph()
#-------------------
edges_G = []
for i1 in range(len(species)):
        for i2 in range(len(species)):
            for j1 in range(len(rates)):
                test_1 = logical_test(model.getReaction(j1), species[i1])
                test_2 = rates[j1].find(species[i2]) != -1
                test_3 = model.getSpecies(i1).getConstant() 
                if test_1 and test_2 and not test_3 and (i2+1, i1+1) not in edges_G and i1 != i2:
                    edges_G.append((i2+1,i1+1))
G.add_edges_from(edges_G)
#-------------------
fig = plt.figure()
fig.canvas.manager.set_window_title(f"Interaction graph: {file_name}")
pos = networkx.kamada_kawai_layout(G)
networkx.draw(G, pos, with_labels=True, arrows=True, node_color='lightblue', node_size=1000, edgecolors='black', linewidths=2)
networkx.draw_networkx_edges(G, pos, arrows=True, arrowsize=25, width=2)
plt.show(block=False) 
#=======================================================================================================================
#=======================================================================================================================
nodes_G = G.nodes()
#---
#My undirected graph
Gun = networkx.Graph()
Gun.add_nodes_from(nodes_G)
for u in nodes_G:
    for v in nodes_G:
        #------------
        if networkx.has_path(G, u, v):
            f_1 =  len(networkx.shortest_path(G, u, v)) - 1
        else:
            f_1 = 0
        #------------
        if networkx.has_path(G, v, u):
            f_2 =  len(networkx.shortest_path(G, v, u)) - 1
        else:
            f_2 = 0
        #------------   
        if u != v and f_1 <= 1 and f_2 <= 1 and f_1*f_2 != 0 and not Gun.has_edge(u,v) and not Gun.has_edge(v,u):
            Gun.add_edge(u, v)
#---------------------------------------------    
fig = plt.figure()
fig.canvas.manager.set_window_title(f"Undirected graph: {file_name}")
pos = networkx.kamada_kawai_layout(Gun)
networkx.draw(Gun, pos, with_labels=True, arrows=True, node_color='lightblue', node_size=1000, edgecolors='black', linewidths=2)
networkx.draw_networkx_edges(Gun, pos, arrows=True, arrowsize=25, width=2)
plt.show(block=False) 
#=======================================================================================================================
#=======================================================================================================================
#rsccs in terms of species ids and names
r_sccs_index = [list(component) for component in networkx.connected_components(Gun)]
r_sccs_species_names = [[species[val-1] for val in sublist] for sublist in r_sccs_index]
#---
print(f"\tThe r-strongly connected components:{r_sccs_species_names}.\n")
print()
print(f"\tThe r-strongly connected components (indices):{r_sccs_index}.\n")
print()
#=======================================================================================================================
#=======================================================================================================================
#Quotient graph
Q = networkx.quotient_graph(G, r_sccs_index)
#---
networkx.relabel_nodes(Q, {f : str(set(f)) for f in Q.nodes()}, copy=False)
#-------------------
fig = plt.figure()
fig.canvas.manager.set_window_title(f"Quotient graph: {file_name}")
#pos = networkx.kamada_kawai_layout(Q)
#pos = networkx.spring_layout(Q)
pos = networkx.circular_layout(Q)
networkx.draw(Q, pos, with_labels=True, arrows=True, node_color='lightblue', node_size=3000, edgecolors='black', linewidths=2)
networkx.draw_networkx_edges(Q, pos, arrows=True, arrowsize=25, width=2)
plt.show(block=False)  
#=======================================================================================================================
#=======================================================================================================================
Q1 = Q.copy()
#---
networkx.set_edge_attributes(Q1, -1, 'weight')
done = False
#---
while not done:
    done = True
    cycles = list(networkx.simple_cycles(Q1))
    #---
    for cycle in cycles:
        cycle_set = [(cycle[-1], cycle[0])] + [(cycle[i], cycle[i + 1]) for i in range(len(cycle) - 1)]
        #---
        if sum(Q1[source][target]['weight'] for source, target in cycle_set) < 0: 
            #---
            done = False
            #---
            for source, target in cycle_set:
                Q1[source][target]['weight'] = - Q1[source][target]['weight']
            #---
            for source, target in cycle_set:
                weight = Q1.get_edge_data(source, target)['weight']
                Q1.remove_edge(source, target)
                Q1.add_edge(target, source)
                Q1[target][source]['weight'] = weight
            #---
            break
#---
negative_edges = []
positive_edges = []
#---
for u, v in Q1.edges():
    if Q1.get_edge_data(u,v)['weight'] == -1:
        negative_edges.append((u,v))
    else:
        positive_edges.append((u,v))
#---
dag = networkx.DiGraph()
dag.add_edges_from(negative_edges)
#---
H = networkx.DiGraph()
H.add_edges_from(positive_edges)
#---
for u, v in H.edges():
    H[u][v]['weight'] = 1 
#---    
agonies = {}
hierarchy_ranks = {node : 0 for node in Q1.nodes()}
#---
networkx.set_node_attributes(Q1, hierarchy_ranks, 'rank') # set every node's rank to zero
weights = networkx.get_edge_attributes(Q1, 'weight')
#---
done = False
while not done:
    for u, v in Q1.edges():
        if Q1.nodes[v]['rank'] < Q1.nodes[u]['rank'] - weights[u, v]: 
            done = False
            Q1.nodes[v]['rank'] = Q1.nodes[u]['rank'] - weights[u, v]
            break
        else: done = True
    else:
        done = True
#---
r_sccs_and_scores = {n : Q1.nodes[n]['rank'] for n in Q1.nodes}
#---
for u, v in dag.edges():
    agonies[u,v] = 0
#---    
for u, v in H.edges():
    agonies[u,v] = Q1.nodes[u]['rank'] - Q1.nodes[v]['rank'] + 1
#---
a_Q = sum(max(r_sccs_and_scores[u] - r_sccs_and_scores[v] + 1, 0) for u, v in Q1.edges)
#---
if Q.edges():
    hierarchy = 1 - a_Q / len(Q.edges)
else:
    hierarchy = 1; 
#---
#Rank the scores and the blocks based on their scores
r_sccs_scores = list([value for value in r_sccs_and_scores.values()])
#---
print(f"\tScores of r-strongly connected components:{r_sccs_scores}\n")  
print()
#---
aut_pairs_index = []
aut_pairs_names = []
#--------
aut_pairs_i = []
aut_pairs_ii = []
for i in range(min(r_sccs_scores), max(r_sccs_scores)+1):
    indices = [index for index, value in enumerate(r_sccs_scores) if value == i]
    for item in indices:
        aut_pairs_i.extend(r_sccs_index[item])
        aut_pairs_ii.extend(r_sccs_species_names[item])
    #--------
    aut_pairs_index.append(aut_pairs_i.copy())
    aut_pairs_names.append(aut_pairs_ii.copy())
#--------    
print(f"\tHiearachical levels:{aut_pairs_names}\n")  
print()
#=======================================================================================================================
#=======================================================================================================================
input()
