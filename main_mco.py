#
# The main filament detection routine
#
# Reads input/input.txt to determine parameters and input files to use.
# Instead of running python Main_mco.py, you could optionally run python mst_ui.py 
# to provide parameters and input files through a user interface.
#
import random
import sys
import math
import os.path
import shutil
from datetime import datetime
import networkx as nx


######################################################################################################################################################
# FILAMENT & GROUP
######################################################################################################################################################


def dist(a, b): 
        return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)

## The tree
def spanning_tree_group(filename):
        try:
                """Creates a spanning tree of group centers"""
                G=nx.Graph()
                file = open(filename) # File with group center
                # Create nodes in G where each node is a galaxy w/ all associated information
                i = 0

                lines = file.readlines()
                while (lines[0].split()[0].isdigit()==False):
                        del lines[0]
                while (lines[(len(lines)-1)] in ['\n', '\r\n']):
                        del lines[(len(lines)-1)]
                        
                for line in lines:
                        words = line.split()
                        grpID = words[0]
                        RA = words[1]
                        Dec = words[2]
                        zRS = words[3]
                        G.add_node(i, grpID = int(grpID), RA = float(RA), Dec = float(Dec), zRS = float(zRS),
                                   x = (C*float(zRS)/100*math.cos(D*float(Dec))*math.cos(D*float(RA))),
                                   y = (C*float(zRS)/100*math.cos(D*float(Dec))*math.sin(D*float(RA))),
                                   z = (C*float(zRS)/100*math.sin(D*float(Dec))))               
                        i += 1
                        
                file.close
                
                for node1 in G: ## add weight to use MST
                        for node2 in range(node1+1, len(G)):
                                distance = dist((G.node[node1]['x'],G.node[node1]['y'],G.node[node1]['z']),(G.node[node2]['x'],G.node[node2]['y'],G.node[node2]['z']))
                                G.add_edge(node1, node2, weight = distance)
                                
                print ("\tSPANNING_TREE_GROUP: (filament_MST_uncut) N b MST:", G.number_of_nodes())
                print ("\tSPANNING_TREE_GROUP: (filament_MST_uncut) E b MST:", G.number_of_edges())
                
                return nx.minimum_spanning_tree(G) ##MST using Kruskal's algorithm
        
        except ValueError:
                print ("\tSPANNING_TREE_GROUP: (file): wrong format!"+
                       "\nINSTRUCTION:\n"+
                       "1. NOT RECCOMMEND table header, NOT RECCOMMEND EMPTY LINE AT THE END -> just the data\n"+
                       "2. SHOULD ONLY contain number\n"+
                       "3. NOTE the columns include: grpID, RA, Dec, zRS (NOTE: It's z, not C*z)\n"+
                       "4. FOLLOW the following format: (4 columns, seperated by tab)\n"+
                       "\t123\t123\t123\t123\n\t...\t...\t...\t...\n\t...\t...\t...\t...\n\t...\t...\t...\t...")              
                sys.exit()
                
        except IOError:
                print (str(filename) +": not found")
                sys.exit()


## decide which filament is deserved to be considered, other group center are taken out
def process_tree(tree, cutoffLength, maxCut):
        '''
        ## tree: the MST which contains all nodes and edges
        ## cutoffLength: edges longer than this lenght are eliminated
        ## maxCut: filaments with less nodes than this are taken out
        '''
        print ("\tPROCESS_TREE: (filament_MST_uncut) N b cut: "), tree.number_of_nodes()
        print ("\tPROCESS_TREE: (filament_MST_uncut) E b trim: "), tree.number_of_edges()
        
        for edge in tree.edges(): # cut long edges
                if tree[edge[0]][edge[1]]['weight'] > cutoffLength:
                        tree.remove_edge(edge[0], edge[1])

        print ("\tPROCESS_TREE: (filament_MST_uncut) E a trim: "), tree.number_of_edges() 
                        
        treeFilamentUntouchedAfterCut = tree.copy() # original copy of the spanning tree after cut off                      
        treeTemp = tree.copy() # this will be torn apart
        filamentsNodesAndEdges = tree.copy() # this will be used to save infos of nodes and edges
        grpCentLonely = nx.Graph() # store groups taken out 
        filamentList = [] # only needed if we develop walk() and make_branch for each tree
        
        i=0
        while treeTemp.nodes(): # when some nodes are still untouched, scanning goes on
                node = random.choice(treeTemp.nodes())
                nbunch = nx.node_connected_component(treeTemp, node)
                temp = treeTemp.subgraph(nbunch)
                if len(temp.nodes()) <= maxCut:
                        filamentsNodesAndEdges.remove_edges_from(temp.edges())
                        for node in nbunch:
                                grpCentLonely.add_node(i, grpID = int(treeTemp.node[node]['grpID']),
                                                       RA = float(treeTemp.node[node]['RA']),
                                                       Dec = float(treeTemp.node[node]['Dec']),
                                                       zRS = float(treeTemp.node[node]['zRS']),
                                                       x = float(treeTemp.node[node]['x']),
                                                       y = float(treeTemp.node[node]['y']),
                                                       z = float(treeTemp.node[node]['z']))
                                i+=1
                        filamentsNodesAndEdges.remove_nodes_from(nbunch)
                        treeTemp.remove_nodes_from(nbunch)
                else:
                        filamentList.append(temp)
                        treeTemp.remove_nodes_from(nbunch)
                        
        print ("\tPROCESS_TREE: (filament) E a cut:"), filamentsNodesAndEdges.number_of_edges()
        print ("\tPROCESS_TREE: (filament) N a cut:"), filamentsNodesAndEdges.number_of_nodes()
        print ("\tPROCESS_TREE: (group_cent_non_filament) N a cut:"), grpCentLonely.number_of_nodes()                      
                        
        return (filamentsNodesAndEdges, grpCentLonely, filamentList, treeFilamentUntouchedAfterCut) 


## find average weight of tree                                
def avg_tree_weight(tree):
        sum=0
        for edge in tree.edges():
                sum+=(tree[edge[0]][edge[1]]['weight'])
        return sum / tree.number_of_edges()


## add galaxies info from text file of all individual galaxies
def add_galaxy(filename):
        try:
                G=nx.Graph()
                file = open(filename) # File with galaxy data

                lines = file.readlines()
                while (lines[0].split()[0].isdigit()==False):
                        del lines[0]
                while (lines[(len(lines)-1)] in ['\n', '\r\n']):
                        del lines[(len(lines)-1)]
                        
                for line in lines:
                        words = line.split()
                        ID = words[0]
                        RA = words[1]
                        Dec = words[2]
                        zRS = words[3]

                        # Assign node name as galaxy ID makes searching in tree faster!
                        G.add_node(int(ID), ID = int(ID), RA = float(RA), Dec = float(Dec), zRS = float(zRS),
                                   x = (C*float(zRS)/100*math.cos(D*float(Dec))*math.cos(D*float(RA))),
                                   y = (C*float(zRS)/100*math.cos(D*float(Dec))*math.sin(D*float(RA))),
                                   z = (C*float(zRS)/100*math.sin(D*float(Dec))))
                file.close
                return G ## galaxies list
        
        except ValueError:
                print ("\tADD_GALAXY: (file): wrong format!"+
                       "\nINSTRUCTION:\n"+
                       "1. NOT RECCOMMEND table header, NOT RECCOMMEND EMPTY LINE AT THE END -> just the data\n"+
                       "2. SHOULD ONLY contain number\n"+
                       "3. NOTE the columns include: ID, RA, Dec, zRS \n"+
                       "4. FOLLOW the following format: (4 columns, seperated by tab)\n"+
                       "\t123\t123\t123\t123\n\t...\t...\t...\t...\n\t...\t...\t...\t...\n\t...\t...\t...\t...")              
                sys.exit()
                
        except IOError:
                print (str(filename) +": not found")
                sys.exit()


## take all the galaxies in group out of the list
def process_galaxy(galL, grpL):
        
        '''
        ## galL: list of galaxy
        ## grpL: list of galaxies in groups
        '''

        print ("\tPROCESS_GALAXY: (galaxy) N b process:"), (galL.number_of_nodes())
        print ("\tPROCESS_GALAXY: (group_full) N:"), (grpL.number_of_nodes())
        
        list = [node for node in galL if node in grpL]
        galL.remove_nodes_from(list)
        
        print ("\tPROCESS_GALAXY: (galaxy) N a process:"), (galL.number_of_nodes())


## add galaxies in group from text file
def add_group(filename):
        try:
                G=nx.Graph()
                file = open(filename) # File with data for galaxies in groups
                
                lines = file.readlines()
                while (lines[0].split()[0].isdigit()==False):
                        del lines[0]
                while (lines[(len(lines)-1)] in ['\n', '\r\n']):
                        del lines[(len(lines)-1)]
                        
                for line in lines:
                        words = line.split()
                        grpID = words[0]
                        ID = words[1]
                        RA = words[2]
                        Dec = words[3]
                        zRS = words[4]
                        G.add_node(int(ID), grpID = int(grpID), ID = int(ID), RA = float(RA), Dec = float(Dec), zRS = float(zRS),
                                   x = (C*float(zRS)/100*math.cos(D*float(Dec))*math.cos(D*float(RA))),
                                   y = (C*float(zRS)/100*math.cos(D*float(Dec))*math.sin(D*float(RA))),                                   
                                   z = (C*float(zRS)/100*math.sin(D*float(Dec))))
                file.close
                return G
        
        except ValueError:
                print ("\tADD_GROUP: (file): wrong format!"+
                       "\nINSTRUCTION:\n"+
                       "1. NOT RECCOMMEND table header, NOT RECCOMMEND EMPTY LINE AT THE END -> just the data\n"+
                       "2. SHOULD ONLY contain number\n"+
                       "3. NOTE the columns include: grpID, ID, RA, Dec, zRS (NOTE: It's z, not C*z)\n"+
                       "4. FOLLOW the following format: (5 columns, seperated by tab)\n"+
                       "\t123\t123\t123\t123\t123\n\t...\t...\t...\t...\t...\n\t...\t...\t...\t...\t...\n\t...\t...\t...\t...\t...")              
                sys.exit()

        except IOError:
                print (str(filename) +": not found")
                sys.exit()


## return list of galaxies in group but not in filament
def process_group(filament, grp):

        '''
        ## filament: filament (which contains all nodes which are grp centers)
        ## grp: list of all galaxies in group
        ## grpMemLonely: list of all galaxies in group but not in filament
        '''
        print ("\tPROCESS_GROUP: (group_full) N: "), (grp.number_of_nodes())
        
        grpMemLonely = grp.copy()
        grpMemFilament = nx.Graph()
        for node in filament:
                for node1 in grpMemLonely.nodes():
                        if grpMemLonely.node[node1]['grpID']==filament.node[node]['grpID']:
                                grpMemFilament.add_node(int(grpMemLonely.node[node1]['ID']),
                                                        grpID = int(grpMemLonely.node[node1]['grpID']),
                                                        ID = grpMemLonely.node[node1]['ID'],
                                                        RA = float(grpMemLonely.node[node1]['RA']),
                                                        Dec = float(grpMemLonely.node[node1]['Dec']),
                                                        zRS = float(grpMemLonely.node[node1]['zRS']),
                                                        x = float(grpMemLonely.node[node1]['x']),
                                                        y = float(grpMemLonely.node[node1]['y']),
                                                        z = float(grpMemLonely.node[node1]['z']))
                                grpMemLonely.remove_node(node1)
                                
        print ("\tPROCESS_GROUP: (group_mem_non_filament) N:"), (grpMemLonely.number_of_nodes())
        print ("\tPROCESS_GROUP: (group_mem_filament) N:"), (grpMemFilament.number_of_nodes())
   
        return (grpMemLonely,grpMemFilament)


######################################################################################################################################################
# SCOOPER
######################################################################################################################################################


## Find distance from edge for scooper
def shortest_dist_from_edge(a,b,c): #Point: c; Edge: a,b (tuples 2D)
        ba = ((b['x']-a['x']), (b['y']-a['y']), (b['z']-a['z']))
        ab = ((a['x']-b['x']), (a['y']-b['y']), (a['z']-b['z']))
        ca = ((c['x']-a['x']), (c['y']-a['y']), (c['z']-a['z']))
        cb = ((c['x']-b['x']), (c['y']-b['y']), (c['z']-b['z']))

        ## Check sweeping angle to make sure the foot of the altitude from the point to the edge does not lie outside of the edge!
        checkB = math.acos((ab[0]*cb[0]+ab[1]*cb[1]+ab[2]*cb[2])/math.sqrt((ab[0]**2+ab[1]**2+ab[2]**2)*(cb[0]**2+cb[1]**2+cb[2]**2)))
        checkA = math.acos((ba[0]*ca[0]+ba[1]*ca[1]+ba[2]*ca[2])/math.sqrt((ba[0]**2+ba[1]**2+ba[2]**2)*(ca[0]**2+ca[1]**2+ca[2]**2)))
        if(checkB>(D*90) or checkA>(D*90)):
              return -1        
        ##sometimes due to floating point problem, checkA and checkB could differ
        return math.sqrt(cb[0]**2+cb[1]**2+cb[2]**2)*(math.sin(checkB)) ##cosine formula for dot product vectors

## Scooper for lonely groups
def scooper_group(grp, G2, G3, r):
        
        '''
        ## grp: Group Centers for Lonely group
        ## G2: list of galaxies
        ## G3: list of galaxies near filaments (gnf)
        '''

        print ("\tSCOOPER_GROUP: (galaxy) N b scoopR:"), G2.number_of_nodes()
        print ("\tSCOOPER_GROUP: (gnf) N b scoopR:"), G3.number_of_nodes()
        
        for node in grp:
                for node1 in G2.nodes():
                        if 0<=dist((grp.node[node]['x'],grp.node[node]['y'],grp.node[node]['z']),
                                (G2.node[node1]['x'],G2.node[node1]['y'],G2.node[node1]['z']))<= r:
                                G3.add_node(node1, ID = G2.node[node1]['ID'], RA = G2.node[node1]['RA'], Dec = G2.node[node1]['Dec'],
                                            zRS = G2.node[node1]['zRS'], x = G2.node[node1]['x'], y = G2.node[node1]['y'], z = G2.node[node1]['z'])
                                G2.remove_node(node1)
                                
        print ("\tSCOOPER_GROUP: (galaxy) N a scoopR:"), G2.number_of_nodes()
        print ("\tSCOOPER_GROUP: (gnf) N a scoopR:"), G3.number_of_nodes()
                                
## Scooper
def scooper(G1, G2 ,r ,d):
        
        '''
        ## G1: filaments (nodes and egdes)
        ## G2: list of galaxies
        ## r: scooper radius around node
        ## d: scooper radius around edge
        
        ## G3: gnf
        '''

        print ("\tSCOOPER: (galaxy) N b scoopR&D:"), G2.number_of_nodes()
        
        G3 = nx.Graph()
        
        # A note on this if both r and d = 0 -> error or r = 0 and d < 0 -> error.
        # Because the same point yields a vector with magnitude zero so the formula becomes acos(infinity) -> domain error
        # This has been overcome by running
              
        for node in G1: # Scan around nodes
                for node1 in G2.nodes():                
                        if 0 <= dist((G1.node[node]['x'],G1.node[node]['y'],G1.node[node]['z']),
                                (G2.node[node1]['x'],G2.node[node1]['y'],G2.node[node1]['z']))<= r:
                                G3.add_node(node1, ID = G2.node[node1]['ID'], RA = G2.node[node1]['RA'], Dec = G2.node[node1]['Dec'],
                                            zRS = G2.node[node1]['zRS'], x = G2.node[node1]['x'], y = G2.node[node1]['y'], z = G2.node[node1]['z'])
                                G2.remove_node(node1)

        print ("\tSCOOPER: (galaxy) N a scoopR:"), G2.number_of_nodes()
        print ("\tSCOOPER: (gnf) N a scoopR:"), G3.number_of_nodes()
        
        for edge in G1.edges(): # Scan around edges
                for node1 in G2.nodes():
                        if 0 <= shortest_dist_from_edge(G1.node[edge[0]], G1.node[edge[1]],G2.node[node1]) <= d:
                                G3.add_node(node1, ID = G2.node[node1]['ID'], RA = G2.node[node1]['RA'], Dec = G2.node[node1]['Dec'],
                                            zRS = G2.node[node1]['zRS'], x = G2.node[node1]['x'], y = G2.node[node1]['y'], z = G2.node[node1]['z'])
                                G2.remove_node(node1)

        print ("\tSCOOPER: (galaxy) N a scoopD:"), G2.number_of_nodes()                                
        print ("\tSCOOPER: (gnf) N a scoopD:"), G3.number_of_nodes()
        
        return G3 ## gnf

# Scooper from Main_Length_Specs.py
def Scoop(G1, G2):
        
        '''
        ## G1: filaments (nodes and egdes)
        ## G2: list of galaxies
        
        ## G3: galaxies with min distance from tree
        '''
        
        G3 = G2.copy()
        nx.set_node_attributes(G3, 'd', -1)
        
        ## A note on this if bot h r and d = 0 -> error or r = 0 and d < 0 -> error.
        ## Because the same point yield vector with magnitude zero so the formula becomes acos(infinity) -> domain error
        ## This has been overcome by running

        for node in G3:
                minD = 10000;
                for node1 in G1:
                        if 0<=dist((G1.node[node1]['x'],G1.node[node1]['y'],G1.node[node1]['z']),
                                (G3.node[node]['x'],G3.node[node]['y'],G3.node[node]['z']))<= minD:
                                minD = dist((G1.node[node1]['x'],G1.node[node1]['y'],G1.node[node1]['z']),
                                (G3.node[node]['x'],G3.node[node]['y'],G3.node[node]['z']))
                for edge in G1.edges():
                        if 0<=shortest_dist_from_edge(G1.node[edge[0]], G1.node[edge[1]],G3.node[node])<=minD:
                                minD = shortest_dist_from_edge(G1.node[edge[0]], G1.node[edge[1]],G3.node[node])
                G3.node[node]['d'] = minD
                
        return G3 ## gnf

######################################################################################################################################################
# TENDRILS & VOID
######################################################################################################################################################

        
## run MST (Kruskal) on the remaing isolated galaxies to find tendrils
def spanning_tree_tendril(filename, maxLength):

        '''
        ## maxLength: edges longer than this are cut
        '''

        ## this precut is to lower the load for system.
        ## for < 20000 galaxy, it should be favorable but if it jumps to 30000 to 60000
        ## this takes 2-6 hours to finish! (Kruskal BigO(N^2logN)

        G=nx.Graph()
        isolatedGal = open(filename)
        
        i=0
        for line in isolatedGal:
                words = line.split()
                ID = words[0]
                RA = words[1]
                Dec = words[2]
                zRS = words[3]
                G.add_node(i, ID = int(ID), RA = float(RA), Dec = float(Dec), zRS = float(zRS),
                           x = (C*float(zRS)/100*math.cos(D*float(Dec))*math.cos(D*float(RA))),
                           y = (C*float(zRS)/100*math.cos(D*float(Dec))*math.sin(D*float(RA))),
                           z = (C*float(zRS)/100*math.sin(D*float(Dec))))
                i+=1
                
        isolatedGal.close                       

        for node1 in G: ## add weight to use MST
                for node2 in range(node1+1, len(G)):
                        distance = dist((G.node[node1]['x'],G.node[node1]['y'],G.node[node1]['z']),(G.node[node2]['x'],G.node[node2]['y'],G.node[node2]['z']))
                        if distance <=maxLength:
                                G.add_edge(node1, node2, weight = distance)
        return nx.minimum_spanning_tree(G)


## I just copy the code from process_tree() because these two are similar
def process_tendril(tree, cutoffLength, maxCut):
        
        '''
        ## tree: the MST which contains all nodes and edges
        ## cutoffLength: edges longer than this lenght are eliminated
        ## maxCut: filaments with less nodes than this are taken out
        '''
        
        print ("\tPROCESS_TENDRIL: (tendril) N b cut:"), tree.number_of_nodes()
        print ("\tPROCESS_TENDRIL: (tendril) E b trim: "), tree.number_of_edges()
        
        for edge in tree.edges(): 
                if tree[edge[0]][edge[1]]['weight'] > cutoffLength:
                        tree.remove_edge(edge[0], edge[1])
                        
        print ("\tPROCESS_TENDRIL: (tendril) E a trim: "), tree.number_of_edges() 
        
        treeFilamentUntouched = tree.copy()                        
        treeTemp = tree.copy()
        filamentsNodesAndEdges = tree.copy()
        grpCentLonely = nx.Graph()
        
        i=1        
        while treeTemp.nodes():
                node = random.choice(treeTemp.nodes())
                nbunch = nx.node_connected_component(treeTemp, node)
                temp = treeTemp.subgraph(nbunch)
                if len(temp.nodes()) <= maxCut:
                        filamentsNodesAndEdges.remove_edges_from(temp.edges())
                        for node in nbunch:
                                grpCentLonely.add_node(i, ID = int(treeTemp.node[node]['ID']),
                                                       RA = float(treeTemp.node[node]['RA']),
                                                       Dec = float(treeTemp.node[node]['Dec']),
                                                       zRS = float(treeTemp.node[node]['zRS']),
                                                       x = float(treeTemp.node[node]['x']),
                                                       y = float(treeTemp.node[node]['y']),
                                                       z = float(treeTemp.node[node]['z']))
                                i+=1
                        filamentsNodesAndEdges.remove_nodes_from(nbunch)
                        treeTemp.remove_nodes_from(nbunch)
                else:
                        treeTemp.remove_nodes_from(nbunch)

        print ("\tPROCESS_TENDRIL: (tendril) E a cut:"), filamentsNodesAndEdges.number_of_edges()
        print ("\tPROCESS_TENDRIL: (tendril) N a cut:"), filamentsNodesAndEdges.number_of_nodes()
        print ("\tPROCESS_TENDRIL: (void) N a cut:"), grpCentLonely.number_of_nodes()
                        
        return (filamentsNodesAndEdges, grpCentLonely, treeFilamentUntouched ) # [0]=tendrils w edges, [1]=void, [2]=remaining galaxies untouched


## add galaxies from text file for tendril and void
def add_galaxy_tendril_and_void(filename):
        
        G=nx.Graph()
        file = open(filename) # File with galaxy data
        # Create nodes in G where each node is a galaxy w/ all associated information
        for line in file:
                words = line.split()
                ID = words[0]
                RA = words[1]
                Dec = words[2]
                zRS = words[3]
                G.add_node(int(ID), ID = int(ID), RA = float(RA), Dec = float(Dec), zRS = float(zRS),
                           x = (C*float(zRS)/100*math.cos(D*float(Dec))*math.cos(D*float(RA))),
                           y = (C*float(zRS)/100*math.cos(D*float(Dec))*math.sin(D*float(RA))),
                           z = (C*float(zRS)/100*math.sin(D*float(Dec))))
        file.close
        return G


######################################################################################################################################################
# OUTPUT
######################################################################################################################################################


## FINAL EXPORT TO LIST WITH ENVIRONMENT ATTRIBUTE
def output_final_list(galaxy_full, inputList):
        
        H = galaxy_full.copy()
       
        nx.set_node_attributes(H, 'ENVIRONMENT', 'n')
#        print "no of nodes: ", H.number_of_nodes()
#        print "nodes(H) ", nodes(H)
#        print "info(H, n=306197) ", info(H,n=306197)
#        print "H.nodes ", H.nodes() #gives just the list of nodes 
#        print "H.nodes data=true ", H.nodes(data=TRUE) #gives all the data in each node
        for G in inputList:
            #print "G= ", G 
            #print "G0 =", G[0]
            #print "G" + ".nodes ", G[0].nodes()
            for node in G[0]:
#                       print "G[0]= ", G[0]
#                        print "G[1]= ", G[1]
#                        print "node=", node
                        H.node[node]['ENVIRONMENT'] = G[1]
                #str(G[0].node[node]['ID']
                # if input node name instead of string into that G.node[] box error is returned!
        return H

## FINAL EXPORT TO LIST WITH ENVIRONMENT ATTRIBUTE SORTED
def output_final_list_sorted(inputList):
        H = nx.Graph()
        
        for G in inputList:
                nx.set_node_attributes(G[0], 'ENVIRONMENT', G[1])
                H.add_nodes_from(G[0].nodes(data=True))
        return H


## export a graph to text file
def output_text_final(graph, filename, folder_name):
        
        dir_path = os.path.join('output', folder_name)  # will return 'output/folder_name'
        if not os.path.isdir(dir_path):
                os.makedirs(dir_path) # create directory [current_path]/output/folder_name
        name_of_file = filename
        file = open(os.path.join(dir_path, name_of_file + ".txt"), "w")
        
        i=0 # make counter to avoid end of file empty line
        countMax = graph.number_of_nodes()
        file.write("# ID  ")
        file.write("RA  ")
        file.write("Dec  ")
        file.write("zRS  ")
        file.write("ENVIRONMENT  \n")        
        for node in graph:
                file.write(str(graph.node[node]['ID']) + "  ")
                file.write(str(graph.node[node]['RA']) + "  ")
                file.write(str(graph.node[node]['Dec']) + "  ")
                file.write(str(graph.node[node]['zRS']) + "  ")
                file.write(str(graph.node[node]['ENVIRONMENT']) + "  ")
                
                i+=1
                if i!= countMax: #check end of file
                        file.write("\n")
        file.close()

# output_text_final from Main_Length_Specs.py
def output_text_main_length_specs(graph, filename): 
    file = open(os.path.join(filename + ".txt"), "w")

    i=0 # make counter to avoid end of file empty line
    countMax = graph.number_of_nodes()
    file.write("# ID  ")
    file.write("RA  ")
    file.write("Dec  ")
    file.write("zRS  ")
    file.write("Distance \n")        
    for node in graph:
            file.write(str(graph.node[node]['ID']) + "  ")
            file.write(str(graph.node[node]['RA']) + "  ")
            file.write(str(graph.node[node]['Dec']) + "  ")
            file.write(str(graph.node[node]['zRS']) + "  ") 
            file.write(str(graph.node[node]['d']) + "  ")
            
            i+=1
            if i!= countMax: #check end of file
                file.write("\n")
    file.close()

## export a graph to text file
def output_text(graph, filename, folder_name):
        
        dir_path = os.path.join('output', folder_name, 'galaxy')  # will return 'output/folder_name'
        if not os.path.isdir(dir_path):
                os.makedirs(dir_path) # create directory [current_path]/output/folder_name
        name_of_file = filename
        file = open(os.path.join(dir_path, name_of_file + ".txt"), "w")
        
        i=0 # make counter to avoid end of file empty line
        countMax = graph.number_of_nodes()
        file.write("# ID  ")
        file.write("RA  ")
        file.write("Dec  ")
        file.write("zRS  \n")
        
        for node in graph:
                file.write(str(graph.node[node]['ID']) + "  ")
                file.write(str(graph.node[node]['RA']) + "  ")
                file.write(str(graph.node[node]['Dec']) + "  ")
                file.write(str(graph.node[node]['zRS']) + "  ")
                
                i+=1
                if i!= countMax: #check end of file
                        file.write("\n")
        file.close()


## create output file to feed data into topcat
def output_text_group_cent(graph, filename, folder_name):
        
        dir_path = os.path.join('output', folder_name, 'group_cent')
        if not os.path.isdir(dir_path):
                os.makedirs(dir_path)
        name_of_file = filename
        file = open(os.path.join(dir_path, name_of_file + ".txt"), "w")
        
        i=0 # make counter to avoid end of file empty line
        countMax = graph.number_of_nodes()
        file.write("# grpID  ")
        file.write("RA  ")
        file.write("Dec  ")
        file.write("zRS  \n")
        
        for node in graph:
                file.write(str(graph.node[node]['grpID']) + "  ")
                file.write(str(graph.node[node]['RA']) + "  ")
                file.write(str(graph.node[node]['Dec']) + "  ")
                file.write(str(graph.node[node]['zRS']) + "  ")
                
                i+=1
                if i!= countMax: #check end of file
                        file.write("\n")
        file.close()


def copy_input(folderName):
        if not os.path.isdir(os.path.join('output', folderName, 'input')):
                os.makedirs(os.path.join('output', folderName, 'input'))
        for files in os.listdir('input'):
                if files.endswith(".txt"):
                        shutil.copy(os.path.join('input', files), os.path.join('output', folderName, 'input'))
                        

######################################################################################################################################################
# LOG
######################################################################################################################################################


## method to copy output to text file
class Log(object):
        def __init__(self, *files):
                self.files = files
        def write(self, obj):
                for f in self.files:
                        f.write(obj)
                        f.flush() # If you want the output to be visible immediately
        def flush(self) :
                for f in self.files:
                        f.flush()


def log_start(folder_name):
        dir_path = os.path.join('output', folder_name)
        if not os.path.isdir(dir_path):
                os.makedirs(dir_path)
        file = open(os.path.join(dir_path, "log.txt"), "w")
        sys.stdout = Log(sys.stdout, file)
        return file


def log_end(logFile):
        logFile.close()


def logSec(section_name):
        print ("\n==============================\n" + section_name +"\n==============================\n")

      
######################################################################################################################################################
# TEMP
######################################################################################################################################################


def output_text_temp(graph, filename):
        if not os.path.isdir('temp'):
                os.makedirs('temp')
        file = open(os.path.join('temp', filename + ".txt"), "w")
        
        i=0
        countMax = graph.number_of_nodes()
        
        for node in graph:
                file.write(str(graph.node[node]['ID']) + "  ")
                file.write(str(graph.node[node]['RA']) + "  ")
                file.write(str(graph.node[node]['Dec']) + "  ")
                file.write(str(graph.node[node]['zRS']) + "  ")
                i+=1
                if i!= countMax:
                        file.write("\n")
        file.close()

def make_path_temp(filename):
        if not os.path.isdir('temp'):
                os.makedirs('temp')
        return os.path.join('temp', filename)

def write_temp(graph, filename):
        if not os.path.isdir('temp'):
                os.makedirs('temp')
        nx.write_gexf(graph,os.path.join('temp', filename + ".gexf"))

def read_temp(filename):
        if not os.path.isdir('temp'):
                os.makedirs('temp')
        return nx.read_gexf(os.path.join('temp', filename + ".gexf"))

def del_temp():
        global temp_choice
        if not os.path.isdir('temp'):
                os.makedirs('temp')
        if os.listdir('temp'):
                choice = raw_input("Remove temp folder? y = YES, other = NO\n")
                if choice == "y":
                        shutil.rmtree('temp')
                        temp_choice = "YES"
                        return
                temp_choice ="NO"
                return
        temp_choice ="EMPTY"


######################################################################################################################################################
#TIME
######################################################################################################################################################

def time_on():
        global startTime
        startTime = datetime.now()
        print (str(startTime) + "\n")
        global timeA
        timeA = datetime.now()
         
def time(process_name):
        global timeB
        global timeA
        timeB = datetime.now()
        print (str(process_name) + " done in: " + str(timeB - timeA))
        timeA = datetime.now()

def time_off():
        global startTime
        print ("TOTAL DURATION: " + str(datetime.now() - startTime))
        print ("\n" + str(datetime.now()) + "\n")


######################################################################################################################################################
#INPUT
######################################################################################################################################################
   

def input_data():
        if not os.path.isdir('input'):
                os.makedirs('input')
                print ("INPUT: (folder): not found")
                return
        try:
                global fileInput
                fileInput = os.path.join('input', "input" + ".txt")
                f = open(fileInput, "r")
        except IOError:
                print ("INPUT: (input.txt): not found")
                sys.exit()
                
        list = []
        for line in f:
                list.append(line.split()[0])
                if len(list) > 11:
                        f.close
                        break
                f.close
                
        if len(list) < 11:
                print("INPUT: (file): wrong format!"+
                      "\nNeed EXACTLY 11 parameters, each parameter typed as the first token on one line.\n"+
                      "Do not start with blank so all parameters lie on the first 11 lines. See example of input.txt below:\n"+
                      "\n==============================\nBEGGINING OF FILE\n==============================\n"+
                      "test_1\nsdss_atlas_90552.txt\nsdss_atlas_90552_grp.txt\nsdss_atlas_90552_cent.txt\n299792.458\n1.29415100275\n"+
                      "3\n3.1\n3.1\n1.1\n2\n## anything on and below this line will not affect the input file:\n"+
                      "PARAMETER by line:\n1. output folder name: anything (parsed to str)\n2. galaxy list: anything.txt (parsed to str)\n"+
                      "3. galaxy in group list: anything.txt (parsed to str)\n4. group center list: anything.txt (parsed to str)\n"+
                      "5. speed of light (km/s): float\n6. trim length - in process_tree() (Mpc/h): float\n"+
                      "7. minimum number of member required to be considered filament: int\n"+
                      "8. radius of scooper from node of filament (Mpc/h): float\n"+
                      "9. distance of scooper from edge of filament (Mpc/h): float\n10. trim length - in process_tendril(): float\n"+
                      "11. minimum number of member required to be considered tendril: int)"+
                      "\n==============================\nEND OF FILE\n==============================\n")             
                sys.exit()

        try:
                global folderName
                folderName = str(list[0])
                global fileAll
                fileAll = os.path.join('input', str(list[1]))
                global fileGrp
                fileGrp = os.path.join('input', str(list[2]))
                global fileGrpCent
                fileGrpCent = os.path.join('input', str(list[3]))
                global C
                C = float(list[4])
                global trimTree
                trimTree = float(list[5])
                global cutTree
                cutTree = int(list[6])
                global scoopR
                scoopR = float(list[7])
                global scoopD
                scoopD = float(list[8])
                global trimTendril
                trimTendril = float(list[9])
                global cutTendril
                cutTendril = int(list[10])
                
        except ValueError:
                print("INPUT: (file): wrong format!"+
                      "\nNeed EXACTLY 11 parameters, each parameter typed as the first token on one line.\n"+
                      "Do not start with blank so all parameters lie on the first 11 lines. See example of input.txt below:\n"+
                      "\n==============================\nBEGGINING OF FILE\n==============================\n"+
                      "test_1\nsdss_atlas_90552.txt\nsdss_atlas_90552_grp.txt\nsdss_atlas_90552_cent.txt\n299792.458\n1.29415100275\n"+
                      "3\n3.1\n3.1\n1.1\n2\n## anything on and below this line will not affect the input file:\n"+
                      "PARAMETER by line:\n1. output folder name: anything (parsed to str)\n2. galaxy list: anything.txt (parsed to str)\n"+
                      "3. galaxy in group list: anything.txt (parsed to str)\n4. group center list: anything.txt (parsed to str)\n"+
                      "5. speed of light (km/s): float\n6. trim length - in process_tree() (Mpc/h): float\n"+
                      "7. minimum number of member required to be considered filament: int\n"+
                      "8. radius of scooper from node of filament (Mpc/h): float\n"+
                      "9. distance of scooper from edge of filament (Mpc/h): float\n10. trim length - in process_tendril(): float\n"+
                      "11. minimum number of member required to be considered tendril: int)"+
                      "\n==============================\nEND OF FILE\n==============================\n")                   
                sys.exit()

## check if there is another output folder with the same name in order to prevent data loss
def input_check_duplicate(folderName):
        if not os.path.isdir('output'):
                os.makedirs('output')
                return
        for files in os.listdir('output'):
                if str(files) == str(folderName):
                        choice = raw_input("Folder " + str(folderName) + " existed. Overwrite? y = YES, other = NO\n")
                        if choice == "y":
                                shutil.rmtree(os.path.join('output', folderName))
                                return
                        else:
                                print ("Please change the folder output name in your input.txt")
                                sys.exit()
        return

######################################################################################################################################################
#VARIABLES & CONSTANTS
######################################################################################################################################################

## Variables

global startTime
global timeA
global timeB
global temp_choice

global fileInput # input/input.txt

global folderName # folder name to store in output
global fileAll # text file contains all galaxy info (note that z is z in this one, not C*z)
global fileGrp # text file contains all galaxy in group info (note that z is C*z in this one, not z)
global fileGrpCent # text file contains all group center info (note that z is C*z in this one, not z)
global C # speed of light (km/s) - use for d = C*z/100 (Mpc/h)
global trimTree # trim length (Mpc/h) - in process_tree()
global cutTree # minimum number of member required to be considered filament - in process_tree()
global scoopR # radius of scooper from node of filament (Mpc/h) - in scooper()
global scoopD # distance of scooper from edge of filament (Mpc/h) - in scooper()
global trimTendril # trim length (Mpc/h) - in process_tendril()
global cutTendril # minimum number of member required to be considered tendril - in process_tendril()


## constants 
D = math.acos(-1.0) / 180 # acos(-1) = Pi


######################################################################################################################################################
#MAIN
######################################################################################################################################################

input_data()
input_check_duplicate(folderName)
del_temp()
logFile = log_start(folderName)
time_on()
print ("DEL_TEMP: " + temp_choice)

## PARAMETER
logSec("PARAMETER")

print ("SOURCE: " + str(fileGrp)) #put the input file here
print ("OUTPUT: output/" + str(folderName)) 

print ("CONSTANT: ")
print ("\tSPEED_OF_LIGHT: (km/s) " + str(C))

print ("IMPORT/ADD: ")
print ("\tGALAXY_LIST: " + str(fileAll))
print ("\tGROUP_LIST: " + str(fileGrp))
print ("\tGROUP_CENT_LIST: " + str(fileGrpCent))

print ("SPANNING_TREE_GROUP:")
print ("\tTRIM_LENGTH: (Mpc/h) " + str(trimTree))
print ("\tMINIMUM_NUMBER_OF_MEMBER: " + str(cutTree))

print ("SCOOPER:")
print ("\tSCOOP_RADIUS: (Mpc/h) " + str(scoopR))
print ("\tSCOOP_DISTANCE: (Mpc/h) " + str(scoopD))

print ("SPANNING_TREE_TENDRIL: ")
print ("\tTRIM_LENGTH: (Mpc/h) " + str(trimTendril))
print ("\tMINIMUM_NUMBER_OF_MEMBER: " + str(cutTendril))



## Import/add
logSec("IMPORT/ADD")
try: # Read MST group cent file
        filament_MST_uncut = read_temp("filament_MST_uncut")
        print ("IMPORT: (filament_MST_uncut) gexf: loaded")
        print ("\tIMPORT: (filament_MST_uncut) N:", filament_MST_uncut.number_of_nodes())
        print ("\tIMPORT: (filament_MST_uncut) E:", filament_MST_uncut.number_of_edges())
except IOError:
        print ("ADD: (filament_MST_uncut) gexf: not found")
        filament_MST_uncut = spanning_tree_group(fileGrpCent) # Takes the spanning tree of the group centers
        print ("finished with spanning_tree_group")
        write_temp(filament_MST_uncut, "filament_MST_uncut")        
        print ("ADD: (filament_MST_uncut) gexf: created")
        print ("\tADD: (filament_MST_uncut) N a MST:", filament_MST_uncut.number_of_nodes())
        print ("\tADD: (filament_MST_uncut) E a MST:", filament_MST_uncut.number_of_edges())
time("IMPORT/ADD: (filament_MST_uncut)")

try: # Read Galaxy list file
        galaxy = read_temp("galaxy")
        galaxy_full = read_temp("galaxy_full")
        print ("IMPORT: (galaxy) gexf: loaded")
except IOError:
        print ("ADD: (galaxy) gexf: not found.")
        galaxy = add_galaxy(fileAll)
        write_temp(galaxy, "galaxy_full")
        galaxy_full = read_temp("galaxy_full")
        print ("ADD: (galaxy_full) gexf: created")
time("IMPORT/ADD: (galaxy)")

try: # Read group file and process it (take out grouped galaxies from galaxy list)
        group_full = read_temp("group_full")
        print ("IMPORT: (group_full) gexf: loaded")
except IOError:
        print ("ADD: (group_full) gexf: not found")
        group_full = add_group(fileGrp)
        write_temp(group_full, "group_full")
        print ("ADD: (group_full) gexf: created")
        process_galaxy(galaxy, group_full)
        write_temp(galaxy, "galaxy") # list of galaxy without ones which belong to groups
        print ("PROCESS_GALAXY: (galaxy) gexf: created")
time("IMPORT/ADD: (group_full)")
#print ("group_full =", group_full.nodes())



## Filament
logSec("FILAMENT")
print ("PROCESS_TREE: (filament_MST_uncut) with trim: " + str(trimTree) + " (Mpc/h), cut: " + str(cutTree) + " (members)")
print ("\tAVG_TREE_WEIGHT: (filament_MST_uncut) length: " + str(avg_tree_weight(filament_MST_uncut)))
resultTree = process_tree(filament_MST_uncut, trimTree, cutTree)

filament = resultTree[0] ## center of filament group (with nodes and edges)
group_cent_non_filament = resultTree[1] # center of non-filament group
##filamentListWithoutNodesAttribute = resultTree[2] # use for walk and makebranch
filament_MST_cut = resultTree[3] # like filament_MST_uncut but already trimmed

write_temp(filament, "filament")
write_temp(group_cent_non_filament, "group_cent_non_filament")
write_temp(filament_MST_cut, "filament_MST_cut")
time("PROCESS_TREE: (filament_MST_uncut)")



## Group
logSec("GROUP")
resultGroup = process_group(filament, group_full)

group_mem_non_filament = resultGroup[0] # galaxy belonged to non-filament group
group_mem_filament = resultGroup[1] # galaxy belonged to filament group

write_temp(group_mem_non_filament, "group_mem_non_filament")
write_temp(group_mem_filament, "group_mem_filament")

time("PROCESS_GROUP: (group_full)")



## Scooper -> Galaxy Near Filament (GNF)
logSec("GALAXIES NEAR FILAMENT(GNF)")
print ("SCOOPER: (galaxy) with radius: " + str(scoopR) + " (Mpc/h), distance: " + str(scoopD) + " (Mpc/h)")
gnf = scooper(filament, galaxy, scoopR, scoopD) # galaxy near filament

##print "SCOOPER_GROUP: (galaxy) with radius: " + str(scoopR) + "(Mpc/h)"
##gnf = scooper_group(group_cent_non_filament, galaxy, gnf, 1) # use if we want to scoop from non-filament group

write_temp(gnf, "gnf")
time("SCOOPER: (galaxy)")


## Tendril and Void
logSec("TENDRIL & VOID")

remaining_isolated_galaxy = galaxy.copy()
write_temp(remaining_isolated_galaxy, "remaining_isolated_galaxy") # to change node ID in order to iterate faster
output_text_temp(remaining_isolated_galaxy, "remaining_isolated_galaxy")

## The following steps might seem redundant but in fact crucial to the script
## as re-importing the graph helps to speed up edge addition in finding MST of tendril'''

tendril_and_void = spanning_tree_tendril(make_path_temp("remaining_isolated_galaxy.txt"), 2)

time("SPANNING_TREE_TENDRIL: (filament_MST_uncut)")
print ("PROCESS_TENDRIL: (tendril_and_void) with trim: " + str(trimTendril) + " (Mpc/h), cut: " + str(cutTendril) + " (members)"
)
resultTendril = process_tendril(tendril_and_void, 1, 2)
output_text_temp(resultTendril[0], "tendril")
output_text_temp(resultTendril[1], "void")
write_temp(add_galaxy_tendril_and_void(make_path_temp("void.txt")),"void") # make void.gexf
write_temp(add_galaxy_tendril_and_void(make_path_temp("tendril.txt")),"tendril") # make tendril.gexf

tendril = read_temp("tendril")
void = read_temp("void")
time("PROCESS_TENDRIL: (filament_MST_uncut)")


## OUTPUT
## This is very important because here we input with textfiles so all the attributes are in float or int
## but for read_gexf all the atrribute are in STRING. this is why at this stage, galaxy_full and group_full and gnf
## all have attributes in string because they all origiate from some gexf file, but not void and tendril so there will
## be probelem with void and tendril galaxies. So we need the line below to read the void and tendril from GEXF
## this must be noted that the galaxy_full, group_full and gnf must also be loaded from gexf file and not text file!

logSec("OUTPUT")
inputList = [(read_temp("group_full"), '1'), (read_temp("gnf"), '2'), (read_temp("tendril"), '3'), (read_temp("void"), '4')]
print "OUTPUT_FINAL_LIST: (notation):"
print "\tGROUP: ", inputList[0][1]
print "\tGNF: ", inputList[1][1]
print "\tTENDRIL: ", inputList[2][1]
print "\tVOID: ", inputList[3][1]

galaxy_full_final = output_final_list(galaxy_full, inputList)
write_temp(galaxy_full_final,"galaxy_full_final")


## OUTPUT TO TEXT
output_text_final(galaxy_full_final, 'galaxy_full_final', folderName)

output_text_group_cent(filament_MST_uncut, 'filament_MST_uncut', folderName)
output_text_group_cent(filament, 'filament', folderName)
output_text_group_cent(group_cent_non_filament, 'group_cent_non_filament', folderName)

output_text(group_full, 'group_full', folderName)
output_text(group_mem_non_filament, 'group_mem_non_filament', folderName)
output_text(group_mem_filament, 'group_mem_filament', folderName)
output_text(gnf, 'gnf', folderName)
output_text(remaining_isolated_galaxy, 'remaining_isolated_galaxy', folderName)
output_text(tendril, 'tendril', folderName)
output_text(void, 'void', folderName)

copy_input(folderName)

time("OUTPUT:")

## STAT
logSec("STATS")

print "GALAXY_COUNT:\n\tALL: ", galaxy_full.number_of_nodes()
print "\tGROUP: ", group_full.number_of_nodes()
print "\tGROUP_MEMBER_IN_FILAMENT: ", group_mem_filament.number_of_nodes()
print "\tGROUP_MEMBER_NOT_IN_FILAMENT: ", group_mem_non_filament.number_of_nodes()
print "\tGALAXIES_NEAR_FILAMENT: ", gnf.number_of_nodes()
print "\tTENDRIL: ", tendril.number_of_nodes()
print "\tVOID: ", void.number_of_nodes()

print "GALAXY_RATIO:\n\tALL: " + "%.2f" % (float(group_full.number_of_nodes()+gnf.number_of_nodes()+tendril.number_of_nodes()+void.number_of_nodes())*100/galaxy_full.number_of_nodes())
print "\tGROUP: " + "%.2f" %  (float(group_full.number_of_nodes()*100)/galaxy_full.number_of_nodes())
print "\tGROUP_MEMBER_IN_FILAMENT: " + "%.2f" % (float(group_mem_filament.number_of_nodes()*100)/galaxy_full.number_of_nodes())
print "\tGROUP_MEMBER_NOT_IN_FILAMENT: " + "%.2f" % (float(group_mem_non_filament.number_of_nodes()*100)/galaxy_full.number_of_nodes())
print "\tGALAXIES_NEAR_FILAMENT: " + "%.2f" % (float(gnf.number_of_nodes()*100)/galaxy_full.number_of_nodes())
print "\tTENDRIL: " + "%.2f" % (float(tendril.number_of_nodes()*100)/galaxy_full.number_of_nodes())
print "\tVOID: " + "%.2f" % (float(void.number_of_nodes()*100)/galaxy_full.number_of_nodes())

# Main_Length_Specs.py
# Calculate the distance from galaxies to the filament backbone and output to galaxyWithD.txt
logSec("Main_Length_Specs")
galaxy = read_temp("galaxy")
filament = read_temp("filament")
group = read_temp("group_mem_non_filament")

gals_processed = galaxy.number_of_nodes() + group.number_of_nodes()

print "Processed " + str(gals_processed) + " galaxies."

galaxyWithD = nx.compose(Scoop(filament, galaxy), Scoop(filament, group))
write_temp(galaxyWithD, "galaxyWithD")
output_text_main_length_specs(galaxyWithD, 'galaxyWithD')

print "Main_Length_Specs:: done in: " + str(datetime.now() - startTime) + "\n"

time_off()
log_end(logFile)