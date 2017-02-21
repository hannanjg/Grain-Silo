import random
import sys
import math
import os.path
import shutil
from datetime import datetime
import networkx as nx

######################################################################################################################################################
# OUTPUT
######################################################################################################################################################

## export a graph to text file
def output_text_final(graph, filename):
        
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
                             
######################################################################################################################################################
# TEMP
######################################################################################################################################################

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

######################################################################################################################################################
# SCOOP
######################################################################################################################################################

def dist(a, b): ##tuples (x,y) of coordinate 2D
        return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)

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
                                
## Scooper
def Scoop(G1, G2):
        
        '''
        ## G1: filaments (nodes and egdes)
        ## G2: list of galaxies
        
        ## G3: galaxies with min distance from tree
        '''

        print "SCOOP for MIN DISTANCE:", G2.number_of_nodes()
        
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


## Variable

startTime = datetime.now()

C  = 299792.458 # speed of light (km/s) - use for d = C*z/100 (Mpc/h)
D = math.acos(-1.0) / 180 # acos(-1) = Pi


######################################################################################################################################################
#MAIN
######################################################################################################################################################

galaxy = read_temp("galaxy")
filament = read_temp("filament")
group = read_temp("group_mem_non_filament")

galaxyWithD = nx.compose(Scoop(filament, galaxy), Scoop(filament, group))
write_temp(galaxyWithD, "galaxyWithD")
output_text_final(galaxyWithD, 'galaxyWithD')

print "TOTAL DURATION: " + str(datetime.now() - startTime)
print "\n" + str(datetime.now()) + "\n"
