from graphspace_python.api.client import GraphSpace
graphspace = GraphSpace('compbio@reed.edu', 'compbio')
from graphspace_python.graphs.classes.gsgraph import GSGraph
import sys 

from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform

import hgraph_utils
from viz_utils import *
    
PAIRS = [['Signaling-by-MST1','Signaling-by-MET',4],
    ['Signaling-by-Activin','Signaling-by-TGF-beta-Receptor-Complex',3],
    ['Signaling-by-BMP','Signaling-by-TGF-beta-Receptor-Complex',3]
    ]

def main(prefix,brelaxdir,outfile):
    names = name_dict()
    H, identifier2id, id2identifier = hgraph_utils.make_hypergraph(prefix)
    
    for source,target,k in PAIRS:
        
        bdist,visited_hedges,source_members,target_members = read_bdist(source,target,brelaxdir,k)
        print('%s -> %s (k=%d)\n  source has %d members and target has %d members' % 
            (source,target,k,len(source_members),len(target_members)))

        subH = get_subhypergraph(H,visited_hedges,id2identifier)
        print('  Induced sub-hypergraph has %d hyperedges and %d nodes' % 
        (stats.number_of_hyperedges(subH),stats.number_of_nodes(subH)))
        
        G = gradient_overlap_survey(subH, source_members, target_members, bdist,k)
        
        G.set_name('%s TO %s k=%d' % (source,target,k))
        G.set_tags(['GLBio2019'])
        print('Posting',source,target,k)
        post(G,graphspace)
        
    return
def get_subhypergraph(H,visited_hedges,id2identifier):
    hids = H.get_hyperedge_id_set()
    for hid in hids:
        if id2identifier[hid] not in visited_hedges:
            H.remove_hyperedge(hid)
    nodes = H.get_node_set()
    for n in nodes:
        if len(H.get_forward_star(n)) == 0 and len(H.get_backward_star(n)) == 0:
            H.remove_node(n)
    return H

def gs_add_hyperedge(G, H, hedge):
    tail = H.get_hyperedge_tail(hedge)
    head = H.get_hyperedge_head(hedge)
    if len(tail) == 1 and len(head) == 1:
        t_node = tail.pop()
        tail.add(t_node)
        h_node = head.pop()
        head.add(t_node)
        G.add_edge(t_node, h_node, directed=True,popup='%s %s' % (hedge,H.get_hyperedge_attribute(hedge,'identifier')))
        G.add_edge_style(t_node,h_node,directed=True)
    else:
        attrs = H.get_hyperedge_attributes(hedge)
        G.add_node(hedge, popup='%s %s' % (hedge,H.get_hyperedge_attribute(hedge,'identifier')))
        G.add_node_style(hedge, height=10, width=10, color='white')
        for n in tail:
            G.add_edge(n,hedge)
        for n in head:
            G.add_edge(hedge, n, directed=True)
            G.add_edge_style(hedge,n,directed=True)
    return

def name_dict():
    namefile = '..//BioPAXSTREAM/output/reactome_limit30_filtered.txt.names'
    names = {}
    with open(namefile) as fin:
        for line in fin:
            row = line.strip().split('\t')
            names[row[0]] = row[1]
    return names

def scalar_mult(scalar, vector):
    return[x*scalar for x in vector]

def vector_add(vector1, vector2):
    if len(vector1) != len(vector2):
        raise ValueError("Error in attempted vector addition: the two vectors need to be the same length")
    ls = []
    for i in range(len(vector1)):
        ls.append(vector1[i] + vector2[i])
    return ls
def assign_gradient(color1, color2, current_step, total_steps):
    term1 = scalar_mult(1 - (current_step / float(total_steps)), color1)
    term2 = scalar_mult(current_step/float(total_steps), color2)
    final_color = [int(x) for x in vector_add(term1, term2)]
    return '#%02x%02x%02x' % tuple(final_color)

# SUPER hacky code to get the names to be pretty.
def parse_name(node_name):
    orig_name = node_name
    for candidate in node_name.split(';'):
        if len(candidate) < len(node_name):
            node_name = candidate
    splitname = node_name.split('-')
    if ':' not in node_name and len(splitname)>2 and any([len(a)>8 for a in splitname]):
        #print(node_name)
        ind = int(len(splitname)/2)
        node_name = '-'.join(splitname[:ind])+'\n'+'-'.join(splitname[ind:])
        #print(node_name)
    node_name = node_name.replace(':','\n')
    if 'STAT3-upregulated' in node_name:
        node_name = node_name.replace(' ','\n')

    width = max([len(l) for l in node_name.split('\n')])*15
    height = len(node_name.split('\n'))*25
    return node_name,orig_name,height,width


def gradient_overlap_survey(H, source_members, target_members, bdist, k):
    names = name_dict()
    G = GSGraph()
    for node in H.node_iterator():
        name = names.get(node,'NONE')
        name,orig_name,height,width = parse_name(name)

        if node in target_members:
            shape='star'
            height=80
            width=80
        else:
            shape='ellipse'
            height=80
            width=80

        if node not in bdist:
            layer_color = '#aab3c1'
            shape='ellipse'
            height=40
            width=40
        else:
            layer_color = assign_gradient([100,100,255],[0,255,0], bdist[node], k)
        popup = '%s<br>%s<br>InSource? %s<br>InTarget? %s<br>InConnectedSet? %s' % (name,node,node in source_members,node in target_members,node in bdist.keys())
        G.add_node(node, label=name, popup = popup)
        G.add_node_style(node, height=height, width=width, color=layer_color, shape=shape)

    for hedge_id in H.hyperedge_id_iterator():
        #print('Adding hyperege ' + hedge_id)
        #print('  %d Tails:' % len(H.get_hyperedge_tail(hedge_id)),[names[n] for n in H.get_hyperedge_tail(hedge_id)])
        #print('  %d Heads:'% len(H.get_hyperedge_head(hedge_id)),[names[n] for n in H.get_hyperedge_head(hedge_id)])
        gs_add_hyperedge(G, H, hedge_id)

    return G

def read_bdist(source,target,brelaxdir,k_limit):
    bfile = '%s/small_molecule_filter_%s_b_relax.txt' % (brelaxdir,source)
    print(bfile)
    bdist = {}
    traversed = set()
    with open(bfile) as fin:
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split()
            if row[0] == '-1':
                source_members = set(row[2].split(';'))
            else:
                k = int(row[0])
                if k > k_limit:
                    break
                for member in row[2].split(';'):
                    bdist[member]=k
                traversed.update(row[3].split(';'))
    ## get target members from k=-1 line.
    bfile = '%s/small_molecule_filter_%s_b_relax.txt' % (brelaxdir,target)
    print(bfile)
    with open(bfile) as fin:
        for line in fin:
            row = line.strip().split()
            if row[0] == '-1':
                target_members = set(row[2].split(';'))
                break

    return bdist,traversed,source_members,target_members



##########################
def post(G,graphspace_client):
  '''
  Post a graph to graphspace. 
  Inputs: Graph (GSGraph Object), GraphSpace client.
  Outputs: the graph ID (int)
  '''
  try:
    # try updating the graph. If the graph does not exist, 
    # this will throw an error.
    graph = graphspace_client.update_graph(G)
  except:
    # catch the error and try posting a new graph.
    graph = graphspace_client.post_graph(G)
  return graph.id

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('USAGE: python post_to_graphspace <PREFIX> <BRELEAX_DIR> <outfile> ')
        sys.exit()
    main(sys.argv[1],sys.argv[2],sys.argv[3])