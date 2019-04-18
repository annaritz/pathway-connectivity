## Modified from /Users/aritz/Research/VT/vt-svn/src/python/SignalingHypergraphs/minHyperpathILP.py
import hgraph_utils as hgraph_utils
from halp import directed_hypergraph
from halp.algorithms import directed_paths as hpaths
from halp.utilities import directed_statistics as stats
from halp.utilities import directed_graph_transformations as transform
import cplex 

EPSILON=0.0001
CONSTANT=1000000
outprefix = 'tmp'
n_map = {}
he_map = {}

def populate_maps(H):
	i = 0
	for node in H.get_node_set():
		n_map[node] = 'a_%d' % (i)
		i+=1
	j = 0
	for hid in H.get_hyperedge_id_set():
		he_map[hid] = 'e_%d' % (j)
		j +=1
	return 

def runILP(H,s,t):
	## H = hypergraph
	## s = source node
	## t = target node
	node_set = H.get_node_set()

	if s not in node_set:
		sys.exit('ERROR: source %s is not a hypernode in H.' % (s))
	if t not in node_set:
		sys.exit('ERROR: target %s is not a hypernode in H.' % (t))

	print('Original Signaling Hypergraph has %d nodes and %d hyperedges' % (stats.number_of_nodes(H),stats.number_of_hyperedges(H)))
	populate_maps(H)

	# write the .lp file - cplex files have an '.lp' suffix.
	print('\nWriting ILP...')
	out = open(outprefix+'.lp','w')
	out = writeObjective(H,out)
	out.write('SUBJECT TO\n')
	c = 0
	out,c = writeActivityVariableConstraints(H,s,out,c)

	## add "t=1" to tofix dict
	tofix = {t:1}
	
	out,c = writeFixedValueConstraints(H,tofix,out,c)
	out,c = writeOrderVariableConstraints(H,out,c)
	out = writeActivityVariableTypes(H,out)
	out.write('END\n')
	out.close()
	print('%d Constraints Written.' % (c))
	print('Done Writing ILP.')

	# solve the ILP for a single, optimal solution.
	objective = solveILP(outprefix) 
	print('Done: objective is %d' % (objective))
	
	## return variables (first solution indexed at 0)
	return objective


######################################
## Objective Function:
######################################
def writeObjective(H,out):
	out.write('MINIMIZE\n')
	out.write('obj: ')
	# Minimize Sum of Active Edges
	for e in H.get_hyperedge_id_set():
		out.write(' + 1 %s' % (he_map[e]))
	out.write('\n')
	return out

######################################
## Activity Variable Constraints
######################################
def writeActivityVariableConstraints(H,s,out,c):

	######################################
	## Active Edges: Constraint for each of T(e),H(e)
	for e in H.get_hyperedge_id_set():
		tail = H.get_hyperedge_tail(e)
		if tail != None:
			## format: \sum_{v \in T(e)} \alpha_v >= |T(e)| \alpha_e
			##       = \sum_{v \in T(e)} \alpha_v - |T(e)| \alpha_e >= 0
			out.write('hyperedge_tails_%s_c%d:' % (he_map[e],c))
			for v in tail:
				out.write(' + %s' % (n_map[v]))
			out.write(' - %d %s >= 0\n' % (len(tail),he_map[e]))
			c+=1 #increment constraint counter

		head = H.get_hyperedge_head(e)
		if head != None:
			## format: \sum_{v \in H(e)} \alpha_v >= |H(e)| \alpha_e
			##       = \sum_{v \in H(e)} \alpha_v - |H(e)| \alpha_e >= 0
			out.write('hyperedge_heads_%s_c%d:' % (he_map[e],c))
			for v in head:
				out.write(' + %s' % (n_map[v]))
			out.write(' - %d %s >= 0\n' % (len(head),he_map[e]))
			c+=1 #increment constraint counter         
   
	######################################
	## Incoming Edges: at least one edge in BS(v) must be active 
	## (except for the src)
	for node in H.get_node_set(): 
		## don't write constraint for the source.
		if node == s:
			continue

		## get backwards star:
		bstar = H.get_backward_star(node)
		if len(bstar)==0:
			## format if empty bstar: \alpha_v = 0
			print('WARNING: hnode %s has an empty backward-star. Alpha var must be 0.' % (node))
			out.write('incoming_%s_c%d: %s = 0\n' % (n_map[node],c,n_map[node]))
		else:
			## format for non-empty bstar: 
			##   \sum_{e \in BS(v)} \alpha_e >= \alpha_v
			## = \sum_{e \in BS(v)} \alpha_e - \alpha_v >= 0
			out.write('incoming_%s_c%d:' % (n_map[node],c))
			for e in bstar:
				out.write(' + %s' % (he_map[e]))
			out.write(' - %s >= 0\n' % (n_map[node]))    
		c+=1 #increment constraint counter

	return out,c

######################################
## Fixed Value Constraints
######################################
def writeFixedValueConstraints(H,tofix,out,c):
	# requre that at least one hypernode containing t is active.
	## tofix is a dictionary of {hypernode: <0 or 1>}

	for t in tofix:
		out.write('fixed_%s_c%d: %s = %d\n' % (n_map[t],c,n_map[t],tofix[t]))  
		c+=1 #increment constraint counter
	return out,c

######################################
## Order Variable Constraints
######################################
def writeOrderVariableConstraints(H,out,c):

	######################################
	## Order Bounds: 0 <= o_v <= \alpha_v \forall v \in V
	for node in H.get_node_set():
		## format: o_v >= 0
		out.write('order_lb_%s_c%d: o_%s >= 0\n' % (n_map[node],c,n_map[node]))
		c+=1 #increment constraint counter

		## format: o_v - \alpha_v <= 0
		out.write('order_ub_%s_%d: o_%s - %s <= 0\n' % (n_map[node],c,n_map[node],n_map[node]))
		c+=1 #increment constraint counter

	######################################
	## Order Constraints: o_u <= o_v - \epsilon + C (1-\alpha_e) \forall u,v in (T(e),H(e)) \forall e \in E 
	for e in H.get_hyperedge_id_set(): 
		for u in H.get_hyperedge_tail(e):
			for v in H.get_hyperedge_head(e):
				## format: o_u <= o_v - \epsilon + C (1-\alpha_e) 
				##       = o_u - o_v + C \alpha_e <= C - \epsilon
				out.write('order_%s_c%d: o_%s - o_%s + %d %s <= %f\n' % \
							  (he_map[e],c,n_map[u],n_map[v],CONSTANT,he_map[e],CONSTANT-EPSILON))
				c+=1 #increment constraint counter
		   
	return out,c

######################################
## Activity Variable Types
######################################
def writeActivityVariableTypes(H,out):
	out.write('BINARY\n')
	for node in H.get_node_set():
		out.write(' %s\n' % (n_map[node]))
	for e in H.get_hyperedge_id_set():
		out.write(' %s\n' % (he_map[e]))
	return out

################################
def solveILP(outprefix):
	print('\nSolving ILP...')

	print('\n' + '-'*20 + 'Cplex Output Start' + '-'*20)
	ilp = cplex.Cplex()
	ilp.read('%s.lp' % (outprefix))

	numsolsfound = 1
	numoptobjective = 0
	maxobj = None
	allvars = []
	
	## Solve ILP
	print('-'*10 + 'Looking for Solution' )
	ilp.solve()
	print('-'*10 + 'Done Looking for Solution')
		
	if ilp.solution.pool.get_num()>0:
		print('Solution Found.'      )
		objective = ilp.solution.pool.get_objective_value(0)
	else:
		print('Infeasible Solution. quitting.')
		sys.exit()
		
	print('-'*20 + 'Cplex Output End' + '-'*20 + '\n')
	return objective

