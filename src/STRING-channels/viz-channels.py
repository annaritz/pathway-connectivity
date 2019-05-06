import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
import glob
import sys
import operator 

def main(infix):
	files = glob.glob('outfiles/%s-*-positive_sets.txt' % (infix))
	#files = ['outfiles/small_molecule_filter-database-positive_sets.txt']
	print('FILES:',files)
	for f in files:
		print('FILE %s' % (f))
		name = f.replace('outfiles/%s-' % (infix),'').replace('-positive_sets.txt','')
		print('NAME %s ' % (name))

		interactions = {}
		any_pathway = set()
		same_pathway = set()
		connected = {}
		with open(f) as fin:
			for line in fin:
				if line[0] == '#':
					continue
				row = line.strip().split()
				n1 = row[0]
				n2 = row[1]
				edge = (n1,n2)
				val = int(row[2])
				interactions[edge] = val

				if int(row[3]) == 1:	
					any_pathway.add(edge)
				if int(row[4]) == 1:
					same_pathway.add(edge)
				if int(row[5]) == 1:
					connected[edge] = int(row[6])

		file_infix = '%s-%s' % (infix,name)
		connected_set = set(connected.keys())
		bconn_set = set([key for key in connected.keys() if connected[key]==0])
		print('%d total; %d any pathway; %d same pathway; %d connected' % (len(interactions),len(any_pathway),len(same_pathway),len(connected_set)))
		venn(any_pathway,same_pathway,connected_set,['Pair in Any Pathway','Pair in Same Pathway','Pair Connected'],file_infix,name)
		precision(interactions,[any_pathway,same_pathway,connected_set,bconn_set],['Pair in Any Pathway','Pair in Same Pathway','Pair Connected','Pair B-Connected'],file_infix,name)
		print('DONE WITH %s\n' % (name))
	return

def venn(pos_1_set,pos_2_set,pos_3_set,pos_names,infix,name):
	plt.figure(figsize=(8,6))
	venn3([pos_1_set,pos_2_set,pos_3_set],tuple(['Pair in Any Pathway','Pair in Same Pathway','Pair Connected']))
	plt.title('STRING\'s "%s" channel' % (name))
	filename = 'outfiles/%s-venn.png' % (infix)
	plt.savefig(filename)
	print('wrote to '+filename)
	filename = filename.replace('png','pdf')
	plt.savefig(filename)
	print('wrote to '+filename)

def precision(interactions,pos_sets,pos_names,infix,name):
	plt.figure(figsize=(10,5))
	resolution=0.0001 ## set to 0 to keep all data points.
	
	counts = []
	for p in pos_sets:
		counts.append([])
	for node,val in sorted(interactions.items(), key=operator.itemgetter(1), reverse=True):
		for i in range(len(pos_sets)):
			if node in pos_sets[i]:
				counts[i].append(1)
			else:
				counts[i].append(0)
	print('Calculating terms...')
	xs_recall =[0]*len(pos_sets)
	xs_fpr =[0]*len(pos_sets)
	ys_prec =[0]*len(pos_sets)
	ys_tpr = [0]*len(pos_sets)
	for i in range(len(pos_sets)):
		xs_recall[i] = []
		xs_fpr[i] = []
		ys_prec[i] = []
		ys_tpr[i] =[]
		TP = 0
		P = len(pos_sets[i])
		N = len(interactions)-len(pos_sets[i])
		for j in range(len(counts[i])):
			TP+=counts[i][j]
			FP = j-TP
			FN = P-TP
			TN = N-FP

			if len(xs_recall[i]) == 0 or not ( abs(xs_recall[i][-1]-TP/(TP+FN)) <= resolution and abs(ys_prec[i][-1]-TP/(TP+FP)) <= resolution ):
				xs_recall[i].append(TP/(TP+FN))
				if TP+FP == 0:
					ys_prec[i].append(0)
				else:	
					ys_prec[i].append(TP/(TP+FP))
			if len(xs_fpr[i]) == 0 or not ( abs(xs_fpr[i][-1]-FP/(FP+TN)) <= resolution and abs(ys_tpr[i][-1]-TP/(TP+FN)) <= resolution ):
				xs_fpr[i].append(FP/(FP+TN))
				ys_tpr[i].append(TP/(TP+FN))

	print('done: dimensions are %d and %d' % (len(xs_recall[i]),len(xs_fpr[i])))

	ax = plt.subplot(1,2,1)
	for i in range(len(pos_sets)):
		ax.plot(xs_recall[i],ys_prec[i],lw=2,label='%s\n(|P|=%d)' % (pos_names[i],len(pos_sets[i])))
	ax.set_xlabel('Recall')
	ax.set_ylabel('Precision')
	ax.set_title('Precision and Recall of\nSTRING\'s "%s" channel\n(%d interactions)' % (name,len(interactions)))
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)
	ax.legend()

	ax = plt.subplot(1,2,2)
	ax.plot([0,0],[1,1],'--k',label='__nolabel__')
	for i in range(len(pos_sets)):
		ax.plot(xs_fpr[i],ys_tpr[i],lw=2,label='%s\n(|P|=%d)' % (pos_names[i],len(pos_sets[i])))
	ax.set_xlabel('False Positive Rate')
	ax.set_ylabel('True Positive Rate')
	ax.set_title('ROC of\nSTRING\'s "%s" channel\n(%d interactions)' % (name,len(interactions)))
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)
	ax.legend()

	plt.tight_layout()
	filename = 'outfiles/%s-precision.png' % (infix)
	plt.savefig(filename)
	print('wrote to '+filename)
	filename = filename.replace('png','pdf')
	plt.savefig(filename)
	print('wrote to '+filename)

# def get_min_dist(n1_nodes,dist_dict_n1,n2_nodes,dist_dict_n2):
# 	val = LARGE_VAL
# 	for n in n2_nodes:
# 		if dist_dict_n1[n] != None:
# 			val = min(dist_dict_n1[n],val)
# 	for n in n1_nodes:
# 		if dist_dict_n2[n] != None:
# 			val = min(dist_dict_n2[n],val)
# 	return val

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print('USAGE: python3 viz-channels.py <INFIX>')
	main(sys.argv[1])