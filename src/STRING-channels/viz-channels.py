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
		viz(interactions,[any_pathway,same_pathway,connected_set,bconn_set],\
			['Pair in\nAny Pathway','Pair in\nSame Pathway','Pair\nConnected','Pair\nB-Connected'],\
			['Pair in Any Pathway','Pair in Same Pathway','Pair Connected','Pair B-Connected'],file_infix,name,brelax=connected)
		print('DONE WITH %s\n' % (name))
		#sys.exit()
	return

def viz(interactions,pos_sets,short_names,pos_names,infix,name,brelax=None):
	plt.figure(figsize=(14,5))
	resolution=0.001 ## set to 0 to keep all data points.

	ax = plt.subplot(1,3,1)
	vd = venn3(pos_sets[:3],set_labels=tuple(short_names[:3]))
	# reposition labels
	
	lbl = vd.get_label_by_id('A')
	x, y = lbl.get_position()
	lbl.set_position((x, y-0.1))
	lbl = vd.get_label_by_id('B')
	x, y = lbl.get_position()
	lbl.set_position((x-0.1, y-0.1))
	plt.title('STRING\'s "%s" channel\n(%d interactions)\n' % (name,len(interactions)))
	
	items = sorted(interactions.items(), key=operator.itemgetter(1), reverse=True)
	node_list = [i[0] for i in items]
	tie_list = [0]*len(node_list)
	for i in range(len(node_list)):
		if i > 0 and items[i][1] == items[i-1][1]:
			tie_list[i] = 1
	print('%d (%.4f) ties' % (sum(tie_list),sum(tie_list)/len(tie_list)))
	print([i[1] for i in items][:10])
	print(tie_list[:10])

	counts = []
	for p in pos_sets:
		counts.append([])
	for node,val in items:
		for i in range(len(pos_sets)):
			if node in pos_sets[i]:
				counts[i].append(1)
			else:
				counts[i].append(0)

	if brelax: # if not None, draw brelax curves. This is a dictionary.
		brelax_thresholds = sorted(list(set(list(brelax.values()))))
		## drop first (bconn) and last (connected) values - these are already considered pos sets.
		brelax_thresholds = brelax_thresholds[1:len(brelax_thresholds)-1]
		print('BRLEAX_thresholds:',brelax_thresholds)
		brelax_counts = []
		brelax_nums = [0]*len(brelax_thresholds)
		for b in brelax_thresholds:
			brelax_counts.append([])
		for node,val in items:
			for i in range(len(brelax_thresholds)):
				b = brelax_thresholds[i]
				if node in brelax and brelax[node] <= b:
					brelax_counts[i].append(1)
					brelax_nums[i]+=1
				else:
					brelax_counts[i].append(0)

	print('Calculating terms...')
	xs_recall =[0]*len(pos_sets)
	xs_fpr =[0]*len(pos_sets)
	ys_prec =[0]*len(pos_sets)
	ys_tpr = [0]*len(pos_sets)
	if brelax: # show intermediate lines for ROC only.
		b_tpr = [0]*len(brelax_thresholds)
		b_fpr = [0]*len(brelax_thresholds)

	for i in range(len(pos_sets)):
		xs_recall[i] = []
		xs_fpr[i] = []
		ys_prec[i] = []
		ys_tpr[i] =[]
		TP = 0
		P = len(pos_sets[i])
		N = len(interactions)-P
		for j in range(len(counts[i])):
			TP+=counts[i][j]
			if tie_list[j] and j != len(counts[i])-1: # skip if there's a tie
				continue

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

		print('  %d: dimensions are %d and %d' % (i,len(xs_recall[i]),len(xs_fpr[i])))

	if brelax:
		print('brelax')
		for i in range(len(brelax_thresholds)):
			b_tpr[i] = []
			b_fpr[i] = []
			TP = 0
			P = brelax_nums[i]
			N = len(interactions)-P
			for j in range(len(brelax_counts[i])):
				TP+=brelax_counts[i][j]
				if tie_list[j] and j != len(brelax_counts[i])-1: # skip if there's a tie
					continue
				FP = j-TP
				FN = P-TP
				TN = N-FP

				if len(b_fpr[i]) == 0 or not ( abs(b_fpr[i][-1]-FP/(FP+TN)) <= resolution and abs(b_tpr[i][-1]-TP/(TP+FN)) <= resolution ):
					b_fpr[i].append(FP/(FP+TN))
					b_tpr[i].append(TP/(TP+FN))
			print('  %d: dimensions are %d' % (i,len(b_tpr[i])))


	ax = plt.subplot(1,3,2)
	for i in range(len(pos_sets)):
		ax.plot(xs_recall[i],ys_prec[i],lw=3,label='%s\n(|P|=%d)' % (pos_names[i],len(pos_sets[i])))
	ax.set_xlabel('Recall')
	ax.set_ylabel('Precision')
	ax.set_title('Precision and Recall of\nSTRING\'s "%s" channel\n(%d interactions)' % (name,len(interactions)))
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)
	ax.legend()

	ax = plt.subplot(1,3,3)
	ax.plot([0,0],[1,1],'--k',label='__nolabel__')
	if brelax:
		for i in range(len(brelax_thresholds)):
			ax.plot(b_fpr[i],b_tpr[i],color=[.8,.8,.8],lw=2,label='__nolabel__')#label='BDist$\leq$%d\n(|P|=%d)' % (brelax_thresholds[i],brelax_nums[i]))
	for i in range(len(pos_sets)):
		ax.plot(xs_fpr[i],ys_tpr[i],lw=2,label='%s\n(|P|=%d)' % (pos_names[i],len(pos_sets[i])))
	if brelax:
		if len(brelax_thresholds) < 10:
			ax.plot(0,0,color=[.8,.8,.8],label='B-Relaxation Thresholds\nd=[%s]' % (','.join([str(s) for s in brelax_thresholds])))
		else:
			ax.plot(0,0,color=[.8,.8,.8],label='B-Relaxation Thresholds\nd=%d...%d' % (brelax_thresholds[0],brelax_thresholds[-1]))
	
	ax.set_xlabel('False Positive Rate')
	ax.set_ylabel('True Positive Rate')
	ax.set_title('ROC of\nSTRING\'s "%s" channel\n(%d interactions)' % (name,len(interactions)))
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)
	ax.legend()

	plt.tight_layout()
	filename = 'outfiles/%s-viz.png' % (infix)
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