import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
import glob
import sys
import operator 
import numpy as np
from scipy.stats import kruskal

## to draw rectangles
# https://matplotlib.org/examples/shapes_and_collections/artist_reference.html
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def main(infix):
	files = glob.glob('outfiles/%s-*-positive_sets.txt' % (infix))
	files = ['outfiles/small_molecule_filter_allpathways-cooccurence-positive_sets.txt']
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

		file_infix = 'outfiles/%s-%s-viz' % (infix,name)
		connected_set = set(connected.keys())
		bconn_set = set([key for key in connected.keys() if connected[key]==0])
		print('%d total; %d any pathway; %d same pathway; %d connected' % (len(interactions),len(any_pathway),len(same_pathway),len(connected_set)))
		viz(interactions,[any_pathway,same_pathway,connected_set,bconn_set],\
			['Pair in\nAny Pathway','Pair in\nSame Pathway','Pair\nConnected','Pair\nB-Connected'],\
			['Pair in Any Pathway','Pair in Same Pathway','Pair Connected','Pair B-Connected'],file_infix,name,brelax=connected)
		print('DONE WITH %s\n' % (name))
		#sys.exit()
	return


def viz(interactions,pos_sets,short_names,pos_names,prefix,name,brelax=None):
	plt.figure(figsize=(10,10))
	#resolution=0.001 ## set to 0 to keep all data points.
	resolution=0.0

	ax = plt.subplot(2,2,1)
	vd = venn3(pos_sets[:3],set_labels=tuple(short_names[:3]))
	# # reposition labels	
	# lbl = vd.get_label_by_id('A')
	# x, y = lbl.get_position()
	# lbl.set_position((x, y-0.1))
	# lbl = vd.get_label_by_id('B')
	# x, y = lbl.get_position()
	# lbl.set_position((x-0.1, y-0.1))
	plt.title('STRING\'s "%s" channel\n(%d interactions map to Reactome)\n' % (name,len(interactions)))
	

	### NEW: Make the "Pair in Any Pathway" the universe
	print('NEW: Make the "Pair in Any Pathway" the universe for plotting.')
	interactions = {e:interactions[e] for e in pos_sets[0]}
	pos_sets = pos_sets[1:]
	pos_names = pos_names[1:]

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
		# tmp_file = 'tmpfiles/'+infix+'-%d.txt' % (i)
		# tmp_out = open(tmp_file,'w')
		# tmp_out.write('#%s\n' % (pos_names[i]))
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

			FP = (j+1)-TP
			FN = P-TP
			TN = N-FP
			
			#tmp_out.write('j=%d: P=%d, N=%d, TP=%d, FP=%d, FN=%d, TN=%d\n' % (j,P,N,TP,FP,FN,TN))

			if len(xs_recall[i]) == 0 or not ( abs(xs_recall[i][-1]-TP/(TP+FN)) <= resolution and abs(ys_prec[i][-1]-TP/(TP+FP)) <= resolution ):
				xs_recall[i].append(TP/(TP+FN))
				if TP+FP == 0:
					ys_prec[i].append(0)
				else:	
					ys_prec[i].append(TP/(TP+FP))
				#tmp_out.write('     REC=%.2f, PREC=%.2f\n' % (xs_recall[i][-1],ys_prec[i][-1]))
			if len(xs_fpr[i]) == 0 or not ( abs(xs_fpr[i][-1]-FP/(FP+TN)) <= resolution and abs(ys_tpr[i][-1]-TP/(TP+FN)) <= resolution ):
				xs_fpr[i].append(FP/(FP+TN))
				ys_tpr[i].append(TP/(TP+FN))
				
				#tmp_out.write('     FPR=%.2f, TPR=%.2f\n' % (xs_fpr[i][-1],ys_tpr[i][-1]))

			assert ys_tpr[i][-1] == xs_recall[i][-1]
		#tmp_out.close()
		#print('wrote to %s' % (tmp_file))
		print('  %d: dimensions are %d and %d' % (i,len(xs_recall[i]),len(xs_fpr[i])))

	if brelax:
		print('brelax...')
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
				FP = (j+1)-TP
				FN = P-TP
				TN = N-FP

				if len(b_fpr[i]) == 0 or not ( abs(b_fpr[i][-1]-FP/(FP+TN)) <= resolution and abs(b_tpr[i][-1]-TP/(TP+FN)) <= resolution ):
					b_fpr[i].append(FP/(FP+TN))
					b_tpr[i].append(TP/(TP+FN))
			#print('  %d: dimensions are %d' % (i,len(b_tpr[i])))


	ax = plt.subplot(2,2,2)
	for i in range(len(pos_sets)):
		ax.plot(xs_recall[i],ys_prec[i],lw=2,label='%s\n(|P|=%d; |N|=%d)' % (pos_names[i],len(pos_sets[i]),len(interactions)-len(pos_sets[i])))
	ax.set_xlabel('Recall')
	ax.set_ylabel('Precision')
	ax.set_title('Precision and Recall of\nSTRING\'s "%s" channel\n(%d interactions in any Reactome pathway)' % (name,len(interactions)))
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)
	ax.legend()

	ax = plt.subplot(2,2,3)
	
	patches = []
	prop_cycle = plt.rcParams['axes.prop_cycle']
	colors = prop_cycle.by_key()['color']
	for i in range(len(pos_sets)):
		P = len(pos_sets[i])
		N = len(interactions)-len(pos_sets[i])
		b_x = xs_fpr[i][0]*N
		t_x = xs_fpr[i][-1]*N
		b_y = ys_tpr[i][0]*P
		t_y = ys_tpr[i][-1]*P
		ax.add_patch(mpatches.Rectangle([b_x,b_y],t_x-b_x, t_y-b_y,alpha=0.2,color=colors[i]))


	#if brelax:
	#	for i in range(len(brelax_thresholds)):
	#		ax.plot(b_fpr[i],b_tpr[i],color=[.8,.8,.8],lw=2,label='__nolabel__')#label='BDist$\leq$%d\n(|P|=%d)' % (brelax_thresholds[i],brelax_nums[i]))
	for i in range(len(pos_sets)):
		P = len(pos_sets[i])
		N = len(interactions)-len(pos_sets[i])
		ax.plot([x*N for x in xs_fpr[i]],[y*P for y in ys_tpr[i]],lw=2,label='%s\n(|P|=%d; |N|=%d)' % (pos_names[i],P,N))
	#if brelax:
	#	if len(brelax_thresholds) < 10:
	#		ax.plot(0,0,color=[.8,.8,.8],label='B-Relaxation Thresholds\nd=[%s]' % (','.join([str(s) for s in brelax_thresholds])))
	#	else:
	#		ax.plot(0,0,color=[.8,.8,.8],label='B-Relaxation Thresholds\nd=%d...%d' % (brelax_thresholds[0],brelax_thresholds[-1]))
	ax.set_xlabel('# of False Positives')
	ax.set_ylabel('# of True Positives')
	ax.set_title('False Positives and True Positives of\nSTRING\'s "%s" channel\n(%d interactions in any Reactome pathway)' % (name,len(interactions)))
	ax.legend()

	ax = plt.subplot(2,2,4)
	#ax.plot([0,1],[0,1],'--k',label='__nolabel__')
	if brelax:
		for i in range(len(brelax_thresholds)):
			ax.plot(b_fpr[i],b_tpr[i],color=[.8,.8,.8],lw=2,label='__nolabel__')#label='BDist$\leq$%d\n(|P|=%d)' % (brelax_thresholds[i],brelax_nums[i]))
	for i in range(len(pos_sets)):
		P = len(pos_sets[i])
		N = len(interactions)-len(pos_sets[i])
		ax.plot(xs_fpr[i],ys_tpr[i],lw=2,label='%s\n(|P|=%d; |N|=%d)' % (pos_names[i],P,N))
	if brelax:
		if len(brelax_thresholds) < 10:
			ax.plot(0,0,color=[.8,.8,.8],label='B-Relaxation Thresholds\nd=[%s]' % (','.join([str(s) for s in brelax_thresholds])))
		else:
			ax.plot(0,0,color=[.8,.8,.8],label='B-Relaxation Thresholds\nd=%d...%d' % (brelax_thresholds[0],brelax_thresholds[-1]))
	ax.set_xlabel('False Positive Rate')
	ax.set_ylabel('True Positive Rate')
	ax.set_title('ROC of\nSTRING\'s "%s" channel\n(%d interactions in any Reactome pathway)' % (name,len(interactions)))
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)
	ax.legend()

	plt.tight_layout()
	filename = '%s.png' % (prefix)
	plt.savefig(filename)
	print('wrote to '+filename)
	filename = filename.replace('png','pdf')
	plt.savefig(filename)
	print('wrote to '+filename)
	return

def viz_box_plot(interactions,pos_sets,short_names,pos_names,prefix,name,brelax=None):
	
	if brelax:

		print('brelax...')
		### NEW: Make the "Pathway Connected" the universe
		print('NEW: Make the "Pathway Connected" the universe for plotting.')
		interactions = {e:interactions[e] for e in pos_sets[2]}
		pos_sets = pos_sets[2:]
		pos_names = pos_names[2:]
		brelax_thresholds = sorted(list(set(list(brelax.values()))))
		## drop first (bconn) and last (connected) values - these are already considered pos sets.
		#brelax_thresholds = brelax_thresholds[1:len(brelax_thresholds)-1]
		print('BRLEAX_thresholds:',brelax_thresholds)
		brelax_counts = []
		for b in brelax_thresholds:
			brelax_counts.append([])
		for node,val in interactions.items():
			for i in range(len(brelax_thresholds)):
				b = brelax_thresholds[i]
				if node in brelax and brelax[node] <= b:
					brelax_counts[i].append(node)

		pos_sets = brelax_counts
		pos_names = ['Pair B-Connected'] + ['B-Dist=%d' % (d) for d in brelax_thresholds[1:-1]] + ['Pair Connected']

		## drop every other if it's too many
		while len(pos_sets) > 15:
			if len(pos_sets) > 30:
				jump = 2
			else:
				jump = 3

			pos_sets = [pos_sets[0]] + [pos_sets[i] for i in range(1,len(pos_sets)-1,jump)] + [pos_sets[-1]]
			pos_names = [pos_names[0]] + [pos_names[i] for i in range(1,len(pos_names)-1,jump)] + [pos_names[-1]]
	else:
		print('NEW: Make the "In Anyt Pathway" the universe for plotting.')
		interactions = {e:interactions[e] for e in pos_sets[0]}
		pos_sets = pos_sets[1:]
		pos_names = pos_names[1:]
	print('PLOTTING %d TOTAL COLS' % (len(pos_sets)))
	items = sorted(interactions.items(), key=operator.itemgetter(1), reverse=True)
	vals = []
	for i in range(len(pos_sets)):
		vals.append([])
	print('assigning to positive sets')
	for i in range(len(pos_sets)):
		this_set = set(pos_sets[i])
		vals[i] = [val for edge,val in items if edge in this_set]

	if not brelax:
		vals = [[i[1] for i in items]] + vals
		pos_names = ['Pair in Any Pathway'] + pos_names
	print('done gathering data.')

	if brelax:
		plt.figure(figsize=(len(pos_sets)*1.2,4))
	else:
		plt.figure(figsize=(8,4))
	#ax1 = plt.subplot(1,2,2)
	ax1 = plt.subplot(1,1,1)
	parts = ax1.violinplot(vals,points=1000,showmeans=False, showmedians=False,showextrema=False)
	labels = ['%s\n($n=%d$)' % (pos_names[i],len(vals[i])) for i in range(len(vals))]
	if brelax:
		colors = ['#541684'] + ['#168426']*(len(pos_names)-2) + ['#541684']
	else:
		colors = ['#494949','#541684','#541684','#541684']
	if len(vals) != len(colors):
		sys.exit('ERROR: len vals is %d and len colors is %d' % (len(vals),len(colors)))
	format_violin(ax1,parts,vals,labels,colors=colors)

	if not brelax:
		print('computing KW test')
		## compute kruskal-wallis test
		# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kruskal.html
		pairs = [[1,2],[1,3],[2,3]]
		shift = max([max(v) for v in vals])/15
		start_point = max([max(v) for v in vals])+shift
		for p1,p2 in pairs:
			statistic,pval = kruskal(vals[p1],vals[p2])
			#print(p1,p2,statistic,pval)
			line_x = [p1+1,p1+1,p2+1,p2+1]
			line_y = [start_point,start_point+shift,start_point+shift,start_point]
			start_point += shift*2
			if pval < 0.01:
				ax1.plot(line_x,line_y,'-k')
			else:
				ax1.plot(line_x,line_y,':k')
			if pval != 0:
				ax1.text((p1+p2+2)/2,start_point-shift*2,'$p=$%.1e' % (pval),horizontalalignment='center')
			else:
				ax1.text((p1+p2+2)/2,start_point-shift*2,'$p<$%.1e' % (sys.float_info.min),horizontalalignment='center')
		
	ax1.set_ylabel('Scores')
	#ax1.set_title('STRING\'s "%s" channel\n(%d interactions map to Reactome)\n' % (name,len(interactions)))
	ax1.set_title('STRING\'s "%s" channel' % (name),fontsize=14)
	print('done with the other subplot')

	plt.tight_layout()
	filename = '%s.png' % (prefix)
	plt.savefig(filename)
	print('wrote to '+filename)
	filename = filename.replace('png','pdf')
	plt.savefig(filename)
	print('wrote to '+filename)
	plt.close()
	
	return

def format_violin(ax,parts,vals,labels,colors=None,alpha=.2):
	for i in range(len(vals)):
		vals[i] = sorted(vals[i])
	i=0
	for pc in parts['bodies']:
		if not colors:
			pc.set_facecolor('#541684')
		else:
			pc.set_facecolor(colors[i])
		i+=1
		pc.set_edgecolor('black')
		pc.set_alpha(alpha)

	quartile1 = [0]*len(vals)
	medians = [0]*len(vals)
	quartile3 = [0]*len(vals)
	for i in range(len(vals)):
		quartile1[i], medians[i], quartile3[i] = np.percentile(vals[i], [25, 50, 75])	
	#quartile1, medians, quartile3 = np.percentile(vals, [25, 50, 75],axis=1)
	whiskers = np.array([
		adjacent_values(sorted_array, q1, q3)
		for sorted_array, q1, q3 in zip(vals, quartile1, quartile3)])
	whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]

	
	for i in range(len(medians)):
		inds = [(i+1)-.1,(i+1)+.1]
		meds = [medians[i],medians[i]]
		ax.plot(inds, meds, color='#ba0000', linestyle='-', lw=5,zorder=3)
		print('  median:',labels[i],medians[i])
	inds = np.arange(1, len(medians) + 1)
	ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
	ax.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)
	
	set_axis_style(ax, labels)
	return

# https://matplotlib.org/gallery/statistics/customized_violin.html
def set_axis_style(ax, labels):
	ax.get_xaxis().set_tick_params(direction='out')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(np.arange(1, len(labels) + 1))
	ax.set_xticklabels(labels)
	ax.set_xlim(0.25, len(labels) + 0.75)

def adjacent_values(vals, q1, q3):
	upper_adjacent_value = q3 + (q3 - q1) * 1.5
	upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

	lower_adjacent_value = q1 - (q3 - q1) * 1.5
	lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
	return lower_adjacent_value, upper_adjacent_value

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print('USAGE: python3 viz-channels.py <INFIX>')
	main(sys.argv[1])