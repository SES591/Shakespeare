########################################################################################################################################
###################################################     Build_Shakespeare Network     ##################################################
##################    Input Raw Movie data, output Networks   -    Tucker Ely, 24 Feb 2016    ##########################################
########################################################################################################################################
import os
import re
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import *
from collections import *
from compiler.ast import flatten


"""	EDIT BELOW	"""

NODE_FILE = '/Users/tuckerely/Google Drive/Ship/9_Complex_systems/My_Network_Attempt/character_list.dat'         # change address for file locations locally
EDGE_FILE_RAW = '/Users/tuckerely/Google Drive/Ship/9_Complex_systems/My_Network_Attempt/Raw_Film_Data.dat'         # change address for file locations locally

output_name = 'MADN_'

weight_treshold = 0 	#	Only look at weights above this value


### Triggers	0 = no, 1 = yes
compute_MI = 0
compute_TE = 1
print_graph = 0
graph_TE = 0
graph_random = 0
graph_casual_edges = 1



### Edge weights
w_g = 1
w_4 = 1
w_3 = 1
w_2 = 1
w_1 = 1
w_ext = 0.5


"""	EDIT ABOVE	"""


#####################    Main    ########################
def main():
	
	f = open(EDGE_FILE_RAW, 'r')
	lines = f.readlines()
	f.close()	


	scenes = count_acts(lines)
	characters = count_characters(NODE_FILE)					#	For calling the indicies of teh stat_grid
	
	### Save the weights here per scene
	stat_grid = np.zeros((len(characters) + 1, scenes + 1)) 	# 	(rows, columns) the +1's are for teh stats rows for MI measures


	### After this is done, i can simply copy the array with the stipulation that the weights be above some specific threshold.


	x = 0
	scene = 0


	### Create blank H for additions after scene 1
	A = build_nodes_list(NODE_FILE)
	H = A.copy() 	#	total_play copy

	global_edge_list = [] 														#	for all character pairs which appear
	character_pairs_global = list(permutations(range(len(characters)-1), 2))

	while x < len(lines):
		
		###	When a new act is encountered, that is not the first act, print the previous one. No += x	
		if re.findall('^Act', lines[x]) and scene > 0:
			### Draw Graph 
			pos=nx.circular_layout(G)
			nx.draw_networkx_nodes(G,pos,node_size=700)
			nx.draw_networkx_labels(G,pos,font_size=10,font_family='sans-serif')
			step_edges = [(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>0]
			weights = [G[u][v]['weight'] for u, v in step_edges]
			# node_size = [(u,v) for (u,v,d) in G.nodes(data=True) if d['weight']>0]
			nx.draw_networkx_edges(G, pos, edgelist=step_edges, width=weights)		#	Error here.		
			plt.axis('off')
			if print_graph == 1:
				plt.title(scene)
				plt.show('')				
				plt.savefig(output_name + '_' + str(scene) + '_.png')
				plt.close()

			### Append stat_grid
			stat = G.degree(weight='weight').items() 
			for item in stat:
				stat_grid[characters.index(item[0])][scene-1] = item[1]



		###	When the end of the file is encountered. print previous, DO += x
		if re.findall('^END', lines[x]): 
			### Draw Graph 
			pos=nx.circular_layout(G)

			nx.draw_networkx_nodes(G,pos,node_size=700)
			nx.draw_networkx_labels(G,pos,font_size=10,font_family='sans-serif')
			step_edges=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>0]
			weights = [G[u][v]['weight'] for u, v in step_edges]
			nx.draw_networkx_edges(G, pos, edgelist=step_edges, width=weights)		#	Error here.		
			plt.axis('off')
			if print_graph == 1:
				plt.title(scene)
				plt.show('')				
				plt.savefig(output_name + '_' + str(scene) + '_.png')
				plt.close()

			### Append stat_grid
			stat = G.degree(weight='weight').items() 
			for item in stat:
				stat_grid[characters.index(item[0])][scene-1] = item[1]



			### Call the node colors
			build_node_color(H)
			step_nodes = [(u) for (u,d) in H.nodes(data=True) if d['total_weight']>0]
			sizes = [(H.node[u]['total_weight']*10) for u in H.nodes()]


			plt.close()
			pos=nx.circular_layout(H)
			nx.draw_networkx_nodes(H,pos, node_size = sizes)
			# nx.draw_networkx_labels(H,pos,font_size=10,font_family='sans-serif')
			step_edges=[(u,v) for (u,v,d) in H.edges(data=True) if d['weight']>0]
			weights = [(H[u][v]['weight']*.5) for u, v in step_edges]
			nx.draw_networkx_edges(H, pos, edgelist=step_edges, width=weights)		#	Error here.		
			plt.axis('off')
			if print_graph == 1:
				plt.title(scene)
				plt.show('')				
				plt.savefig(output_name + '_All.png')
				plt.close()



			x += 1			
			break

		### Intra-act behavior. Building the weights into the current Act.
		scene += 1	

		G = A.copy()					#	build new blank graph with NODE FILE  (char. list), need copy here to preserve the order of the characters in teh circle
		x += 1
		
		while not re.findall('^Act', lines[x]) and not re.findall('^END', lines[x]): 								#	within scene, add weight to edges.
			if re.findall('^g', lines[x]):
				weight_assignment(lines, x, 'g', w_g, G, H, global_edge_list)
				x += 1
			if re.findall('^1', lines[x]):
				weight_assignment(lines, x, '1', w_1, G, H, global_edge_list)
				x += 1
			if re.findall('^2', lines[x]):
				weight_assignment(lines, x, '2', w_2, G, H, global_edge_list)
				x += 1
			if re.findall('^3', lines[x]):			
				weight_assignment(lines, x, '3', w_3, G, H, global_edge_list)
				x += 1		
			if re.findall('^4', lines[x]):
				weight_assignment(lines, x, '4', w_4, G, H, global_edge_list)
				x += 1



	if compute_MI == 1: 														#	Compute mutual info *****Abandonded as unrealistic for this project*****
		MI_data_table = get_MI_data_table(stat_grid, characters)

	if compute_TE == 1: 														#	Compute the transfer entropy
		TE_grid_global, one_ratio = TE_data_1(stat_grid, characters)

	if graph_TE == 1: 															#	TE(Rank order) plot for given k
		TE_plots = graph_TE_attributes(TE_grid_global, 'TE_3_thresh_0.png')

	if graph_random ==1:
		TE_grid_random = random_graph(characters, scenes, one_ratio) 			#	Build randomily generated TE for comparison
		
		# random_plots = graph_TE_random(TE_grid_random, 'TE_random_all.png')
		
	if graph_casual_edges ==1:
		causal_edges = build_causal_edges(global_edge_list, TE_grid_global, characters, character_pairs_global)

		# causal_edges_random = build_causal_edges(global_edge_list, TE_grid_global, characters, character_pairs_global)


########################################################
################    Functions     ######################
class weight_assignment(object):
 	kind = 'variable'
 	
 	def __init__(self, lines, x, group_size, w, G, H, global_edge_list): 		#	
		self.lines = lines
		self.x = x
		self.group_size = group_size
		self.w = w 
		self.G = G									
		self.global_edge_list = global_edge_list
		a = str(lines[x]) 													#	Clean up the data in line (below)
		b = re.split(r'\t+', a.rstrip('\t'))
		b[0] = b[0].replace(group_size,"")
		b[0] = b[0].replace("(","")
		b[0] = b[0].replace(")","")
		b[0] = b[0].replace("\n","")
		c = b[0].split(', ')
		d = list(combinations(c,2)) 										#	All pairwise interactions of cahacters in the line 
		if len(c) > 1:
			for e in d: 													#	add edge to local(scene) as well as global 
				G[e[0]][e[1]]['weight'] += w
				H[e[0]][e[1]]['weight'] += w 								
				global_edge_list.append((e[0],e[1]))

		else:																#	Deal with monologs
			G[c[0]][c[0]]['weight'] += w
			H[c[0]][c[0]]['weight'] += w
			global_edge_list.append((c[0],c[0]))

def build_nodes_list(NODE_FILE): 			#	Load Character list
    nodes_list = []
    for line in open(NODE_FILE, 'r').readlines():
        items = [x.strip() for x in line.rstrip().split('\t')]
        if line[0] == '#' or line=='':
            continue
        nodes_list.append(items[0])
        #some_attribute_list.append(items[1]) # this will let me put attributes in a second tabed column to identify characters with.
    G = nx.Graph()
    
    
    for item in nodes_list:
    	G.add_node(str(item), total_weight='total_weight')

	a = list(combinations_with_replacement(nodes_list, 2)) 	#	Non-directional(orde doesnt matter)
    for e in a:
    	G.add_edge(e[0] ,e[1] ,weight=0) 					#	I set up all potential edges with weight 0 here, so later on i can update my edge weights without worrying if they are present yet
    return G


def build_node_color(H):
	### sum the total incoming weigts of the total_play copy, set as color
	
	y = list(product(range(len(H.nodes())), range(len(H.nodes()))))


	x = 0
	while x < len(y):
		if H.nodes()[y[x][0]] == H.degree(weight='weight').items()[y[x][1]][0]:
			H.node[H.nodes()[y[x][0]]]['total_weight'] = H.degree(weight='weight').items()[y[x][1]][1]
			x += 1
		else:
			x += 1
	return H


def count_acts(lines):
	x = 0
	scenes = 0
	while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
		if re.findall('^Act', lines[x]): 
			scenes += 1
			x += 1
		else:
			x += 1
	return scenes


def count_characters(NODE_FILE):
	f = open(NODE_FILE, 'r')
	lines = f.readlines()
	f.close()
	characters = [line.rstrip() for line in lines]
	return characters


def get_MI_data_table(stat_grid, characters):
	stat_grid_threshold = np.copy(stat_grid)
	for row in stat_grid_threshold:
		for (i, item) in enumerate(row):
			if item <= weight_treshold:
			 	row[i] = 0
			else:
				row[i] = 1

	### Scenes summed to bottom (add final row) This is the (global value) P(x,y,, . . )
	trans_array = np.transpose(stat_grid_threshold)
	for row in trans_array: 		
		a = [str(item) for item in row]
		b = [item.rstrip('0').rstrip('.') for item in a]
		c = b[:-1]
		d = "".join(c)
		f = int(str(d),2)
		row[-1] = f	
	stat_grid_threshold = np.transpose(trans_array)
	
	### characters summed to right (fill in P col)this is the P-value for the individual characters (P(x), P(y) ect.
	for row in stat_grid_threshold: 
		a = float(sum(row)) / float(len(row))		
		row[-1] = a
	stat_grid_threshold.astype(int)

	return stat_grid_threshold


def TE_data_1(stat_grid, characters): 	### My first attmpt at solving the various p-values

	### make binary character grid
	stat_grid_threshold = np.copy(stat_grid)  								#	Binary grid
	for row in stat_grid_threshold:
		for (i, item) in enumerate(row):
			if item <= weight_treshold:
			 	row[i] = 0
			else:
				row[i] = 1


	one_ratio = float(np.count_nonzero(stat_grid_threshold)) / float((len(stat_grid_threshold)-1)*(len(np.transpose(stat_grid_threshold))-1))



	character_pairs = list(permutations(range(len(stat_grid_threshold)-2), 2))	

	k_choices = 8															#	I am stacking the k_choices into the length of the TE_grid, this affects my c count
	
	TE_grid_global = [] 													#	All combinad k grids (TE_grid_local) 

	k = 3
	while k <= 3:#k_choices: 													#	iterate over all k

		TE_grid_local = np.zeros((len(character_pairs), 4)) 				# 	(rows, columns)
		
		c = 0
		while c < len(character_pairs): 									# 	iterate over all character pairs			
			local_TE = 0
			n_steps = (len(np.transpose(stat_grid_threshold))-1) - (k+1) +1
			step_size = 1
			r = []
			
			n = 1
			while n <= n_steps: 												
				
				# x_n = (n*step_size - 1)					
				x_n = (k+n) -1 							#	(Xn+1)  -1 for index of list adjustment
				y_n = x_n-1								#	(Yn)
				k_local = range(x_n - k, x_n) 			#	(Xn^k) indicies of stat_grid_threshold

				### Corretly apply character indices to grid
				char1 = character_pairs[c][0]
				char2 = character_pairs[c][1]
				all_combo_number = 2**k

				### r_1 = all combos. 
				r.append(flatten([list(stat_grid_threshold[char1][k_local]), stat_grid_threshold[char1][x_n], stat_grid_threshold[char2][y_n]]))
				r_1 = []
				for i in r:
					f = [str(item) for item in i]
					g = [item.replace('.0','') for item in f]
					h = ''.join(g)
					r_1.append(h)

				n += 1
			
			
			set_r_1 = list(set(r_1)) 								#	find unique states
			i = 0
			while i < len(set_r_1):									#	find probabilities needed in TI
				
				### r_p = p(Xnk, Xn+1, Yn)
				r_check = [item for item in r_1 if item == set_r_1[i]]
				r_count = len(r_check)
				r_p = float(r_count) / float(n_steps)
				
				### s_p = p(Xn+1 | Xnk, Yn)
				s_exp = len([item for item in r_1 if item[:k] == set_r_1[i][:k] and item[-1] == set_r_1[i][-1]]) # this divisor is limited to theose that met the condition
				s_p = float(r_count) / float(s_exp) 

				### t_p = p(Xn+1 | Xnk)
				t_check = [item for item in r_1 if item[:(k+1)] == set_r_1[i][:(k+1)]]
				t_count = len(t_check)
				t_exp = len([item for item in r_1 if item[:k] == set_r_1[i][:k]])
				t_p = float(t_count) / float(t_exp)

				local_TE = local_TE + (r_p*np.log2(s_p / t_p)) 		#	sum local TE of all combos

				i += 1



			TE_grid_local[c][0] = k
			TE_grid_local[c][1] = character_pairs[c][0] 									#	character 1 id
			TE_grid_local[c][2] = character_pairs[c][1] 									#	character 2 id
			TE_grid_local[c][3] = local_TE 													# 	One value for every character pair(c)
		
			c += 1
		TE_grid_local = TE_grid_local[TE_grid_local[:,3].argsort()[::-1]]							#	sort grid by 4th col. TE info.

		TE_grid_global.append(TE_grid_local)
		
		name = 'TE_grid_' + str(k) + '.csv'
		np.savetxt(name, TE_grid_local, delimiter=',')

		k += 1

	### grab one_ratio
	# one_ratio = float(np.count_nonzero(TE_grid_global)) / float((len(TE_grid_local)-1)*(len(np.transpose(TE_grid_local))-1))

	return TE_grid_global, one_ratio


def graph_TE_attributes(TE_grid_global, output_name):
	
	plt.close()
	for item in TE_grid_global:
		x = range(len(item))
		y = list(np.transpose(item)[-1])
		plt.plot(x,y)
	
	plt.xlabel('Rank')
	plt.ylabel('TE (bits)')

	# plt.text(20,0.005, 'k = ' + str(item[0][0]))
	plt.title('Much Ado About Nothing')
	plt.show('')
	# output_name = 'TE_graph_all.png' 
	# output_name = 'TE_graph_' + str(item[0][0])	 + '.png'			
	plt.savefig(output_name)
	plt.close()


def random_graph(characters, scenes, one_ratio):			#	 build random graphs a bunch of times, then average them.
	
	### for a random graph based ont eh correct one_ratio
	random_graph_out = np.zeros((len(characters) + 1, scenes + 1))
	y = 0
	while y < len(random_graph_out)-1: 	# iterate over rows
		x = 0
		while x < len(np.transpose(random_graph_out)-1):	#	iterate within rows
			dice_roll = np.random.random()
			if dice_roll < one_ratio:
				random_graph_out[y][x] = 1.0			
				x += 1
			else:
				random_graph_out[y][x] = 0.0
				x += 1
		y += 1		
	

	character_pairs = list(permutations(range(len(random_graph_out)-2), 2))	

	k_choices = 8															#	I am stacking the k_choices into the length of the TE_grid, this affects my c count
	
	TE_grid_global_random = [] 													#	All combinad k grids (TE_grid_local) 
	iterations = 100
	# iter_grid = np.array([[]]*iterations)
	iter_grid = np.zeros((iterations, len(character_pairs)))

	e = 0
	k = 3
	while e < iterations:
		### for a random graph based ont eh correct one_ratio
		random_graph_out = np.zeros((len(characters) + 1, scenes + 1))
		y = 0
		while y < len(random_graph_out)-1: 	# iterate over rows
			x = 0
			while x < len(np.transpose(random_graph_out)-1):	#	iterate within rows
				dice_roll = np.random.random()
				if dice_roll < one_ratio:
					random_graph_out[y][x] = 1.0			
					x += 1
				else:
					random_graph_out[y][x] = 0.0
					x += 1
			y += 1	



		TE_grid_local = [] #np.zeros((len(character_pairs), 1)) 				# 	(rows, columns)
		
		c = 0
		while c < len(character_pairs): 									# 	iterate over all character pairs			
			local_TE = 0
			n_steps = (len(np.transpose(random_graph_out))-1) - (k+1) +1
			step_size = 1
			r = []
			
			n = 1
			while n <= n_steps: 												
				
				# x_n = (n*step_size - 1)					
				x_n = (k+n) -1 							#	(Xn+1)  -1 for index of list adjustment
				y_n = x_n-1								#	(Yn)
				k_local = range(x_n - k, x_n) 			#	(Xn^k) indicies of random_graph_out

				### Corretly apply character indices to grid
				char1 = character_pairs[c][0]
				char2 = character_pairs[c][1]
				all_combo_number = 2**k

				### r_1 = all combos. 
				r.append(flatten([list(random_graph_out[char1][k_local]), random_graph_out[char1][x_n], random_graph_out[char2][y_n]]))
				r_1 = []
				for i in r:
					f = [str(item) for item in i]
					g = [item.replace('.0','') for item in f]
					h = ''.join(g)
					r_1.append(h)

				n += 1
			
			
			set_r_1 = list(set(r_1)) 								#	find unique states
			i = 0
			while i < len(set_r_1):									#	find probabilities needed in TI
				
				### r_p = p(Xnk, Xn+1, Yn)
				r_check = [item for item in r_1 if item == set_r_1[i]]
				r_count = len(r_check)
				r_p = float(r_count) / float(n_steps)
				
				### s_p = p(Xn+1 | Xnk, Yn)
				s_exp = len([item for item in r_1 if item[:k] == set_r_1[i][:k] and item[-1] == set_r_1[i][-1]]) # this divisor is limited to theose that met the condition
				s_p = float(r_count) / float(s_exp) 

				### t_p = p(Xn+1 | Xnk)
				t_check = [item for item in r_1 if item[:(k+1)] == set_r_1[i][:(k+1)]]
				t_count = len(t_check)
				t_exp = len([item for item in r_1 if item[:k] == set_r_1[i][:k]])
				t_p = float(t_count) / float(t_exp)

				local_TE = local_TE + (r_p*np.log2(s_p / t_p)) 		#	sum local TE of all combos

				i += 1

			TE_grid_local.append(local_TE)													# 	One value for every character pair(c)
		
			c += 1
		
		TE_grid_local.sort(reverse=True)

		iter_grid[e] = TE_grid_local
		
		e += 1

	random_mean = np.mean(iter_grid, axis = 0)


	random_std = np.std(iter_grid, axis = 0)
	print random_std
	random_high = random_mean + random_std 
	random_low = random_mean - random_std

	name = 'TE_random_mean_' + str(k) + '.csv'
	np.savetxt(name, random_mean, delimiter=',')

	name = 'TE_random_sd_high_' + str(k) + '.csv'
	np.savetxt(name, random_high, delimiter=',')

	name = 'TE_random_sd_low_' + str(k) + '.csv'
	np.savetxt(name, random_low, delimiter=',')


def build_causal_edges(global_edge_list, TE_grid_global, characters, character_pairs_global):

	x = 0
	while x < len(global_edge_list): 					#	Convert global edge list to indicies
		global_edge_list[x] = (characters.index(global_edge_list[x][0]),characters.index(global_edge_list[x][1]))
		# global_edge_list[x][1] = characters.index(global_edge_list[x][1])
		x += 1


	global_edge_list = set(global_edge_list) 			#	Get rid of duplicates in global edge list

	### 	make TE > 0 List
	global_TE_list = []				
	x = 0
	while x < len(TE_grid_global[0]):
		if TE_grid_global[0][x][3] > 0:
			global_TE_list.append((TE_grid_global[0][x][1],TE_grid_global[0][x][2]))
		x += 1


	### has edge, has TE > 0
	# in global_edge_list, value in column 4 of TE_grid_global > 0:
	edge_with_TE = 0
	for item in global_edge_list:
		if item in global_TE_list:
			edge_with_TE += 1

	### has edge, has TE = 0
	# in global_edge_list, value in column 4 of TE_grid_global == 0:
	edge_without_TE = 0
	for item in global_edge_list:
		if item not in global_TE_list:
			edge_without_TE += 1

	### no edge, has TE > 0
	# not in global_edge_list, value in column 4 of TE_grid_global > 0:
	no_edge_with_TE = 0
	for item in global_TE_list:
		if item not in global_edge_list:
			no_edge_with_TE += 1

	### no edge, has TE = 0
	# not in global_edge_list, value in column 4 of TE_grid_global == 0:
	no_edge_without_TE = float(len(character_pairs_global)) - (edge_with_TE + edge_without_TE + no_edge_with_TE)

	p_edge_with_TE = float(edge_with_TE) / float(len(character_pairs_global))
	p_edge_without_TE = float(edge_without_TE) / float(len(character_pairs_global))
	p_no_edge_with_TE = float(no_edge_with_TE) / float(len(character_pairs_global))
	p_no_edge_without_TE = float(no_edge_without_TE) / float(len(character_pairs_global))

	print p_edge_with_TE
	print p_edge_without_TE
	print p_no_edge_with_TE
	print p_no_edge_without_TE

###
# def TE_causal_edge_graph():

# >>> import numpy
# >>> a = numpy.random.random(siaze=100) * 100 
# >>> numpy.histogram(a, bins=(0.0, 7.3, 22.4, 55.5, 77, 79, 98, 100))
# (array([ 8, 14, 34, 31,  0, 12,  1]), 
#  array([   0. ,    7.3,   22.4,   55.5,   77. ,   79. ,   98. ,  100. ]))


#################################################################
#################################################################


if __name__=='__main__':
    main()








####################################################################################################################
###############      				Visualizations left to draft into model    					  ##################
# Change color of token based on some attribute gotten from second column of Node File, such as character attributes.
# Add connections for subject matter of conversations
# Add title which includes Act_Name
# Add a running weighted sum to token which considers the total weight from all scences, put this into the color.
# incorporate word count info somewhere
# learn to combine multiple graphs on top of each other.



####################################################################################################################
###############      APPLY IMFORMATION THEORY ! ! !      FOCUS ON ONE INFORMATION EXPERIMENT       #################
#  Potential experiments:
#  (1) Can I tease out a de-facto control Kernal?
#  (2) WHO processes the most information, and how do i Measure this?
#  (3) Can i look at how much individual pairwise interactions predict future network behavior?
#  (4) which edges are corrolated? Simple starting point maybe.
#  (5) could i corrolate my architect with how far ahead/behind the audience is? (Get Phil to map this)
#  (6) Sara suggested point-wise information measures might fit my system better.
#  (7) where is info stroage relative to processing
#  (8) What are the causal structurs? (What does that even mean?)
#  (9) Instead of focusing on edges, which if i focus on states, what are the states of my characters? can I say that they are happy or sad
#  (10) What is active information anyway?
#  (11) Change point detection (cole has software for this), think ENRON Eamil network
#  (12) I can think of % of time a character had an edge, which is kind of the same as summing the wieghts
#
# (13) Generated P-distros of characters, how often a character is talked about? ie, if one character is talked about, how well does that predict others.
# (14) I can simply look at how the sum of interactions predicts future interactions.  Sum-edges()
#
#  Integrated info, look at partitions of network,a nd ask how much info they have related.Treat the network as a whole as a state.
# What is the info propagation inthe network over time. The whole network evolves diferent than the inidvdual nodes, This is the observer (whole network)
# 
# precentage of calusal (edge). generate all pairwise interactions, Average TE(% causal edges)


####################################################################################################################
###############     					 General thoughts on the play      						   #################
# What is the goal of a play (maximize conflict resolution int he final act), <-- ask Phil about this
# how much information do I gain about the observe state from pairwise interactions.
# could coarse grain by number of participants . .  ie the observe gains more info from 3-ways than 2-ways.
# 




