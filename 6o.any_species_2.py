########################################################################################################################################
###################################################     6oany_speces_grab    ################################################################
##################    to grab 6o output phases as needed    Tucker Ely, 19 October 2015    ###########################################
########################################################################################################################################
# Currently works for any aq and s species. Gas not yet built in.
# Erros generated form 3-digit exponents will be dealt with using a try: statment on the float conversion
# 
# Need to add aq_element total capability (this is done in the search file)
# Need to add trc capability.
#



import re, sys, os, csv
from collections import defaultdict
import fileinput
import shutil
import numpy as np
import pandas as pd
from string import digits
from operator import itemgetter
from collections import OrderedDict
import itertools


""" EDIT BELOW """
run  = 'PAL_350_1-1'
directory = '/Users/tuckerely/Google Drive/Ship/5_OPERATIONS/3_Gale_Seg/3_Analyze_6o'
solid_sp_template = '/Users/tuckerely/Google Drive/Ship/5_OPERATIONS/3_Gale_Seg/wtb_solid_stochiometry.csv'
info_file = '/Users/tuckerely/Google Drive/Ship/5_OPERATIONS/3_Gale_Seg/3_Analyze_6o/Dalton_Segment_catalog_total.csv'



print_final_output = 1 				# CSV desired ?   1 = yes, 0 = no

# For listed names: metacharacters: if this are in a name, they must be precedded with a \   . ^ $ * + ? { } [ ] \ | ( )
# Remember:  + is a metacharacter and  - is not

aq_grab = 1 						# 1 = yes   0 = no
aq_sp = ['SO4-2','CO,aq','CO2,aq','HCO3-','CO3-2','O2,aq','HS-','H2S,aq','H2,aq','METHANE,AQ','FORMATE,AQ','ACETATE,AQ','MTHEANOL,AQ','ACETIC-ACID,AQ','FORMALDEHYDE,AQ','ACETALDEHYDE,AQ','ETHANE,AQ','SiO2,aq','HPO4-2','Fe\+2','Fe\+3','Mn\+2','MnSO4,aq'] 	# List desired
aq_choice = 1 						# data type:	mols = 1  grams = 2  conc = 3  log conc = 4  log g = 5  log act = 6
aq_unit_print = 'molal'

aq_99_grab = 0 						# 1 = yes   0 = no
aq_99_of_basis = ['SO4-2'] 			# Only one basis species at a time
aq_99_of_basis_choice = 1 			# data type:	molal conc = 1   percent = 2
aq_99_of_basis_unit_print = 'molal conc'

s_grab = 1 							# 1 = yes   0 = no
s_sp = ['PYRITE', 'QUARTZ','HEMATITE','EPIDOTE'] 		# List desired
s_choice  = 2 						#	data type: 		log mols = 1 	mols = 2 	grams = 3 	volume = 4
s_unit_print = 'mols'

g_grab = 0 							# 1 = yes   0 = no
g_sp = ['H2S,g', 'h2\(g\)','CH4,g','CO2,g','CO,g','o2\(g\)','S2,g','SO2,g'] 		# List desired
g_choice = 1 						# data type: 	log f = 1 	f = 2
g_unit_print = 'log(f)'



"""Solid Element Grab"""
s_element_grab = 1 					# 1 = yes   0 = no

s_element_unit_print = 'mols'
rnt_form = 1 						# 0 = mineral, 1 = oxide, 2 = aq phase
rock_name = 'basalt'				# reactant name if oxide is selected


grab_fe = 1
grab_mg = 0
grab_si = 0
grab_al = 0
grab_ca = 0
grab_na = 0
grab_k = 0
grab_ti = 0
grab_mn = 0
grab_s_6 = 0
grab_s_neg2 = 0
grab_s_neg1 = 0


""" EDIT ABOVE """
#fitting coefficients, eq# denotes order

# 2-350K and 500bar
neutral_ph_eq0 = 7.34334 
neutral_ph_eq1 = - 0.0193815
neutral_ph_eq2 = 0.0000796789 
neutral_ph_eq3 = - 1.75343E-7
neutral_ph_eq4 = 1.83503E-10

fmq_eq0 = -87.9119 
fmq_eq1 = 0.335374
fmq_eq2 = - 0.000959826
fmq_eq3 = 1.84084E-6
fmq_eq4 = - 1.57724E-9

#############################################################################################################
##########################################
"""Files"""
file_name = []
file_list = []
for root, dirs, files in os.walk(directory):  # this works. 
    for file in files:
        if file.endswith(".6o"):
        	file_name.append(file)
        	file_list.append(os.path.join(root, file))


##########################################
"""Functions"""

def grab_t():
	t = []
	x = 0
	while x < len(lines): 		# builds t and zi in master file,   per .6o
		if re.findall('^ {21}temperature    =', lines[x]):
			a = str(lines[x])
			b = a.split()
			t.append(float(b[2]))
			x += 1
		else:
			x += 1      
	return t;

def grab_fo2():
	fo2 = []
	x = 0
	while x < len(lines):
		if re.findall('^ {10,10}log oxygen fugacity = ', lines[x]): 
			a = str(lines[x])
			b = a.split()
			fo2.append(float(b[4])) 	# 	log value 	
			x += 1
		else:
			x += 1  

	fmq = [(fmq_eq0 + fmq_eq1*i + fmq_eq2*i**2 + fmq_eq3*i**3 + fmq_eq4*i**4) for i in  t]
	modified_fo2 = [0]*len(fmq)
	x = 0
	while x < len(fmq):
		modified_fo2[x] = fo2[x] - fmq[x] 			#	This will mean that - values are more acidic thn neutral, and + more basic
		x += 1 

	return modified_fo2, fmq, fo2;

def grab_ph_eh_pe():
	ph = []
	eh = []
	pe = []
	x = 0
	while x < len(lines): #	Modified nbs ph scale
		if re.findall('^ {5,5}modified nbs', lines[x]): 
			a = str(lines[x])
			b = a.split()
			ph.append(float(b[4])) 
			eh.append(float(b[5]))
			pe.append(float(b[6]))
			x += 1
		else:
			x += 1
	
	# incorporat neutral ph change. modified then becomes the pH as if 7 were always neutral
	neutral_ph = [(neutral_ph_eq0 + neutral_ph_eq1*i + neutral_ph_eq2*i**2 + neutral_ph_eq3*i**3 + neutral_ph_eq4*i**4) for i in  t]
	modified_ph = [0]*len(neutral_ph)
	x = 0
	while x < len(neutral_ph):
		modified_ph[x] = ph[x] - neutral_ph[x] 			#	This will mean that - values are more acidic thn neutral, and + more basic
		x += 1 



	return ph, modified_ph, neutral_ph, eh, pe;

def grab_aq():
	zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
	z = 0
	aq_array = np.zeros((len(aq_sp), len(t))) 					#	Array to contain all Fe_2 output values
	while z < len(aq_sp): 											#	Loop through all Fe contian members of fe_2 list		
		zeros = [0]*len(t)
		summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
		name = '^ {3}'  +  aq_sp[z]  +  ' {'  +  str(21-len(aq_sp[z]))  +  '}'#	Name  search string
		count = -1 														#	Must beign at -1, since the first count must move the list to position 0
		x = 0
		while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
			if re.findall('^ {21}temperature    =', lines[x]): 
				count += 1
				x += 1
			elif re.findall('^ {3}species {16}moles {8}grams {10}conc {8}log conc {6}log g {6}log act', lines[x]): 
				summary_of_prduct_trigger += 1
				x += 1
			elif re.findall('^ {11}--- major aqueous species contributing to mass balances ---', lines[x]): 
				summary_of_prduct_trigger -= 1
				x += 1
			elif re.findall(name, lines[x]) and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
				a = str(lines[x])
				b = a.split()
				try:
					zeros[count] = float(b[aq_choice])
					x += 1
				except:
					zeros[count] = 0
					x += 1
			else:
				x += 1
		aq_array[z] = zeros 										#	Add to build array,
		z += 1

	return aq_array

def grab_aq_99():
	aq_99_of_basis_expanded_list = [] 									#	Build species list for basis chosen
	summary_of_prduct_trigger = 0 										# 	This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
	name = '^ {3}\D.{21}\d'
	count = -1 															#	Must beign at -1, since the first count must move the list to position 0
	x = 0
	while x < len(lines): 												#	find fe_2_sp. and put it in appropriate step (count) position.
		if re.findall('^ {21}temperature    =', lines[x]): 
			count += 1
			x += 1                                                                                                                                     
		elif re.findall('^ aqueous species accounting for 99% or more of '  +  aq_99_of_basis[0], lines[x]): 
			summary_of_prduct_trigger += 1
			x += 1
		elif re.findall('^ aqueous species accounting for 99% or more of', lines[x]) and summary_of_prduct_trigger == 1:  # terminate at next line like this
			summary_of_prduct_trigger -= 1
			x += 1
		elif re.findall(name, lines[x]) and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
			a = str(lines[x])
			b = a.split()
			if b[0] not in aq_99_of_basis_expanded_list:
				aq_99_of_basis_expanded_list.extend([b[0]])
				x += 1
			else:
				x += 1
		else:
			x += 1

	# need to make aq_99_of_basis_expanded_list compatible with    . ^ $ * + ? { } [ ] \ | ( )
	x = 0
	while x < len(aq_99_of_basis_expanded_list):
		a = aq_99_of_basis_expanded_list[x]
		aq_99_of_basis_expanded_list[x] = a.replace(".", "\.")
		aq_99_of_basis_expanded_list[x] = a.replace("(", "\(")
		aq_99_of_basis_expanded_list[x] = a.replace(")", "\)")
		aq_99_of_basis_expanded_list[x] = a.replace("+", "\+")
 		x += 1

	
	zeros = [0]*len(t) 														#	add a list of 0's for eventual use by some species.
	z = 0
	aq_99_array = np.zeros((len(aq_99_of_basis_expanded_list), len(t))) 	#	Array to contain all Fe_2 output values
	while z < len(aq_99_of_basis_expanded_list): 													#	Loop through all Fe contian members of fe_2 list		
		zeros = [0]*len(t)
		summary_of_prduct_trigger = 0 										# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
		name = '^ {3}'  +  aq_99_of_basis_expanded_list[z]  +  ' {'  +  str(22-len(aq_99_of_basis_expanded_list[z]))  +  '}'#	Name  search string
		count = -1 															#	Must beign at -1, since the first count must move the list to position 0
		x = 0
		while x < len(lines): 												#	find fe_2_sp. and put it in appropriate step (count) position.
			if re.findall('^ {21}temperature    =', lines[x]): 
				count += 1
				x += 1
			elif re.findall('^ aqueous species accounting for 99% or more of '  +  aq_99_of_basis[0], lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
			elif re.findall('^ aqueous species accounting for 99% or more of', lines[x]) and summary_of_prduct_trigger == 1:
					summary_of_prduct_trigger -= 1
					x += 1
			elif re.findall(name, lines[x]) and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
				a = str(lines[x])
				b = a.split()
				try:
					zeros[count] = float(b[aq_99_of_basis_choice])
					x += 1
				except:
					zeros[count] = 0
					x += 1
			else:
				x += 1
		aq_99_array[z] = zeros 										#	Add to build array,
		z += 1
	
	#indexed_df1 = pandas.DataFrame.set_index([aq_99_array])
	#print indexed_df1
	# for x in 
	# 	d = {aq_99_of_basis_expanded_list[x]:aq_99_array[x]}

	# df = pandas.DataFrame(data=d)
	# print df
	#aq_99_array = np.hstack((aq_99_of_basis_expanded_list, aq_99_array))
	return aq_99_array, aq_99_of_basis_expanded_list

def grab_s():
	zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
	z = 0
	s_array = np.zeros((len(s_sp), len(t))) 					#	Array to contain all Fe_2 output values
	while z < len(s_sp): 											#	Loop through all Fe contian members of fe_2 list		
		zeros = [0]*len(t)
		summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
		name = '^ {3}'  +  s_sp[z]  +  ' {'  +  str(25-len(s_sp[z]))  +  '}'#	Name  search string
		count = -1 														#	Must beign at -1, since the first count must move the list to position 0
		x = 0
		while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
			if re.findall('^ {21}temperature    =', lines[x]): 
				count += 1
				x += 1
			elif re.findall('^ {2}product {20}log moles {8}moles {8}grams {8}volume, cc', lines[x]): 
				summary_of_prduct_trigger += 1
				x += 1
			elif re.findall('^ {27}mass, grams {8}volume, cc', lines[x]): 
				summary_of_prduct_trigger -= 1
				x += 1
			elif re.findall(name, lines[x]) and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
				a = str(lines[x])
				b = a.split()
				try:
					zeros[count] = float(b[s_choice])
					x += 1
				except:
					zeros[count] = 0
					x += 1
			else:
				x += 1
		s_array[z] = zeros 										#	Add to build array,
		z += 1
	return s_array

def grab_g():
	zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
	z = 0
	g_array = np.zeros((len(g_sp), len(t))) 					
	while z < len(g_sp): 		
		summary_of_prduct_trigger = 0 							
		zeros = [0]*len(t)								
		name = '^ {3}'  +  g_sp[z]  +  ' {'  +  str(27-len(g_sp[z]))  +  '}'
		count = -1 														
		x = 0
		while x < len(lines): 												#	find fe_2_sp. and put it in appropriate step (count) position.
			if re.findall('^ {21}temperature    =', lines[x]): 
				count += 1
				x += 1
			elif re.findall('^ {11}--- summary of gas species ---', lines[x]): 
				summary_of_prduct_trigger += 1
				x += 1
			elif re.findall('^ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -', lines[x]) and summary_of_prduct_trigger == 1:
				summary_of_prduct_trigger -= 1
				x += 1
			elif re.findall(name, lines[x]) and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
				a = str(lines[x])
				b = a.split()
				try:
					zeros[count] = float(b[g_choice])
					x += 1
				except:
					zeros[count] = 0
					x += 1
			else:
				x += 1
		g_array[z] = zeros 										#	Add to build array,
		z += 1
	return g_array

def grab_solid_element_tot():
	####  Molecular weights od the oxides  ####
	sio2_mw = 60.0843
	tio2_mw = 79.87880
	al2o3_mw = 101.96128
	fe2o3_mw = 159.6922
	cr2o3_mw = 151.9904
	feo_mw = 71.84640
	mno_mw = 70.93025
	mgo_mw = 40.30440
	nio_mw = 74.68940
	cao_mw = 56.07740
	na2o_mw = 61.97894
	k2o_mw = 94.19
	p2o5_mw = 141.94424

	"""Determine sp. of interest from solid stochometry csv file"""

	with open(solid_sp_template, 'rU') as file:  # Load file
	    contents = csv.reader(file)
	    matrix = list()
	    for row in contents:
	    	matrix.append(row)
	fe_0_list = []
	fe_0_coe = []

	fe_2_list = []
	fe_2_coe = []

	fe_3_list = []
	fe_3_coe = []

	mg_2_list = []
	mg_2_coe = []

	si_4_list = []
	si_4_coe = []

	al_3_list = []
	al_3_coe = []

	ca_2_list = []
	ca_2_coe = []

	na_1_list = []
	na_1_coe = []

	k_1_list = []
	k_1_coe = []

	ti_4_list = []
	ti_4_coe = []

	mn_2_list = []
	mn_2_coe = []

	s_6_list = []
	s_6_coe = []

	s_neg2_list = []
	s_neg2_coe = []

	s_neg1_list = []
	s_neg1_coe = []

	x = 1
	while x < len(matrix[0]): 						#	search matrix for sp. and add to sp. lists
		if float(matrix[6][x]) != 0.0:		#	fe_3
			fe_3_list.extend([matrix[0][x]])
			fe_3_coe.extend([matrix[6][x]])
		if float(matrix[37][x]) != 0.0:
			fe_3_list.extend([matrix[0][x]])
			fe_3_coe.extend([matrix[37][x]])

		if float(matrix[7][x]) != 0.0:		#	fe_2
			fe_2_list.extend([matrix[0][x]])
			fe_2_coe.extend([matrix[7][x]])
		if float(matrix[36][x]) != 0.0:
			fe_2_list.extend([matrix[0][x]])
			fe_2_coe.extend([matrix[36][x]])

		if float(matrix[35][x]) != 0.0:		#	fe_0
			fe_0_list.extend([matrix[0][x]])
			fe_0_coe.extend([matrix[35][x]])

		if float(matrix[9][x]) != 0.0: 		#	mg_4
			mg_2_list.extend([matrix[0][x]])
			mg_2_coe.extend([matrix[9][x]])

		if float(matrix[3][x]) != 0.0:  	#	si_4
			si_4_list.extend([matrix[0][x]])
			si_4_coe.extend([matrix[3][x]])

		if float(matrix[5][x]) != 0.0: 		#	al_3 
			al_3_list.extend([matrix[0][x]])
			al_3_coe.extend([matrix[5][x]])

		if float(matrix[11][x]) != 0.0: 		#	ca_2 
			ca_2_list.extend([matrix[0][x]])
			ca_2_coe.extend([matrix[11][x]])
		if float(matrix[33][x]) != 0.0: 		
			ca_2_list.extend([matrix[0][x]])
			ca_2_coe.extend([matrix[33][x]])

		if float(matrix[12][x]) != 0.0:  	#	na_1
			na_1_list.extend([matrix[0][x]])
			na_1_coe.extend([matrix[12][x]])
		if float(matrix[43][x]) != 0.0:  	
			na_1_list.extend([matrix[0][x]])
			na_1_coe.extend([matrix[43][x]])

		if float(matrix[13][x]) != 0.0: 		#	k_1 
			k_1_list.extend([matrix[0][x]])
			k_1_coe.extend([matrix[13][x]])
		if float(matrix[46][x]) != 0.0: 	
			k_1_list.extend([matrix[0][x]])
			k_1_coe.extend([matrix[46][x]])

		if float(matrix[4][x]) != 0.0: 		#	ti_4 
			ti_4_list.extend([matrix[0][x]])
			ti_4_coe.extend([matrix[4][x]])

		if float(matrix[8][x]) != 0.0:  	#	mn_2
			mn_2_list.extend([matrix[0][x]])
			mn_2_coe.extend([matrix[8][x]])
		if float(matrix[28][x]) != 0.0: 
			mn_2_list.extend([matrix[0][x]])
			mn_2_coe.extend([matrix[28][x]])

		if float(matrix[16][x]) != 0.0: 		#	s_6
			s_6_list.extend([matrix[0][x]])
			s_6_coe.extend([matrix[16][x]])

		if float(matrix[25][x]) != 0.0: 		#	s_neg2
			s_neg2_list.extend([matrix[0][x]])
			s_neg2_coe.extend([matrix[25][x]])

		if float(matrix[26][x]) != 0.0: 		#	s_neg1
			s_neg1_list.extend([matrix[0][x]])
			s_neg1_coe.extend([matrix[26][x]])

		x += 1



	x = 05
	t = []
	zi = []
	w_r = []
	rnt_mols = []
	h2o_mols = []
	while x < len(lines): 		# builds t and zi in master file,   per .6o
		if re.findall('^ {21}temperature    =', lines[x]):
			a = str(lines[x])
			b = a.split()
			t.append(float(b[2]))
			x += 1
		if re.findall('^ {21}reaction progress {8}=', lines[x]):
			a = str(lines[x])
			b = a.split()
			zi.append(float(b[3]))
			x += 1
		if re.findall('^ {20}moles of solvent h2o =', lines[x]):
			a = str(lines[x])
			b = a.split()
			h2o_mols.append(float(b[5]))
			count = +1
			x += 1
		else:
			x += 1


	#################################################
	# Build rnt list, and then find all of the delta mols values for all rnt's and sum them into a list call rnt_mols
	# this is used to determine w_r

	if rnt_form == 0: 			#	minerals
		rnt_list = [] 										#	build in all reactants
		x = 50 												#	Look only in the first 50 lines. For files with excessive notes at the this, this may need to be changed.
		while x < 200:	
			if re.findall('^ {3}reactant= ', lines[x]):
				a = str(lines[x])
				b = a.split()
				rnt_list.append(b[1])
				x += 1
			else:
				x += 1

		build_array = np.zeros((len(rnt_list), len(t)))		#	build in all reactants
		
		z = 0
		while z < len(rnt_list): 
			zeros = [0]*len(t)
			name = '^ {2}'  +  str(rnt_list[z])
			trigger = 0
			count = -1
			x = 0
			while x < len(lines): 
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				if re.findall(name, lines[x]) and trigger == 1: 		
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				elif re.findall('^ {5}reactant {18}moles {5}delta moles {6}grams {5}delta grams', lines[x]): 
					trigger += 1
					x += 1
				elif re.findall('^ {5}reactant {17}affinity   rel. rate', lines[x]): 
					trigger -= 1
					x += 1
				else:
					x += 1	
			build_array[z] = zeros 
			z += 1

		rnt_mols = build_array.sum(axis = 0) 	 	#	Sum the mols of all reactants repsent
		h2o_mols = np.divide(h2o_mols, 55.51) 		# 	I believe peter mentioning that the fluid phase is difined as one mol / kg
		rnt_mols[0] = rnt_mols[1]
		w_r = np.divide(h2o_mols, rnt_mols) 		

	if rnt_form == 1: 			#	Oxides
		rnt_list = []
		rnt_mol_sp_tot = [] 									#	build in all reactants
		x = 50 												#	Look only in the first 50 lines. For files with excessive notes at the this, this may need to be changed.
		while x < 200:	
			if re.findall('^ {4}\D. {11}\d.\d{15}E.\d{2}', lines[x]):
				a = str(lines[x])
				b = a.split()
				rnt_list.append(b[0])
				rnt_mol_sp_tot.append(float(b[1]))
				x += 1
			else:
				x += 1
		
		rnt_mol_sp_tot = np.asarray(rnt_mol_sp_tot) 		# 	Convert mol_sp_ totals to array so it can be summed
		rnt_mol_sum = rnt_mol_sp_tot.sum(axis = 0)
		rnt_mol_percent = np.divide(rnt_mol_sp_tot, float(rnt_mol_sum)) 


		build_array = np.zeros((len(rnt_list), len(t)))		#	build in all reactants
		z = 0
		while z < len(rnt_list): 
			zeros = [0]*len(t)
			name = '^ {2}'  +  str(rock_name)
			trigger = 0
			count = -1
			x = 0
			while x < len(lines): 
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				if re.findall(name, lines[x]) and trigger == 1: 		
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])*rnt_mol_percent[z] 		#	Mol_percent must be multiplied here, as the only reactant tracked at each step is the rock itself
					x += 1
				elif re.findall('^ {5}reactant {18}moles {5}delta moles {6}grams {5}delta grams', lines[x]): 
					trigger += 1
					x += 1
				elif re.findall('^ {5}reactant {17}affinity   rel. rate', lines[x]): 
					trigger -= 1
					x += 1
				else:
					x += 1	
			build_array[z] = zeros 
			z += 1
	

		rnt_mols = build_array.sum(axis = 0) 	 	#	Sum the mols of all reactants repsent
		h2o_mols = np.divide(h2o_mols, 55.51) 		# 	I believe peter mentioning that the fluid phase is difined as one mol / kg
		rnt_mols[0] = rnt_mols[1]
		w_r = np.divide(h2o_mols, rnt_mols) 


		#################################################
		# determine oxide weigth percent

		rnt_sio2_wt = float(rnt_mol_sp_tot[0] * sio2_mw)/10		
		rnt_tio2_wt = float(rnt_mol_sp_tot[7] * tio2_mw)/10
		rnt_al2o3_wt = float(rnt_mol_sp_tot[1] * al2o3_mw)/10
		rnt_feo_wt = float(rnt_mol_sp_tot[2] * feo_mw)/10
		rnt_mno_wt = float(rnt_mol_sp_tot[10] * mno_mw)/10
		rnt_mgo_wt = float(rnt_mol_sp_tot[3] * mgo_mw)/10
		rnt_cao_wt = float(rnt_mol_sp_tot[4] * cao_mw)/10
		rnt_na2o_wt = float(rnt_mol_sp_tot[5] * na2o_mw)/10
		rnt_k2o_wt = float(rnt_mol_sp_tot[6] * k2o_mw)/10
		rnt_p2o5_wt = float(rnt_mol_sp_tot[9] * p2o5_mw)/10

		#################################################
		# determine rnt mg_number

		rnt_mg_number = 100*(float(rnt_mol_sp_tot[3]) / float(rnt_mol_sp_tot[3] + rnt_mol_sp_tot[2]))


		#################################################
		# determine rnt Fe Conc 
		rnt_fe_conc = 0


	#build master output array and name list
	master_output_array = [0]*len(t)
	master_output_names = []

	#################################################
	# IRON
	if grab_fe == 1:


		#################################################
		# Find all Fe_2_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(fe_2_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(fe_2_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  fe_2_list[z]  +  ' {'  +  str(22-len(fe_2_list[z]))  +  '}' 	#	Name search string
			coe = [float(fe_2_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('^ {11}--- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		fe_2_total = build_array.sum(axis = 0) 	
		

		#################################################
		# Find all Fe_3_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(fe_3_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(fe_3_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  fe_3_list[z]  +  ' {'  +  str(22-len(fe_3_list[z]))  +  '}' 	#	Name search string
			coe = [float(fe_3_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		fe_3_total = build_array.sum(axis = 0) 	



		#################################################
		# Find all Fe_3_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(fe_0_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(fe_0_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  fe_0_list[z]  +  ' {'  +  str(22-len(fe_0_list[z]))  +  '}' 	#	Name search string
			coe = [float(fe_0_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		fe_0_total = build_array.sum(axis = 0) 	


		#################################################
		# Generate Ceff
		
		build_array = np.zeros((3, len(t))) 

		build_array[0] = fe_0_total
		build_array[1] = fe_2_total
		build_array[2] = fe_3_total

		tot_fe = build_array.sum(axis = 0)
		old_err_state = np.seterr(divide='ignore')
		ceff = np.divide(fe_3_total, tot_fe)

		# Add Fe sp. to master output array
		master_output_array = np.vstack((master_output_array, fe_2_total, fe_3_total, fe_0_total, ceff))
		master_output_names.extend(["Fe\+2","Fe\+3","Fe","Fe\+3/Fe-tot"])


	#################################################
	# MAGNESIUM 
	if grab_mg == 1: 


		#################################################
		# Find all Mg_0_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(mg_2_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(mg_2_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  mg_2_list[z]  +  ' {'  +  str(22-len(mg_2_list[z]))  +  '}' 	#	Name search string
			coe = [float(mg_2_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		mg_2_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, mg_2_total))
		master_output_names.extend(["Mg"])

	#################################################
	# SILICON 
	if grab_si == 1: 


		#################################################
		# Find all si_4_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(si_4_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(si_4_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  si_4_list[z]  +  ' {'  +  str(22-len(si_4_list[z]))  +  '}' 	#	Name search string
			coe = [float(si_4_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		si_4_total = build_array.sum(axis = 0) 	
		master_output_array = np.vstack((master_output_array, si_4_total))
		master_output_names.extend(["Si"])

	#################################################
	# ALUMINIUM 
	if grab_al == 1: 


		#################################################
		# Find all al_3_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(al_3_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(al_3_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  al_3_list[z]  +  ' {'  +  str(22-len(al_3_list[z]))  +  '}' 	#	Name search string
			coe = [float(al_3_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		al_3_total = build_array.sum(axis = 0) 	
		
		master_output_array = np.vstack((master_output_array, al_3_total))
		master_output_names.extend(["Al"])


	#################################################
	# CALCIUM
	if grab_ca == 1: 


		#################################################
		# Find all ca_2_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(ca_2_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(ca_2_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  ca_2_list[z]  +  ' {'  +  str(22-len(ca_2_list[z]))  +  '}' 	#	Name search string
			coe = [float(ca_2_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		ca_2_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, ca_2_total))
		master_output_names.extend(["Ca"])

	#################################################
	# SODIUM
	if grab_na == 1: 


		#################################################
		# Find all na_1_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(na_1_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(na_1_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  na_1_list[z]  +  ' {'  +  str(22-len(na_1_list[z]))  +  '}' 	#	Name search string
			coe = [float(na_1_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		na_1_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, na_1_total))
		master_output_names.extend(["Na"])

	#################################################
	# POTASSIUM
	if grab_k == 1: 


		#################################################
		# Find all k_4_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(k_1_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(k_1_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  k_1_list[z]  +  ' {'  +  str(22-len(k_1_list[z]))  +  '}' 	#	Name search string
			coe = [float(k_1_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		k_1_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, k_1_total))
		master_output_names.extend(["K"])

	#################################################
	# TITANIUM
	if grab_ti == 1: 


		#################################################
		# Find all ti_4_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(ti_4_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(ti_4_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  ti_4_list[z]  +  ' {'  +  str(22-len(ti_4_list[z]))  +  '}' 	#	Name search string
			coe = [float(ti_4_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		ti_4_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, ti_4_total))
		master_output_names.extend(["Ti"])

	#################################################
	# MANGANESE
	if grab_mn == 1: 


		#################################################
		# Find all mn_2_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(mn_2_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(mn_2_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  mn_2_list[z]  +  ' {'  +  str(22-len(mn_2_list[z]))  +  '}' 	#	Name search string
			coe = [float(mn_2_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		mn_2_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, mn_2_total))
		master_output_names.extend(["Mn"])

	#################################################
	# SULFUR 6
	if grab_s_6 == 1: 


		#################################################
		# Find all s_6_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(s_6_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(s_6_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  s_6_list[z]  +  ' {'  +  str(22-len(s_6_list[z]))  +  '}' 	#	Name search string
			coe = [float(s_6_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		s_6_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, s_6_total))
		master_output_names.extend(["S6"])

	#################################################
	# SULFUR NEG2
	if grab_s_neg2 == 1: 


		#################################################
		# Find all s_neg2_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(s_neg2_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(s_neg2_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  s_neg2_list[z]  +  ' {'  +  str(22-len(s_neg2_list[z]))  +  '}' 	#	Name search string
			coe = [float(s_neg2_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		s_neg2_total = build_array.sum(axis = 0) 	

		master_output_array = np.vstack((master_output_array, s_neg2_total))
		master_output_names.extend(["S-2"])

	#################################################
	# SULFUR NEG1
	if grab_s_neg1 == 1: 


		#################################################
		# Find all s_neg1_ list 
		zeros = [0]*len(t) 	#	add a list of 0's for eventual use by some species.
		z = 0
		build_array = np.zeros((len(s_neg1_list), len(t))) 					#	Array to contain all Fe_2 output values
		while z < len(s_neg1_list): 											#	Loop through all Fe contian members of fe_2 list		
			zeros = [0]*len(t)
			summary_of_prduct_trigger = 0 									# This allows me to "turn on" the retrieval of info between the solid product summary and the grand summary
			name = '^ {3}'  +  s_neg1_list[z]  +  ' {'  +  str(22-len(s_neg1_list[z]))  +  '}' 	#	Name search string
			coe = [float(s_neg1_coe[z])]*len(t)								#	determin minerl --> fe mols conversion CoE
			count = -1 														#	Must beign at -1, since the first count must move the list to position 0
			x = 0
			while x < len(lines): 											#	find fe_2_sp. and put it in appropriate step (count) position.
				if re.findall('^ {21}temperature    =', lines[x]): 
					count += 1
					x += 1
				elif re.findall('--- summary of solid product phases---', lines[x]): 
					summary_of_prduct_trigger += 1
					x += 1
				elif re.findall('           --- mineral saturation state summary ---', lines[x]): 
					summary_of_prduct_trigger -= 1
					x += 1
				elif re.findall(name, lines[x]) and zeros[count] == 0 and summary_of_prduct_trigger == 1: 		#	I dont technically need ==0 here, as the second occurace of annite will cover the first, given that the count doesnt change between them, but this is more correct, so maybe it will avoid errors later.
					a = str(lines[x])
					b = a.split()
					zeros[count] = float(b[2])
					x += 1
				else:
					x += 1
			build_array[z] = zeros 										#	Add to build array,
			build_array[z] = build_array[z]*coe							#	convert to mols Fe
			z += 1
		s_neg1_total = build_array.sum(axis = 0) 

		master_output_array = np.vstack((master_output_array, s_neg1_total))
		master_output_names.extend(["S-1"])



	master_output_array = np.delete(master_output_array, (0), axis=0)
	return master_output_array, master_output_names, w_r 

def grab_seg_value():
	extra_data = []
	extra_data_names = []
	with open(info_file, 'rU') as csvfile:
		csvfile.seek(0)
		reader = csv.reader(csvfile)
		for row in reader:
			name = str(file_name[w][0:-3])
			if row[0] == name:
				a = row[1::]
				x = 0
				while x < len(a):
					extra_data.append(float(a[x]))
					x += 1
			if row[0] == 'SegmentName':
				a = row[1::]
				x = 0
				while x < len(a):
					extra_data_names.append(a[x])
					x += 1

	extra_array = np.zeros((len(t), len(extra_data)))   #	build array for of identicle values for all steps
	z = 0
	while z < len(t):
		x = 0
		while x < len(extra_data):
			extra_array[z][x] = extra_data[x]
			x += 1
		z += 1


	extra_array = np.transpose(extra_array)

	return extra_array, extra_data_names

##########################################
"""Execute"""
w = 0
while w < 1:#len(file_name):
	print file_name[w] 							# Output tracking.
	print w_r


	f = open(file_list[w], "r") 				# Read current file
	lines = f.readlines()
	f.close()
	t = grab_t()
	modified_fo2, fmq, fo2 = grab_fo2()
	ph, modified_ph, neutral_ph, eh, pe, = grab_ph_eh_pe()

	output_array = [0]*len(t)  					# build blank row 1. This is left empty for any additional data about teh reactants i may want to stuff in.
	output_array = np.vstack((output_array, t, modified_fo2, fmq, fo2, ph, modified_ph, neutral_ph, eh, pe))

	all_names = []
	all_names.extend([file_name[w][0:-3]]) 		# insert file name 
	all_names.extend(['T','modified_fo2','fmq','fO2','pH','modified_ph','neutral_ph','eh','pe'])

	unit_list = ['','C','','','','','','','',''] 						# default info row and temperature row
	phase_list = ['','','','','','','','','',''] 						# Default info row and temperautre row

	#Grab all pertinant data, and their names, units, and phases
	if aq_grab == 1:
		aq_array = grab_aq()
		all_names.extend(aq_sp)
		output_array = np.vstack((output_array, aq_array))
		unit_list.extend([aq_unit_print]*len(aq_sp))
		phase_list.extend(['aq']*len(aq_sp))

	if aq_99_grab == 1:
		aq_99_array, aq_99_of_basis_expanded_list = grab_aq_99()
		all_names.extend(aq_99_of_basis_expanded_list)
		output_array = np.vstack((output_array, aq_99_array))
		unit_list.extend([aq_99_of_basis_unit_print]*len(aq_99_of_basis_expanded_list))
		phase_list.extend(['aq']*len(aq_99_of_basis_expanded_list))

	if s_grab == 1:
		s_array = grab_s()
		all_names.extend(s_sp)
		output_array = np.vstack((output_array, s_array))
		unit_list.extend([s_unit_print]*len(s_sp))
		phase_list.extend(['s']*len(s_sp))

	if g_grab == 1:
		g_array = grab_g()
		all_names.extend(g_sp)
		output_array = np.vstack((output_array, g_array))
		unit_list.extend([g_unit_print]*len(g_sp))
		phase_list.extend(['g']*len(g_sp))

	if s_element_grab ==1:
		master_output_array, master_output_names, w_r = grab_solid_element_tot()
		all_names.extend(master_output_names)
		output_array = np.vstack((output_array, master_output_array))
		unit_list.extend([s_element_unit_print]*(len(master_output_names)-1))
		unit_list.extend([''])
		phase_list.extend(['s_ele']*len(master_output_names))

	
	# Postfix Segment info:

	extra_array, extra_data_names = grab_seg_value()
	all_names.extend(extra_data_names)
	output_array = np.vstack((output_array, extra_array))
	unit_list.extend(['','','','','','','','','','','','','','','','','','','','','','',''])
	phase_list.extend(['','','','','','','','','','','','','','','','','','','','','','',''])






	# Remove backslashes from names before printing
	x = 0
	while x < len(all_names):
		a = all_names[x]
		all_names[x] = a.replace("\\", "") 							# Remove \ from names used to regex
		if ',g' in all_names[x]:										# Remove commas from sp_names and add (g) to gasses missing this
 			all_names[x] = all_names[x][:all_names[x].find(',')] + '(g)'
 		if ',' in all_names[x]:										
 			all_names[x] = all_names[x][:all_names[x].find(',')]
 		x += 1


	# Combine names and data using pandas
	x = 0
	while x < len(all_names):
		if x == 0:
			df1 = pd.DataFrame({all_names[x]:output_array[x]})
			x += 1
		else:
			df2 = pd.DataFrame({all_names[x]:output_array[x]})
			df1 = pd.concat([df1, df2], axis=1)
			x += 1

	df1 = pd.DataFrame.transpose(df1)

	df1.insert(0,0,phase_list,1) 		#	Insert units list as 2nd array column (first dataframe col)
	df1.insert(0,1,unit_list,1) 		#	Insert units list as 'New' 2nd array column (first dataframe col)
	

	short_file_name = file_name[w][0:-3]  +  '_'  +  run  +  '_any.csv'
	df1.to_csv(short_file_name,na_rep="NaN") 	#export to csv
	w += 1
















