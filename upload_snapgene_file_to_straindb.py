import psycopg2
import pandas as pd 
import csv
import argparse
import subprocess
import shutil
import glob
import sys
import numpy as np 
import scipy.spatial.distance
import random
from fp_seq import * 
from snapgene_reader_w_primers import * 

DB = "***"
USERNAME = "***"
PASSWORD = "***"
HOST = "***"
PORT = "***"

valid_insert_type_list= ['FP', 'LKR', 'TAG', 'POI', 'CLV']

features_types_that_should_not_overlap = ['ori', 'promoter', 'terminator', 'antibiotic', 'inducer']

asterisk_line = '******************************************************************************'

plasmid_num_fname_dict = {}

def get_ori(snap_gene_info):

	ori_output = [] 

	keys = list(snap_gene_info.keys())
	keys = remove_items_from_list(keys, ['plasmid_name', 'plasmid_num', 'vector', 'species_type']) 

	for feature in keys:
		if 'type' in snap_gene_info[feature]:
			if snap_gene_info[feature]['type'] == 'rep_origin':
				ori = feature.split('_')[0] #.replace('high copy','')
				ori_output.append(ori)

	if len(ori_output) == 0:
		return None
	else:
		return ori_output[0]



def get_antibiotic(snap_gene_info):


	antibiotic_output = [] 

	keys = list(snap_gene_info.keys())
	keys = remove_items_from_list(keys, ['plasmid_name', 'plasmid_num', 'vector', 'species_type']) 

	for feature in keys:
		if 'type' in snap_gene_info[feature]:
			if snap_gene_info[feature]['type'] == 'antibiotic':
				antibiotic_output.append(feature.split('_')[0]) #remove _0 at end

	if len(antibiotic_output) == 0:
		return None,None,None
	else:
		if len(antibiotic_output) > 3:
			print(asterisk_line)
			print('ERROR: more than 3 antibiotics present')
			print(asterisk_line)
			sys.exit() 
		else:
			if len(antibiotic_output) == 1:
				return antibiotic_output[0],None,None
			elif len(antibiotic_output) == 2:
				return antibiotic_output[0],antibiotic_output[1],None
			else:
				return antibiotic_output[0],antibiotic_output[1],antibiotic_output[2] 

def get_inducer(snap_gene_info):

	keys = list(snap_gene_info.keys())
	keys = remove_items_from_list(keys, ['plasmid_name', 'plasmid_num', 'vector', 'species_type']) 

	for feature in keys:
		if 'type' in snap_gene_info[feature]:
			if snap_gene_info[feature]['type'] == 'inducer':
				return(feature.split('_')[0])

	return None



def insert_plasmid_wide_table(snap_gene_info, promoter_insert_terminator_featurename_dict):

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database")

	cur = conn.cursor() 

	plasmid_name = 	snap_gene_info['plasmid_name']		
	vector = snap_gene_info['vector']
	species_type = snap_gene_info['species_type'] 
	plasmid_num = snap_gene_info['plasmid_num']

	#based on how many promoter/insert/terminator present, we need different insert statements

	pit1_present = False
	pit2_present = False
	pit3_present = False

	promoter1 = None
	insert1 = None
	terminator1 = None

	promoter2 = None
	insert2 = None
	terminator2 = None

	promoter3 = None
	insert3 = None
	terminator3 = None

	if 1 in promoter_insert_terminator_featurename_dict:
		promoter1 = promoter_insert_terminator_featurename_dict[1]['promoter']
		insert1 = promoter_insert_terminator_featurename_dict[1]['insert']
		if 'terminator' in promoter_insert_terminator_featurename_dict[1]: #in case terminator does not exist (see plasmid 539)
			terminator1 = promoter_insert_terminator_featurename_dict[1]['terminator']
		pit1_present = True 

	if 2 in promoter_insert_terminator_featurename_dict:
		promoter2 = promoter_insert_terminator_featurename_dict[2]['promoter']
		insert2 = promoter_insert_terminator_featurename_dict[2]['insert']
		if 'terminator' in promoter_insert_terminator_featurename_dict[2]:
			terminator2 = promoter_insert_terminator_featurename_dict[2]['terminator']
		pit2_present = True 

	if 3 in promoter_insert_terminator_featurename_dict:
		promoter3 = promoter_insert_terminator_featurename_dict[3]['promoter']
		insert3 = promoter_insert_terminator_featurename_dict[3]['insert']
		if 'terminator' in promoter_insert_terminator_featurename_dict[3]:
			terminator3 = promoter_insert_terminator_featurename_dict[3]['terminator']
		pit3_present = True


	ori = get_ori(snap_gene_info)
	antibiotic1,antibiotic2,antibiotic3 = get_antibiotic(snap_gene_info)
	inducer = get_inducer(snap_gene_info)
	nucleotide_sequence = snap_gene_info['seq']

	print(asterisk_line)
	print('antibiotics found:')
	print(antibiotic1)
	print(antibiotic2)
	print(antibiotic3)
	print(asterisk_line)

	cur.execute('SELECT plasmid_num FROM public.plasmid_wide WHERE nucleotide_sequence = %s;', [nucleotide_sequence,])
	plasmid_num_query_output = cur.fetchall()

	if len(plasmid_num_query_output) > 0: #nucleotide_sequence already exists in DB
		conflict_plasmid_num = plasmid_num_query_output[0][0]
		if int(conflict_plasmid_num) != int(plasmid_num):
			print(asterisk_line)
			print('ERROR: nucleotide sequence already exists in DB as a different plasmid_num %s - cannot overwrite row with plasmid_num %s. can only overwrite with custom_plasmid_num if nucleotide sequence does not exist in DB' % (conflict_plasmid_num,plasmid_num))
			print(asterisk_line)
			sys.exit() 

	#on plasmid_num conflict, update row
	##its possible the nucleotide seq changed, so update that
	insert_query = 'INSERT INTO public.plasmid_wide (plasmid_num, plasmid_name, vector, promoter1, insert1, terminator1, promoter2, insert2, terminator2, promoter3, insert3, terminator3, ori, antibiotic1, antibiotic2, antibiotic3, inducer, species_type, nucleotide_sequence) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) ON CONFLICT (plasmid_num) DO UPDATE SET (plasmid_name, vector, promoter1, insert1, terminator1, promoter2, insert2, terminator2, promoter3, insert3, terminator3, ori, antibiotic1, antibiotic2, antibiotic3, inducer, species_type, nucleotide_sequence) = ROW(EXCLUDED.plasmid_name, EXCLUDED.vector, EXCLUDED.promoter1, EXCLUDED.insert1, EXCLUDED.terminator1, EXCLUDED.promoter2, EXCLUDED.insert2, EXCLUDED.terminator2, EXCLUDED.promoter3, EXCLUDED.insert3, EXCLUDED.terminator3, EXCLUDED.ori, EXCLUDED.antibiotic1, EXCLUDED.antibiotic2, EXCLUDED.antibiotic3, EXCLUDED.inducer, EXCLUDED.species_type, EXCLUDED.nucleotide_sequence);'
	insert_data = (plasmid_num, plasmid_name, vector, promoter1, insert1, terminator1, promoter2, insert2, terminator2, promoter3, insert3, terminator3, ori, antibiotic1, antibiotic2, antibiotic3, inducer, species_type, nucleotide_sequence)

	print(insert_data)
	cur.execute(insert_query, insert_data)

	conn.commit() 
	conn.close()
	cur.close()

	return plasmid_num


def insert_plasmid_long_table(plasmid_long_data, snap_gene_info):

	nucleotide_sequence = snap_gene_info['seq']
	plasmid_num = snap_gene_info['plasmid_num']

	start_list_circular_order = plasmid_long_data[0]
	end_list_circular_order = plasmid_long_data[1]
	feature_list_circular_order = plasmid_long_data[2]
	type_list_circular_order = plasmid_long_data[3]
	aa_list_circular_order = plasmid_long_data[4]

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database")

	cur = conn.cursor() 

	for i in range(0,len(start_list_circular_order)):

		start_nucleotide = int(start_list_circular_order[i])
		end_nucleotide = int(end_list_circular_order[i])
		feature = feature_list_circular_order[i]
		type = type_list_circular_order[i]
		aa_sequence = aa_list_circular_order[i]
	
		insert_query = 'INSERT INTO public.plasmid_long (plasmid_num, feature, type, start_nucleotide, end_nucleotide, aa_sequence) VALUES (%s, %s, %s, %s, %s, %s) ON CONFLICT (plasmid_num, feature, type) DO UPDATE SET (start_nucleotide, end_nucleotide, aa_sequence) = ROW(EXCLUDED.start_nucleotide, EXCLUDED.end_nucleotide, EXCLUDED.aa_sequence);'
		insert_data = (plasmid_num, feature, type, start_nucleotide, end_nucleotide, aa_sequence)
		#print(insert_data)
		cur.execute(insert_query, insert_data)

	conn.commit() 
	conn.close()
	cur.close()




def get_feature_map(snapgene_file):


	plasmid_name = (snapgene_file.split('.dna')[0]).split('/')[-1] #remove .dna and path info
	if plasmid_name.count('_') < 2:
		print(asterisk_line)
		print('ERROR: %s is not formatted correctly' % plasmid_name)
		print("Should be of the format pLLXXX_vector_M/B_XYZ_123")
		sys.exit()
	species_type = plasmid_name.split('_')[2]


	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database")

	cur = conn.cursor() 

	feature_type_map = {} 

	if species_type in ['M','I']:
		cur.execute('SELECT * FROM public.feature_type_mapping_mammalian;')
		out = cur.fetchall()
		for i in range(0,len(out)):
			feature = out[i][1]
			type = out[i][2]
			feature_type_map[feature] = type
	else:
		cur.execute('SELECT * FROM public.feature_type_mapping_bacterial;')
		out = cur.fetchall()
		for i in range(0,len(out)):
			feature = out[i][1]
			type = out[i][2]
			feature_type_map[feature] = type


	return(feature_type_map)

	conn.commit() 
	conn.close()
	cur.close()



#check to see if FP matches list of expected FPs and sequence matches expected sequence 
def apply_FP_check(snap_gene_info):

	print(fp_seq_map)

	keys = list(snap_gene_info.keys())
	keys = remove_items_from_list(keys, ['plasmid_name', 'plasmid_num', 'vector', 'species_type']) 

	for feature in keys:
		
		if 'type' in snap_gene_info[feature]:

			type = snap_gene_info[feature]['type']

			if type == 'FP':

				translation = snap_gene_info[feature]['translation'].replace(',','') #the replace statement is needed
				translation = translation.replace('*','')
				
				if translation[0] == 'M':
					M_beg = True
					comparison_start_idx = 0
				else:
					M_beg = False 
					comparison_start_idx = 1 

				fp = feature.split('_')[0]

				if fp not in fp_seq_map:
					print(asterisk_line)
					print('ERROR: FP_%s does not match list of expected:' % fp)
					print('Feature: %s' % feature)				
					print(asterisk_line)
					sys.exit()
				elif fp_seq_map[fp][comparison_start_idx:] != translation: #ignore M at beginning 
					print(asterisk_line)
					print('ERROR: FP_%s translation does not match expected' % fp)
					print('Feature: %s' % feature)
					print('Expected Sequence: %s' % fp_seq_map[fp][comparison_start_idx:])
					print('Actual Sequence: %s' % translation)
					print(asterisk_line)
					user_input = input("Do you still want to upload? Enter either Y or N: ").strip().lower()
					if user_input == 'n':
						sys.exit() 

				print(asterisk_line)
				print('PASS: FP %s does not have any errors' % feature.split('_')[0])
				print(asterisk_line)


#check to see if POI naming convention is correct 
def apply_POI_check(snap_gene_info):

	keys = list(snap_gene_info.keys())
	keys = remove_items_from_list(keys, ['plasmid_name', 'plasmid_num', 'vector', 'species_type']) 

	naming_convention_exceptions = ['turboid']

	for feature in keys:

		if 'type' in snap_gene_info[feature]:

			type = snap_gene_info[feature]['type']

			if type == 'POI':

				if any(exception in feature.lower() for exception in naming_convention_exceptions):
					continue
				else:
					if feature.count('-') < 2:
						print(asterisk_line)
						print('WARNING: POI %s should be of format genusname-genename-WT/mutantinfo' % feature)
						user_input = input("Do you still want to upload? Enter either Y or N: ").strip().lower()
						if user_input == 'n':
							sys.exit() 


				print(asterisk_line)
				print('PASS: POI %s does not have any errors' % feature)
				print(asterisk_line)


def remove_items_from_list(l, items_to_remove):
    return [i for i in l if not i in items_to_remove]

#check to see if any feature have overlap in nucleotide sequence. this only applies to features with a defined type 

def apply_overlap_check(snap_gene_info):

	type_startend_dict = {} 
	
	startend_dict = {}
	
	keys = list(snap_gene_info.keys())
	keys = remove_items_from_list(keys, ['plasmid_name', 'plasmid_num', 'vector', 'species_type'])


	#first check that end is greater than start 
	for i in range(0,len(keys)):
		
		feature1 = keys[i]

		if 'type' in snap_gene_info[feature1]:
		
			start1 = snap_gene_info[feature1]['start']
			end1 = snap_gene_info[feature1]['end']
			type1 = snap_gene_info[feature1]['type']

			if end1 <= start1:
				print('ERROR: end is less than start for:')
				print('feature1: %s' % feature1)
				print(snap_gene_info[feature1])
				print('To the adjust the starting position of the plasmid in SnapGene, select the corresponding origin of replication feature and click `View>Set Origin` to set that point as base 1.')


	for i in range(0,len(keys)):
		
		feature1 = keys[i]

		if 'type' in snap_gene_info[feature1]:
		
			start1 = snap_gene_info[feature1]['start']
			end1 = snap_gene_info[feature1]['end']
			type1 = snap_gene_info[feature1]['type']

			for j in range(0,len(keys)):

				if i != j:
				
					feature2 = keys[j]

					if 'type' in snap_gene_info[feature2]:

						start2 = snap_gene_info[feature2]['start']
						end2 = snap_gene_info[feature2]['end']
						type2 = snap_gene_info[feature2]['type']
			
						if type1 is not None and type2 is not None:

							if type1 != 'misc_feature' or type2 != 'misc_feature':

							#if type1 in features_types_that_should_not_overlap and type2 in features_types_that_should_not_overlap:
		
								if (start2 >= start1 and start2 <= end1):						
									print('ERROR: overlap exists bewteen:')
									print('feature1: %s' % feature1)
									print(snap_gene_info[feature1])
									print('feature2: %s' % feature2)
									print(snap_gene_info[feature2])
									sys.exit()
			


	print(asterisk_line)
	print('PASS: No 2 features overlap in terms of start/end')
	print(asterisk_line)


#stop codon should exist at end of insert and there should be only one stop codon  
def apply_stop_codon_check(snap_gene_info, plasmid_long_data, promoter_insert_terminator_featurename_dict, promoter_insert_terminator_featureidx_dict):

	start_list_circular_order = plasmid_long_data[0]
	end_list_circular_order = plasmid_long_data[1]
	feature_list_circular_order = plasmid_long_data[2]
	type_list_circular_order = plasmid_long_data[3]
	aa_list_circular_order = plasmid_long_data[4]
	strand_list_circular_order = plasmid_long_data[5]

	for key in promoter_insert_terminator_featureidx_dict:

		insert_idx_list = promoter_insert_terminator_featureidx_dict[key]['insert']
		insert_aa = '' 

		for insert_idx in insert_idx_list:
			insert_aa += aa_list_circular_order[insert_idx]

		if insert_aa[-1] != '*':
			print(asterisk_line)
			print('WARNING: insert sequence %s does not have stop codon at end' % insert_aa)
			print('this insert corresponds to the following promoter/insert/terminator:')
			print(promoter_insert_terminator_featurename_dict[key])
			print(asterisk_line)
			user_input = input("Do you still want to upload? Enter either Y or N: ").strip().lower()
			if user_input == 'n':
				sys.exit() 
		if insert_aa.count('*') > 1:
			print(asterisk_line)
			print('ERROR: insert sequence %s has more than one stop codon' % insert_aa)
			print('this insert corresponds to the following promoter/insert/terminator:')
			print(promoter_insert_terminator_featurename_dict[key])
			print(asterisk_line)
			user_input = input("Do you still want to upload? Enter either Y or N: ").strip().lower()
			if user_input == 'y':
				sys.exit() 

		print(asterisk_line)
		print('PASS: insert sequence %s does not have any errors' % insert_aa)
		print('this insert corresponds to the following promoter/insert/terminator:')
		print(promoter_insert_terminator_featurename_dict[key])
		print(asterisk_line)



#need to reconstruct circular plasmid and then you can infer promoter_insert_terminator 
def promoter_insert_terminator(snap_gene_info, primer_info, ignore_pit):

	nearest_type_dict = {} 

	start_list = []
	end_list = [] 
	feature_list = [] 
	type_list = [] 
	aa_list = [] 
	strand_list = []

	#create list consisting of start and end values
	#get pairwise distance start-end and end-start

	keys = list(snap_gene_info.keys())
	keys = remove_items_from_list(keys, ['plasmid_name', 'plasmid_num', 'vector', 'species_type']) 

	for feature in keys:

		if 'type' in snap_gene_info[feature]:

			type = snap_gene_info[feature]['type']

			if type is not None:

				start = snap_gene_info[feature]['start']
				end = snap_gene_info[feature]['end']

				start_list.append(start)
				end_list.append(end)
				feature_list.append(feature.split('_')[0])
				type_list.append(type)

				if 'strand' in snap_gene_info[feature]:
					strand_list.append(snap_gene_info[feature]['strand'])
				else:
					strand_list.append('')

				if 'translation' in snap_gene_info[feature]:
					aa_list.append(snap_gene_info[feature]['translation'].replace(',',''))
				else:
					aa_list.append('')


	print(start_list)
	print(end_list)
	print(feature_list)
	print(type_list)
	print(asterisk_line)

	#find the smallest start value and the largest end value
	min_start = min(start_list)
	max_end = max(end_list)
	num_features = len(start_list)


	min_start_idx = np.argmin(np.array(start_list))
	max_end_idx = np.argmax(np.array(end_list))

	start_list = np.array(start_list)
	end_list = np.array(end_list)

	nearest_feature_start_end = {} #maps to corresponding index
	nearest_feature_end_start = {} #maps to corresponding index

	dist_matrix_input = np.zeros((len(start_list)+len(end_list),))

	dist_matrix_input[0:len(start_list)] = start_list
	dist_matrix_input[len(start_list):] = end_list

	pairwise_dist = (dist_matrix_input[:, None] - dist_matrix_input)
	pairwise_dist_start_end = (pairwise_dist[0:len(start_list),len(start_list):]) #pairwise distance between each start and end index

	#for each start find nearest end that is to the LEFT (counterclockwise) --> you are allowed to look for any feature whose end is less than current start OR IF at first feature, then you can look at the last feature 
	for i in range(0,pairwise_dist_start_end.shape[1]):
		pairwise_dist_start_end_copy = np.copy(pairwise_dist_start_end)
		invalid_idx = (pairwise_dist_start_end_copy[i,:] <= 0)
		if i == min_start_idx:
			invalid_idx[max_end_idx] = False #last feature is valid to look at if we are at first feature 
		invalid_idx[i] = True #you can't include yourself for finding nearest 
		pairwise_dist_start_end_copy[i,invalid_idx] = np.inf 
		nearest_end_idx = np.argmin(np.abs(pairwise_dist_start_end_copy[i,]))
		nearest_feature_start_end[i] = nearest_end_idx

	#this is reverse order to make the signs work 
	pairwise_dist = (dist_matrix_input - dist_matrix_input[:, None])
	pairwise_dist_end_start = (pairwise_dist[len(start_list):,0:len(start_list)]) #pairwise distance between each end and start index

	#for each end find nearest start that is to the RIGHT (clockwise) --> you are allowed to look for any feature whose start is greater than current end  OR IF at last feature, then you can look at the first feature e
	for i in range(0,pairwise_dist_end_start.shape[1]):
		pairwise_dist_end_start_copy = np.copy(pairwise_dist_end_start)
		invalid_idx = pairwise_dist_end_start_copy[i,:] <= 0 
		if i == max_end_idx:
			invalid_idx[min_start_idx] = False #first feature is valid to look at if we are at last feature 
		invalid_idx[i] = True #you can't include yourself for finding nearest 
		pairwise_dist_end_start_copy[i,invalid_idx] = np.inf 
		nearest_start_idx = np.argmin(np.abs(pairwise_dist_end_start_copy[i,]))
		nearest_feature_end_start[i] = nearest_start_idx
	
	
	#we only need one of nearest_feature_start_end/nearest_feature_end_start
	#however they should be consistent --> we verified that they are consistent 

	print(nearest_feature_start_end)
	print(nearest_feature_end_start)

	if (len(nearest_feature_start_end) != num_features) or (len(nearest_feature_end_start) != num_features):
		print('ERROR: overlapping features exist')
		sys.exit()

	######constructing the correctly ordered circular plasmid

	#we will start with the ori 
	if ignore_pit == False:
		rep_origin_idx = -1 
		for i in range(0,len(type_list)):
			if type_list[i] == 'rep_origin':
				rep_origin_idx = i 

		if rep_origin_idx == -1:
			print(asterisk_line)
			print('ERROR: no ori exists in plasmid')
			print(asterisk_line)
			sys.exit() 
	else:
		rep_origin_idx = 0 



	##create circular plasmid for start/end and end/start 

	#starting idx should be rep_origin 
	starting_type_idx = rep_origin_idx

	curr_idx = starting_type_idx
	circular_plasmid_idx_list = [starting_type_idx]
	circular_plasmid_idx_str = str(starting_type_idx)
	plasmid_construction_incomplete = True

	num_iter = 0 

	while plasmid_construction_incomplete:

		if nearest_feature_start_end[curr_idx] == starting_type_idx:
			plasmid_construction_incomplete = False #succesfully mapped plasmid 
		else:
			circular_plasmid_idx_str = circular_plasmid_idx_str + '_' + str(nearest_feature_start_end[curr_idx])
			circular_plasmid_idx_list.append(nearest_feature_start_end[curr_idx])
			curr_idx = nearest_feature_start_end[curr_idx]

		num_iter += 1

		if num_iter > num_features:
			print("ERROR: cannot construct circular plasmid")
			feature_list_circular_order = ([feature_list[i] for i in circular_plasmid_idx_list])
			print(feature_list_circular_order)
			print(circular_plasmid_idx_list)
			sys.exit()


	circular_plasmid_idx_str_start_end = circular_plasmid_idx_str
	circular_plasmid_idx_list_start_end = circular_plasmid_idx_list.copy()

	#starting idx should be rep_origin 
	starting_type_idx = rep_origin_idx

	curr_idx = starting_type_idx
	circular_plasmid_idx_list = [starting_type_idx]
	circular_plasmid_idx_str = str(starting_type_idx)
	plasmid_construction_incomplete = True

	num_iter = 0 

	while plasmid_construction_incomplete:

		if nearest_feature_end_start[curr_idx] == starting_type_idx:
			plasmid_construction_incomplete = False #succesfully mapped plasmid 
		else:
			circular_plasmid_idx_str = circular_plasmid_idx_str + '_' + str(nearest_feature_end_start[curr_idx])
			circular_plasmid_idx_list.append(nearest_feature_end_start[curr_idx])
			curr_idx = nearest_feature_end_start[curr_idx]

		num_iter += 1

		if num_iter > num_features:
			print("ERROR: cannot construct circular plasmid")
			feature_list_circular_order = ([feature_list[i] for i in circular_plasmid_idx_list])
			print(feature_list_circular_order)
			print(circular_plasmid_idx_list)
			sys.exit()


	rev_list = circular_plasmid_idx_list_start_end[1:][::-1]
	if circular_plasmid_idx_list[1:] != rev_list:
		print('ERROR: circular plasmid inconsistent between start/end and end/start')
		print(circular_plasmid_idx_list[1:])
		print(rev_list)
		sys.exit()



	#print(circular_plasmid_idx_str)

	start_list_circular_order = ([start_list[i] for i in circular_plasmid_idx_list])	
	end_list_circular_order = ([end_list[i] for i in circular_plasmid_idx_list])	
	feature_list_circular_order = ([feature_list[i] for i in circular_plasmid_idx_list])
	type_list_circular_order = ([type_list[i] for i in circular_plasmid_idx_list])
	aa_list_circular_order = ([aa_list[i] for i in circular_plasmid_idx_list])
	strand_list_circular_order = ([strand_list[i] for i in circular_plasmid_idx_list])


	plasmid_long_data = [start_list_circular_order, end_list_circular_order, feature_list_circular_order, type_list_circular_order, aa_list_circular_order, strand_list_circular_order]

	if ignore_pit:
		return plasmid_long_data, {}, {} 

	if len(primer_info) > 0:
		for i in range(0,len(primer_info)):
			plasmid_long_data[0].append(primer_info[i]['start'])
			plasmid_long_data[1].append(primer_info[i]['end'])
			plasmid_long_data[2].append(primer_info[i]['name'])
			plasmid_long_data[3].append(primer_info[i]['type'])
			plasmid_long_data[4].append(primer_info[i]['seq'])


	print("CIRCULAR ORDER")
			
	print(feature_list_circular_order)
	print(type_list_circular_order)
	print(circular_plasmid_idx_list)


	#assemble promoter/insert/terminator 
	#note we can have multiple promoters/inserts/terminators (up to 3)

	#for each promoter assess if there is a FP/LKR/TAG/POI/CLV next to it --> note: NEED to use circular index for this 
	#if this exists, keep appending to insert until a terminator is reached 

	promoter_insert_terminator_featurename_dict = {} #this stores a mapping between p/i/t and feature name 
	promoter_insert_terminator_featureidx_dict = {} #this stores a mapping between p/i/t and feature index (circular index)
	promoter_insert_terminator_dict_idx = 1 

	for i in range(0,len(type_list_circular_order)):

		if type_list_circular_order[i] == 'promoter':

			print('at promoter:')
			print(i)
			print(feature_list_circular_order[i])
			print(feature_list_circular_order)
			print(feature_list_circular_order[(i+1)%len(type_list_circular_order)])
			print(type_list_circular_order[(i+1)%len(type_list_circular_order)])
			
			
			#example: if i = 3 and there are 8 elements in type_list_circular_order:
			##then circular_idx = [4,5,6,7,0,1,2]
			circular_idx = []
			for j in range(i+1,len(type_list_circular_order)):
				circular_idx.append(j)
			for j in range(0,i):
				circular_idx.append(j)
			print(circular_idx)

			if (promoter_insert_terminator_dict_idx not in promoter_insert_terminator_featurename_dict) and (promoter_insert_terminator_dict_idx <= 3): #can't have more than 3 promoter/insert/terminator

				promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx] = {}
				promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx] = {}

				promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]['promoter'] = feature_list_circular_order[i]
				promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx]['promoter'] = [i] 


			for j in circular_idx: #use circular index 

				print(j)

				if type_list_circular_order[j] == 'promoter': 

					if 'insert' not in promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx]: #if insert hasn't been found, don't use current promoter, use a different one
						del promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx]
						del promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]
						break #breaks out of inner for loop 

				if type_list_circular_order[j] in valid_insert_type_list:

					print('insert')

					if 'insert' not in promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]:
						promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]['insert'] = ''
						promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx]['insert'] = [] 

					if j not in promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx]['insert']: #don't append duplicates. this can happen if a terminator doesn't exist, but there are two promoters. see plasmid 539 for example
						promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]['insert'] += feature_list_circular_order[j] + '_'
						promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx]['insert'].append(j)


				elif type_list_circular_order[j] == 'terminator': #this should work unless terminator is in between insert which shouldn't happen....
					
					print('terminator')

					if 'insert' in promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]: #only add terminator if insert exists
						promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]['insert'] = promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]['insert'][0:-1] #remove the last '_' that is appended to insert 
						promoter_insert_terminator_featurename_dict[promoter_insert_terminator_dict_idx]['terminator'] = feature_list_circular_order[j]
						promoter_insert_terminator_featureidx_dict[promoter_insert_terminator_dict_idx]['terminator'] = [j] 
						
						promoter_insert_terminator_dict_idx += 1  #only increment if you success
						break 



	print(promoter_insert_terminator_featurename_dict)
	print(promoter_insert_terminator_featureidx_dict)


	keys_to_delete = []
	if len(promoter_insert_terminator_featurename_dict.keys()) == 0:
		print('ERROR: No Promoter Found')
		sys.exit()
	else:
		for key in promoter_insert_terminator_featurename_dict:

			if 'insert' not in promoter_insert_terminator_featurename_dict[key]:
				print('Ignoring the following promoter:') 
				print(promoter_insert_terminator_featurename_dict[key])
				keys_to_delete.append(key)
			else:

				if promoter_insert_terminator_featurename_dict[key]['insert'][-1] == '_':
					promoter_insert_terminator_featurename_dict[key]['insert'] = promoter_insert_terminator_featurename_dict[key]['insert'][0:-1]

				if 'terminator' not in promoter_insert_terminator_featurename_dict[key]:
					print('WARNING: No Terminator Found In:') 
					print(promoter_insert_terminator_featurename_dict[key])
					user_input = input("Do you still want to upload? Enter either Y or N: ").strip().lower()
					if user_input == 'n':
						sys.exit() 


	for key in keys_to_delete:
		del promoter_insert_terminator_featureidx_dict[key]
		del promoter_insert_terminator_featurename_dict[key]

	return(plasmid_long_data, promoter_insert_terminator_featurename_dict, promoter_insert_terminator_featureidx_dict)




def process_snapgene_file(snapgene_file, feature_type_map):

	verbose = True

	snap_gene_info = {} 

	plasmid_name = (snapgene_file.split('.dna')[0]).split('/')[-1] #remove .dna and path info
	if plasmid_name.count('_') < 2:
		print(asterisk_line)
		print('ERROR: %s is not formatted correctly' % plasmid_name)
		print("Should be of the format pLLXXX_vector_M/B_XYZ_123")
		sys.exit()

	snap_gene_info['plasmid_name'] = plasmid_name
	snap_gene_info['plasmid_num'] = int(plasmid_name.split('_')[0].replace('pLL',''))
	snap_gene_info['vector'] = plasmid_name.split('_')[1]
	snap_gene_info['species_type'] = plasmid_name.split('_')[2]

	print(asterisk_line)
	print(asterisk_line)
	print('Plasmid Name Information:')
	print(snap_gene_info)
	print(asterisk_line)
	print(asterisk_line)
	print('\n')


	if snap_gene_info['species_type'] not in ['M','B','I','NA']:
		print(asterisk_line)
		print('ERROR: %s is not formatted correctly' % snap_gene_info['plasmid_name'])
		print("Filename is missing species name (M or B or I)")
		print(asterisk_line)
		sys.exit() 		

	snapgene_dict = snapgene_file_to_dict(snapgene_file)
	seqrecord = snapgene_file_to_seqrecord(snapgene_file)

	snap_gene_info['seq'] = snapgene_dict['seq'].upper()

	'''if verbose:
		for key in snapgene_dict:
			print("KEY: %s" % key)
			print(snapgene_dict[key])'''

	if 'features' not in snapgene_dict:
		print(asterisk_line)
		print("ERROR: features not present in snapgene file")
		print(asterisk_line)
		sys.exit() 

	features_info = snapgene_dict['features']
	primer_info = snapgene_dict['primers']

	feature_name_count = {}

	for i in range(0,len(features_info)):

		if verbose:
			print(features_info[i])
			print('\n')

		curr_feature_data = features_info[i]
		feature_name = curr_feature_data['name']

		strand = curr_feature_data['strand']
		start_idx = curr_feature_data['start']
		end_idx = curr_feature_data['end']

		start = start_idx+1 #0-idx
		end = end_idx

		if ('FP_' in feature_name) or ('LKR_' in feature_name) or ('POI_' in feature_name) or ('TAG_' in feature_name) or ('CLV_' in feature_name): #feature corresponds to insert

			if feature_name in feature_name_count:
				feature_name_count[feature_name] += 1
			else:
				feature_name_count[feature_name] = 0 


			if feature_name.count('_') != 1:
				print(asterisk_line)
				print('ERROR: %s is formatted incorrectly. should be formatted as FP_fpname, LKR_linkername, TAG_tagname, POI_genus-genename-WT, POI_genus-genename-mutantinfo, CLV_cleavagesite' % feature_name) 
				print(asterisk_line)
				sys.exit()
			else:
				print(feature_name)
				print(feature_name_count)
				feature_name_insert_subset = feature_name.split('_')[1]
				type = feature_name.split('_')[0]
				feature_name_insert_subset_w_num = feature_name_insert_subset + '_' + str(feature_name_count[feature_name])

				snap_gene_info[feature_name_insert_subset_w_num] = {}
				snap_gene_info[feature_name_insert_subset_w_num]['start'] = start
				snap_gene_info[feature_name_insert_subset_w_num]['end'] = end	
				snap_gene_info[feature_name_insert_subset_w_num]['type'] = type

				if 'qualifiers' in curr_feature_data:
					if 'translation' in curr_feature_data['qualifiers']:
						translation = features_info[i]['qualifiers']['translation']
						snap_gene_info[feature_name_insert_subset_w_num]['translation'] = translation
					else:
						print(asterisk_line)
						print('ERROR: %s does not have protein translation. In SnapGene, make sure translate this feature is CHECKED' % feature_name) 
						print(asterisk_line)
						sys.exit()

		else: #feature is a common feature 

			feature_name_wo_space = feature_name.replace('_',' ') #snapgene_reader python library replaces space with _ --> in original snapgene file, there are no underscores (only spaces)

			feature_name_wo_space = feature_name_wo_space.replace('Promoter','promoter')
			feature_name_wo_space = feature_name_wo_space.replace('Enhancer','enhancer')
			feature_name_wo_space = feature_name_wo_space.replace('Terminator','terminator')
			feature_name_wo_space = feature_name_wo_space.strip()

			if feature_name_wo_space in feature_name_count:
				feature_name_count[feature_name_wo_space] += 1
			else:
				feature_name_count[feature_name_wo_space] = 0 

			#this is to account for duplicates [i.e I can have the same promoter twice]
			feature_name_wo_space_w_num = feature_name_wo_space + '_' + str(feature_name_count[feature_name_wo_space])

			snap_gene_info[feature_name_wo_space_w_num] = {}
			snap_gene_info[feature_name_wo_space_w_num]['start'] = start
			snap_gene_info[feature_name_wo_space_w_num]['end'] = end
			
			if feature_name_wo_space in feature_type_map:
				type = feature_type_map[feature_name_wo_space]
				snap_gene_info[feature_name_wo_space_w_num]['type'] = type
				print("PRESENT in feature_type mapping")
				print(feature_name_wo_space)
			else:
				print("NOT PRESENT in feature_type mapping")
				print(feature_name_wo_space)


			if 'qualifiers' in curr_feature_data:
				if 'translation' in curr_feature_data['qualifiers']:
					translation = features_info[i]['qualifiers']['translation']
					snap_gene_info[feature_name_wo_space_w_num]['translation'] = translation


	return snap_gene_info, primer_info 




if __name__=="__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument("--snapgene_filepath", required=True, help="can either be a .DNA file or a folder containing a list of snapgene files")
	#parser.add_argument("--ignore_pit", action='store_true', help="if promoter/insert/terminator not present, ignore these fields when uploading")
	
	args = parser.parse_args()

	snapgene_file_list = [] 

	if '.dna' in args.snapgene_filepath:
		snapgene_file_list.append(args.snapgene_filepath)  
	else:
		for snapgene_fn in glob.glob("%s/*.dna" % args.snapgene_filepath):
			snapgene_file_list.append(snapgene_fn) 

	print(asterisk_line)
	print("These are the following snapgene files that will be uploaded to DB:")
	print(snapgene_file_list)
	print(asterisk_line)
	

	for snapgene_file in sorted(snapgene_file_list):

		feature_type_map = get_feature_map(snapgene_file)
		
		print('\n ')
		print(asterisk_line)
		print(asterisk_line)
		print("PROCESING SNAPGENE FILE %s:" % snapgene_file)
		print(asterisk_line)
		print(asterisk_line)
		print('\n ')
		snap_gene_info, primer_info = process_snapgene_file(snapgene_file, feature_type_map)
		print("PRIMER INFO:")
		print(primer_info)
		print('\n ')
		print(asterisk_line)
		print("PROCESSED SNAPGENE FILE:")
		print(snap_gene_info)
		print(asterisk_line)
		print('\n ')

		if snap_gene_info['vector'] == 'NA' or snap_gene_info['species_type'] == 'NA':
			ignore_pit = True 
		else:
			ignore_pit = False 

		print('\n ')
		print(asterisk_line)
		print("APPLY overlap CHECK (check to see if type has multiple features with same nucleotide sequence range):")
		apply_overlap_check(snap_gene_info)
		print(asterisk_line)
		print('\n ')


		print('\n ')
		print(asterisk_line)
		print("DETECTING PROMOTER/INSERT/TERMINATOR REGIONS")
		print(asterisk_line)
		print('\n ')

		plasmid_long_data, promoter_insert_terminator_featurename_dict, promoter_insert_terminator_featureidx_dict = promoter_insert_terminator(snap_gene_info, primer_info, ignore_pit)


		if ignore_pit == False:
			print('\n ')
			print(asterisk_line)
			print("DETECTED PROMOTER/INSERT/TERMINATOR REGIONS:")
			print(promoter_insert_terminator_featurename_dict)
			print(asterisk_line)
			print('\n ')


			print('\n ')
			print(asterisk_line)
			print("APPLY FP CHECK (if FP exists, checks to see if sequence matches expected):")
			apply_FP_check(snap_gene_info)
			print(asterisk_line)
			print('\n ')


			print('\n ')
			print(asterisk_line)
			print("APPLY stop codon CHECK (stop codon should be present at end of each insert and there should only be one stop codon in insert):")
			apply_stop_codon_check(snap_gene_info, plasmid_long_data, promoter_insert_terminator_featurename_dict, promoter_insert_terminator_featureidx_dict)		
			print(asterisk_line)
			print('\n ')

			print('\n ')
			print(asterisk_line)
			print("APPLY POI CHECK (if POI exists, checks to see if it is formatted correctly):")
			apply_POI_check(snap_gene_info)
			print(asterisk_line)
			print('\n ')


		print('\n ')
		print(asterisk_line)
		print("INSERTING INTO PLASMID WIDE:")
		print(asterisk_line)
		print('\n ')

		plasmid_num = insert_plasmid_wide_table(snap_gene_info, promoter_insert_terminator_featurename_dict)

		print('\n ')
		print(asterisk_line)
		print("INSERTED PLASMID NUM %d INTO PLASMID WIDE:" % plasmid_num)
		print(asterisk_line)
		print('\n ')

		print('\n ')
		print(asterisk_line)
		print("INSERTING INTO PLASMID LONG:")
		print(asterisk_line)
		print('\n ')

		insert_plasmid_long_table(plasmid_long_data, snap_gene_info)

		print('\n ')
		print(asterisk_line)
		print("INSERTED PLASMID NUM %d INTO PLASMID LONG:" % plasmid_num)
		print("COMPLETED UPLOAD OF PLASMID NUM %d" % plasmid_num)
		print(asterisk_line)
		print('\n ')

		plasmid_num_fname_dict[plasmid_num] = snapgene_file

		print(asterisk_line)
		print(asterisk_line)
		print(asterisk_line)
		print(asterisk_line)
		print('\n')

	
	print(plasmid_num_fname_dict)






