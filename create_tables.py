import psycopg2

DB = "***"
USERNAME = "***"
PASSWORD = "***"
HOST = "***"
PORT = "***" 


def create_plasmid_wide_table(table_name):

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database") 

	cur = conn.cursor()
	query = "CREATE TABLE IF NOT EXISTS %s (id serial PRIMARY KEY, plasmid_num int UNIQUE, plasmid_name varchar, vector varchar, promoter1 varchar, insert1 varchar, terminator1 varchar, promoter2 varchar, insert2 varchar, terminator2 varchar, promoter3 varchar, insert3 varchar, terminator3 varchar, ori varchar, antibiotic1 varchar, antibiotic2 varchar, antibiotic3 varchar, inducer varchar, species_type varchar, construction_method varchar, cloning_strain varchar, constructed_by varchar, construction_date date, snapgene_file_link varchar, sequencing_results_link varchar, nucleotide_sequence varchar, cloning_strain_num int);" % table_name
	cur.execute(query)

	query = "CREATE UNIQUE INDEX ON %s (md5(nucleotide_sequence));" % table_name
	cur.execute(query)

	conn.commit()
	conn.close()
	cur.close()

#create_plasmid_wide_table('public.plasmid_wide_m')



def create_plasmid_long_table(table_name):

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database") 

	cur = conn.cursor()
	query = "CREATE TABLE IF NOT EXISTS %s (id serial PRIMARY KEY, plasmid_num int REFERENCES public.plasmid_wide(plasmid_num), feature varchar, type varchar, start_nucleotide int, end_nucleotide int, aa_sequence varchar, UNIQUE(plasmid_num, feature, type));" % table_name
	cur.execute(query)

	conn.commit()
	conn.close()
	cur.close()

#create_plasmid_long_table('public.plasmid_long')




def create_strain_table(table_name):

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database") 

	cur = conn.cursor()
	
	query = "CREATE TABLE IF NOT EXISTS %s (strain_num int UNIQUE, plasmid_num int, species_name varchar, strain_name varchar, comments varchar, frozen_date date, frozen_by varchar, freezer_number int, shelf_number int, rack_number int, box_number int, box_index int, aliquot_num int, antibiotic varchar, inducer varchar, PRIMARY KEY (strain_num));" % table_name
	cur.execute(query)

	conn.commit()
	conn.close()
	cur.close()

create_strain_table('public.strain')


####To initialize table: right click on the table name, click import data, toggle header to on, and select the appropriate columns, 
#data file is 'snapgene_feature_list.csv'

def create_typenames_table(table_name):

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database") 

	cur = conn.cursor()
	query = "CREATE TABLE IF NOT EXISTS %s (id serial PRIMARY KEY, type varchar UNIQUE, source varchar);" % table_name
	cur.execute(query)

	conn.commit()
	conn.close()
	cur.close()

#create_typenames_table('public.type_names')



####To initialize table: right click on the table name, click import data, toggle header to on, and select the appropriate columns, 
#data file is 'snapgene_feature_mapping - Feature_Type_Mapping_B.csv' and 'snapgene_feature_mapping - Feature_Type_Mapping_M.csv

def create_featuretypemapping_table(table_name):

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database") 

	cur = conn.cursor()
	query = "CREATE TABLE IF NOT EXISTS %s (id serial PRIMARY KEY, feature varchar UNIQUE, type varchar, source varchar, annotator varchar);" % table_name
	cur.execute(query)

	conn.commit()
	conn.close()
	cur.close()

#create_featuretypemapping_table('public.feature_type_mapping_bacterial')
#create_featuretypemapping_table('public.feature_type_mapping_mammalian')








