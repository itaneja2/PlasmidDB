import psycopg2

DB = "***"
USERNAME = "***"
PASSWORD = "***"
HOST = "***"
PORT = "***"


def create_strain_view():

	try:
	    conn = psycopg2.connect(database = DB, user = USERNAME, password = PASSWORD, host = HOST, port = PORT)
	except:
	    print("unable to connect to the database") 

	cur = conn.cursor()
	query = "CREATE VIEW public.strainplasmid AS SELECT strain.strain_num, strain.plasmid_num, plasmid_wide.plasmid_name, plasmid_wide.vector, plasmid_wide.promoter1, plasmid_wide.insert1, plasmid_wide.terminator1, strain.species_name, strain.strain_name, strain.comments, strain.frozen_date, strain.frozen_by, strain.freezer_number, strain.shelf_number, strain.rack_number, strain.box_number, strain.box_index, strain.aliquot_num FROM strain LEFT JOIN plasmid_wide ON strain.plasmid_num = plasmid_wide.plasmid_num;" 
	cur.execute(query)

	conn.commit()
	conn.close()
	cur.close()

create_strain_view()
