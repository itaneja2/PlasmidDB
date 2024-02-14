import csv
import argparse
import subprocess
import shutil
import glob
import pandas as pd


snapgene_file_list = []

for snapgene_fn in glob.glob("%s/*.dna" % './SnapGene_StandardFeatures'):
	snapgene_file_list.append(snapgene_fn) 

snapgene_file_list = sorted(snapgene_file_list, key=str.casefold)

#print(snapgene_file_list)

snapgene_feature_list_cleaned = [] 

for snapgene_fn in snapgene_file_list:
	feature_name = (snapgene_fn.split('/')[-1].split('.'))[0]

	last_chars = feature_name.split(' ')[-1]
	if last_chars[0] == '(' and last_chars[-1] == ')': #duplicate, don't include
		pass
	elif 'ABE(' in feature_name:
		pass
	else:
		if feature_name not in snapgene_feature_list_cleaned:
			snapgene_feature_list_cleaned.append(feature_name)


snapgene_feature_list_cleaned = sorted(snapgene_feature_list_cleaned, key=str.casefold)
df = pd.DataFrame(snapgene_feature_list_cleaned, columns=['element'])
csv_data = df.to_csv(index=False)
df.to_csv('./snapgene_feature_list.csv', index=False)
