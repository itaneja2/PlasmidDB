How to format SnapGene files so they can be uploaded to DB
==============================


### File naming convention

1. The general format is:

		pLLPLASMIDNUM_VECTOR_SPECIESTYPE
		Example: pLL323_pBADHisA_B.dna

	* `pLLPLASMIDNUM`: pLL323		
	* `VECTOR`: Name of backbone vector 
	* `SPECIESTYPE`: B (bacterial) or M (mammalian)
	* Note 1: if a vector or species type is not present/known/applicable, use `NA`. For instance: `pLL323_pBADHisA_NA.dna`
	* Note 2: if you want to provide more information about your plasmid, you are totally free to do that as long as the first three components (plasmid number, vector, species type) are present. For instance: `pLL323_pBADHisA_B_XYZ_123.dna`

	
### Creating SnapGene files

1. Add common features. 
2. If you have any custom features in your plasmid that you want to be documented in the database, make sure it is separated by **spaces**, not `_`.     
3. Make sure no common features (promoters, terminators, antibiotics, inducers, oris, etc.) overlap in nucleotide sequence. If there are, delete them. Even if you don't remember to do this, the upload program will catch it.  
	* Note: It is OK if features overlap if they do not have a corresponding type annotated (either explicitly via FP/TAG/CLV/LKR or through the `feature_type_mapping` table). 
	* What if you want to keep overlapping features? Let's say you have a feature that is of type `protein_bind` and it overlaps with a feature that is of type `promoter`. If you change that feature type to `misc_feature` in the feature_type_mapping table, then you can upload the file with overlapping features. 
4. Annotate your insert by specifying a unique type for each fluorescent protein (FP), tag (TAG), cleavage (CLV), linker (LKR), or protein of interest (POI) followed by `_FEATURENAME`. So for instance, `FP_mCherry`. Example annotations include: 

	* `FP_mCherry`: mCherry, eYFP, sfGFP, mMaple3, etc.. Make sure the FP name follows the naming convention from [fpbase](https://fpbase.org)
	* `TAG `: his tag, solubility tag, etc...
	* `CLV`: cleavage site
	* `LKR`: GS, PAG, etc.
	* `POI`: genusname-geneName-WT or genusname-geneName-mutantinfo. Examples: Make sure you use the proper geneName specified by Uniprot. Note genus name is specified as first initial followed by the remaining word (e.g Ecoli). 
5. Make sure the last feature in your insert (e.g POI) contains a stop codon in its sequence. In other words, the feature in the snapgene file should contain the stop codon. **In order for SnapGene to detect the stop codon, make sure the last feature includes the stop codon. This may require you to edit the feature and extend the location of the feature by three nucleotides.**
6. For each feature in the insert, edit the feature and 1) specify its type as CDS and 2) click the checkbox `Translate this feature in Sequence view`.
7. [TROUBLESHOOTING] Make sure the following properties hold:
	* For each insert, there should be a promoter numerically BEFORE it and a terminator numerically AFTER it. More precisely, the end base pair of the promoter should be less than the start base pair of the insert and the start base pair of the terminator should be greater than the end base pair of the insert. 
		* One potential solution is to flip the sequence `View>Flip Sequence>Save`. 
	
	* The program will complain if for a given feature its end nucleotide number is greater than its start nucleotide number. To the adjust the starting position of the plasmid in SnapGene, select the corresponding origin of replication feature and click `View>Set Origin` to set that point as base 1.  
	* No common features overlap 


	
### Making sure your features are detected by the upload program 

If you have a feature that should be documented in the database, you will need to explicitly specify the type of that feature **in the database**. The good news is you only have to do this once for a given feature, and there is a good chance someone may have already annotated it. This bad news is if you forget to do this, the program may complain saying that it can't find your insert because it can't detect your promoter or terminator. (For reference, the type of each feature in the insert is already implicitly annotated because of the prefix such as FP or LKR). 


* To specifcy the type of each feature go to `feature_type_mapping` table. The list of valid types is in the `type_names` table. 
	* **IF YOUR FEATURE IS NOT PRESENT**: You can add a new feature via the `Add Row` (top-left). Make sure any custom feature you add is separated by **spaces**, not `_` because thats how snapgene does its default features. 
	* To specify the type, go to the `type_names` table, find the corresponding type, and copy paste the **exact type** into `type` column in `feature_type_mapping`
		* If your type is not present in `type_names`, add a new `Type` and specify the source as `Custom` in the `type_names` table. 
	* To finish, write your initials under the `annotator` column if you are adding a new feature-type mapping. 

