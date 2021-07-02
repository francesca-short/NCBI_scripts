#!/usr/bin/env python
from Bio import SeqIO
from Bio import SeqFeature
from Bio import Entrez
import argparse
import csv

#Prepare input arguments
parser = argparse.ArgumentParser(description=
                                 '''This script takes a list of non-redundant Refseq protein IDs 
                                 (for example, hits from a BLAST search against RefSeq select), 
                                 looks up the identical proteins in NCBI and writes the flanking region
                                  for one representative of each identical protein group to a genbank file. 
                                  Useful for finding genomic context of genes similar to a gene of interest. 
                                  Requires list of protein IDs and the length of surrounding DNA to retrieve
                                   (default = 5kb). Also your email address for NCBI download management. Note 
                                   this script doesn't do any sanity checking, so it will fail if you try to 
                                   search for something that isn't an ID.
            						''')
parser.add_argument('email', help = 'Your email address')
parser.add_argument('prot_ID_list', help="Protein ID list .txt file")
parser.add_argument('-n', '--number_bases', help = 'number of bases upstream and downstream to retrieve', type = int)
parser.add_argument('-o', '--output', help="Output csv file prefix")

args = parser.parse_args()
if args.number_bases == None:
	args.number_bases = 5000
if args.output == None:
	args.output = args.prot_ID_list

#Define function to get identical protein info from NCBI for a given refseq protein ID
def get_ipg_search_info(ID):
	print(ID)
	Entrez.email=args.email
    # Downloading...
	net_handle = Entrez.efetch(db="protein", id=ID, rettype = "ipg", retmode = "xml")
	IPG_info = Entrez.read(net_handle, validate = True)
	return(IPG_info)

#Read in list of protein IDs
file_handle = open(args.prot_ID_list, 'r')
prot_list = set()
for protein in file_handle.readlines():
	prot_list.add(protein.rstrip())

print("Will search for identical protein groups for:\n", prot_list)

#Make list to hold identical protein group info
prot_info_list = []

#Loop over list of protein IDs and retrieve info from NCBI
for ID in prot_list:
	prot_info_list.append(get_ipg_search_info(ID))

#Function to get relevant info on identical proteins
def extract_protein_info(item):
	if "IPGReport" not in item.keys():
		print("Could not find info for ", item)
	else:
		print("Info for protein ", item['IPGReport']['Product'].attributes["accver"])
		protein_to_get = item['IPGReport']['ProteinList'][0]['CDSList'][0]
		protein_to_get_info = protein_to_get.attributes
		protein_to_get_info['Query_prot_ID'] = item['IPGReport']['Product'].attributes["accver"]
		protein_to_get_info['Source'] = item['IPGReport']['ProteinList'][0].attributes["source"]
		#{'Query_prot_ID':item['IPGReport']['Product'].attributes["accver"], 'Start_nt':protein_to_get.attributes.get("start"), 'End_nt':protein_to_get.attributes.get("stop"), 'Accession':protein_to_get.attributes.get("accver"), '}
		#print("Found search info for identical protein in assembly " + protein_to_get.attributes.get("assembly") +", " +protein_to_get.attributes.get("org"), protein_to_get.attributes.get("strain"))
		print(protein_to_get_info)
		return(protein_to_get_info)

#Get identical protein info for next search
new_search_info = []
for item in prot_info_list:
	new_search_info.append(extract_protein_info(item))

#Write useful stuff to csv
filename = (args.output + '.csv')
keys = new_search_info[0].keys()
with open(filename, 'w') as outfile:
	dict_writer = csv.DictWriter(outfile, keys)
	dict_writer.writeheader()
	dict_writer.writerows(new_search_info)

#Batch query mode
accession_list = []
for item in new_search_info:
	accession_list.append(item['accver'])
	
#Get genbank info for identical protein - single query mode, only downloads defined start and end bits
for item in new_search_info:
	print("Retrieving ", item['accver'], item["Query_prot_ID"], " flanking region from NCBI.")
	try:
		start = int(item['start']) - args.number_bases
		stop = int(item['stop']) + args.number_bases
		handle = Entrez.efetch(db="nucleotide", id = item['accver'], rettype = 'gbwithparts', retmode='text', seq_start=start, seq_stop=stop)
		record = SeqIO.read(handle, "genbank")
		handle.close()
		filename = (record.id + "_" + str(start) + "-" + str(stop) + ".gb")
		SeqIO.write(record, filename, 'gb')
		print("Flanking region written to ", filename)	
	except:
		print("An error occurred.")

