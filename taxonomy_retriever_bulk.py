#!/usr/bin/env python
from Bio import SeqIO
from Bio import SeqFeature
from Bio import Entrez
import argparse
import csv

#Prepare input arguments
parser = argparse.ArgumentParser(description=
                                 '''This script takes a list of taxIDs and reports their lineage info as csv.
            						''')
parser.add_argument('email', help = 'Your email address')
parser.add_argument('tax_ID_list', help="Taxonomy ID list .txt file")
parser.add_argument('-o', '--output', help="Output csv file prefix")

args = parser.parse_args()
if args.output == None:
	args.output = args.tax_ID_list

#Read in list of taxonomy IDs
file_handle = open(args.tax_ID_list, 'r')
taxID_list = []
for taxID in file_handle.readlines():
	taxID_list.append(taxID.rstrip())
print("Searching for lineage information for:\n", taxID_list)

#Alternative search for whole taxID list at once
Entrez.email = args.email
search_results = Entrez.read(Entrez.epost("taxonomy", id=",".join(taxID_list)))
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
fetch_handle = Entrez.efetch(db = "taxonomy", rettype = "xml", webenv=webenv, query_key=query_key)
tax_info = Entrez.read(fetch_handle, validate=True)

taxonomy_list = []

for item in tax_info:
	taxonomy_list.append([item['TaxId'], item['LineageEx']])

#Extract useful info to dictionary
for item in taxonomy_list:
	print(item)

#Define function to get scientific names for lineage
def search(rank, lineage):
	for item in lineage:
		if item['Rank'] == rank:
			return(item['ScientificName'])

#get info
taxonomy_info = []
for item in taxonomy_list:
	TaxID = item[0]
	Phylum = search('phylum', item[1])
	Class = search('class', item[1])
	Order = search('order', item[1])
	Family = search('family', item[1])
	Genus = search('genus', item[1])
	Species = search('species', item[1])
	info_to_keep = [TaxID, Phylum, Class, Order, Family, Genus, Species]
	taxonomy_info.append(info_to_keep)

#Write useful stuff to csv
filename = (args.output + '.csv')
header = ("TaxID","Phylum","Class","Order","Family","Genus","Species")
with open(filename, "w", newline="") as f:
	writer = csv.writer(f)
	writer.writerow(header)
	writer.writerows(taxonomy_info)