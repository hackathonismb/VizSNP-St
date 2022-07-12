#!/usr/local/bin/python3
from argparse import ArgumentParser
from collections import defaultdict
import sys
from pysam import VariantFile
import requests
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime

def get_protein_id(gene):
    """ get SWISSPROT ID"""

    URL = 'https://rest.uniprot.org/idmapping'

    params = {
    'from': 'Ensembl',
    'to': 'UniProtKB-Swiss-Prot',
    'ids': gene
    }
    response = requests.post(f'{URL}/run', params)
    job_id = response.json()['jobId']
    job_status = requests.get(f'{URL}/status/{job_id}')
    d = job_status.json()
    # List to store the mapped ID
    #SwissProt_IDs = []

    # Make three attemps to get the results
    # for i in range(3):
    if d.get("job_status") == 'FINISHED' or d.get('results'):
        job_results = requests.get(f'{URL}/results/{job_id}')
        results = job_results.json()
        for obj in results['results']:
            Swisprot_id=(obj["to"])
            #print(Swisprot_id)
            return Swisprot_id
        #break

    time.sleep(1)


def get_gene_info(Swiss_id, info):
    """get gene chr, coordinate and protein id info from ensmble"""
    # genes = gene.split(",")
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/" + Swiss_id + "?expand=1"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    info[Swiss_id]["id"] = decoded["Transcript"][0]["Translation"]["id"]
    info[Swiss_id]["chr"] = decoded["Transcript"][0]["seq_region_name"]
    info[Swiss_id]["start"] = decoded["Transcript"][0]["Translation"]["start"]
    info[Swiss_id]["end"] = decoded["Transcript"][0]["Translation"]["end"]
    print(info)


def extract_vcf(vcff, info):
    """ extract vcf file variants that match the given gene"""
    vcf_reader = VariantFile(vcff)
    data = ""
    for key in info.keys():
        for record in vcf_reader.fetch(info[key]["chr"], info[key]["start"], info[key]["end"]):
            line = str(record)
            dlist = line.split()
            data += "\"" + " ".join(dlist) + "\" ,"
    return data

def muttaster_output(data):
    """ Run Mutation Taster and capture """
    url = 'https://www.genecascade.org/MT2021/MT_API102.cgi?variants='

    variants = data.replace("\"","").split(",")
    for i in variants:
        var=i.split(" ")
        if(len(var)>4):
            print(var[0],var[1],var[3],var[4])
            variant = str(var[0]) + ":" + str(var[1]) + var[3] + '%3E' + var[4]
            snpurl = url + variant
            req = urllib.request.Request(snpurl)
            with urllib.request.urlopen(req) as f:
                response = f.read()
                print(response.decode('utf-8'))

def vep_output(data, iCn3D_sift, iCn3D_polyphen):
    """ Run VEP with the identified variants and capture sift and polyphen scores"""

    server = "https://rest.ensembl.org"
    ext = "/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    vcf_lines = "{\"variants\" : [" + data[:-1] + "]}"
    r = requests.post(server + ext, headers=headers, data=vcf_lines)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    for i in decoded:
        for key,val in i.items():
             if isinstance(i[key], list):
                 for j in i[key]:
                     if "sift_prediction" in j or "polyphen_prediction" in j:
                            alt_aa = j.get("amino_acids").split("/")[1]
                            iCn3D_sift[j.get("gene_id")][j.get("protein_start")][alt_aa] = {
                                     "sift_prediction": j.get("sift_prediction"), "sift_score": j["sift_score"]}
                            iCn3D_polyphen[j.get("gene_id")][j.get("protein_start")][alt_aa] = {
                                     "polyphen_prediction": j.get("polyphen_prediction"), "polyphen_score": j["polyphen_score"]}

def get_iCn3D_path(sift, polyphen,Swisprot_id):
    """ generates the iCn3D path based on the variants"""
    date = datetime.now()
    url='https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?'

    # print("\n")
    # for SwissProt_ID in pid:
    print("UniProt Primary Accession:",Swisprot_id)
    iCn3Durl =  url + 'afid=' + Swisprot_id + '&date=' + date.strftime("%Y%m%d") + '&v=3.11.5&command=view annotations; set annotation cdd; '
    sift_str = variant_string(sift)
    polyphen_str = variant_string(polyphen)
    if(sift_str):
        iCn3Durl += 'add track | chainid ' + Swisprot_id + '_A | title SIFT_predict | text ' + sift_str
    if(polyphen_str):
        iCn3Durl += '; add track | chainid ' + Swisprot_id + '_A | title PolyPhen_predict | text ' + polyphen_str

    print("Here is your iCn3D link:")
    print(iCn3Durl, "\n")

def variant_string(predict):
    """ extract the variants that are deleterious from sift and polyphen dicts and returns a combined string per prediction"""
    variants = ""
    for gene in predict.keys():
        for pos in predict[gene].keys():
            for aa in predict[gene][pos].keys():
                if 'sift_prediction' in predict[gene][pos][aa]:
                    if predict[gene][pos][aa]['sift_prediction'] == 'deleterious':
                        variants += str(pos) + ' ' + aa + ','
                elif 'polyphen_prediction' in predict[gene][pos][aa]:
                    if predict[gene][pos][aa]['polyphen_prediction'] == 'benign':
                        variants += str(pos) + ' ' + aa + ','

    return variants

def cli():

    desc = ("PredVEP2iCn3D.py: Run VEP on certain gene variants extracted from a VCF file \
             and generates an iCn3D link")
    parser = ArgumentParser(description=desc)
    parser.add_argument('-g', metavar='GENE', help="Ensemble Gene or Transcript IDs seperated by comma")
    parser.add_argument('-v', metavar='VCF', help="VCF file to extract the variants; gziped and tabix")
    parser.add_argument('-o', metavar='OUTPUT', help="Output file with iCn3D link")
    return parser

def main(args):
    info = defaultdict(dict)
    iCn3D_sift = defaultdict(lambda: defaultdict(dict))  # dict to hold sift scores per variants
    iCn3D_polyphen = defaultdict(lambda: defaultdict(dict))  # dict to hold sift scores per variants
    Ensembl_ids=args.g.split(",")

    for gen in Ensembl_ids:
        pid = get_protein_id(gen)
        print("Swiss prot id is " + pid)
        get_gene_info(gen, info)
        #print(inf)
        #print(pid)
        # get_iCn3D_path(iCn3D_sift, iCn3D_polyphen,pid)
        #
        # get_gene_info(args.g, info)
    #vep_output(data, iCn3D_sift, iCn3D_polyphen)

    # loop over the gene names
    #print(len(args.g.split(",")))
    # ids=args.g.split(",")
    # for gen in ids:
    #     pid = get_protein_id(gen)
    #     get_iCn3D_path(iCn3D_sift, iCn3D_polyphen,pid)

if __name__ == '__main__':
    main(cli().parse_args())
