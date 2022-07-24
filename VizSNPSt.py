#!/usr/local/bin/python3
from argparse import ArgumentParser
from collections import defaultdict
import sys
from pysam import VariantFile
import requests
import json
import time
from datetime import datetime
import webbrowser

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

    # Make three attemps to get the results
    for i in range(3):
        if d.get("job_status") == 'FINISHED' or d.get('results'):
            job_results = requests.get(f'{URL}/results/{job_id}')
            results = job_results.json()
            for obj in results['results']:
                SwissProt_ID = obj["to"]
            break
        time.sleep(1)
    return SwissProt_ID

def get_gene_info(gene, info):

    """get gene chr, coordinate and protein id info from ensmble"""

    server = "https://rest.ensembl.org"
    ext = "/lookup/id/" + gene + "?expand=1"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    info[gene]["id"] = decoded["Transcript"][0]["Translation"]["id"]
    info[gene]["chr"] = decoded["Transcript"][0]["seq_region_name"]
    info[gene]["start"] = decoded["Transcript"][0]["Translation"]["start"]
    info[gene]["end"] = decoded["Transcript"][0]["Translation"]["end"]

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
        for key in i:
             if isinstance(i[key], list):
                 for j in i[key]:
                     if "sift_prediction" in j or "polyphen_prediction" in j:
                            alt_aa = j.get("amino_acids").split("/")[1]
                            iCn3D_sift[j.get("gene_id")][j.get("protein_start")][alt_aa] = {
                                     "sift_prediction": j.get("sift_prediction"), "sift_score": j["sift_score"]}
                            iCn3D_polyphen[j.get("gene_id")][j.get("protein_start")][alt_aa] = {
                                     "polyphen_prediction": j.get("polyphen_prediction"), "polyphen_score": j["polyphen_score"]}

def variant_string(predict):
    """ extract the variants that are deleterious from sift and polyphen dicts and returns a combined string per prediction"""

    variants = ""
    for gene in predict.keys():
        for pos in predict[gene]:
            for aa in predict[gene][pos]:
                if 'sift_prediction' in predict[gene][pos][aa]:
                    if predict[gene][pos][aa]['sift_prediction'] == 'deleterious':
                        if variants == '':
                            variants += str(pos) + ' ' + aa
                        else:
                            variants += str(pos) + ' ' + aa + ','
                elif 'polyphen_prediction' in predict[gene][pos][aa]:
                    if predict[gene][pos][aa]['polyphen_prediction'] == 'deleterious':
                        if variants == '':
                            variants += str(pos) + ' ' + aa
                        else:
                            variants += str(pos) + ' ' + aa + ','
    
    return variants


def get_iCn3D_path(sift, polyphen, pid):
    """ generates the iCn3D path based on the variants"""
    date = datetime.now()
    url='https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?'

    print("\nUniProt Primary Accession:", pid, "\n")
    iCn3Durl =  url + 'afid=' + pid + '&date=' + date.strftime("%Y%m%d") + '&v=3.11.5&command=view annotations; set annotation cdd; set view detailed view; '
    sift_str = variant_string(sift)
    polyphen_str = variant_string(polyphen)
    if(sift_str):
        scap_str = 'scap interaction ' + pid + '_A_' + sift_str.replace(" ", "_")
        iCn3Durl += 'add track | chainid ' + pid + '_A | title SIFT_predict | text ' + sift_str + ";" + scap_str
    if(polyphen_str):
        # only need to add scap_str once, need to check if sift_str exists already
        if not scap_str:
            iCn3Durl += '; add track | chainid ' + pid + '_A | title PolyPhen_predict | text ' + polyphen_str
        else:
            scap_str = 'scap interaction ' + pid + '_A_' + sift_str.replace(" ", "_")
            iCn3Durl += '; add track | chainid ' + pid + '_A | title PolyPhen_predict | text ' + polyphen_str + ";" + scap_str

    print("Here is your iCn3D link:")
    print(iCn3Durl, "\n")
    webbrowser.open(iCn3Durl)

def cli():

    desc = ("PredVEP2iCn3D.py: Run VEP on certain gene variants extracted from a VCF file \
             and generates an iCn3D link")
    parser = ArgumentParser(description=desc)
    parser.add_argument('-g', metavar='GENE', help="Ensemble Gene ID of interest")
    parser.add_argument('-v', metavar='VCF', help="VCF file to extract the variants; gziped and tabix")
    return parser

def main(args):
    info = defaultdict(dict)
    iCn3D_sift = defaultdict(lambda: defaultdict(dict))  # dict to hold sift scores per variants
    iCn3D_polyphen = defaultdict(lambda: defaultdict(dict))  # dict to hold sift scores per variants

    get_gene_info(args.g, info)
    data = extract_vcf(args.v, info)
    vep_output(data, iCn3D_sift, iCn3D_polyphen)
    pid = get_protein_id(args.g)
    get_iCn3D_path(iCn3D_sift, iCn3D_polyphen, pid)

if __name__ == '__main__':
    main(cli().parse_args())    
