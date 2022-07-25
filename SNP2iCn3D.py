#!/usr/local/bin/python3

'''
SNP2iCn3D.py
Reads in a VCF file, submits variants to VEP server, generates an iCn3D link that 
shows the variants in the sequence track and highlights the first mutant in the 3D viewer

Original script: Shashi Ratnayake, NCI
Modifications by: Michael Sierk, NCI, Manoj Wagle, University of Grenoble Alpes

'''
from argparse import ArgumentParser
from collections import defaultdict
import sys
from pysam import VariantFile
import requests
import json
import time
from datetime import datetime
import webbrowser
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError

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

    SwissProt_ID = None
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

class EnsemblRestClient(object):
    # from https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
           
        return data

    def get_gene_ids(self, coords):

        genes_list = defaultdict(dict)

        for c in coords:
            loc = c[0]+":"+c[1]
            #print("loc: ", loc)
            #ext = "/overlap/region/human/7:140424943-140624564?feature=gene"
            ext = "/overlap/region/human/" + c[0] + ":" + c[1] + "-" + c[1]
            gene = self.perform_rest_action(
                endpoint=ext, 
                params={'feature': 'gene'}
            )
            #print(type(gene))
            if gene:
                for g in gene:
                    print("retrieved gene id: ", g["gene_id"])
                    # get the SwissProt ID here (genes_list[loc]["ens_id"] and genes_list[loc]["sp_id"])
                    genes_list[loc]["ens_id"] = g["gene_id"]
                    genes_list[loc]["sp_id"] = get_protein_id(g["gene_id"])
        return genes_list

    def get_gene_info(self, genes, genes_info):

        """get gene chr, coordinate and protein id info from Ensembl"""
        for g in genes:
            print("g: ", g)
            ext = "/lookup/id/" + g
            g_info = self.perform_rest_action(    
                endpoint=ext, 
                params={'expand': '1'}
            )
            genes_info[g]["id"] = g_info["Transcript"][0]["Translation"]["id"] # not used
            genes_info[g]["chr"] = g_info["Transcript"][0]["seq_region_name"]
            genes_info[g]["start"] = g_info["Transcript"][0]["Translation"]["start"]
            genes_info[g]["end"] = g_info["Transcript"][0]["Translation"]["end"]
    
    def get_vep(self, variants):
        ext = "/vep/homo_sapiens/region"
        hdrs = {}
        hdrs['Accept'] = "application/json"

        vcf_lines = "{\"variants\" : [" + variants + "]}"
        print(vcf_lines)

        r = self.perform_rest_action(
            endpoint=ext,
            hdrs=hdrs,
            params = {'data': vcf_lines} # doesn't work with urlencode, spaces get converted to +
        )
        print(json.dumps(r, indent=2))

        return r.json()

def extract_vcf(vcf, gene_ids):
    """ extract vcf file variants that match the given gene"""
    #data = [] #defaultdict(dict)
    variants = ""
    for l in gene_ids.keys():
        # if gene_ids[l]['ens_id'] == gene:
        for c in vcf:
            loc = c[0]+":"+c[1]
            if l == loc:
                variants += "\"" + " ".join(c[0:5]) + " . . ." + "\" ," 
                #print(variants) # "19 15256965 . T G . . ." ,
        #data = variants
    #print(variants)
    return variants[:-1]

def vep_output(variants, iCn3D_sift, iCn3D_polyphen):
    from difflib import SequenceMatcher
    """ Run VEP with the identified variants and capture sift and polyphen scores"""
    #client = EnsemblRestClient()
    print("\nRetrieving VEP results...\n")
    
    server = "https://rest.ensembl.org"
    ext = "/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    vcf_lines = "{\"variants\" : [" + variants + "]}"
    #print(vcf_lines)
    #test = '{"variants" : ["1 930963 . G A . . ." ,"1 3625748 rs200029021 T G . . ." ,"1 8951432 rs201002216 C G . . ." ,"1 11858106 rs115986244 C T . . ." ,"1 12725716 rs774509385 T C . . ." ,"1 12859409 rs200583219 C T . . ." ,"1 34905012 . C A . . ." ,"1 56744276 . A C . . ." ,"1 85090058 . C A . . ." ,"1 102996069 rs764459522 C T . . ." ]}'
    #print(test)
    #if test == vcf_lines:
    #    print("test matched vcf_lines")
    #else:
    #    print('match: ', SequenceMatcher(a=vcf_lines, b=test).ratio())
    r = requests.post(server + ext, headers=headers, data=vcf_lines)
    #r = requests.post(server + ext, headers=headers, data='{"variants" : ["1 930963 . G A . . .", "1 3625748 rs200029021 T G . . .", "1 8951432 rs201002216 C G . . .","1 11858106 rs115986244 C T . . .","1 12725716 rs774509385 T C . . .","1 12859409 rs200583219 C T . . .","1 34905012 . C A . . .","1 56744276 . A C . . .","1 85090058 . C A . . ." ,"1 102996069 rs764459522 C T . . ."]}')
    #r = requests.post(server+ext, headers=headers, data='{ "variants" : ["21  26960070  rs116645811 G A . . .", "21  26965148  rs1135638 G A . . ." ] }')
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    #decoded = client.get_vep(variants)
    decoded = r.json()
    #print(json.dumps(decoded, indent=2))
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
    #print("iCn3D_sift: ", iCn3D_sift)
    #print("iCn3D_polyphen: ", iCn3D_polyphen)

def get_iCn3D_path(sift, polyphen, gene, pid):
    """ generates the iCn3D path based on the variants"""
    date = datetime.now()
    url='https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?'

    print("UniProt Primary Accession:", pid, "\n")
    iCn3Durl =  url + 'afid=' + pid + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command=view annotations; set annotation cdd; set view detailed view;'

    sift_str = variant_string(sift[gene]) # need to modify to produce BED file to show SIFT/Polyphen score
    polyphen_str = variant_string(polyphen[gene])

    if(sift_str):
        scap_str = 'scap interaction ' + pid + '_A_' + sift_str.replace(" ", "_")
        iCn3Durl += 'add track | chainid ' + pid + '_A | title SIFT_predict | text ' + sift_str + ";" + scap_str
    
    if(polyphen_str):
        # only need to add scap_str once, need to check if scap_str exists already
        if scap_str:
            iCn3Durl += '; add track | chainid ' + pid + '_A | title PolyPhen_predict | text ' + polyphen_str
        else:
            scap_str = 'scap interaction ' + pid + '_A_' + polyphen_str.replace(" ", "_")
            iCn3Durl += '; add track | chainid ' + pid + '_A | title PolyPhen_predict | text ' + polyphen_str + ";" + scap_str
    
    if (sift_str) or (polyphen_str):
        print("Here is your iCn3D link:")
        print(iCn3Durl, "\n")
    else:
        print ("No deleterious mutations found!")


    #webbrowser.open(iCn3Durl)

def variant_string(predict):
    """ extract the variants that are deleterious from sift and polyphen dicts and returns a combined string per prediction"""
    # this is putting all the variants into one string, regardless of gene
    variants = ""
    #for gene in predict.keys():
    for pos in predict:
        for aa in predict[pos]:
            if 'sift_prediction' in predict[pos][aa]:
                print('sift pred:', predict[pos][aa]['sift_prediction'])
                if predict[pos][aa]['sift_prediction'] == 'deleterious':
                    if variants == '':
                        variants += str(pos) + ' ' + aa
                    else:
                        variants += ',' + str(pos) + ' ' + aa
            elif 'polyphen_prediction' in predict[pos][aa]:
                print('polyp pred:', predict[pos][aa]['polyphen_prediction'])
                if predict[pos][aa]['polyphen_prediction'] == 'probably_damaging':
                    if variants == '':
                        variants += str(pos) + ' ' + aa 
                    else:
                        variants += ',' + str(pos) + ' ' + aa

    return variants

def cli():

    desc = ("PredVEP2iCn3D.py: Runs VEP on certain gene variants extracted from a VCF file \
             and generates iCn3D links for predicted deleterious mutations.  \
             Two modes: \
                1) provide specific Ensembl genes with the -g flag, \
                only locations matching those genes will be extracted; \
                2) do not provide the -g flag, all variants \
                in the VCF will be run through VEP.")
    parser = ArgumentParser(description=desc)
    parser.add_argument('-g', metavar='GENE', help="Ensemble Gene IDs of interest")
    parser.add_argument('-v', metavar='VCF', help="VCF file to extract the variants; must be compressed with bgzip \
        and the tabix .tbi index file must be present.")
    return parser

def get_vcf(vcff):
    # read in the VCF file
    vcf = []
    vcf_reader = VariantFile(vcff)

    for record in vcf_reader.fetch():
        line = str(record)
        vcf_line = line.split()
        # only get variants
        #print(vcf_line)
        if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
            vcf_line[0] = vcf_line[0].replace("chr","")
            vcf.append(vcf_line[0:6])

    return(vcf)

def main(args):
    # read in the whole VCF file, use client.get_gene_ids to fill genes dictionary: genes[loc]["ens_id"] and genes[loc]["sp_id"]
    # if args.g, select subset of genes
    # go through genes, get coords from vcf, get VEP, generate iCn3D link for each gene

    gene_info = defaultdict(dict)
    iCn3D_sift = defaultdict(lambda: defaultdict(dict))  # nested dict to hold sift/polyphen scores per variants
    iCn3D_polyphen = defaultdict(lambda: defaultdict(dict))  # need to be able to access 'gene' if not in dict

    vcf = get_vcf(args.v)
    client = EnsemblRestClient()
    gene_ids = client.get_gene_ids(vcf) # extract Ensemble gene and SwissProt ids from the vcf coordinates
 
    #client.get_gene_info(gene_id_list, gene_info)    # fills gene info dictionary
    gene_ids_select = defaultdict(dict) # subset gene_ids if using -g flag
    if args.g:
        gene_id_list = [args.g]
        for l in gene_ids.keys():
            if gene_ids[l]['ens_id'] in gene_id_list:
                gene_ids_select[l] = gene_ids[l] 
    else:
        gene_id_list = [gene_ids[x]['ens_id'] for x in gene_ids.keys()]
        gene_ids_select = gene_ids
    
    variants = extract_vcf(vcf, gene_ids_select) 
    vep_output(variants, iCn3D_sift, iCn3D_polyphen) # also want to be able to read in vep results if in VCF file (e.g. from TCGA)
    
    # create dict to look up pid from ens_id
    gene_to_pid = defaultdict(dict)
    for loc in gene_ids.keys():
        if not gene_to_pid[gene_ids[loc]['ens_id']]:
            gene_to_pid[gene_ids[loc]['ens_id']] = gene_ids[loc]['sp_id']

    for gene in gene_id_list:
        print("=========================\nGetting link for ", gene)
        if (gene_to_pid[gene]):
            get_iCn3D_path(iCn3D_sift, iCn3D_polyphen, gene, gene_to_pid[gene])
        else:
            print("No SwissProt ID!")
if __name__ == '__main__':
    main(cli().parse_args())    
