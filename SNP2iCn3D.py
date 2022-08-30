#!/usr/local/bin/python3

'''
SNP2iCn3D.py
Reads in a VCF file, submits variants to VEP server, generates an iCn3D link that 
shows the variants in the sequence track and highlights the first mutant in the 3D viewer

Original script: Shashi Ratnayake, NCI
Modifications by: Michael Sierk, NCI, Manoj Wagle, University of Grenoble Alpes

TODO (8/8/22):
    - html output
        - get variant string
    - .tsv output for importing into a spreadsheet
    - retrieve PDB if available
    - load SIFT/PolyPhen scores into iCn3D
        - requires BED file, has to be loaded manually into iCn3D?

'''
from argparse import ArgumentParser
from collections import defaultdict
import sys
import re
from pysam import VariantFile
import requests
import json
import time
from datetime import datetime
import webbrowser
from urllib.parse import quote, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError


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
        ''' get the Ensembl gene ID for each coordinate location from the VCF file'''
        genes_list = defaultdict(dict)
        print("Getting Ensembl gene IDs and SwissProt IDs...")
        for c in coords:
            loc = c[0]+":"+c[1] # need to account for same chr, same coord, different base (e.g. A->T vs A->C)
            #ext = "/overlap/region/human/7:140424943-140624564?feature=gene"
            ext = "/overlap/region/human/" + c[0] + ":" + c[1] + "-" + c[1]
            gene = self.perform_rest_action(
                endpoint=ext, 
                params={'feature': 'gene'}
            )
            if gene:
                #print(gene)
                for g in gene:
                    # get the Uniprot ID here 
                    genes_list[loc]["ens_id"] = g["gene_id"]
                    genes_list[loc]["sp_id"] = get_protein_id(g["gene_id"])
                    #print(genes_list[loc]["sp_id"], "\n")
            else:
                print("no gene found for " + loc)
        return genes_list

    def get_gene_info(self, genes, genes_info):

        """
        NOT USED
        get gene chr, coordinate and protein id info from Ensembl
            - only used if pulling specific coords for a gene from VCF file
        """
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
        '''attempt to get VEP output via EnsemblRestClient class
         - problem with encoding when submitting variants string'''

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
        #print(json.dumps(r, indent=2))

        return r.json()

def extract_vcf(vcf, gene_ids):
    """ extract vcf variants that match a list of given genes"""
    variants = ""
    for l in gene_ids.keys():
        for c in vcf:
            loc = c[0]+":"+c[1]
            if l == loc:
                variants += "\"" + " ".join(c[0:5]) + " . . ." + "\" ," 
                #print(variants) # "19 15256965 . T G . . ." ,
    return variants[:-1]

def get_protein_id(gene):
    """ get Uniprot ID using Ensembl gene ID
    https://www.biostars.org/p/9529129/#9529154
    """
    #print("getting Uniprot ID for", gene)

    URL = 'https://rest.uniprot.org/idmapping'

    params = {
        'from': 'Ensembl',
        #'to': 'UniProtKB',
        'to': 'UniProtKB-Swiss-Prot',
        'ids': gene
    }

    response = requests.post(f'{URL}/run', params)
    #print(response)

    job_id = response.json()['jobId']
    #print('job id:', job_id)
    job_status = requests.get(f'{URL}/status/{job_id}')
    d = job_status.json()

    SwissProt_ID = None
    # Make three attemps to get the results
    for i in range(3):
        #print(d.get('jobStatus'))
        if d.get("jobStatus") == 'FINISHED' or d.get('results'):
            job_results = requests.get(f'{URL}/results/{job_id}')
            results = job_results.json()
            #print(json.dumps(results, indent=2))
            for obj in results['results']:
                SwissProt_ID = obj["to"]
            break
        time.sleep(1)
    return SwissProt_ID

def vep_output(variants, sift, polyphen):
    #from difflib import SequenceMatcher
    """ Run VEP with the identified variants and capture sift and polyphen scores"""
    #client = EnsemblRestClient()
    print("\nRetrieving VEP results...\n")
    
    server = "https://rest.ensembl.org"
    ext = "/vep/homo_sapiens/region?uniprot=1"
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
                            aa = j.get("amino_acids").split("/")[0]
                            aa = aa.strip('\s+')
                            aanum = str(j.get("protein_start"))
                            alt_aa = j.get("amino_acids").split("/")[1]
                            alt_aa = alt_aa.strip()
                            mutaa = aa + aanum + alt_aa # e.g. T355C
                            sift[j.get("gene_id")][mutaa] = {
                                     "sift_prediction": j.get("sift_prediction"), "sift_score": j["sift_score"]}
                            polyphen[j.get("gene_id")][mutaa] = {
                                     "polyphen_prediction": j.get("polyphen_prediction"), "polyphen_score": j["polyphen_score"]}

def get_pdb_id(pid):
    print("get_pdb_id for " + pid)
    url = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/" + pid 
    r = requests.get(url)
    #print(r.json())
    if r.status_code == 200:
        return(r.json())
    else:
        return("None")

def get_iCn3D_path(sift, polyphen, gene, pid):
    """ generates the iCn3D path based on the variants"""
    date = datetime.now()
    url_path='https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?'
    #url_query=''
    #url_command=''
    print("UniProt Primary Accession:", pid, "\n")

    # see if PDB entry exists for pid
    # problem: may need different PDB IDs for different locations in uniprotID
    #  - need the Uniprot residue #, which comes from VEP output
    #  - VEP output is stored in sift, polyphen dicts
    #  - variant_string returns all variants for a gene
    # 1) All variants in 1 PDB
    # 2) Some variants in 1 PDB, rest in Uniprot ID (alphafold)
    # 3) Different variants in different PDBs

    sift_str = variant_string(sift[gene]) # need to modify to produce BED file to show SIFT/Polyphen score
    polyphen_str = variant_string(polyphen[gene])

    pdbid_list = get_pdb_id(pid) 
    print("pdbid_list:\n", pdbid_list)
    #{"P07550":[{"end":64,"entity_id":1,"chain_id":"A","pdb_id":"6ps2","start":54,"unp_end":40,"coverage":1.0,"unp_start":30,"resolution":2.4,"experimental_method":"X-ray diffraction","tax_id":9606},
   
    # fill all the locations in struct_id
    struct_id = dict() 
    if (sift_str): 
        sift_split = sift_str.split(", ") 
        for s in sift_split:
            print('s:', s)
            #s = s.strip()
            aanum = s.split()[0]
            print("aanum: ", aanum)
            struct_id[aanum] = ''

    if (polyphen_str):
        polyph_split = polyphen_str.split(", ")
        for p in polyph_split:
            print('p:', p)
            #p = p.strip()
            aanum = p.split()[0]
            print("aanum: ", aanum, sep='')
            struct_id[aanum] = ''

    # determine PDB/AFID(default) for each aanum
    if (struct_id):
        for n in struct_id:
            print('n:', n)
            struct_id[n] = pid
            if (pdbid_list != 'None'):
                for m in pdbid_list:
                    print('m:', m)
                    #print('start:', m.get('unp_start'), 'end:', m.get('unp_end'))
            
    if (pdbid_list != 'None'):
        for m in pdbid_list:
            print('m:', m)
            for j in pdbid_list[m]:
                #if isinstance(m[id], list):
                #for j in m[id]:
                    #if "sift_prediction" in j or "polyphen_prediction" in j:
                id = j["pdb_id"]
                chain = j["chain_id"]
                resolution = j["resolution"]
                start = j["unp_start"]
                end = j['unp_end']
                print('id:', id, 'chain:', chain, 'res:', resolution, 'start:', start, 'end:', end)
                #print('start:', m.get('unp_start'), 'end:', m.get('unp_end'))


    url_query =  'afid=' + pid + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command='
    url_command='view annotations; set annotation cdd; set view detailed view;'    

    if(sift_str):
        # Note: scap interaction only displays 1st snp in list
        scap_str = 'scap interaction ' + pid + '_A_' + sift_str.replace(" ", "_")
        url_command += 'add track | chainid ' + pid + '_A | title SIFT_predict | text ' + sift_str + ";" + scap_str
    
    if(polyphen_str):
        # only need to add scap_str once, need to check if scap_str exists already
        if scap_str:
            url_command += '; add track | chainid ' + pid + '_A | title PolyPhen_predict | text ' + polyphen_str
        else:
            scap_str = 'scap interaction ' + pid + '_A_' + polyphen_str.replace(" ", "_")
            url_command += '; add track | chainid ' + pid + '_A | title PolyPhen_predict | text ' + polyphen_str + ";" + scap_str
    
    url_command = quote(url_command) # encode the spaces
    iCn3Durl = url_path + url_query + url_command

    if (sift_str) or (polyphen_str):
        print("Here is your iCn3D link:")
        print(iCn3Durl, "\n")
        #iCn3Durl = quote(iCn3Durl, safe=)
    else:
        iCn3Durl = "No deleterious mutations found for " + pid
    
    return(iCn3Durl) #, sift_str, polyphen_str)
    #webbrowser.open(iCn3Durl)

def variant_string(predict):
    """ extract the variants that are deleterious from sift and polyphen dicts and returns 
    a combined string per prediction for the iCn3D url command
    - this includes all the variants from one gene 
    """
    variants = ""
    for mutaa in predict:
        s = re.split(r'(\d+)', mutaa) # need the number and new aa
        if 'sift_prediction' in predict[mutaa]:
            #print('sift pred:', predict[mutaa]['sift_prediction'])
            if predict[mutaa]['sift_prediction'] == 'deleterious':
                if variants == '':
                    variants += s[1] + ' ' + s[2]
                else:
                    variants += ',' + s[1] + ' ' + s[2]
        elif 'polyphen_prediction' in predict[mutaa]:
            #print('polyp pred:', predict[mutaa]['polyphen_prediction'])
            if predict[mutaa]['polyphen_prediction'] == 'probably_damaging':
                if variants == '':
                    variants += s[1] + ' ' + s[2]
                else:
                    variants += ',' + s[1] + ' ' + s[2]
    return variants

def variant_string2(predict):
    """ returns a variant string for html output
        - puts all variants for one gene into the string
    """
    variants = ''
    for mutaa in predict:
        if 'sift_prediction' in predict[mutaa]:
            #print('sift pred:', predict[pos][aa]['sift_prediction'])
            str_add = mutaa + ' ' + predict[mutaa]['sift_prediction'] + \
                        ' ' + str(predict[mutaa]['sift_score'])
            if variants == '':
                variants += str_add
            else:
                variants += "\n" + str_add
        elif 'polyphen_prediction' in predict[mutaa]:
            #print('polyp pred:', predict[pos][aa]['polyphen_prediction'])
            str_add = mutaa + ' ' + predict[mutaa]['polyphen_prediction'] + \
                    ' ' + str(predict[mutaa]['polyphen_score'])
            if variants == '':
                variants += str_add
            else:
                variants += "\n" + str_add
            
    return variants

def print_html(vcf_file, url_list, sift, polyphen, gene_to_pid):

    # to open/create a new html file in the write mode
    f = open('iCn3D_links.html', 'w')
  
    # html code 
    html = """<html>
    <head>
    <title>iCn3D links</title>
    </head>
    <body>
    <h2>Click on a link to open the variant in iCn3D</h2>
    <table border='1'>"""
    html += "Input file: " + vcf_file

    # output table
    html += "<tr><th>Gene</th><th>UniprotID</th><th>SIFT</th><th>PolyPhen</th><th>iCn3D link</th></tr>"
    for gene in gene_to_pid.keys():

        if re.match(r"^https", url_list[gene]):
            url = "<a href=" + ''.join(url_list[gene]) + " target=\"_blank\">iCn3D link</a>"
        else:
            url = url_list[gene] # text output, no link

        sif = variant_string2(sift[gene])
        pol = variant_string2(polyphen[gene])

        unid = "NA"
        if gene_to_pid[gene]:
            unid = gene_to_pid[gene]
        html += "<tr><td>" + gene + "</td><td>" + unid + "</td><td>" + sif + "</td><td>" + \
            pol + "</td><td>" + url + "</td></tr>"

    html += "</table>"
    html += "</body>"
    html += "</html>"

    f.write(html)
    f.close()
    
    webbrowser.open('iCn3D_links.html')

def get_vcf(vcff):
    # read in the VCF file
    vcf = []
    vcf_reader = VariantFile(vcff)

    for record in vcf_reader.fetch():
        line = str(record)
        vcf_line = line.split()
        # only get variants
        if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
            vcf_line[0] = vcf_line[0].replace("chr","")
            vcf.append(vcf_line[0:6])

    return(vcf)

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

def main(args):
    import dill as pickle

    # read in the whole VCF file, use client.get_gene_ids to fill genes dictionary: genes[loc]["ens_id"] and genes[loc]["sp_id"]
    # if args.g, select subset of genes
    # go through genes, get coords from vcf, get VEP, generate iCn3D link for each gene
    # TODO: put in a cli arg for html vs. .tsv output

    #gene_info = defaultdict(dict) # used when picking out coordinates of a known gene from VCF
    sift = defaultdict(lambda: defaultdict(dict))      # nested dict to hold sift/polyphen scores per variants
    polyphen = defaultdict(lambda: defaultdict(dict))  

    test = 0 # set = 1 to avoid regenerating everything if just testing html output
    if test == 0:
        vcf = get_vcf(args.v)
        client = EnsemblRestClient()
        gene_ids = client.get_gene_ids(vcf) # extract Ensemble gene and SwissProt ids from the vcf coordinates
 
        # change from Ensembl gene ID to Uniprot ID?
        # rearrange so args.g is submitted to get_gene_ids so not parsing vcf twice
        gene_ids_select = defaultdict(dict) # subset gene_ids if using -g flag
        if args.g:
            gene_id_list = [args.g]
            for l in gene_ids.keys():
                if gene_ids[l]['ens_id'] in gene_id_list:
                    gene_ids_select[l] = gene_ids[l] 
        else:
            gene_id_list = [gene_ids[x]['ens_id'] for x in gene_ids.keys()]
            gene_ids_select = gene_ids
    
        #client.get_gene_info(gene_id_list, gene_info)    # fills gene info dictionary

        # get list of variants (as a string) to submit to VEP
        variants = extract_vcf(vcf, gene_ids_select) 
        vep_output(variants, sift, polyphen) # also want to be able to read in vep results if in VCF file (e.g. from TCGA)
    
        #create dict to look up pid from ens_id
        gene_to_pid = defaultdict(dict)
        for loc in gene_ids.keys():
            if not gene_to_pid[gene_ids[loc]['ens_id']]:
                gene_to_pid[gene_ids[loc]['ens_id']] = gene_ids[loc]['sp_id']
                #print(gene_ids[loc]['ens_id'], gene_to_pid[gene_ids[loc]['ens_id']])


        # save data so we don't have to rerun
        pfile = open(r'iCn3D_data.pkl', 'wb')
        pickle.dump(gene_ids_select, pfile)
        pickle.dump(sift, pfile)
        pickle.dump(polyphen, pfile)
        pickle.dump(gene_to_pid, pfile)
        #pickle.dump(url_list, pfile)
        pfile.close()
    # end if test 

    else: 
        #reload object from file
        file = open(r'iCn3D_data.pkl', 'rb')
        gene_ids_select = pickle.load(file)
        #url_list = pickle.load(file)
        sift = pickle.load(file)
        polyphen = pickle.load(file)
        gene_to_pid = pickle.load(file)
        file.close()

    #gene_id_list = [gene_ids_select[x]['ens_id'] for x in gene_ids_select.keys()]
    url_list = defaultdict(dict)
    for gene in gene_id_list:
        print("=========================\nGetting link for ", gene)
        if (gene_to_pid[gene]):
            #print("exists")
            url_list[gene] = get_iCn3D_path(sift, polyphen, gene, gene_to_pid[gene])
        else:
            url_list[gene] = "No SwissProt ID for " + gene + "!"

    #generate an html page with all the iCn3D links
    print_html(args.v, url_list, sift, polyphen, gene_to_pid)

if __name__ == '__main__':
    main(cli().parse_args())    
