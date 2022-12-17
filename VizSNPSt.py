#!/usr/local/bin/python3

'''
SNP2iCn3D.py
Reads in a VCF file, submits variants to VEP server, generates an iCn3D link that 
shows the variants in the sequence track and highlights the first mutant in the 3D viewer

Original script: Shashi Ratnayake, CGBB, CBIIT, NCI
Modifications by: Michael Sierk, NCI 
                  Manoj M Wagle, Université Grenoble Alpes; Manipal Academy of Higher Education

TODO (12/1/22):
    - load SIFT/PolyPhen scores into iCn3D
        - requires BED file, has to be loaded manually into iCn3D?
    - option to include SIFT/Polyphen score cutoff
    - add other options for prediction, such as open cravat
'''
from argparse import ArgumentParser, HelpFormatter
import textwrap
from collections import defaultdict
import sys
import re
from pysam import VariantFile
import requests
import json
import time
from datetime import datetime
from urllib.parse import quote, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError
import pandas as pd
from math import nan


    
def get_protein_id(gene, rest):
    """ get Uniprot ID using Ensembl gene ID
    https://www.biostars.org/p/9529129/#9529154
    """
    #print("getting Uniprot ID for", gene)
    SwissProt_ID = None

    # use REST API (default; slower, sometimes is down)
    if rest:
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
    else: 
        # set up a local file to retrieve Ensembl->SwissProt mapping
        # can use mapping file: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
        # Use following shell commands (can put into a shell script and run it):
        #  grep Ensembl HUMAN_9606_idmapping.dat| grep -v "_" > tmp.dat
        #  cut -f3 tmp.dat | cut -d"." -f1 > ens
        #  cut -f1 tmp.dat > sp
        #  paste ens sp > HUMAN_9606_idmapping_Ensembl.dat
        #  rm tmp.dat ens sp

        # creates a two column file:
        #  % head HUMAN_9606_idmapping_Ensembl.dat
        #  ENSG00000126653	A0A024QZ33
        #  ENSG00000249915	A0A024QZ42
        #  ENSG00000170312	A0A024QZP7
        # TODO: Need to deal with multiple UniprotIDs for 1 Ensembl ID.
        with open('HUMAN_9606_idmapping_Ensembl.dat') as f:
            ens_to_swissprot = dict([line.split() for line in f])
        SwissProt_ID = ens_to_swissprot[gene]
        #ENSG00000131686	P23280
        #ENSG00000131686	Q8N4G4

    return SwissProt_ID

def vep_output(variants, args):
    """ Run VEP with the identified variants and capture sift and polyphen scores"""
    
    species = args.s
    print("species:", species)
    server = "https://rest.ensembl.org"
    ext = "/vep/" + species + "/region?uniprot=1"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # combine all variants in a string to submit to VEP
    vcf_lines = "{\"variants\" : [" + variants + "]}"
    #print(vcf_lines)

    r = requests.post(server + ext, headers=headers, data=vcf_lines) 
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    #with open('data.json', 'w') as f:
    #    json.dump(decoded, f)
    #print(json.dumps(decoded, indent=2))
    '''
    {
    "input": "1 3625748 rs200029021 T G .",
    "transcript_consequences": [
      {
        "transcript_id": "ENST00000344579",
        "consequence_terms": [
          "intron_variant"
        ],
      },
      {
        "codons": "gTg/gGg",
        "biotype": "protein_coding",
        "impact": "MODERATE",
        "polyphen_prediction": "probably_damaging",
        "cds_end": 329,
        "consequence_terms": [
          "missense_variant"
        ],
        "strand": 1,
        "gene_symbol_source": "HGNC",
        "variant_allele": "G",
        "cds_start": 329,
        "hgnc_id": "HGNC:27007",
        "cdna_start": 387,
        "protein_start": 110,
        "amino_acids": "V/G",
        "sift_score": 0,
        "protein_end": 110,
        "gene_id": "ENSG00000158109",
        "cdna_end": 387,
        "polyphen_score": 0.982,
        "uniparc": [
          "UPI000014067B"
        ],
        "swissprot": [
          "Q5T0D9.129"
        ],
        "transcript_id": "ENST00000378344",
        "gene_symbol": "TPRG1L",
        "sift_prediction": "deleterious"
      },
    '''

    colnames = ["variant", "EnsID", "SPID", "PDBID", "mutaa", "SIpred", "SIscore", "PPpred", "PPscore"]
    results = pd.DataFrame(columns=colnames)
    for i in decoded:
        var = i.get("input")
        #print("input:", var, "\n")
        for key in i:
            if key == "transcript_consequences":
                for j in i[key]:
                    # get all predictions, not just deleterious
                    if "sift_prediction" in j or "polyphen_prediction" in j:
                        gene_id = j.get("gene_id")

                        if "swissprot" in j:
                            sp = j.get("swissprot")[0].split('.')[0]

                            # deal with trembl-only entries later
                            #elif "trembl" in j:
                            #    sp = j.get("trembl")[0].split('.')[0]
                            #else:
                            #    sp = 'NA'

                            # get the amino acid mutation
                            aa = j.get("amino_acids").split("/")[0]  # A/C
                            aa = aa.strip('\s+')
                            aanum = str(j.get("protein_start"))
                            alt_aa = j.get("amino_acids").split("/")[1]
                            alt_aa = alt_aa.strip()
                            mutaa = aa + aanum + alt_aa  # e.g. T355C

                            row = pd.Series({'variant': var,
                                            'EnsID': gene_id,
                                            'SPID': sp,
                                            'PDBID': '',
                                            'mutaa': mutaa,
                                            'SIpred': j.get("sift_prediction"),
                                            'SIscore': j["sift_score"],
                                            'PPpred': j.get("polyphen_prediction"),
                                            'PPscore': j["polyphen_score"]})

                            #print("var:", var, "mutaa:", mutaa)
                            results = pd.concat([results, row.to_frame().T])

                            break # just take the first hit (avoid isoforms)
    results.set_index('variant', inplace=True)

    return(results)

def get_pdb_id(results):
    '''
    Retrieve the list of PDB IDs for each Uniprot ID
    For each amino acid position:
        1. identify if there is an x-ray or EM structure
        2. if so, get the PDB ID for the highest resolution structure
        3. if not, get the PDB ID of any NMR structures available
    '''
    spid_list = list(set(results.SPID)) # unique set of SwissProt ids

    for spid in spid_list:
        print("spid:", spid)
        url = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/" + spid 
        r = requests.get(url)
        #print(r.json())
        if r.status_code == 200:
            pdbid_list = r.json()
        else:
            pdbid_list = "None"

        if (pdbid_list != 'None'):
            aa_list = results.loc[results['SPID']==spid,"mutaa"]
            for aa in aa_list:
                #print("aa:", aa)
                aanum = aa[1:-1] 
                for m in pdbid_list:
                    #print('pdbid:', m)
                    maxres = 10
                    for j in pdbid_list[m]:
                        #print(aanum)
                        #print('pdbid:', j['pdb_id'],'chain:',j['chain_id'],'resolution:',j['resolution'],'start:', j["unp_start"], 'end:', j['unp_end'])
                        # if n in range of pdb, use pdb id instead of uniprot id
                        if ((int(aanum) >= j['unp_start']) & (int(aanum) <= j['unp_end'])):
                            #print(type(j['resolution']),' resolution:',j['resolution'],"|",sep='',)
                            if j["resolution"] is None:
                                res = 0
                            else:
                                res = j["resolution"]
                            #print("res:", res, "maxres:", maxres)
                            if (res == 0) & (maxres < 10):
                                # NMR, but already have xray
                                break
                            elif ((res == 0) & (maxres == 10) | (0 < res < maxres)):
                                # either NMR, no xray, or xray with better resolution
                                # note: replaces existing PDBID if NMR only
                                pdbid = j["pdb_id"].upper() + "_" + j["chain_id"]
                                results.loc[(results['SPID']==spid) & (results['mutaa']==aa),'PDBID'] = pdbid
                                if res > 0:
                                    maxres = res # if xray, reset maxres

def variant_string(mut_list):
    ''' 
    create a string that includes all the variants from one PDB ID or Uniprot ID for iCn3D
    '''
    variant_str = ""
    for mutaa in mut_list:
        s = re.split(r'(\d+)', mutaa) # need the number and new aa
        if variant_str == '':
            variant_str += s[1] + ' ' + s[2]
        else:
            variant_str += ',' + s[1] + ' ' + s[2]
 
    return variant_str

def get_iCn3D_path(gene_res):
    '''
    generates the iCn3D path(s) based on the variants for a given gene
    TODO: need to modify to produce BED file to show SIFT/Polyphen score
    need a separate URL for each PDB/AFID
    '''
    
    date = datetime.now()
    url_path='https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?'

    url_list = dict()

    spid = gene_res.SPID.unique().tolist()[0] # should be only 1 swissprotID per gene
    print("SwissProt ID:", spid)
    sift_str = ''
    poly_str = ''
    scap_str = ''

    afid = 0
    url_list[spid] = "No deleterious mutations found for " + spid
    # check to see if we are using AlphaFold structure
    for row in gene_res.itertuples():
        if (row.PDBID == ''):
            if ((row.SIpred == 'deleterious') | (row.PPpred == 'probably_damaging')):
                afid = 1
    
    # all mutations in the alphafold structure in the same URL
    if afid == 1:
        print("Getting alphafold url...")
        for row in gene_res.itertuples():
            mutaa = row.mutaa
            sift = 0
            # check if any deleterious mutations
            if (row.SIpred == 'deleterious'):
                s = re.split(r'(\d+)', mutaa) # need the number and new aa
                sift_str += ',' +  s[1] + ' ' + s[2]
                scap_str += ',' + spid + '_A' + '_' + s[1] + "_" + s[2] # e.g. P16860_A_113_Y
                sift = 1
            if (row.PPpred == 'probably_damaging'):
                p = re.split(r'(\d+)', mutaa) 
                poly_str += ',' +  p[1] + ' ' + p[2]
                if sift == 0: # avoid repeating scap
                    scap_str += ',' + spid + '_A' + '_' + p[1] + "_" + p[2]

        if scap_str != '': 
            sift_str = sift_str.lstrip(',')
            poly_str = poly_str.lstrip(',')
            scap_str = scap_str.lstrip(',')
            
            url_query =  'afid=' + spid + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command='
            url_command = 'view annotations; set annotation cdd; set view detailed view;  set thickness | stickrad 0.2'    
            url_command += '; add track | chainid ' + spid + '_A' + ' | title SIFT_predict | text ' + sift_str
            url_command += '; add track | chainid ' + spid + '_A' + ' | title PolyPhen_predict | text ' + poly_str
            url_command += '; scap interaction ' + scap_str

            iCn3Durl = url_path + url_query + url_command

            url_command = quote(url_command) # encode the spaces for URL
            iCn3Durl = url_path + url_query + url_command
            url_list[spid] = iCn3Durl

    # go through the PDB IDs
    print('getting PDB urls...')
    for structid in gene_res.PDBID.unique().tolist():
        if structid == '':
            continue
        #print('structid:', structid)
        sift_str = ''
        poly_str = ''
        scap_str = ''
        url_query = ''
        url_command = ''
        iCn3Durl = ''

        # have to check for offset between Uniprot -> PDB residue number mapping
        pdb = structid.split("_")[0]
        chainid = structid.split("_")[1]
        url = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/" + pdb
        r = requests.get(url)
        if r.status_code == 200:
            uniprot_map = r.json()
        else:
            uniprot_map = "None"
        
        for id in uniprot_map:
            if id.upper() == pdb:
                mappings = uniprot_map[id]["UniProt"][spid]["mappings"]
                if mappings[0]['chain_id'] == chainid:
                    pdb_start_num = mappings[0]['start']['author_residue_number']
                    uniprot_start_num = mappings[0]['unp_start']
        
        if pdb_start_num != None:
            diff = uniprot_start_num - int(pdb_start_num)
        else:
            diff = 0

        # put all mutations for a given PDB in same url
        for row in gene_res.itertuples():
            mutaa = row.mutaa
            s = re.split(r'(\d+)', mutaa) # need the number and new aa
            pdb_mutaa_num = int(s[1]) - diff # correct for uniprot->pdb offset 
            pdbid = row.PDBID
            sift = 0
            if pdbid == structid:
                # check if any deleterious mutations
                if (row.SIpred == 'deleterious'):
                    sift_str += ',' +  s[1] + ' ' + s[2]
                    scap_str += ',' + structid + '_' + str(pdb_mutaa_num) + "_" + s[2] # e.g. 1HLZ_A_113_Y
                    sift = 1
                if (row.PPpred == 'probably_damaging'):
                    poly_str += ',' +  s[1] + ' ' + s[2]
                    if sift == 0: # avoid repeating scap
                        scap_str += ',' + structid + '_' + str(pdb_mutaa_num) + "_" + p[2]

        if scap_str != '':
            sift_str = sift_str.lstrip(',')
            poly_str = poly_str.lstrip(',')
            scap_str = scap_str.lstrip(',')
            url_query =  'pdbid=' + structid.split("_")[0] + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command='
            url_command = 'view annotations; set annotation cdd; set view detailed view;  set thickness | stickrad 0.2'    
            url_command += '; add track | chainid ' + structid + ' | title SIFT_predict | text ' + sift_str
            url_command += '; add track | chainid ' + structid + ' | title PolyPhen_predict | text ' + poly_str
            url_command += '; scap interaction ' + scap_str

            iCn3Durl = url_path + url_query + url_command
            url_command = quote(url_command) # encode the spaces for URL
            iCn3Durl = url_path + url_query + url_command
            url_list[structid] = iCn3Durl 
        else:
            url_list[structid] = "No deleterious mutations found for " + structid
    
    return(url_list) 

def print_html(args, url_list, results):
    '''
    Write out the results dataframe to an HTML file with iCn3D links
    '''
    fout = args.v + "_output.html"
    f = open(fout, 'w')
  
    # html code 
    html = """<html>
    <head>
    <title>iCn3D links</title>
    <style>
    table, th, td {
        border: 1px solid black;
        border-collapse: collapse;
    }
    th, td {
        padding: 10px;
    }
    th {
        background-color: #D3D3D3;
    }
    tr {
        border-bottom: 1px solid #ddd;
    }
    </style>
    </head>
    <body>
    <h2>Click on a link to open the deleterious variants in iCn3D</h2>
     <table border='1'>"""
    
    html += "Input file: " + args.v + "<br>"
    if args.t:
        html += "SIFT & PolyPhen scores taken from TCGA VCF file.<br>"

    # output table
    html += "<tr><th>#</th><th>Variant</th><th>Gene</th><th>UniprotID</th><th>PDB ID</th><th>mutaa</th><th>SIFT</th><th>PolyPhen</th><th>iCn3D link</th></tr>"
    n = 0
    for row in results.itertuples():
       # limit output to 1000 rows
        n += 1
        if n > 1000:
            print("HTML output limited to 1000 rows...")
            break
        url = ''
        if row.PDBID == '':
            url1 = url_list[row.EnsID][row.SPID]
        else:
            url1 = url_list[row.EnsID][row.PDBID]

        # URL is repeated for same PDB/Alphafold ID
        if re.match(r"^https", url1):
            url = "<a href=" + ''.join(url1) + " target=\"_blank\">iCn3D link</a><br>"
        else:
            url = url1 # text output, no link

        html += "<tr><td>" + str(n) + "</td><td>" + str(row.Index) + "</td><td>" + row.EnsID + "</td><td>" + row.SPID + "</td><td>" \
            + row.PDBID + "</td><td>" + row.mutaa + "</td><td>" \
            + row.SIpred + " " + str(row.SIscore) + "</td><td>" \
            + row.PPpred + " " + str(row.PPscore) + "</td><td>" + url + "</td></tr>"

    html += "</table>"
    html += "</body>"
    html += "</html>"

    f.write(html)
    f.close()
    
def print_csv(args, url_list, results):
    '''
    print the results dataframe in a .csv file
    '''
    fout = args.v + "_results.csv"
    for index,row in results.iterrows():

        url = ''
        if row.PDBID == '':
            url1 = url_list[row.EnsID][row.SPID]
        else:
            url1 = url_list[row.EnsID][row.PDBID]
        
        # URL is repeated for same PDB/Alphafold ID
        if re.match(r"^https", url1):
            url = '=HYPERLINK("' + url1 + '","iCn3D link")'
        else:
            url = url1 # text output, no link
        results.loc[index, 'Link'] = url
    
    results.to_csv(fout)

def get_vcf(vcff):
    '''
    read in the VCF file & extract variants
    '''
    vcf = []
    vcf_reader = VariantFile(vcff)
    n = 0
    for record in vcf_reader.fetch():
        n += 1
        line = str(record)
        vcf_line = line.split()
        # only get variants
        if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
            vcf_line[0] = vcf_line[0].replace("chr","") # can't have chr in VEP REST API submission
            vcf.append(vcf_line[0:6])

            # can make limit a command-line argument
            if len(vcf) == 1000:
                print("Stopping after 1000 variants...(out of",n,"vcf lines)")
                break

    return(vcf)

def get_vcf_tcga(vcff):
    ''' 
    Get Ensembl ID, SwissProt ID, SIFT, Polyphen predictions & scores directly from TCGA VCF file
    Note: not recommended, there are lots of discrepancies between TCGA values & VEP/Ensembl REST values
    '''

    # read in the VCF file
    vcf = []
    vcf_reader = VariantFile(vcff)

    colnames = ["variant", "EnsID", "SPID", "PDBID", "mutaa", "SIpred", "SIscore", "PPpred", "PPscore"]
    results = pd.DataFrame(columns=colnames)

    for record in vcf_reader.fetch():
        line = str(record)
        vcf_line = line.split()
        variant, mutaa, spid = '', '', ''
        c,d,e,sift_split,polyph_split = (), (), (), (), ()
        
        # only get variants
        if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
            vcf_line[0] = vcf_line[0].replace("chr","") 
            variant = vcf_line[0:6]
            vcf.append(variant)

            c = record.info["CSQ"] # consequence: multiple INFO strings stored as a tuple
            #Gene: position 5
            #Protein: 15 110/272
            #Amino acids: 16 V/G
            #SwissProt: 31
            #SIFT: 36 deleterious(0)
            #PolyPhen: 37 probably_damaging(0.993)

            for d in c: 
                e = d.split("|")
                if e[1] == "missense_variant":
                    #print("gene id:", e[4], "swissprot:", e[31]) # lots of Trembl IDs, need Swissprot to get PDB IDs
                    if e[31] == '': # lot of these are missing in TCGA files
                        spid = get_protein_id(e[4], 1)
                        print("getting spid from REST API:", spid)
                    else:
                        spid = e[31]

                    mutaa = e[15].split("/")[0] + e[14].split("/")[0] + e[15].split("/")[1]

                    sift_split = ['','',0,'']
                    if e[35] != '':
                        sift_split = re.split(r'(\w+)\(([0-9]+\.?\d*)\)', e[35])

                    polyph_split = ['','',0,''] 
                    if e[36] != '':
                        polyph_split = re.split(r'(\w+)\(([0-9]+\.?\d*)\)', e[36])

                    row = pd.Series({'variant': variant,
                                    'EnsID': e[4],
                                    'SPID': spid,
                                    'PDBID': '',
                                    'mutaa': mutaa,
                                    'SIpred': sift_split[1],
                                    'SIscore': float(sift_split[2]),
                                    'PPpred': polyph_split[1],
                                    'PPscore': float(polyph_split[2])})
                    results = pd.concat([results, row.to_frame().T])

            if len(vcf) == 1000:
                print("Stopping after 1000 variants...")
                break

    results.set_index('variant', inplace=True)

    return(vcf, results)

def cli():

    # format the description
    class RawFormatter(HelpFormatter):
        def _fill_text(self, text, width, indent):
            return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

    desc = f'''
            SNP2iCn3D.py: Runs VEP on single-nucleotide variants extracted from a VCF file 
            and generates iCn3D links for predicted deleterious mutations.  
            If the VCF file is from TCGA, the -t flag extracts the VEP data from the VCF file instead of 
            submitting to the VEP server.

            Two modes: 
                1) provide a comma-separated list of specific Ensembl genes with the -g flag, 
                   only locations matching those genes will be extracted; 
                2) do not provide the -g flag, all variants 
                   in the VCF will be run through VEP. 
                    
            Output:    
                - an html file listing the variants, SIFT & PolyPhen scores, and iCn3D links
                - a .csv file that can be imported into Numbers or Google Sheets (URLs are broken in Excel)
            '''
    parser = ArgumentParser(description=desc, formatter_class=RawFormatter)
    parser.add_argument('-g', metavar='GENE', help="Select only variants from Ensembl Gene IDs of interest")
    parser.add_argument('-v', required=True, metavar='VCF', help="VCF file to extract the variants; must be compressed with bgzip \
        and the tabix .tbi index file must be present.")
    parser.add_argument('-t', action='store_true', help="Extract SIFT & PolyPhen scores from VCF file from TCGA instead of submitting to VEP")
                        #action='store_true' means default is false 
    parser.add_argument('-s', type=str, default='human', help="species (default human) (use a common name from http://rest.ensembl.org/info/species.json)")
    return parser

def main(args):
    '''
    Output:
        Input VCF file: <vcf_file> # could name <vcf_file>_out.csv?
        Coord EnsGene SPid PDBID mutaa SIFT Polyphen iCn3Dlink
          -- mutaa, iCn3dlink grouped by pdbid, spid
    
     -read in the whole VCF file, use client.get_gene_ids to fill gene_ids dictionary: 
        gene_ids[loc]["ens_id"] and gene_ids[loc]["sp_id"]
        -if args.g, select subset of genes
     -go through genes, get VEP results 
        -if args.t, get EnsID, SwissProtID, SIFT, Polyphen results from TCGA VCF file (SwissProt ID may not be current -> replace with get_protein_id)
     -generate iCn3D link for each gene
    
     Data structures
      vcf: list of single nucleotide variants pulled from VCF file "19 15256965 . T G . . ."
      gene_ids: gene_ids[loc]["ens_id"] and gene_ids[loc]["sp_id"]
        -> requires vcf
        -> limited to genes in args.g if present
      gene_id_list: list of Ensembl gene ids
        -> extracted from gene_ids
      variants (subset of vcf): string of variants for submitting to vep "19 15256965 . T G . . ."
        -> requires vcf, gene_ids_select
      sift: dictionary sift["gene_id"][mutaa] = {"sift_prediction": j.get("sift_prediction"), "sift_score": j["sift_score"]}
        -> requires variants
      polyphen: dictonary polyphen["gene_id"][mutaa] = {"polyphen_prediction": j.get("polyphen_prediction"), "polyphen_score": j["polyphen_score"]}
        -> requires variants
      gene_to_pid: dictionary gene_to_pid['ens_id'] = spid
        -> requires gene_ids
     
     Pandas data frame including all of above:
      variant EnsID SPID mutaa PDBID SIpred SIscore PPpred PPscore 

      url_list: dictionary url_list[gene] 
        -> requires sift, polyphen, gene, gene_to_pid[gene]
        -> 1..n urls per gene
        iCn3D link:
            - each gene has a list of URLs (1...n)
            - need afid/pdbid for each aa position
                struct_id[n] = pid or pdbid
            - need SIFT, PolyPhen strings for each afid/pdbid 
                variant_string
            - need scap strings for each afid/pdbid (different from SIFT/PolyPhen strings)

     HTML:
     EnsGene SPid PDBID SIFT PolyPhen iCn3Dlink
      - all mutations for a given PDB/SPid together in same iCn3D link
        variant_string2
     .csv:
     Coord EnsGene SPid PDBID mutaa SIFT Polyphen iCn3Dlink
    '''

    colnames = ["variant", "EnsID", "SPID", "PDBID", "mutaa", "SIpred", "SIscore", "PPpred", "PPscore"]
    results = pd.DataFrame(columns=colnames)

    # Extract SIFT & PolyPhen scores from TCGA file
    if args.t:
        print("Getting VEP values from TCGA file...")
        vcf, results = get_vcf_tcga(args.v) # returns vcf, gene_ids dict, fills sift & polyphen
        get_pdb_id(results) 

    # Get variants from VCF file, SIFT & PolyPhen from VEP REST API (default)
    else:    
        # extract the variants from the vcf file
        vcf = get_vcf(args.v)

        # need to break into chunks of 200 variants - maximum POST size is 200
        def divide_variants_list(list, n):
            for i in range(0, len(list), n):
                yield list[i:i + n]

        variant_list = list(divide_variants_list(vcf, 200)) # returns a list of lists

        # combine all variants in a string to submit to VEP
        print("\nSubmitting variants to VEP server (batch = 200 variants)...")
        n = 0
        for v in variant_list:
            n += 1
            print("processing batch", n)
            variants = ''
            for c in v:
                loc = " ".join(c) 
                variants += "\"" + loc + "\" ," 
            variants = variants[:-1]
            vep_result = vep_output(variants, args)
            results = pd.concat([results, vep_result])
            time.sleep(2)
        print("Done")
            
        # find if mutations are in PDB structures or not
        print("\nGetting PDB IDs...", end='')
        get_pdb_id(results)
        print("Done") 

    # end if args.t

    print("\nGenerating iCn3D URLs...")
    url_list = defaultdict(dict)
    gene_id_list = list(set(results.EnsID)) # unique set of gene ids

    for gene in gene_id_list:
        print("=========================\nGetting link for ", gene)

        gene_res = results[results['EnsID'] == gene]
        url_list[gene] = get_iCn3D_path(gene_res)

    # generate an html page with results dataframe, iCn3D links
    print("\nPrinting html file...")
    print_html(args, url_list, results)

    # generate a .csv file with results dataframe, iCn3D links
    print("Printing .csv file...")
    print_csv(args, url_list, results)

if __name__ == '__main__':   
    main(cli().parse_args())   
