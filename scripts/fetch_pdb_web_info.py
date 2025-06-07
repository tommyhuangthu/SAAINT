#!/usr/bin/env python3
from utils import read_list, read_fasta
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from utils import is_fasta_legal
import os, sys, re, html

abs_path = os.path.dirname(os.path.realpath(__file__))

url_tag = 'https://www.rcsb.org/structure'

def read_url_to_string(url):
    return urlopen(url).read().decode('utf-8')

def get_organisms(content):
    organisms = []
    p1 = re.compile(r'Organism\(s\):&nbsp<\/strong><a (.*?)<\/li>')
    matches = re.findall(p1, content)
    for match in matches:
        p2 = re.compile(r'>(.*?)<\/a>')
        hits = re.findall(p2, match)
        for hit in hits:
            organisms.append(hit)
    return '; '.join(organisms)

def get_mutation_status(content):
    p1 = re.compile(r'Mutation\(s\)\:&nbsp<\/strong>([YesNo]+)&nbsp')
    matches = re.findall(p1, content)
    if matches:
        return matches[0].lower()
    else:
        return 'N.A.'

def get_classification(content):
    #p1 = re.compile(r'Classification:&nbsp<a href="\/search\?q=struct_keywords\.pdbx_keywords:([0-9A-Za-z\/\(\)\- ]+)" onClick')
    p1 = re.compile(r'Classification:&nbsp<a href="\/search\?q=struct_keywords\.pdbx_keywords:(.*?)" onClick')
    matches = re.findall(p1, content)
    if matches:
        return matches[0].lower()
    else:
        return 'N.A.'

def get_deposited_date(content):
    p1 = re.compile(r'Deposited:&nbsp<\/strong>([0-9-]+)&nbsp')
    matches = re.findall(p1, content)
    if matches:
        return matches[0]
    else:
        return 'N.A.'

def get_released_date(content):
    p1 = re.compile(r'Released:&nbsp<\/strong>([0-9-]+)&nbsp')
    matches = re.findall(p1, content)
    if matches:
        return matches[0]
    else:
        return 'N.A.'

def get_exp_method(content):
    p1 = re.compile(r'Method:&nbsp<\/strong>([A-Z- ]+)<\/li>')
    matches = re.findall(p1, content)
    if matches:
        return matches[0].lower()
    else:
        return 'N.A.'


def get_resolution(content):
    p1 = re.compile(r'Resolution:&nbsp<\/strong>([0-9.]+)')
    matches = re.findall(p1, content)
    if matches:
        return matches[0]
    else:
        return 'N.A.'

def get_xray_r_free(content):
    p1 = re.compile(r'R-Value Free:&nbsp<\/strong><div> ([0-9.]+) \(')
    matches = re.findall(p1, content)
    if matches:
        return matches[0]
    else:
        return 'N.A.'

def get_xray_r_work(content):
    p1 = re.compile(r'R-Value Work:&nbsp<\/strong><div> ([0-9.]+) \(')
    matches = re.findall(p1, content)
    if matches:
        return matches[0]
    else:
        return 'N.A.'

def get_em_aggregation_state(content):
    p1 = re.compile(r'Aggregation State:&nbsp<\/strong>([A-Z]+)&nbsp')
    matches = re.findall(p1, content)
    if matches:
        return matches[0].lower()
    else:
        return 'N.A.'

def get_em_reconstruction_method(content):
    p1 = re.compile(r'Reconstruction Method:&nbsp<\/strong>([A-Z ]+)&nbsp')
    matches = re.findall(p1, content)
    if matches:
        return matches[0].lower()
    else:
        return 'N.A.'

def get_pubmed_id(content):
    p1 = re.compile(r'pubmed_id:([0-9]+)')
    matches = re.findall(p1, content)
    if matches:
        return matches[0]
    else:
        return 'N.A.'

def get_doi(content):
    p1 = re.compile(r'<strong>DOI:&nbsp<\/strong><a href="https:\/\/doi.org\/(.*?)" ')
    matches = re.findall(p1, content)
    if matches:
        return matches[0]
    else:
        return 'N.A.'

def get_rcsb_title(content):
    p1 = re.compile(r'<title>RCSB PDB - \w{4}: (.*?)<\/title>')
    matches = re.findall(p1, content)
    if matches:
        return matches[0].lower()
    else:
        return 'N.A.'


def fetch_pdb_web_info(fasta_path, pdbid):
    # truncate pdb name: 1e28-assembly => 1e28
    pdbid = pdbid[0:4]
    url = f'{url_tag}/{pdbid}'
    print(f'fetching pdb information for {pdbid} from {url}')
    
    # quick exit if synchronized fasta file does not exist
    fasta = f'{fasta_path}/{pdbid[1:3]}/{pdbid}.fasta'
    if not os.path.exists(fasta):
        print(f'done fetching ({fasta} does not exist, please skip)')
        exit(0)

    if not is_fasta_legal(fasta):
        print(f'done fetching ({fasta} is not a legal fasta file, please skip)')
        exit(0)
    

    try:
        response = urlopen(url)
        content = response.read().decode('utf-8')
    except HTTPError as e:
        print(f"HTTP Error: {e.code} - {e.reason}")
        exit(0)
    except URLError as e:
        print(f"URL Error: {e.reason}")
        exit(0)
    
    # parse webpage content
    organism_str = get_organisms(content)
    mut_status = get_mutation_status(content)
    class_ = html.unescape(get_classification(content))
    deposit_date = get_deposited_date(content)
    release_date = get_released_date(content)
    method = get_exp_method(content)
    pmid = get_pubmed_id(content)
    doi = html.unescape(get_doi(content))
    title = html.unescape(get_rcsb_title(content))
    resolution = get_resolution(content)
    r_free = get_xray_r_free(content)
    r_work = get_xray_r_work(content)
    aggregation = get_em_aggregation_state(content)
    reconstruction = get_em_reconstruction_method(content)

    return mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title


def write_pdb_info(tsv, mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title):
    with open(tsv, 'w') as f:
        f.write('\t'.join(['PDB_ID', 'Mutation(s)', 'Classification', 'Deposit_date', 'Release_date', 'Method', 'Resolution', 'R_free', 'R_work', 'PMID', 'DOI', 'Title\n']))
        f.write('\t'.join(['"'+pdbid+'"', mut_status, '"'+class_+'"', deposit_date, release_date, '"'+method+'"', resolution, r_free, r_work, pmid, doi, '"'+title+'"'+'\n']))
    return
    

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'Usage: {sys.argv[0]} pdbid')
        print(f'Example: {sys.argv[0]} 5zxv')

    pdbid = sys.argv[1].lower().replace('_', '')
    fasta_path = f'{abs_path}/../database/fasta_divided'
    mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title = fetch_pdb_info_from_web(fasta_path, pdbid)
    write_pdb_info(f'{pdbid}.tsv', mut_status, class_, deposit_date, release_date, method, resolution, r_free, r_work, pmid, doi, title)

