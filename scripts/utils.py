#!/usr/bin/env python3
import os, re, glob
from Bio import PDB, SeqIO, BiopythonWarning
import numpy as np

# global variables
homedir      = '/home/xiaoqiah/turbo/work'
unidesign    = f'{homedir}/UniDesign/UniDesign'
faspr        = f'{homedir}/FASPR/FASPR'
pulchra      = f'{homedir}/Programs/pulchra/pulchra'
abs_path     = os.path.dirname(os.path.realpath(__file__))
abrsa        = f'{abs_path}/../ext_bin/AbRSA/AbRSA'
abrsa_pdb    = f'{abs_path}/../ext_bin/AbRSA_PDB'
abalign      = f'{abs_path}/../ext_bin/Abalign_V2/Abalign_2024.10'
abalign_lib  = f'{abs_path}/../ext_bin/Abalign_V2/lib'
tmalign      = f'{abs_path}/../ext_bin/TMalign_cpp'
struct_ref_vl = f'{abs_path}/../ext_bin/4unu_A.pdb'
struct_ref_vh = f'{abs_path}/../ext_bin/6cnw_A.pdb'


aadict_3to1 = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
        'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'}


def residue_has_ca_atom(residue):
    ca_found = False
    for atom in residue.get_atoms():
        if atom.get_name() == 'CA':
            return True
    return False


class SelectChains(PDB.Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        if chain.get_id() in self.chain_ids: 
            return 1
        return 0

    def accept_atom(self, atom):
        if not atom.get_name().startswith('H'):
            return 1
        return 0


class Entity:
    def __init__(self, ids, name, seq, species):
        self.chain_ids = ids
        self.name = name
        self.fas_seq = seq
        self.fas_ab_res = []
        self.species = species
        self.chain_types = [] # protein, DNA, RNA, or unknown
        self.real_pdb_seqs = []
        self.filled_pdb_seqs = []
        self.abrsa_types = [] # 'heavy', 'light', 'heavy_light', 'nonab'
        self.cdr_res = []
        self.vgene_subgroups = []
        self.paired_chains = []
        self.mean_radii = []

    def get_chain_ids(self):
        return self.chain_ids

    def set_chain_ids(self, new_ids):
        if len(new_ids) == len(self.chain_ids):
            self.chain_ids = new_ids
        return

    def get_chain_id(self, index):
        return self.chain_ids[index]

    def get_name(self):
        return self.name

    def get_fas_seq(self):
        return self.fas_seq

    def get_fas_seq_len(self):
        return len(self.fas_seq)

    def get_fas_ab_res(self):
        return self.fas_ab_res

    def set_fas_ab_res(self, ab_res):
        self.fas_ab_res = ab_res
        return

    def init_chain_types(self):
        self.chain_types = ['' for chain_id in self.chain_ids]
        return

    def get_chain_type(self, index):
        return self.chain_types[index]

    def set_chain_type(self, index, ctype):
        self.chain_types[index] = ctype
        return

    def init_real_pdb_seqs(self):
        self.real_pdb_seqs = ['' for chain_id in self.chain_ids]
        return

    def get_real_pdb_seq(self, index):
        return self.real_pdb_seqs[index]

    def set_real_pdb_seq(self, index, seq):
        self.real_pdb_seqs[index] = seq
        return
    
    def get_real_pdb_len(self, index):
        return len(self.real_pdb_seqs[index])


    def init_filled_pdb_seqs(self):
        self.filled_pdb_seqs = ['' for chain_id in self.chain_ids]
        return

    def get_filled_pdb_seq(self, index):
        return self.filled_pdb_seqs[index]

    def set_filled_pdb_seq(self, index, seq):
        self.filled_pdb_seqs[index] = seq
        return
    
    def get_filled_pdb_len(self, index):
        return len(self.filled_pdb_seqs[index])

    def get_species(self):
        return self.species

    def init_abrsa_types(self):
        self.abrsa_types = ['' for chain_id in self.chain_ids]
        return
    
    def set_abrsa_type(self, index, abrsa_type):
        self.abrsa_types[index] = abrsa_type
        return

    def get_abrsa_type(self, index):
        return self.abrsa_types[index]

    def get_abrsa_types(self):
        return self.abrsa_types

    def merge_entities(self, new_entity):
        for chain_id in new_entity.get_chain_ids():
            self.chain_ids.append(chain_id)
        return

    def get_cdr_res(self):
        return self.cdr_res if self.cdr_res else None

    def set_cdr_res(self, hcdr_res, lcdr_res):
        if not self.cdr_res:
            self.cdr_res.append(hcdr_res)
            self.cdr_res.append(lcdr_res)
        else:
            h_list = self.cdr_res[0]
            h_list = list(set(h_list + hcdr_res))
            h_list.sort()
            self.cdr_res[0] = h_list
            l_list = self.cdr_res[1]
            l_list = list(set(l_list + lcdr_res))
            l_list.sort()
            self.cdr_res[1] = l_list
        return

    def init_vgene_subgroups(self):
        self.vgene_subgroups = [[] for chain_id in self.chain_ids]
        return

    def set_vgene_subgroup(self, index, subgroup):
        self.vgene_subgroups[index] = subgroup
        return

    def get_vgene_subgroup(self, index):
        return self.vgene_subgroups[index]

    def init_paired_chains(self):
        self.paired_chains = [[] for chain_id in self.chain_ids]
        return

    def get_paired_chain(self, index):
        return self.paired_chains[index]

    def init_mean_radii(self):
        self.mean_radii = [0.0 for chain_id in self.chain_ids]
        return

    def set_mean_radius(self, index, radius):
        self.mean_radii[index] = radius
        return

    def get_mean_radius(self, index):
        return self.mean_radii[index]


def read_list(list_file, replace=False):
    lst = []
    with open(list_file) as f:
        for line in f.readlines():
            if line.startswith('#'): continue
            line = line.rstrip()
            if replace == True:
                line = line.replace('.00','')
                line = line.replace('+','')
                line = line.replace('-','')
                line = line.replace('_','')
            lst.append(line)
    return lst


def read_fasta(fasta):
    header_list, seq_dict = [], dict()
    with open(fasta, 'r') as f: 
        blocks = f.read().lstrip('>').split('\n>')
    for block in blocks:
        lines = block.splitlines()
        header = lines[0]
        sequence = ''.join(lines[1:])
        header_list.append(header)
        seq_dict[header] = sequence
    return header_list, seq_dict


def write_fasta(header_list, seq_dict, file):
    with open(file, 'w') as f:
        for i in range(0, len(header_list)): 
            f.write('>' + header_list[i] + '\n' + seq_dict[header_list[i]])
    return


def run_abrsa(abrsa, fasta, abrsa_out):
    # -k, kabat; -c, chothia; -g, imgt
    os.system(f'{abrsa} -c -i {fasta} > {abrsa_out}')
    return


def run_abrsa_pdb(abrsa_pdb, pdbfile, abrsa_out):
    os.system(f'{abrsa_pdb} -S {pdbfile} > {abrsa_out}')
    return


def run_unidesign_find_interface_residues(unidesign, pdb, chains_1, chains_2, unidesign_out):
    os.system(f'{unidesign} --command FindInterfaceRes --pdb {pdb} --split_chains {chains_1},{chains_2} > {unidesign_out}')
    return


def run_faspr_repack_pdb(faspr, pdb):
    pdbid = pdb.replace('.pdb', '')
    os.system(f'{faspr} -i {pdb} -o {pdbid}.faspr.pdb > /dev/null')
    os.system(f'mv {pdbid}.faspr.pdb {pdb}')
    return


def run_pulchra(pulchra, pdb):
    pdbid = pdb.replace('.pdb', '')
    # remove alt_id in the PDB file
    fout = open(f'{pdbid}.pulchra.pdb', 'w')
    with open(pdb, 'r') as f:
        for line in f.readlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                newline = line[:16]+' '+line[17:]
                fout.write(f'{newline}')
            elif line.startswith('TER'):
                fout.write(f'{line}')
    fout.close()
    os.system(f'mv {pdbid}.pulchra.pdb {pdb}')
    os.system(f'{pulchra} {pdb} > /dev/null')
    os.system(f'mv {pdbid}.rebuilt.pdb {pdb}')
    return


def parse_pdb_fasta_header(header):
    strs = header.split('|')
    pdbid_ndx, chain_str, name_str, species_str = strs[0], strs[1], strs[2], strs[3]
    species = re.findall(r'(.*?) \(\d+\)', species_str)
    if species:
        species = species[0].lower()
    elif species_str:
        species = species_str.lower()
    else:
        species = 'N.A.'
    return pdbid_ndx, chain_str, name_str.lower(), species


def reformat_pulchra_rebuilt_pdb(pdb, chain_id):
    pdbid = pdb.replace('.pdb', '')
    fout = open(pdbid+'.reformat.pdb', 'w')
    with open(pdb, 'r') as f:
        for line in f.readlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                newline = line[:21] + chain_id + line[22:54]+'  1.00  0.00           ' + line[13]
                fout.write(f'{newline}\n')
            elif line.startswith('TER'):
                fout.write(f'{line}')
    fout.close()
    os.system(f'mv {pdbid}.reformat.pdb {pdb}')
    return


def run_pdbfixer(pdb, chain_id):
    pdbid = pdb.replace('.pdb', '')
    os.system(f'pdbfixer {pdb} --output={pdbid}.fixer.pdb --add-atoms=heavy')
    os.system(f'mv {pdbid}.fixer.pdb {pdb}')
    return


def reformat_pdbfixer_pdb(pdb):
    pdbid = pdb.replace('.pdb', '')
    fout = open(pdbid+'.fixer.pdb', 'w')
    with open(pdb, 'r') as f:
        for line in f.readlines():
            if line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('TER'):
                fout.write(f'{line}')
    fout.close()
    os.system(f'mv {pdbid}.fixer.pdb {pdb}')
    return

def determine_chain_type_by_pdb_content(model, chain_id):
    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                if residue.get_resname() in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 
                        'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
                        'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                    return 'protein'
                elif residue.get_resname() in ['A', 'C', 'G', 'U']:
                    return 'RNA'
                elif residue.get_resname() in ['DA', 'DC', 'DG', 'DT']:
                    return 'DNA'
    return 'unknown'


def reindex_pdb(pdb, start_index=1):
    f = open(pdb, 'r')
    pdb_txt = f.read()
    f.close()
    
    pdb_txt_reindex = ''
    res_old_index = ''
    res_atoms = []
    for line in pdb_txt.splitlines():
        if line.startswith('TER'):
            pdb_txt_reindex += 'TER\n'
            continue
        elif not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        else:
            if not line[16] in ['A', '1', ' ']: # alternative location identifier
                continue
        res_index = line[22:27] # residue index on the chain (insertion code included)
        res_chain_id = line[21] # chain identifier
        res_atom_name = line[12:16] # residue atom name 

        if not res_old_index: # first residue encountered
            res_old_index = res_index
            res_new_index = int(start_index)
            chain_id = res_chain_id
            res_new_index_str = str(res_new_index)
            res_new_index_str = ' '*(4-len(res_new_index_str)) + res_new_index_str + ' '
        elif res_chain_id != chain_id: # new chain encountered
            continue
        elif res_index != res_old_index: # new residue encountered
            res_atoms = [] # reset to empty
            res_new_index += 1
            res_old_index = res_index
            res_new_index_str = str(res_new_index)
            res_new_index_str = ' '*(4-len(res_new_index_str)) + res_new_index_str + ' '
        # atoms on the old residue
        elif res_atom_name in res_atoms:
            continue
        res_atoms.append(res_atom_name) # save the atom name
        pdb_txt_reindex += line[:16] + ' ' + line[17:22] + res_new_index_str + line[27:] + '\n'

    f = open(pdb, 'w')
    f.write(pdb_txt_reindex)
    f.close()
    return


def reindex_pdb_by_list(pdb, respos_list):
    f = open(pdb, 'r')
    pdb_txt = f.read()
    f.close()
    
    pdb_txt_reindex = ''
    res_old_index = ''
    res_atoms = []
    res_counter = 0
    for line in pdb_txt.splitlines():
        if line.startswith('TER'):
            pdb_txt_reindex += 'TER\n'
            continue
        elif not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        else:
            if not line[16] in ['A', '1', ' ']: # alternative location identifier
                continue
        res_index = line[22:27] # residue index on the chain (insertion code included)
        res_chain_id = line[21] # chain identifier
        res_atom_name = line[12:16] # residue atom name 
        if not res_old_index: # first residue encountered
            res_old_index = res_index
            #res_new_index = int(start_index)
            res_new_index = respos_list[res_counter]
            chain_id = res_chain_id
            res_new_index_str = str(res_new_index)
            res_new_index_str = ' '*(4-len(res_new_index_str)) + res_new_index_str + ' '
        elif res_chain_id != chain_id: # new chain encountered
            continue
        elif res_index != res_old_index: # new residue encountered
            res_atoms = [] # reset to empty
            #res_new_index += 1
            res_counter += 1
            res_new_index = respos_list[res_counter]
            res_old_index = res_index
            res_new_index_str = str(res_new_index)
            res_new_index_str = ' '*(4-len(res_new_index_str)) + res_new_index_str + ' '
        # atoms on the old residue
        elif res_atom_name in res_atoms:
            continue
        res_atoms.append(res_atom_name) # save the atom name
        pdb_txt_reindex += line[:16] + ' ' + line[17:22] + res_new_index_str + line[27:] + '\n'

    f = open(pdb, 'w')
    f.write(pdb_txt_reindex)
    f.close()
    return




def extract_pdb_chains_to_file(in_pdb, chain_ids, out_pdb):
    fo = open(out_pdb, 'w')
    with open(in_pdb, 'r') as f:
        flag = False
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('TER')):
                fo.write(f'{line}')
            elif line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id in chain_ids:
                    flag = True
                    fo.write(f'{line}')
            elif line.startswith('TER'):
                if flag:
                    flag = False
                    fo.write(f'TER\n')
    fo.close()
    return


def get_l2codes(file_path):
    lst = []
    for l2full in glob.glob(f'{file_path}/[0-9a-z][0-9a-z]'):
        l2code = l2full[len(l2full)-2:]
        lst.append(l2code)
    return lst


def get_pdbids_by_l2code(file_path, l2code):
    lst = []
    if not os.path.exists(os.path.join(file_path, l2code)):
        return lst
    for ent in glob.glob(f'{file_path}/{l2code}/[0-9a-z]{l2code}[0-9a-z]*.cif*'):
        slash_index = -1
        for i in range(len(ent)-1, -1, -1):
            if ent[i] == '/':
                slash_index = i
                break
        entcode = ent[slash_index+1:slash_index+5]
        if entcode not in lst:
            lst.append(entcode)
    return lst


def is_fasta_legal(fasta):
    with open(fasta, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                return True
            else:
                return False


def pdb_to_fasta(pdb, fasta):
    fo = open(fasta, 'w')
    fo.write('>pdb_chain_with_CA_atoms\n')
    with open(pdb, 'r') as fi:
        for line in fi.readlines():
            if line.startswith('ATOM'):
                atomname = line[13:15]
                if atomname == 'CA':
                    resname3 = line[17:20]
                    resname1= aadict_3to1[resname3]
                    fo.write(resname1)
    fo.write('\n')
    fo.close()
    return


def calculate_mean_radius(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    coords = []
    for atom in structure.get_atoms():
        coords.append(atom.coord)

    coords = np.array(coords)
    center = np.mean(coords, axis=0)
    square = np.mean(np.sum((coords - center)**2, axis=1))
    radius = np.sqrt(square)

    return radius


def parse_rsync_mmcif_out(rsync_out):
    obsolete_list, update_dict = [], dict()
    with open(rsync_out, 'r') as f:
        for line in f.readlines():
            line = line.rstrip()
            if line.startswith('deleting'):
                slash_index = -1
                for i in range(len(line)-1, -1, -1):
                    if line[i] == '/':
                        slash_index = i
                        break
                pdbid = line[slash_index+1:-7]
                if pdbid not in obsolete_list:
                    obsolete_list.append(pdbid)
            elif line[-7:] == '.cif.gz':
                slash_index = -1
                for i in range(len(line)-1, -1, -1):
                    if line[i] == '/':
                        slash_index = i
                        break
                pdbid = line[slash_index+1:-7]
                l2code = pdbid[1:3]
                if l2code not in update_dict:
                    update_dict[l2code] = [pdbid]
                else:
                    if pdbid not in update_dict[l2code]:
                        update_dict[l2code].append(pdbid)
    return obsolete_list, update_dict



