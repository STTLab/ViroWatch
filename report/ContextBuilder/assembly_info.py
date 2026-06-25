import csv
import hashlib
from Bio import SeqIO

def read_fasta(fasta):
    return [ (seq.id, len(str(seq.seq)), str(seq.seq)) for seq in SeqIO.parse(fasta, 'fasta')]

def hash_sequence(sequence):
    return hashlib.md5(sequence.encode()).hexdigest()

def read_flye_assembly_info(flye_info):
    _data = {}
    with open(flye_info, encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            _data.update({
                row['#seq_name']: {                                     # Contig ID produced by Flye
                    'length': row['length'],                            # Contig length in base pairs
                    'coverage': row['cov.'],                            # Average read coverage for that contig
                    'is_circular': row['circ.'] == 'Y',                 # Whether Flye thinks the contig is circular (Y/N)
                    'is_repeated_region': row['repeat'] == 'Y',         # Whether the contig represents a repeat region
                    'intracontig_copy_number': row['mult.'],            # Estimated copy number in the genome
                    'alternative_contigs_in_group': row['alt_group'],   # Alternative contigs from the same graph branch
                    'graph_path': row['graph_path']                     # Path through the assembly graph nodes
                }
            })
    return _data

def get_length_with_flye_info(fasta, flye_info, flatten=True):
    fasta_data = read_fasta(fasta)
    flye_data = { # Filter only contigs that are left in medeka polish.
        k: v for k, v in read_flye_assembly_info(flye_info).items() 
        if k in [i[0] for i in fasta_data]
    }
    for contig_name, length, sequence in fasta_data:
        flye_data[contig_name].update({
            'sequence': sequence,
            'length': length,
            'sequence_hash': hash_sequence(sequence),
            'hash_algorithm': 'md5'
        })
    if flatten:
        return [ dict({'contig_id': contig_id}, **record ) for contig_id, record in flye_data.items() ]
    return flye_data
