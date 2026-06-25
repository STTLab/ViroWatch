import os
from datetime import datetime
from typing import Literal
import json
import hashlib
from functools import cache
from collections import defaultdict
import copy

class SierraPyResult:
    def __init__(self, hit_records, json_file_path):
        self._json_file_path = json_file_path
        # All records in the output
        self._hit_records = hit_records
        # All query names in the output
        self._query_names = self.get_query_names()
        
        # Only quries that have been validated i.e. exclude CRITICAL which has no PR, RT, and IN genes found, refuse to process
        self._validated_query_names = self._get_validated_query_names(exclude_levels=['CRITICAL'])
        self._validated_hit_records = self._get_validated_records(self._validated_query_names)

    @classmethod
    def read_sierrapy_json(cls, f_path: os.PathLike):
        with open(f_path, 'r', encoding='utf-8') as f:
            data = json.loads(f.read())
        return cls(data, f_path)

    def get_query_names(self):
        return tuple([rec.get('inputSequence').get('header') for rec in self._hit_records])
    
    def _get_validated_query_names(self, exclude_levels:list[Literal['OK','NOTE','WARNING','SEVERE_WARNING', 'CRITICAL']]):
        q_names = list(self._query_names)
        to_remove = []
        for vlevel in exclude_levels:
            for q in self.get_validation_results(level=vlevel).keys():
                if q not in to_remove: to_remove.append(q)
        return tuple([ q for q in q_names if q not in to_remove ])

    def _get_validated_records(self, validated_query_names):
        headers_to_keep = set(validated_query_names)
        records = self._hit_records
        return [
            r for r in records
            if r.get("inputSequence", {}).get("header") in headers_to_keep
        ]

    def remove_hits_with_warnings(self):
        validated_query_names = self._get_validated_query_names(['CRITICAL', 'SEVERE_WARNING', 'WARNING'])
        data = self._get_validated_records(validated_query_names)
        return SierraPyResult(data)

    def get_aligned_genes(self):
        queries = self._query_names
        aligned_genes = zip(queries, [ rec.get('alignedGeneSequences') for rec in self._hit_records ])
        extracted_data = dict()
        for query, gene_data in aligned_genes:
            extracted_data[query] = tuple([ {
                'name': subrec.get('gene').get('name'),
                'length': subrec.get('gene').get('length'),
                'firstAA':subrec.get('firstAA'),
                'lastAA': subrec.get('lastAA')
            } for subrec in gene_data ])
        return extracted_data
    
    def get_validation_results(
            self,
            level: Literal['OK','NOTE','WARNING','SEVERE_WARNING', 'CRITICAL']|None = None
        ):
        queries = self._query_names
        aligned_genes = zip(queries, [ rec.get('validationResults') for rec in self._hit_records ])
        extracted_data = dict()
        for query, data in aligned_genes:
            if not data:
                continue
            if level is None:
                filtered = data
            else:
                filtered = [d for d in data if d.get('level') == level]
            if filtered:
                extracted_data[query] = data
        return extracted_data
    
    def get_mutations(
            self,
            query_names:list[str]|str|None = None,
            gene_filter:list[Literal['PR','RT','IN']]|Literal['PR','RT','IN', None] = None,
            simple_return: bool = True,
            drm_only: bool = False
        ):
        queries = self._query_names
        if query_names:
            if isinstance(query_names, str):
                query_names = [query_names]
            if not set(query_names).issubset(set(queries)):
                raise ValueError("query_names must be a subset of queries")
            
        aligned_genes = zip(queries, [ rec.get('alignedGeneSequences') for rec in self._hit_records ])
        extracted_data = dict()
        for query, gene_data in aligned_genes:
            if query_names and not query in query_names:
                continue
            extracted_data.setdefault(query, {})
            for record in gene_data:
                gene_name = record.get("gene", {}).get("name")
                if gene_filter:
                    if isinstance(gene_filter, str):
                        gene_filter = [gene_filter]
                    if not gene_name in gene_filter:
                        continue
                mutations = record.get("mutations", [])
                extracted_data[query].update({
                    gene_name: (
                        tuple(m for m in mutations if m.get("isSDRM"))
                        if drm_only else mutations
                    )
                })
        if query_names:
            return dict(sorted(extracted_data.items(), key=lambda x: query_names.index(x[0])))
        if simple_return:
            return {k: v for k, v in extracted_data.items() if v}
        return extracted_data

    @classmethod
    def transform_mutations(cls, data, sep=None):
        """
        Reorganize mutation annotations from a contig-centric structure to a
        gene-centric structure while tracking which contigs contain each mutation.

        Input structure:
            {
                contig_id: {
                    gene: [
                        mutation_dict,
                        ...
                    ]
                }
            }

        Where each mutation_dict contains fields such as:
            consensus, position, AAs, isInsertion, isDeletion, isApobecMutation,
            isApobecDRM, isUnusual, isSDRM, hasStop, primaryType, text

        Transformation:
            The function groups mutations by gene and mutation identifier
            (using the 'text' field, e.g. "L10I"). Mutations observed on multiple
            contigs are merged into a single mutation entry and annotated with a
            list of contig IDs where the mutation occurs.

        Output structure:
            {
                gene: [
                    {
                        ...mutation fields...,
                        "on_contigs": [contig_id, ...]
                    },
                    ...
                ]
            }

        Example:
            Input:
                {
                    "contig_3": {"PR": [{"text": "L10I"}]},
                    "contig_4": {"PR": [{"text": "L10I"}]}
                }

            Output:
                {
                    "PR": [
                        {
                            "text": "L10I",
                            "on_contigs": ["contig_3", "contig_4"]
                        }
                    ]
                }

        Returns:
            dict
                Gene-centric dictionary where mutations are aggregated and
                annotated with the contigs on which they were observed.

        Notes:
            - Mutation uniqueness is determined using the 'text' field.
            - The 'on_contigs' list is sorted for deterministic output.
            - Mutation dictionaries are copied to avoid modifying the input data.
            - if sep=None will return a list i.e. ['contig_1', 'contig_2'] if set it will concat to string
        """
        gene_mutations = {}
        contig_sets = defaultdict(lambda: defaultdict(set))

        for contig, genes in data.items():
            for gene, muts in genes.items():
                for m in muts:
                    key = m["text"]  # unique mutation identifier

                    if gene not in gene_mutations:
                        gene_mutations[gene] = {}

                    if key not in gene_mutations[gene]:
                        gene_mutations[gene][key] = copy.deepcopy(m)

                    contig_sets[gene][key].add(contig)

        result = {}

        for gene, muts in gene_mutations.items():
            result[gene] = []
            for key, mut in muts.items():
                if isinstance(sep, str):
                    mut["on_contigs"] = sep.join(sorted(contig_sets[gene][key]))
                else:
                    mut["on_contigs"] = sorted(contig_sets[gene][key])
                result[gene].append(mut)

        return result

    def get_drm_mutations(
        self,
        query_names:list[str]|str|None = None,
        types: list[str]|str|None = None,
        simple_return: bool = True
    ):
        if types:
            if isinstance(types, str):
                types = [types]
            types = set(types)
        data = self.get_mutations(query_names, drm_only=True)
        if types is None:
            return {k: v for k, v in data.items() if v and any(v.values())}
        return_data = {
            contig: {
                gene: tuple(m for m in muts if m.get("primaryType") in types)
                for gene, muts in genes.items()
                if any(m.get("primaryType") in types for m in muts)
            }
            for contig, genes in data.items()
        }
        if simple_return:
            return {k: v for k, v in return_data.items() if v and any(v.values())}
        return return_data
    
    def get_drug_resistance_db_version(self):
        queries = self._query_names
        return list(zip(queries, [ rec.get('drugResistance') for rec in self._hit_records ]))


    def get_drug_resistance_panal(self):
        queries = self._query_names
        return list(zip(queries, [ rec.get('drugResistance') for rec in self._hit_records ]))

    @cache
    def get_result_file_hash(self, algorithm='sha256'):
        """Compute the hash of a file using the specified algorithm."""
        hash_func = hashlib.new(algorithm)
        
        with open(self._json_file_path, 'rb') as file:
            # Read the file in chunks of 8192 bytes
            while chunk := file.read(8192):
                hash_func.update(chunk)
        
        return hash_func.hexdigest()

    @cache
    def get_resutl_file_datetime(self):
        ts = os.path.getctime(self._json_file_path)
        return datetime.fromtimestamp(ts).astimezone().isoformat()
