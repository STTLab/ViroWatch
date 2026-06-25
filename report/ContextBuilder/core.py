from os import PathLike
import json

class BLASTResult:
    def __init__(self, data):
        self._dump = data

    @classmethod
    def read_json(cls, f_path: PathLike):
        with open(f_path, 'r', encoding='utf-8') as f:
            data = json.loads(f.read())
        return cls(data)
    
    def _filter_by_query_name(self, query_name):
        '''
        Return results section of query with that name
        '''
        _filtered = tuple(filter(lambda q: q.get('report').get('results').get('search').get('query_title')==query_name, self._dump.get('BlastOutput2')))
        if len(_filtered) == 0: raise ValueError('query name not found in the result')
        return _filtered[0].get('report').get('results')

    def dump(self):
        return self._dump
    
    def get_query_names(self):
        return tuple([ q.get('report').get('results').get('search').get('query_title') for q in self._dump.get('BlastOutput2') ])

    def get_query_length(self, query_name):
        filtered = self._filter_by_query_name(query_name)
        return filtered.get('search').get('query_len')
    
    def _flatten_hsps(self, hsps_list):
        for hsp in hsps_list:
            pass

    def get_hits(self, query_name, n:int=-1):
        filtered = self._filter_by_query_name(query_name)
        query_len = self.get_query_length(query_name)
        hit_records = []
        for hit in filtered.get('search').get('hits'):
            description = hit.get('description', [{}])[0]
            hsps = hit.get('hsps')
            aligned_len = sum([hsp.get('align_len') for hsp in hsps])
            subject_len = hit.get('len')
            rec = {
                'query_name': query_name,
                'accession': description.get('accession'),
                'subject_title': description.get('title'),
                'taxid': description.get('taxid'),
                'sciname': description.get('sciname'),
                'subject_length': subject_len,
                'aligned_length': aligned_len,
                'piden': round(sum([hsp.get('identity') for hsp in hsps])/aligned_len, 3),
                'query_fraction': round(aligned_len/query_len, 3),
                'subject_fraction': round(aligned_len/subject_len, 3),
                'hsps_query': [(hsp.get('query_strand'), hsp.get('query_from'), hsp.get('query_to')) for hsp in hsps],
                'hsps_subject': [(hsp.get('hit_strand'), hsp.get('hit_from'), hsp.get('hit_to')) for hsp in hsps],
                'evalues_per_hsp': [hsp.get('evalue') for hsp in hsps],
                'bit_score': sum([hsp.get('bit_score') for hsp in hsps])
            }
            rec['hsps_compact'] = tuple(zip(rec['hsps_subject'], rec['hsps_query'], rec['evalues_per_hsp']))
            hit_records.append(rec)
            if n > 0 and len(hit_records) == n:
                break
        return hit_records
