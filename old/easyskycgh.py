# -*- encoding: utf-8 -*-
import shlex

def as_is(value):
    return value

def to_fields_dict(fields, opt_fields=None):
    if opt_fields is None: opt_fields = ()

    def _inner(value):
        values = shlex.split(value)

        num_values = len(values)
        num_fields = len(fields)
        num_opt_fields = len(opt_fields)

        # This work iff the optional fields are not alone
        # on a new continuation line (this is not the our case... for NOW!)
        if num_fields <= num_values <= (num_fields + num_opt_fields):
            return dict(zip(fields + opt_fields, values))
        else:
            return value

    return _inner

fields = {'SkyData': ('description', 'submitter', 'email'),
          'SkyCase': ('casename', 'organism', 'age', 'sex', 'imunoType',
                      'cellSource', 'disorderTumor', 'diseaseStage',
                      'diseaseStatus', 'tissue', 'treatment',
                      'hereditarySyn', 'priv'),
          'SkyCite': (('authors', 'title', 'journal', 'volume',
                       'year', 'issue', 'pages'), ('pubmedId',)),
          'SkyCell': (('chromoCt', 'ploidy', 'cellct', 'imagefile',
                       'isClone'), ('comments',)),
          'NormalChromosome': ('chromosome', 'normalCount', 'actualCount'),
          'SkyChromosome': (('karyotypelocation', 'aberCpyCnt'), ('recurCnt',)),
          'SkyFrag': (('position', 'parentChrom', 'fstart', 'fstop',
                       'isHsr', 'sizeEstimate'), ('gene',))
}

keywords = {'SkyData': to_fields_dict(fields['SkyData']),
            'SkyCase': to_fields_dict(fields['SkyCase']),
            'SkyCite': to_fields_dict(*fields['SkyCite']),
            'diagnosis': as_is,
            'icd_code': as_is,
            'site': as_is,
            'icd_site': as_is,
            'stage': as_is,
            'comments': as_is,
            'cyto_comments': as_is,
            'cell_line_name': as_is,
            'source_cat_no': as_is,
            'culture_med': as_is,
            'SkyCell': to_fields_dict(*fields['SkyCell']),
            'Karyotype': as_is,
            'NormalChromosome': to_fields_dict(fields['NormalChromosome']),
            'SkyChromosome': to_fields_dict(*fields['SkyChromosome']),
            'SkyFrag': to_fields_dict(*fields['SkyFrag']),
            'CGHSample': as_is,
            'CGHBin': as_is,
            'CGHBinFrag': as_is,
            'CGHFrag': as_is
}

def esi_parser(file_path):
    esi_file = open(file_path)

    out = dict()
    data = dict()
    prev_key = None
    for line in esi_file:
        # Pass on comments and blank line
        if not line.strip() or line.startswith('#'): continue

        # Key splitting
        key, value = line.strip().split(' ', 1)

        # New temporary structure
        if key == 'SkyData':
            if data:
                out[data['SkyCase']['casename']] = data
            data = dict()

        # Indented non-blank, non-comment and non-keyword line: continuation line
        if (line.startswith(' ') or line.startswith('\t')) and not key in keywords:
            key, value = prev_key, (data[prev_key] + line).strip()
        prev_key = key

        # Value saving (if the keyword is not present, save the value as is)
        result_value = keywords.get(key, as_is)(value)
        if key == 'NormalChromosome':
            data.setdefault(key, []).append(result_value)
        elif key == 'SkyChromosome':
            result_value['SkyFrag'] = list() # we prepare the list to be filled
            data.setdefault(key, []).append(result_value)
        elif key == 'SkyFrag':
            # Append to the last SkyChromosome
            data['SkyChromosome'][-1]['SkyFrag'].append(result_value)
        else:
            data[key] = result_value

    esi_file.close()
    return out

if __name__ == '__main__':
    import pprint
    out = esi_parser('NCI60_cell_line_panel_Genetics_Branch_I.R.Kirsch.esi')
    pprint.pprint(out['HT-29 (submitter:Roschke)'])

    print out.keys(), len(out)
