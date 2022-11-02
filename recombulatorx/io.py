import numpy
import pandas
from networkx import DiGraph 

def ped2graph(path):
    """
    This function takes a ped file as input and build a graph with the relationships.
    It returns a list of tuples, each composed by the graph, a dictionary (with iid as key and their tab row as value) and the family identifier.
    """

    # read ped file
    # FIXME separator may be any space
    # FIXME na values should be 0 (or -9 for phenotype)
    # FIXME handle headerless peds
    if path.endswith('.tsv'):
        tab = pandas.read_csv(path, sep='\t', dtype=str, na_values=['0'], index_col=False)
        #tab = pandas.read_csv(path, sep='\t', dtype={'SEX': str}, na_values=['0'])
    elif path.endswith('.xlsx'):
        tab = pandas.read_excel(path)
    #else:
    #    tab = pandas.read_csv(path, sep='\t', dtype={'SEX': str})

    # check even number of columns
    assert tab.shape[1] % 2 == 0

    #tab.iloc[:, 6:] = tab.iloc[:, 6:].astype(float)

    hap1_markers = tab.columns[6::2]
    hap2_markers = tab.columns[7::2]
    markers = [col.split('-')[0] for col in hap1_markers] # FIXME: marker names
    
    # recode sex
    possible_sex_codings = [
        ['XY', 'XX'],
        ['1', '2'],
        ['M', 'F'],
        ['MALE', 'FEMALE'],
    ]
    orig_sex = tab['SEX'].str.upper()
    for sex_coding in possible_sex_codings:
        if orig_sex.isin(sex_coding).all():
            tab['SEX'] = orig_sex.replace(sex_coding, ['male', 'female'])
            break
    else:
        values = tab['SEX'].value_counts(sort=True).index.tolist()
        value_str = ', '.join(map(repr, values[:5]))
        if len(values) > 5:
            value_str += '...' 
        raise ValueError(f'Unknown values {value_str} in SEX column of "{path}"')


    # get graph nodes and edges
    fams_list = []
    for fid, fam in tab.groupby('FID'):
        #fam_dict = {}
        G = DiGraph()
        family_iids = set(fam['IID'])

        for _, row in fam.iterrows():
            iid = row.at['IID']
            pat = row.at['PAT']
            mat = row.at['MAT']        
            
            #fam_dict[iid] = row
            if row.at['SEX'] == 'female':
                haps = [
                    row.loc[hap1_markers],
                    row.loc[hap2_markers],
                ]
            else:
                haps = [row.loc[hap1_markers]]

            geno = numpy.array(haps).T
            try:
                geno = geno.astype(float)
                #print('STR markers detected')
            except ValueError:
                #print('SNP markers detected')
                raise NotImplementedError('SNP markers detected')


            G.add_node(iid, sex=row.at['SEX'], geno=geno)
            if not pandas.isnull(pat):
                if pat not in family_iids:
                    raise IOError(f'The father "{pat}" of "{iid}" is not in family "{fid}"')
                G.add_edge(pat, iid)
            if not pandas.isnull(mat):
                if mat not in family_iids:
                    raise IOError(f'The mother "{mat}" of "{iid}" is not in family "{fid}"')
                G.add_edge(mat, iid)
    
        fam_tuple = (fid, G)
        fams_list.append(fam_tuple)

    return fams_list, markers

