y_chr_dataset = {
    'Chagyrskaya2' : ['Chagyrskaya2',            'Neandertal',  ['hg19'], 80000],
    'Mezmaiskaya2' : ['Mezmaiskaya2\n~43 kya',   'Neandertal',  ['hg19'], 43449], # 44965 – 42156 [IntCal20]   43449
    'A00_I10871'   : ['A00 (Shum Laka)\n~8 kya', 'African',     ['hg19', 'A0b', 'T2T'], 7913], # 8008-7839  [IntCal20]      7913
    'A00'          : ['A00',                     'African',     ['hg19']],
    'HG02982'      : ['A0a1',                    'African',     ['hg19', 'A0b', 'T2T']],
    'HG02666'      : ['A1a',                     'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP01029'    : ['A1b1a1a1',                'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP01406'    : ['A1b1b2b',                 'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00931'    : ['B-M181',                  'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00478'    : ['B2a',                     'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00453'    : ['B2b2',                    'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00475'    : ['B2b1a1c',                 'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00984'    : ['B2b1a1a',                 'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00757'    : ['D1b2a2',                  'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP01214'    : ['D1a1a',                   'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP01031'    : ['E2b1a1',                  'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP01200'    : ['E1a2a1b1',                'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00908'    : ['E1b1a1a1c2c3b',           'African',     ['hg19', 'A0b', 'T2T']],
    'HGDP00103'    : ['C2b1c',                   'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP00213'    : ['G2b1',                    'non-African', ['hg19', 'A0b', 'T2T']],
    'Loschbour'    : ['I2 (Loschbour)\n~8 kya',  'non-African', ['hg19'], 8025], # 8171-7936   [IntCal20]     8025
    'HGDP00127'    : ['I2a2a1a2a2',              'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP00057'    : ['J1a2b',                   'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP00056'    : ['J2a1',                    'non-African', ['hg19', 'A0b', 'T2T']],
    'Ust_Ishim'    : ['K (Ust\'Ishim)\n~44 kya', 'non-African', ['hg19'], 44366], # 45930-42904 [IntCal20]     44366
    'HGDP00549'    : ['M1',                      'non-African', ['hg19', 'A0b', 'T2T']],
    'Yana1'        : ['P1 (Yana1)\n~32 kya',     'non-African', ['hg19'], 31850], # 32200 - 31532 [IntCal20]   31850
    'HGDP01009'    : ['Q1a2a1a1',                'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP00136'    : ['R1a1a1b2',                'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP01075'    : ['R1b1a2a1a2b',             'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP00218'    : ['R2a',                     'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP01298'    : ['L1a2',                    'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP01192'    : ['N1',                      'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP01225'    : ['O2a2b1a1a',               'non-African', ['hg19', 'A0b', 'T2T']],
    'HGDP01190'    : ['O1b1a1a1a',               'non-African', ['hg19', 'A0b', 'T2T']],
}

def GetHaplogroup(key):
    return y_chr_dataset[key][0]

def GetPop(key):
    return y_chr_dataset[key][1]

def GetReference(key):
    return y_chr_dataset[key][2]

def GetAge(key):
    return y_chr_dataset[key][3]

samples = list(y_chr_dataset.keys())
ancient_samples = ['Chagyrskaya2', 'Mezmaiskaya2', 'Ust_Ishim', 'Yana1', 'Loschbour', 'A00_I10871']
non_african_samples = [x for x in samples if ((GetPop(x) == 'non-African') & (x not in ancient_samples))]
