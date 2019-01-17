import requests.packages.urllib3
import datetime
from ga4gh.client import client
from itertools import chain
requests.packages.urllib3.disable_warnings()


httpClient = client.HttpClient("https://brcaexchange.org/backend/data/ga4gh/v0.6.0a7/")
chrom = {
    "BRCA1": "chr17",
    "BRCA2": "chr13"
}

annotCols = ['id', 'Pathogenicity_all']

def brca_query(gene, start_pos, end_pos):
    
    query = httpClient.search_variants(reference_name = chrom[gene], variant_set_id = "brca-hg37",
                                        start = int(start_pos), end = int(end_pos)+1) 

    listOutput = []
    for var in query:
        posInfo = [
            var.info['Gene_Symbol'].values[0].string_value,
            var.reference_name,
            var.start,
            var.end,
            str(var.reference_bases),
            str(var.alternate_bases[0])
            ]
        annotInfo = [var.info[x].values[0].string_value if x in var.info.keys() else '' for x in annotCols]
        out = posInfo + annotInfo
        listOutput.append(out)
    
    return(listOutput)
