# working on dir:
# /Users/user/Documents/Doctorado/Analysis/Meta-analysis_COPD/GEO_data

import GEOparse
import padas as pd

gse = ['GSE106899','GSE106986','GSE76925','GSE37768','GSE57148','GSE56342','GSE46903','GSE47460','GSE37147','GSE38974','GSE13896','GSE8608','GSE3510','GSE1786','GSE1650','GSE1122','GSE475','GSE27543','GSE26296','GSE103174','GSE124180','GSE61397','GSE19407']

df = []

for i in gse[:2]:
    path = GEOparse.get_GEO_file(geo= i, destdir="./soft/")
    info = GEOparse.get_GEO(filepath=path[0])
    info = info.metadata
    df = pd.DataFrame(data=[df,info])
