import os

ucsc_dir = os.path.dirname(__file__)
UCSC = {'hg18': dict(), 'hg19': dict()}

for k in UCSC:
    for file in os.listdir(os.path.join(ucsc_dir, k)):
        name, ext = os.path.splitext(file)
        UCSC[k][name] = os.path.join(ucsc_dir, k, file)
