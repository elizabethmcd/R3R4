import os 

# Read in lists
with open('refList.txt', 'r') as f:
    refGenList = f.read().splitlines()

with open('metagenomeList.txt', 'r') as f:
    metaList = f.read().splitlines()

# Make combos
with open('mappingCombos.txt', 'w') as outfile:
    for ref in refGenList:
        for meta in metaList:
            refbase = os.path.splitext(os.path.basename(ref))[0]
            metabase = os.path.splitext(os.path.basename(meta))[0]
            outname = 'mappingResults/{0}-vs-{1}.bam'.format(refbase, metabase)
            outfile.write('{0}\t{1}\t{2}\n'.format(ref,meta,outname))
