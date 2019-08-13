## STRING Interactions

### Mapped STRING interactions
This directory contains processed STRING interactions.  

1. Download detailed STRING links filtered by human from the [STRING Downloads site](https://string-db.org/cgi/download.pl).

```
wget https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
gunzip 9606.protein.links.detailed.v11.0.txt.gz
```

2. Download [Uniprot to STRING identifiers](https://www.uniprot.org/uniprot/?query=database:(type:string)&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22).  Click the link and click `Download.'  Rename this file as `uniprot_to_string.tab`. (this is already committed in the repo).


3. Make a `processed/` directory.

```mkdir processed```

4. Create separate files for each STRING channel.  It uses `9606.protein.links.detailed.v11.0.txt` and `uniprot_to_string.tab`

```
python3 process_string.py 
```

The output is the number of interactions for each channel:

```
neighborhood 396132
fusion 23388
cooccurence 60404
coexpression 6080248
experimental 4574384
database 706440
textmining 9653980
combined_score 11759454
```



##########################################################
old code

awk '($6>0) {print $1,$2,$6}' 9606.protein.links.full.v11.0.txt | head

wc -l 9606.protein.links.full.v11.0.txt
# 11759455 9606.protein.links.full.v11.0.txt
# 11,759,455 (11.76M)


awk '($6>0) {print $1,$2,$6}' 9606.protein.links.full.v11.0.txt > 96066.protein.links.cooccurence.v11.0.txt
wc -l 96066.protein.links.cooccurence.v11.0.txt
# 60405 96066.protein.links.cooccurence.v11.0.txt
# 60,405

##Download mapping files from STRING (human.uniprot_2_string.2018)

## still missing many mappings (994), map 22 proteins by hand. Search for them in UniProt.
## human.uniprot_2_string.MISSING.tsv 

## map STRING Ids to UniPROT IDs.
python3 string2uniprot.py 96066.protein.links.cooccurence.v11.0.txt

60298 written, 106 interactions skipped, and 5 nodes skipped
   9606.ENSP00000433415
   9606.ENSP00000417240
   9606.ENSP00000439534
   9606.ENSP00000464265
   9606.ENSP00000347739


#########