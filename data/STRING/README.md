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

```
mkdir processed
```

4. Create separate files for each STRING channel.  It uses `9606.protein.links.detailed.v11.0.txt` and `uniprot_to_string.tab` and outputs files in the `processed/` directory.

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
