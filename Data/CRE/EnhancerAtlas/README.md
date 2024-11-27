Human keratinocyte enhancer-gene interactions downloaded from the (EnhancerAtlas database)[http://enhanceratlas.org/downloadv2.php] on 2024-06-19
Mouse kidney enhancer-gene interactions downloaded from the (EnhancerAtlas database)[http://enhanceratlas.org/downloadv2.php] on 2024-07-01
Fly BG3-c2 enhancer-gene interactions downloaded from the (EnhancerAtlas database)[http://enhanceratlas.org/downloadv2.php] on 2024-07-03

Data were processed using the following commands:
```
cut -f 1 enhancer_atlas_human_keratinocyte | sort | uniq | cut -d $ -f 1 | cut -d _ -f 2 | sort | uniq -c > human_keratinocyte_enhancer_counts_per_gene.csv
cut -f 1 enhancer_atlas_mouse_kidney | sort | uniq | cut -d $ -f 1 | cut -d _ -f 2 | sort | uniq -c > mouse_kidney_enhancer_counts_per_gene.csv
cut -f 1 enhancer_atlas_fly_bg3c2| sort | uniq | cut -d $ -f 1 | cut -d _ -f 2 | sort | uniq -c > fly_bg3c2_enhancer_counts_per_gene.csv
```