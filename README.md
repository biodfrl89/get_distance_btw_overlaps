# get_distance_btw_overlaps

Rscript to get the distance between the borders of overlaping sequences.

Once the redundancy has been cleaned and the manual annotation BED file as been filtered for the desired genetic feature to analize, this script will produce a plot for the distances between the 5' and 3' ends of the genetic features in the manual annotation file vs. the reduced file.

## Syntax

Rscript get_distance_btw_overlaps.R -q [FILENAME] -s [FILENAME] -f [STRING] -o [FILENAME]

## Options

| Short Option | Long option | Description |
| --- | --- | --- |
| -q | --query | Manual anotation with the genetic feature to be analyzed, as a BED formated file. For cleaner results, filter your BED file for the desired genetic feature. |
| -s | --subject | Reduced GFF file, produced by gff_disambiguation.R script. |
| -f | --feature | The name of the genetic feature to be analized. |
| -o | --outfile | Filename for the plot to be produced. |
