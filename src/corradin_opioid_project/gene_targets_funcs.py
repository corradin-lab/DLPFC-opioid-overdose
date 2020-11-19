import pandas as pd
import pybedtools

def get_gene_targets(peak_list, region_split_char, gene_target_file):
    bed_info = pd.DataFrame(peak_list, columns = ["peak"])
    bed_info[["chromosome", "start", "stop"]] = bed_info.peak.str.split(region_split_char, expand = True)
    bed_info = bed_info.drop(columns = "peak")
    display(bed_info.head())

    bed_file = pybedtools.BedTool.from_dataframe(bed_info)
    intersected = bed_file.intersect(gene_target_file, wao=True).to_dataframe().query("score != -1")
    display(intersected.head())

    intersected_exploded = intersected.assign(thickStart=intersected.thickStart.str.split(";")).explode("thickStart")
    intersected_exploded.drop_duplicates(["chrom", "start", "end", "thickStart"])
    return intersected_exploded
