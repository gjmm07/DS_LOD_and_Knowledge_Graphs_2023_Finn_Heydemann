import time
from bioservices import EUtils
import pandas as pd
from operator import itemgetter
import urllib
from functools import lru_cache


_s = EUtils(email="finn.heydemann@th-koeln.de")


@lru_cache
def _fetch_refseq_nucleotide(refseq_nucleotid_id: str) -> dict[str, str | float]:
    print(refseq_nucleotid_id)
    """
    Based on a ncbi refseq nucleotid ID the encoded protein, gene,
    bacterial strain and NCBI taxonomy ID is returned
    """
    found_data: dict[str, str] = {"gene": float("NaN"),
                                  "protein": float("NaN"),
                                  "genome": float("NaN"),
                                  "organism": float("NaN"),
                                  "taxonomy ID": float("NaN")}
    time.sleep(0.5)  # wait to ensure to not be blocked
    try:
        raw_data = _s.EFetch("nucleotide", refseq_nucleotid_id, retmode="dict", rettype="summary")["GBSet"]["GBSeq"]
    except (AttributeError, TypeError):
        return found_data
    definition = raw_data["GBSeq_definition"]
    if "gene" in definition:
        try:
            gene, protein = definition.split("gene for")
        except ValueError:
            gene, protein = definition.split("gene")
        idx = protein.upper().find(", COMPLETE CDS")
        protein = protein[:idx]
        # protein = protein.upper().replace(", COMPLETE CDS", "")
        organism = raw_data["GBSeq_organism"]
        found_data["gene"] = gene.strip()
        found_data["protein"] = protein.strip()
        found_data["organism"] = organism.strip()
    elif "genome" in definition:
        genome = definition.split(", complete genome")[0]
        genome = genome.split(", whole genome")[0]
        found_data["genome"] = genome.strip()
    # Looking for NCBI taxonomy ID
    items = raw_data["GBSeq_feature-table"]["GBFeature"]
    items = items if isinstance(items, list) else [items]
    for z in items:
        if z["GBFeature_key"] == "source":
            for y in z["GBFeature_quals"]["GBQualifier"]:
                if y["GBQualifier_name"] == "db_xref":
                    found_data["taxonomy ID"] = y["GBQualifier_value"]
    return found_data


@lru_cache
def _fetch_refseq_protein(refseq_protein_id: str) -> tuple[str | float, str | float, str | float]:
    print(refseq_protein_id)
    """
    Based on the ncbi refseq protein ID the parent taxon, the protein name is returned
    If needed the same parent taxon is returned based on another entry
    """
    time.sleep(0.5)  # wait to ensure to not be blocked
    try:
        data = _s.EFetch("protein", refseq_protein_id, retmode="dict", rettype="summary")["GBSet"]["GBSeq"]
        print("hi")
    except (AttributeError, TypeError):
        return float("NaN"), float("NaN"), float("NaN")
    parent_taxon: str or list[str] = data["GBSeq_organism"]
    if isinstance(parent_taxon, list):
        parent_taxon = parent_taxon[-1]
    parent_taxon = parent_taxon.strip()
    definition = data["GBSeq_definition"]
    start_index, stop_index = (n := definition.find("[")), definition.find("]", n)
    same_parent_taxon: str = definition[start_index + 1: stop_index].strip()
    protein = (definition[:start_index] + definition[stop_index + 1:])
    protein: str = protein.replace("MULTISPECIES: ", "").strip()
    return parent_taxon, protein, same_parent_taxon


@lru_cache
def _fetch_genbank_protein(genbank_protein: str) -> tuple[str | float, str | float]:
    print(genbank_protein)
    strain, organism = float("NaN"), float("NaN")
    time.sleep(0.5)
    try:
        data = _s.EFetch("protein", genbank_protein, retmode="dict", rettype="summary")["GBSet"]["GBSeq"]
    except (AttributeError, TypeError):
        return organism, strain
    # GB source
    for t in data["GBSeq_feature-table"]["GBFeature"]:
        if t["GBFeature_key"] == "source":
            for x in t["GBFeature_quals"]["GBQualifier"]:
                if x["GBQualifier_name"] == "organism":
                    organism = x["GBQualifier_value"]
                if x["GBQualifier_name"] == "strain":
                    strain = x["GBQualifier_value"]
    return organism, strain


@lru_cache
def _fetch_genbank_nucleotide(genbank_nucleotide: str) -> tuple[str | float, str | float, str | float]:
    print(genbank_nucleotide)
    strain, organism, taxonomy_id = float("NaN"), float("NaN"), float("NaN")
    while True:
        time.sleep(0.5)
        try:
            data = _s.EFetch("nucleotide", genbank_nucleotide, retmode="dict", rettype="summary")
        except (AttributeError, TypeError):
            return organism, strain, taxonomy_id
        if data is not None:
            data = data.get("GBSet").get("GBSeq").get("GBSeq_feature-table").get("GBFeature")
            break
        else:
            continue
    if not isinstance(data, list):
        data = [data]
    for t in data:
        if t["GBFeature_key"] == "source":
            for x in t["GBFeature_quals"]["GBQualifier"]:
                match x["GBQualifier_name"]:
                    case "organism": organism = x["GBQualifier_value"]
                    case "strain": strain = x["GBQualifier_value"]
                    case "db_xref": taxonomy_id = x["GBQualifier_value"]
    #elif isinstance(data, dict):
    #    if data["GBFeature_key"] == "source":
    #        for x in data["GBFeature_quals"]["GBQualifier"]:
    #            if x["GBQualifier_name"] == "organism":
    #                organism = x["GBQualifier_value"]
    #            if x["GBQualifier_name"] == "strain":
    #                strain = x["GBQualifier_value"]
    return organism, strain, taxonomy_id


"""
################################
################################
#####    Callable Funcs    #####
################################
################################ 
"""


def get_protein_and_parent(df_row: pd.Series) -> pd.Series:
    """
    To a pandas dataframe row this function return a pd Series with the bacteria parent taxon and the protein it encodes
    """
    parent_taxon, protein, same_parent_taxon = _fetch_refseq_protein(df_row.refseq_protein_accession)
    return pd.Series([parent_taxon, protein, same_parent_taxon])


def get_strain_and_gene(df_row: pd.Series) -> pd.Series:
    # _, gene, strain = _fetch_refseq_nucleotid(df_row.refseq_nucleotide_accession)
    return pd.Series(itemgetter("gene",
                                "protein",
                                "genome",
                                "organism",
                                "taxonomy ID")(_fetch_refseq_nucleotide(df_row.refseq_nucleotide_accession)))


def get_organism_strain_via_prot(df_row: pd.Series) -> pd.Series:
    """
    return organism and strain
    """
    return pd.Series(_fetch_genbank_protein(df_row.genbank_protein_accession))


def get_organism_strain_via_nuc(df_row: pd.Series) -> pd.Series:
    """
    return organism and strain
    """
    return pd.Series(_fetch_genbank_nucleotide(df_row.genbank_nucleotide_accession))


if __name__ == "__main__":
    # url = ("https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/"
    #        "AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt")
    # df = pd.read_csv(urllib.request.urlopen(url), delimiter="\t")
    # df[["genbank_organism2", "genbank_strain2", "genbank_tax_id"]] = df.apply(get_organism_strain_via_nuc, axis=1)
    # df.to_csv("resistance_dfD.csv", index=False)
    print(_fetch_genbank_nucleotide("AE004091.2"))
    print(_fetch_genbank_nucleotide("AE004091.2"))

