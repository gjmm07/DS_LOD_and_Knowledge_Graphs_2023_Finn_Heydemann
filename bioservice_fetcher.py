import time
from bioservices import EUtils
import pandas as pd
import urllib
import numpy as np
from operator import itemgetter

_s = EUtils(email="finn.heydemann@th-koeln.de")


def _fetch_refseq_nucleotid(refseq_nucleotid_id: str) -> dict[str, str]:
    """
    Based on a ncbi refseq nucleotid ID the encoded protein, gene and bacterial strain is returned 
    """
    found_data = {"gene": "", "protein": "", "genome": "", "organism": ""}
    try:
        raw_data = _s.EFetch("nucleotide", refseq_nucleotid_id, retmode="dict", rettype="summary")["GBSet"]["GBSeq"]
    except AttributeError:
        return found_data
    except TypeError:
        print(refseq_nucleotid_id)
    time.sleep(0.5)  # wait to ensure to not be blocked
    definition = raw_data["GBSeq_definition"]
    if "gene" in definition:
        gene, protein = definition.split("gene for")
        protein = protein.replace(", complete CDS", "")
        organism = raw_data["GBSeq_organism"]
        found_data["gene"] = gene.strip()
        found_data["protein"] = protein.strip()
        found_data["organism"] = organism.strip()
    elif "genome" in definition:
        genome = definition.split(", complete genome")[0]
        genome = genome.split(", whole genome")[0]
        found_data["genome"] = genome.strip()
    return found_data


def _fetch_refseq_protein(refseq_protein_id: str) -> tuple[str, str, str]:
    """
    Based on the ncbi refseq protein ID the parent taxon, the protein name is returned
    If needed the same parent taxon is returned based on another entry
    """
    try:
        data = _s.EFetch("protein", refseq_protein_id, retmode="dict", rettype="summary")["GBSet"]["GBSeq"]
    except AttributeError:
        return "", "", ""
    time.sleep(0.5)  # wait to ensure to not be blocked
    parent_taxon: str = data["GBSeq_organism"].strip()
    definition = data["GBSeq_definition"]
    start_index, stop_index = (n := definition.find("[")), definition.find("]", n)
    same_parent_taxon: str = definition[start_index + 1: stop_index].strip()
    protein = (definition[:start_index] + definition[stop_index + 1:])
    protein: str = protein.replace("MULTISPECIES: ", "").strip()
    return parent_taxon, protein, same_parent_taxon


def _fetch_genbank_protein(genbank_protein: str) -> tuple[str, str]:
    strain, organism = "", ""
    try:
        data = _s.EFetch("protein", genbank_protein, retmode="dict", rettype="summary")["GBSet"]["GBSeq"]
    except AttributeError:
        return organism, strain
    # GB source
    time.sleep(0.5)
    for t in data["GBSeq_feature-table"]["GBFeature"]:
        if t["GBFeature_key"] == "source":
            for x in t["GBFeature_quals"]["GBQualifier"]:
                if x["GBQualifier_name"] == "organism":
                    organism = x["GBQualifier_value"]
                if x["GBQualifier_name"] == "strain":
                    strain = x["GBQualifier_value"]
    return organism, strain


def _fetch_genbank_nucleotide(genbank_nucleotide: str) -> tuple[str, str]:
    time.sleep(0.5)
    data = _s.EFetch("nucleotide", genbank_nucleotide, retmode="dict", rettype="summary")
    data = data.get("GBSet").get("GBSeq").get("GBSeq_feature-table").get("GBFeature")
    strain, organism = "", ""
    for t in data:
        if t["GBFeature_key"] == "source":
            for x in t["GBFeature_quals"]["GBQualifier"]:
                if x["GBQualifier_name"] == "organism":
                    organism = x["GBQualifier_value"]
                if x["GBQualifier_name"] == "strain":
                    strain = x["GBQualifier_value"]
    return organism, strain


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
    parent_taxon, protein, *_ = _fetch_refseq_protein(df_row.refseq_protein_accession)
    return pd.Series([parent_taxon, protein])


def get_strain_and_gene(df_row: pd.Series) -> pd.Series:
    # _, gene, strain = _fetch_refseq_nucleotid(df_row.refseq_nucleotide_accession)
    return pd.Series(itemgetter("gene",
                                "protein",
                                "genome",
                                "organism")(_fetch_refseq_nucleotid(df_row.refseq_nucleotide_accession)))


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
    url = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"
    df = pd.read_csv(urllib.request.urlopen(url), delimiter="\t")
    #
    # sample_df = df.sample(5, random_state=42)

    print(_fetch_refseq_nucleotid("NG_052278.1"))
    print(df)

    # df[["refseq_gene",
    #     "refseq_protein",
    #     "refseq_genome",
    #     "refseq_organism"]] = df.apply(get_strain_and_gene, axis=1)

    print(0)
    # df[["refseq_parent", "rrefseq_protein"]] = df.apply(get_protein_and_parent, axis=1)

    print(1)

    # df[["genbank_organsim_nuc", "genbank_strain_nuc"]] = df.apply(get_organism_strain_via_nuc, axis=1)

    print(2)

    # df[["genbank_organsim_prot", "genbank_strain_prot"]] = df.apply(get_organism_strain_via_prot, axis=1)

    df.to_csv("samples.csv")
