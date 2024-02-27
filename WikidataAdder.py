import pywikibot
from SPARQLWrapper import SPARQLWrapper, JSON
import pandas as pd
from typing import Optional


sparql = SPARQLWrapper("https://query.wikidata.org/sparql")
site = pywikibot.Site("wikidata")
repo = site.data_repository()


def search_for_exist(wd_key: str):
    query = f""" SELECT ?item WHERE 
                    {{?item rdfs:label "{wd_key}"@en}}
                    """
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert().get("results").get("bindings")
    print(results)
    return results


def _add_string_claim(page: pywikibot.ItemPage, p_id: str, s: str, summary: Optional[str] = ""):
    if p_id not in page.claims:
        stringclaim = pywikibot.Claim(repo, p_id)
        stringclaim.setTarget(s)
        page.addClaim(stringclaim, summary=summary)


def _add_claim(page: pywikibot.ItemPage, p_id: str, q_id: str, summary: Optional[str] = ""):
    if p_id in page.claims:
        existing_claims = [claim.target.getID() for claim in page.claims[p_id]]
        if q_id in existing_claims:
            print("claim already existent")
            return
    claim = pywikibot.Claim(repo, p_id)
    target = pywikibot.ItemPage(repo, q_id)
    claim.setTarget(target)
    page.addClaim(claim, summary=summary)


def _add_identifier(page: pywikibot.ItemPage, p_id: str, identifier, summary: Optional[str] = ""):
    if p_id in page.claims:
        existing_identifiers = [claim.target for claim in page.claims[p_id]]
        if identifier in existing_identifiers:
            return
    new_claim = pywikibot.Claim(repo, p_id)
    new_claim.setTarget(identifier)
    page.addClaim(new_claim, summary=summary)


def _edit_description(page: pywikibot.ItemPage, description: str):
    if len(page.descriptions):
        if page.descriptions["en"] == description:
            print("description already exists")
            return
    page.editDescriptions({"en": description}, summary="Setting new description")


def _edit_aliases(page: pywikibot.ItemPage, aliases: list[str]):
    if len(page.aliases):
        if aliases[0] in page.aliases["en"]:
            print("alias already exists")
            return
    page.editAliases({"en": aliases}, summary="Setting new alias")


class StrainWikidataAdder:

    def __init__(self, df_row: pd.Series, sim: bool):
        self.df_row = df_row
        self.strain_page: pywikibot.ItemPage | None = None
        self.sim = sim  # for debugging purposes

    def create_strain(self):
        if self.sim:
            print("Label: " + self.df_row["genbank_organism2"])
            print("Description: " + "bacterial strain")
            print("Alias: " + self.df_row["genbank_organism2"].split()[-1])
        else:
            results = search_for_exist(self.df_row["genbank_organism2"])
            if len(results) >= 2:
                print("more than one item found")
                self.sim = True
                return
            elif len(results) == 1:
                self.strain_page = pywikibot.ItemPage(repo, results[0].get("item").get("value").split("/")[-1])
            elif not results:
                self.strain_page = pywikibot.ItemPage(site)
                self.strain_page.editLabels({"en": self.df_row["genbank_organism2"]}, summary="Setting new label")
            _edit_description(self.strain_page, "bacterial strain")
            _edit_aliases(self.strain_page, [self.df_row["genbank_organism2"].split()[-1]])
            print(self.strain_page.getID())

    def add_instance_of_strain(self):
        if not self.sim:
            _add_claim(self.strain_page,
                       u"P31",
                       u"Q855769",
                       "add instance of strain")
        else:
            print("Instance of strain")

    def add_taxon_name(self):
        if not self.sim:
            _add_string_claim(self.strain_page,
                              u"P225",
                              self.df_row["genbank_organism2"],
                              "adding taxon name")
        else:
            print("taxon name: " + self.df_row["genbank_organism2"])

    def add_parent_taxon(self):
        if not self.sim:
            _add_claim(self.strain_page,
                       u"P171",
                       self.df_row["parent_taxon"].split("/")[-1],
                       "adding parent taxon")
        else:
            print("parent_taxon: " + self.df_row["parent_taxon"].split("/")[-1])

    def add_ncbi_taxonomy_id(self):
        ncbi_taxonomy_id = self.df_row["genbank_tax_id"].split(":")[-1].strip()
        if not self.sim:
            _add_identifier(self.strain_page,
                            "P685",
                            ncbi_taxonomy_id,
                            "adding taxonomy ID")
        else:
            print("taxon ID: " + ncbi_taxonomy_id)


class GeneWikidataAdder:

    def __init__(self, df_row: pd.Series, sim: bool):
        self.df_row = df_row
        self.sim = sim
        self.gene_page: pywikibot.ItemPage | None = None

    def create_gene(self):
        strain = self.df_row["genbank_organism2"]
        if self.sim:
            print("Label: " + self.df_row["refseq_gene"])
            print("Desciption: " + f"microbial gene found in {strain}")
        else:
            results = search_for_exist(self.df_row["refseq_gene"])
            if len(results) >= 2:
                print("more than one item found")
                self.sim = True
                return
            elif len(results) == 1:
                self.gene_page = pywikibot.ItemPage(repo, results[0].get("item").get("value").split("/")[-1])
            elif not results:
                self.gene_page = pywikibot.ItemPage(site)
                self.gene_page.editLabels({"en": self.df_row["refseq_gene"]}, summary="Setting new label")
            _edit_description(self.gene_page, f"microbial gene found in {strain}")
            # _edit_aliases(self.gene_page, [self.df_row["genbank_organism2"].split()[-1]])
            print(self.gene_page.getID())

    def add_instance_of_gene(self):
        if not self.sim:
            _add_claim(self.gene_page,
                       u"P31",
                       u"Q7187",
                       "add instance of gene")
        else:
            print("Instance of gene")

    def add_subclass_of(self):
        if not self.sim:
            _add_claim(self.gene_page,
                       u"P279",
                       u"Q7187",
                       "add subclass of gene")
            _add_claim(self.gene_page,
                       u"P279",
                       u"Q20747295",
                       "add subclass of protein coding gene")
        else:
            print("subclass of gene and protein-coding-gene")

    def add_found_in_taxon(self):
        if not self.sim:
            _add_claim(self.gene_page,
                       u"P703",
                       self.df_row["strain_wd_id"],
                       "add found in taxon")
        else:
            print("found in taxon: " + self.df_row["strain_wd_id"])

    def add_genbank_start(self):
        # I'm not sure if this makes any sense if the chromosome is not included --> will skip
        pass

    def add_genbank_end(self):
        # same as above
        pass

    def add_entrez_gene_id(self):
        # Property P351
        pass

    def encodes_protein(self):
        # Encodes a protein -- protein is only present in next step
        pass


class ProteinWikidataAdder:

    def __init__(self, df_row: pd.Series, sim: bool):
        self.df_row = df_row
        self.sim = sim
        self.protein_page: pywikibot.ItemPage | None = None

    def create_protein(self):
        strain = self.df_row["genbank_organism2"]
        if self.sim:
            print("Label: " + self.df_row["product_name"])
            print("Desciption: " + f"microbial protein found in {strain}")
        else:
            results = search_for_exist(self.df_row["product_name"])
            if len(results) >= 2:
                print("more than one item found")
                self.sim = True
                return
            elif len(results) == 1:
                self.protein_page = pywikibot.ItemPage(repo, results[0].get("item").get("value").split("/")[-1])
            elif not results:
                self.protein_page = pywikibot.ItemPage(site)
                self.protein_page.editLabels({"en": self.df_row["product_name"]}, summary="Setting new label")
            _edit_description(self.gene_page, f"microbial protein found in {strain}")
            # _edit_aliases(self.gene_page, [self.df_row["genbank_organism2"].split()[-1]])
            print(self.gene_page.getID())

    def instance_of(self):
        if not self.sim:
            pass
        else:
            print("Instance of (P31) Protein (Q8054)")

    def subclass_of(self):
        if not self.sim:
            pass
        else:
            print("Subclass of (P279) Protein (Q8054)")

    def found_in_taxon(self):
        pass

    def encoded_by(self):
        pass

    def add_identifier(self):
        pass

    def add_antibiotic_resistance(self):
        pass





