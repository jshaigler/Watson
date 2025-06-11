from Bio import Entrez, SeqIO
from Bio.SeqUtils import seq3, nt_search
from Bio.SeqUtils import GC123
import RNA
import re
from collections import Counter

class ProteinToCRISPR:
    def __init__(self, email):
        self.email = email
        Entrez.email = email

    def fetch_gene_sequence(self, gene_name):
        """Fetch the DNA sequence for a given gene from NCBI."""
        try:
            search_term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if not record['IdList']:
                return None

            gene_id = record['IdList'][0]
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()

            return seq_record

        except Exception as e:
            print(f"Error fetching gene sequence: {e}")
            return None

    def find_pam_sites(self, sequence, pam="NGG"):
        """Identify all possible PAM sites in the sequence."""
        pam_regex = pam.replace("N", "[ACGT]")
        return [match.start() for match in re.finditer(pam_regex, str(sequence))]

    def design_grna(self, sequence, pam_sites, edit_function):
        """Design a gRNA based on the PAM site and edit function."""
        for pam_site in pam_sites:
            if pam_site < 20:
                continue

            grna = str(sequence[pam_site - 20:pam_site])
            pam = str(sequence[pam_site:pam_site + 3])
            
            # Basic validation of gRNA sequence
            if 'N' not in grna and len(grna) == 20:
                return grna, pam

        return None, None

    def get_crispr_design(self, gene_name, edit_function):
        """Complete pipeline for CRISPR design."""
        try:
            # Fetch gene sequence
            seq_record = self.fetch_gene_sequence(gene_name)
            if not seq_record:
                return {
                    "success": False,
                    "error": f"Could not find gene sequence for: {gene_name}"
                }

            # Find PAM sites
            pam_sites = self.find_pam_sites(seq_record.seq)
            if not pam_sites:
                return {
                    "success": False,
                    "error": "No PAM sites found in the sequence."
                }

            # Design gRNA
            grna, pam = self.design_grna(seq_record.seq, pam_sites, edit_function)
            if not grna:
                return {
                    "success": False,
                    "error": "Could not design suitable gRNA."
                }

            # Predict off-target effects (simplified)
            off_targets = len(nt_search(str(seq_record.seq), grna)[1:])

            # Calculate secondary structure
            sec_struct = RNA.fold(grna)[0]

            # Select appropriate Cas protein
            cas_protein = self.select_cas_protein(edit_function)

            # Suggest delivery systems
            delivery = self.suggest_delivery_system(edit_function)

            return {
                "success": True,
                "Gene": gene_name,
                "gRNA": grna,
                "PAM": pam,
                "Cas Protein": cas_protein,
                "Off-Target Effects": {"Count": off_targets},
                "Secondary Structure": sec_struct,
                "Suggested Delivery Systems": delivery
            }

        except Exception as e:
            print(f"Error in CRISPR design: {e}")
            return {"success": False, "error": str(e)}

    def select_cas_protein(self, edit_function):
        cas_proteins = {
            "delete": "Cas9",
            "modify": "Cas9-HF1",
            "insert": "Cas12a",
            "activate": "dCas9-VP64",
            "interfere": "dCas9-KRAB"
        }
        return cas_proteins.get(edit_function.lower(), "Cas9")

    def suggest_delivery_system(self, edit_function):
        delivery_systems = {
            "delete": ["Plasmid", "RNP Complex"],
            "modify": ["Lentiviral Vector", "RNP Complex"],
            "insert": ["Adeno-Associated Virus (AAV)", "RNP Complex"],
            "activate": ["Lentiviral Vector", "Plasmid"],
            "interfere": ["Lentiviral Vector", "Plasmid"]
        }
        return delivery_systems.get(edit_function.lower(), ["Plasmid"])