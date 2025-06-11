import requests
from Bio import Entrez, SeqIO
import re

class ProteinToSiRNA:
    def __init__(self, email):
        self.email = email
        Entrez.email = email

    def fetch_mrna_from_protein(self, protein_name):
        """
        Uses Entrez to find human mRNA corresponding to a given protein.
        Returns the mRNA sequence (as string), RefSeq ID, and description.
        """
        try:
            search_term = f"{protein_name}[Protein Name] AND Homo sapiens[Organism] AND mRNA[Filter]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmode="xml", retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                return None, None, None

            nucleotide_id = record["IdList"][0]
            handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="text")
            seq_record = SeqIO.read(handle, "genbank")
            handle.close()

            seq = str(seq_record.seq).upper()
            description = seq_record.description
            accession = seq_record.id
            return seq, accession, description

        except Exception as e:
            print(f"Error fetching mRNA: {e}")
            return None, None, None

    def design_sirna(self, mrna_seq):
        """
        Scans real mRNA for siRNA candidates with proper GC content (30-60%),
        returns sense + antisense strands with positions and GC%.
        """
        candidates = []
        for i in range(0, len(mrna_seq) - 21):
            seq = mrna_seq[i:i+21]
            if "N" in seq:
                continue
            gc = (seq.count("G") + seq.count("C")) / len(seq)
            if 0.3 <= gc <= 0.6:
                candidates.append({
                    "target_region": i,
                    "siRNA_sense": seq,
                    "siRNA_antisense": self.reverse_complement_rna(seq),
                    "GC_content": round(gc * 100, 2)
                })
        return candidates[:5]

    def reverse_complement_rna(self, seq):
        """Returns the antisense strand (reverse complement) of the siRNA sense."""
        complement = str.maketrans("AUCG", "UAGC")
        return seq.translate(complement)[::-1]

    def get_sirna_from_protein(self, protein_name):
        mrna_seq, accession, desc = self.fetch_mrna_from_protein(protein_name)
        if not mrna_seq:
            return {"success": False, "error": f"Could not find mRNA for: {protein_name}"}

        sirnas = self.design_sirna(mrna_seq)

        return {
            "success": True,
            "protein_name": protein_name,
            "refseq_id": accession,
            "description": desc,
            "siRNA_candidates": sirnas,
            "num_candidates": len(sirnas)
        }