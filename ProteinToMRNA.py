import requests
import pandas as pd
import re
from io import StringIO  # Use io.StringIO instead of pd.compat.StringIO
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

class ProteinToMRNA:
    def __init__(self, email):
        """
        Initialize the ProteinToMRNA converter.

        Args:
            email (str): Email address for NCBI Entrez queries
        """
        self.email = email
        Entrez.email = email

        # Standard codon table
        self.codon_table = CodonTable.unambiguous_dna_by_name["Standard"]

        # Common codon optimization table for humans (frequency-based)
        # Values represent relative usage frequency in human cells
        self.human_codon_optimization = {
            'A': ['GCC', 'GCT', 'GCA', 'GCG'],  # Alanine (most to least frequent)
            'R': ['CGG', 'AGA', 'AGG', 'CGC', 'CGA', 'CGT'],  # Arginine
            'N': ['AAC', 'AAT'],  # Asparagine
            'D': ['GAC', 'GAT'],  # Aspartic acid
            'C': ['TGC', 'TGT'],  # Cysteine
            'Q': ['CAG', 'CAA'],  # Glutamine
            'E': ['GAG', 'GAA'],  # Glutamic acid
            'G': ['GGC', 'GGG', 'GGA', 'GGT'],  # Glycine
            'H': ['CAC', 'CAT'],  # Histidine
            'I': ['ATC', 'ATT', 'ATA'],  # Isoleucine
            'L': ['CTG', 'CTC', 'TTG', 'CTT', 'TTA', 'CTA'],  # Leucine
            'K': ['AAG', 'AAA'],  # Lysine
            'M': ['ATG'],  # Methionine (start)
            'F': ['TTC', 'TTT'],  # Phenylalanine
            'P': ['CCC', 'CCT', 'CCA', 'CCG'],  # Proline
            'S': ['AGC', 'TCC', 'TCT', 'AGT', 'TCA', 'TCG'],  # Serine
            'T': ['ACC', 'ACT', 'ACA', 'ACG'],  # Threonine
            'W': ['TGG'],  # Tryptophan
            'Y': ['TAC', 'TAT'],  # Tyrosine
            'V': ['GTG', 'GTC', 'GTT', 'GTA'],  # Valine
            '*': ['TGA', 'TAA', 'TAG']  # Stop codons
        }

    def fetch_protein_sequence(self, protein_name):
        """
        Fetch protein sequence from UniProt or NCBI based on protein name.

        Args:
            protein_name (str): Name of the protein

        Returns:
            tuple: (protein_sequence, protein_id, description)
        """
        # Try UniProt first
        try:
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/search?query={protein_name}+AND+reviewed:true&format=tsv&fields=id,protein_name,sequence,organism_name"
            response = requests.get(uniprot_url)

            if response.status_code == 200 and len(response.text) > 0:
                df = pd.read_csv(StringIO(response.text), sep='\t')

                if not df.empty:
                    # Filter for human proteins
                    human_df = df[df['Organism'] == 'Homo sapiens']
                    if not human_df.empty:
                        df = human_df

                    protein_id = df.iloc[0]['Entry']
                    description = f"{df.iloc[0]['Protein names']} [{df.iloc[0]['Organism']}]"
                    sequence = df.iloc[0]['Sequence']
                    return sequence, protein_id, description
        except Exception as e:
            print(f"UniProt search failed: {e}")

        # Fall back to NCBI Protein database
        try:
            search_term = f"{protein_name}[Protein Name] AND Homo sapiens[Organism] AND refseq[Filter]"
            search_handle = Entrez.esearch(db="protein", term=search_term, retmax=5)
            search_record = Entrez.read(search_handle)
            search_handle.close()

            if int(search_record["Count"]) > 0:
                protein_ids = search_record["IdList"]
                fetch_handle = Entrez.efetch(db="protein", id=protein_ids[0], rettype="fasta", retmode="text")
                record = SeqIO.read(fetch_handle, "fasta")
                fetch_handle.close()

                return str(record.seq), record.id, record.description
        except Exception as e:
            print(f"NCBI search failed: {e}")

        return None, None, None

    def back_translate(self, protein_sequence, optimize=True):
        """
        Perform back-translation of a protein sequence to mRNA.

        Args:
            protein_sequence (str): The amino acid sequence
            optimize (bool): Whether to optimize codons for human expression (frequency-based)

        Returns:
            str: The mRNA sequence (using T instead of U for DNA representation)
        """
        if not protein_sequence:
            return None

        mrna = ""

        # Add start codon
        mrna += "ATG"

        for aa in protein_sequence:
            if aa not in self.human_codon_optimization:
                continue  # Skip invalid amino acids

            if optimize:
                # Use the most frequent codon for human expression
                mrna += self.human_codon_optimization[aa][0]
            else:
                # Use the first valid codon from the standard codon table
                codons = [c for c in self.codon_table.forward_table.keys()
                          if self.codon_table.forward_table[c] == aa]
                if codons:
                    mrna += codons[0]

        # Add stop codon
        mrna += "TGA"

        return mrna

    def back_translate_structure_aware(self, protein_sequence, target_gc=0.50):
        """
        Perform back-translation while choosing codons that drive the global GC content
        toward a target value (here a proxy for structure optimization).

        Args:
            protein_sequence (str): The amino acid sequence
            target_gc (float): Desired GC ratio (default 0.50)

        Returns:
            str: The structure-aware optimized mRNA sequence (with T for U)
        """
        if not protein_sequence:
            return None

        mrna = "ATG"  # Start codon
        for aa in protein_sequence:
            if aa not in self.human_codon_optimization:
                continue
            candidates = self.human_codon_optimization[aa]
            best_codon = None
            best_diff = float('inf')
            # Choose the candidate that, when appended, brings the overall GC ratio closest to target_gc
            for codon in candidates:
                test_seq = mrna + codon
                gc_count = test_seq.count('G') + test_seq.count('C')
                gc_ratio = gc_count / len(test_seq)
                diff = abs(gc_ratio - target_gc)
                if diff < best_diff:
                    best_diff = diff
                    best_codon = codon
            if best_codon:
                mrna += best_codon
        mrna += "TGA"  # Stop codon
        return mrna

    def mrna_to_vaccine_format(self, mrna_sequence, pseudouridine=False):
        """
        Convert DNA sequence to proper mRNA format and add vaccine-specific elements.
        Optionally replace uridines with pseudouridine (Ψ).

        Args:
            mrna_sequence (str): DNA sequence (with T)
            pseudouridine (bool): If True, replace U with Ψ in the coding region

        Returns:
            str: Modified mRNA vaccine sequence with vaccine-specific elements
        """
        # Convert T to U for the coding region (or Ψ if pseudouridine is used)
        if pseudouridine:
            coding = mrna_sequence.replace('T', 'Ψ')
        else:
            coding = mrna_sequence.replace('T', 'U')

        # 5' cap representation
        cap = "m7G(5')ppp(5')G"
        # 5' UTR with Kozak sequence
        utr_5 = "GGGAAAUAAGAGAGAAAAGAAGAGUAAGAAGAAAUAUAAGACCCCGGCGCCGCCACC"
        # 3' UTR for stability
        utr_3 = "UGAUACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        # Poly(A) tail
        poly_a = "A" * 100

        complete_mrna = f"{cap}-{utr_5}-{coding}-{utr_3}-{poly_a}"
        return complete_mrna

    def get_optimized_mrna(self, protein_name, use_structure_aware=False, use_pseudouridine=False):
        """
        Complete pipeline to go from protein name to optimized mRNA sequence.

        Args:
            protein_name (str): Name of the protein
            use_structure_aware (bool): If True, also generate a structure-aware mRNA version
            use_pseudouridine (bool): If True, generate vaccine format mRNA with Ψ modifications

        Returns:
            dict: Results including mRNA sequences and metadata
        """
        # Step 1: Fetch protein sequence
        protein_seq, protein_id, description = self.fetch_protein_sequence(protein_name)

        if not protein_seq:
            return {
                "success": False,
                "error": f"Could not find protein: {protein_name}"
            }

        # Step 2: Back-translate to mRNA (basic, non-optimized)
        basic_mrna = self.back_translate(protein_seq, optimize=False)

        # Step 3: Back-translate to mRNA (optimized, frequency-based)
        optimized_mrna = self.back_translate(protein_seq, optimize=True)

        # Step 4: Structure-aware back-translation (if requested)
        structure_aware_mrna = None
        if use_structure_aware:
            structure_aware_mrna = self.back_translate_structure_aware(protein_seq)

        # Step 5: Convert optimized mRNA to vaccine format (with optional pseudouridine modification)
        vaccine_mrna = self.mrna_to_vaccine_format(optimized_mrna, pseudouridine=use_pseudouridine)

        # Step 6: Calculate GC content (important for stability)
        gc_content = ((optimized_mrna.count('G') + optimized_mrna.count('C')) /
                      len(optimized_mrna) * 100)

        result = {
            "success": True,
            "protein_name": protein_name,
            "protein_id": protein_id,
            "description": description,
            "protein_sequence": protein_seq,
            "protein_length": len(protein_seq),
            "basic_mrna": basic_mrna,
            "optimized_mrna": optimized_mrna,
            "vaccine_format_mrna": vaccine_mrna,
            "mrna_length": len(optimized_mrna),
            "gc_content": round(gc_content, 2)
        }

        if use_structure_aware and structure_aware_mrna:
            gc_content_sa = ((structure_aware_mrna.count('G') + structure_aware_mrna.count('C')) /
                             len(structure_aware_mrna) * 100)
            result["structure_aware_mrna"] = structure_aware_mrna
            result["structure_aware_gc_content"] = round(gc_content_sa, 2)

        return result