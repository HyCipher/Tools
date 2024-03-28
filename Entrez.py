from Bio import Entrez
from Bio import SeqIO

# Set your email address
Entrez.email = 'your_email@example.com'

def get_gene_id_from_ncbi(protein_accession):
    try:
        # Use Entrez's efetch to get GenBank record
        handle = Entrez.efetch(db="protein", id=protein_accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")

        # Extract gene ID from GenBank record
        # Note: Not all GenBank records contain gene IDs, you may need to check if the field is present in record
        gene_id = record.annotations.get('gene')

        return gene_id
    except Exception as e:
        print(f"Error retrieving gene ID for NCBI protein accession {protein_accession}: {e}")
        return None
    finally:
        handle.close()

# Example usage
protein_accession = 'NP_724342.1'  # Replace with the protein accession number you want to query
gene_id = get_gene_id_from_ncbi(protein_accession)
print(f"Gene ID for NCBI protein accession {protein_accession}: {gene_id}")
