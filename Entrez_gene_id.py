from Bio import Entrez

# Verification with email
Entrez.email = "your.email@example.com"

# gene_ID
gene_id = "LOC6737448"

def fetch_gene_info(gene_id):
    """
    Use Entrez to look up gene information and extract gene names and species information.
    
    :param gene_id: The gene ID to look up
    :return: None
    """
    # Search gene_ID with Entrez.esearch
    with Entrez.esearch(db="gene", term=gene_id, retmode="xml") as handle:
        record = Entrez.read(handle)
        
        # Acquire the ID list from searching result (Usually only one ID, but may be more than one)
        gene_ids = record["IdList"]
        
        # If gene_ID was found, get gene information with Entrez.efetch.
        if gene_ids:
            gene_id = gene_ids[0]  # Take the first ID, if there are multiple IDs, you may need to iterate through them
            with Entrez.efetch(db="gene", id=gene_id, rettype="xml", retmode="xml") as handle:
                gene_record = Entrez.read(handle)
                
                # 打印返回的记录，查看实际结构
                # print(gene_record)
                
                # Extract gene name and species name
                try:
                    for gene in gene_record:
                        # gene name
                        gene_name = gene.get('Entrezgene_gene', {}).get('Gene-ref', {}).get('Gene-ref_locus', 'Unknown')
                        gene_desc = gene.get('Entrezgene_gene', {}).get('Gene-ref', {}).get('Gene-ref_desc', 'No Description')
                        
                        # species name
                        organism = gene.get('Entrezgene_source', {}).get('BioSource', {}).get('BioSource_org', {}).get('Org-ref', {}).get('Org-ref_taxname', 'Unknown')
                        
                        print(f"Gene Name: {gene_name}")
                        print(f"Gene Description: {gene_desc}")
                        print(f"Organism: {organism}")
                except KeyError as e:
                    print(f"KeyError: {e}. Please check the structure of the gene record.")
        else:
            print(f"No gene found for ID {gene_id}")

fetch_gene_info(gene_id)
