from Bio import Entrez

def get_gene_name(gene_id):
    Entrez.email = "kang.haoyu@.skku.edu"
    handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="text")
    record = handle.read().strip().split("\n")
    gene_name = None
    for line in record:
        if "Symbol" in line:
            gene_name = line.split("\t")[1]
            break
    handle.close()
    return gene_name

# 用法示例
gene_id = "LOC113564130"  # 你的基因ID
gene_name = get_gene_name(gene_id)
if gene_name:
    print(f"The gene name for ID {gene_id} is: {gene_name}")
else:
    print(f"No gene name found for ID {gene_id}")
