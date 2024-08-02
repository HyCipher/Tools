from Bio import Entrez

# 设置NCBI的邮箱，用于身份验证
Entrez.email = "your.email@example.com"

# 基因ID
gene_id = "LOC108112063"

def fetch_gene_info(gene_id):
    """
    使用 Entrez 查找基因信息，并提取基因名称和物种信息。
    
    :param gene_id: 要查找的基因ID
    :return: None
    """
    # 使用 Entrez.esearch 搜索基因ID
    with Entrez.esearch(db="gene", term=gene_id, retmode="xml") as handle:
        record = Entrez.read(handle)
        
        # 获取搜索结果中的 ID 列表（通常只有一个 ID，但也可能有多个）
        gene_ids = record["IdList"]
        
        # 如果找到了基因ID，使用 Entrez.efetch 获取基因信息
        if gene_ids:
            gene_id = gene_ids[0]  # 取第一个ID，如果有多个ID，你可能需要遍历它们
            with Entrez.efetch(db="gene", id=gene_id, rettype="xml", retmode="xml") as handle:
                gene_record = Entrez.read(handle)
                
                # 打印返回的记录，查看实际结构
                print(gene_record)
                
                # 提取基因名称和物种名称
                try:
                    for gene in gene_record:
                        # 提取基因名称
                        gene_name = gene.get('Entrezgene_gene', {}).get('Gene-ref', {}).get('Gene-ref_locus', 'Unknown')
                        gene_desc = gene.get('Entrezgene_gene', {}).get('Gene-ref', {}).get('Gene-ref_desc', 'No Description')
                        
                        # 提取物种名称
                        organism = gene.get('Entrezgene_source', {}).get('BioSource', {}).get('BioSource_org', {}).get('Org-ref', {}).get('Org-ref_taxname', 'Unknown')
                        
                        print(f"Gene Name: {gene_name}")
                        print(f"Gene Description: {gene_desc}")
                        print(f"Organism: {organism}")
                except KeyError as e:
                    print(f"KeyError: {e}. Please check the structure of the gene record.")
        else:
            print(f"No gene found for ID {gene_id}")

# 执行函数
fetch_gene_info(gene_id)
