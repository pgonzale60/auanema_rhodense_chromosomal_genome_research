import csv
 
def generate_table_s5():
    """
    Generates Table S5: Genes present in eliminated DNA of Auanema rhodense.
    """
    input_file = "analyses/genes/PDE_characterization/all_eliminated_genes_annotated.tsv"
    
    print("# Table S5. Genes present in eliminated DNA of Auanema rhodense")
    print("\nThis table lists the 360 gene models identified within the germline-restricted (eliminated) regions. Functional annotations, expression levels (TPM), and association with transposable elements (TEs) are provided.\n")
    
    try:
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            headers = reader.fieldnames
            
            # Print Markdown Header
            print("| " + " | ".join(headers) + " |")
            print("| " + " | ".join(["---"] * len(headers)) + " |")
            
            for row in reader:
                print("| " + " | ".join(row[h] for h in headers) + " |")
                
    except FileNotFoundError:
        print(f"Error: {input_file} not found.")
 
if __name__ == "__main__":
    generate_table_s5()
