import os
import sys
from Bio import Entrez, SeqIO

Entrez.email = "researcher@example.com"

# Dictionary of species and their GenBank accessions
# For C. auriculariae, we have the extracted sequence from nxCaeAuri1.1.primary.fa
species_info = {
    "Auanema_rhodense": {
        "18S": "EU196004",
        "28S": "EU195960"
    },
    "Oscheius_tipulae": {
        "18S": "EU196009",
        "28S": "EU195969"
    },
    "Mesorhabditis_belari": {
        "18S": "MH983024",
        "28S": "EF417149"
    },
    "Parascaris_univalens": {
        "18S": "U94378",
        "28S": "DQ145716"
    },
    "Ascaris_suum": {
        "18S": "PX659182",
        "28S": "DQ145715"
    },
    "Ascaris_lumbricoides": {
        "18S": "U94366",
        "28S": "AY210806"
    },
    "Toxocara_canis": {
        "18S": "JN256977",
        "28S": "FJ418783"
    },
    "Caenorhabditis_elegans": {
        "18S": "NR_000054",
        "28S": "MN519140"
    },
    "Caenorhabditis_monodelphis": {
        "18S": "JN636070",
        "28S": "JN636070"
    },
    "Caenorhabditis_parvicauda": {
        "18S": "MH800325",
        "28S": "MH800325"
    },
    "Caenorhabditis_briggsae": {
        "18S": "MN519141",
        "28S": "MN519141"
    },
    "Trichinella_spiralis": {
        "18S": "AY009111",
        "28S": "AY009111"
    }
}

os.makedirs("phylogenetic_analysis/data", exist_ok=True)

seqs_18S = {}
seqs_28S = {}

# Load C. auriculariae
rec18 = SeqIO.read("phylogenetic_analysis/c_auriculariae_18S.fa", "fasta")
seqs_18S["Caenorhabditis_auriculariae"] = str(rec18.seq)

rec28 = SeqIO.read("phylogenetic_analysis/c_auriculariae_28S.fa", "fasta")
seqs_28S["Caenorhabditis_auriculariae"] = str(rec28.seq)

print("Fetching GenBank sequences for remaining 10 species...")
for sp, accs in species_info.items():
    print(f"Fetching {sp}...")
    
    # 18S
    acc18 = accs["18S"]
    h18 = Entrez.efetch(db="nucleotide", id=acc18, rettype="gb", retmode="text")
    rec_gb_18 = SeqIO.read(h18, "genbank")
    
    # Check if full sequence or rRNA feature
    if acc18 in ["JN636070", "MH800325"]:
        # Contains multiple features
        extracted = False
        for f in rec_gb_18.features:
            if f.type == "rRNA" and ("18S" in f.qualifiers.get("product", [""])[0] or "small" in f.qualifiers.get("product", [""])[0]):
                seqs_18S[sp] = str(f.extract(rec_gb_18.seq))
                extracted = True
                break
        if not extracted:
            # If feature not explicitly tagged as 18S, grab 1..1737 for JN636070, or 1..1698 for MH800325
            if acc18 == "JN636070":
                seqs_18S[sp] = str(rec_gb_18.seq[:1737])
            elif acc18 == "MH800325":
                seqs_18S[sp] = str(rec_gb_18.seq[:1700])
    else:
        seqs_18S[sp] = str(rec_gb_18.seq)
        
    # 28S
    acc28 = accs["28S"]
    if acc28 == acc18:
        rec_gb_28 = rec_gb_18
    else:
        h28 = Entrez.efetch(db="nucleotide", id=acc28, rettype="gb", retmode="text")
        rec_gb_28 = SeqIO.read(h28, "genbank")
        
    if acc28 in ["MN519140", "JN636070", "MH800325", "AY210806"]:
        extracted = False
        for f in rec_gb_28.features:
            if f.type == "rRNA" and any(k in f.qualifiers.get("product", [""])[0] for k in ["28S", "26S", "large"]):
                seqs_28S[sp] = str(f.extract(rec_gb_28.seq))
                extracted = True
                break
        if not extracted:
            if acc28 == "JN636070":
                seqs_28S[sp] = str(rec_gb_28.seq[2673:6173])
            elif acc28 == "MH800325":
                seqs_28S[sp] = str(rec_gb_28.seq[2700:])
            elif acc28 == "AY210806":
                seqs_28S[sp] = str(rec_gb_28.seq[323:3979])
    else:
        seqs_28S[sp] = str(rec_gb_28.seq)

print("\n--- Sequence Summary ---")
for sp in sorted(seqs_18S.keys()):
    print(f"{sp:30s} | 18S: {len(seqs_18S[sp])} bp | 28S: {len(seqs_28S[sp])} bp")

with open("phylogenetic_analysis/data/all_18S_13taxa.fasta", "w") as f:
    for sp, s in seqs_18S.items():
        f.write(f">{sp}\n{s}\n")

with open("phylogenetic_analysis/data/all_28S_13taxa.fasta", "w") as f:
    for sp, s in seqs_28S.items():
        f.write(f">{sp}\n{s}\n")

print("\nSaved all_18S_13taxa.fasta and all_28S_13taxa.fasta successfully!")
