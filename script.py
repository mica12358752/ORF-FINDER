from Bio import SeqIO

def find_orfs(seq, frame):
    """
    Identifica Open Reading Frames (ORFs) en una secuencia de ADN.
    Retorna una lista de tuplas con la posición de inicio y la secuencia del ORF.
    """
    seq = str(seq)
    orfs = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    for i in range(frame - 1, len(seq), 3):
        if seq[i:i+3] == start_codon:
            for j in range(i, len(seq), 3):
                codon = seq[j:j+3]
                if codon in stop_codons:
                    orf = seq[i:j+3]
                    orfs.append((i + 1, orf)) # Posición basada en 1
                    break
    return orfs

def main():
    file_path = "dna2.fasta"
    records = list(SeqIO.parse(file_path, "fasta"))
    print(f"Total de registros procesados: {len(records)}")

    # Análisis de ORFs en diferentes marcos de lectura
    longest_f2 = 0
    longest_f3 = ("", 0)

    for record in records:
        # Análisis en Frame 2
        orfs2 = find_orfs(record.seq, 2)
        for pos, orf in orfs2:
            if len(orf) > longest_f2:
                longest_f2 = len(orf)

        # Análisis en Frame 3
        orfs3 = find_orfs(record.seq, 3)
        for pos, orf in orfs3:
            if len(orf) > len(longest_f3[0]):
                longest_f3 = (orf, pos)

    print(f"Longitud máxima en Frame 2: {longest_f2}")
    print(f"Posición de inicio del ORF más largo en Frame 3: {longest_f3[1]}")

if __name__ == "__main__":
    main()
