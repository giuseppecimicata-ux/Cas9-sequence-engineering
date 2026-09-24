import pandas as pd


# ============================================================================
# USER INPUT
# ============================================================================

INPUT_FILE = "/path/to/variant_sequences.xlsx"
OUTPUT_FILE = "cas9_variant_dna_sequences.xlsx"

# Column containing the amino-acid variant sequences.
SEQUENCE_COLUMN = "aa_seq_nuova"

# Reference protein and DNA sequences.
REFERENCE_AA_SEQUENCE = (
    "GALLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFH"
)

REFERENCE_DNA_SEQUENCE = (
    "ggagccctgctgttcgacagcggcgaaacagccgaggccgctcggctgtctagagaacgtagaagaagat"
    "acgaaagacggcgtaaccggttaaaatatctggaagagatcttcgctgatgagatggccaaggaagacgac"
    "agcttcttccac"
)

# Predefined codon assignment used for engineered residues.
CODON_TABLE = {
    "A": "GCT",
    "R": "CGT",
    "N": "AAT",
    "D": "GAT",
    "C": "TGT",
    "Q": "CAA",
    "E": "GAA",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "L": "TTA",
    "K": "AAA",
    "M": "ATG",
    "F": "TTT",
    "P": "CCT",
    "S": "TCT",
    "T": "ACT",
    "W": "TGG",
    "Y": "TAT",
    "V": "GTT",
    "*": "TAA",
}


# ============================================================================
# DNA DESIGN
# ============================================================================

def generate_modified_dna(
    reference_aa,
    reference_dna,
    variant_aa,
    codon_table,
):
    """
    Generate a DNA sequence corresponding to an amino-acid variant.

    Codons from the reference DNA sequence are preserved at unchanged
    amino-acid positions. Predefined codons are assigned to engineered
    residues.
    """
    if len(reference_aa) * 3 != len(reference_dna):
        raise ValueError(
            "Reference DNA length does not match the reference "
            "protein sequence."
        )

    if len(variant_aa) != len(reference_aa):
        raise ValueError(
            "Variant protein sequence must have the same length "
            "as the reference sequence."
        )

    modified_dna = []

    for position, (reference_residue, variant_residue) in enumerate(
        zip(reference_aa, variant_aa)
    ):

        if variant_residue not in codon_table:
            raise ValueError(
                f"Residue '{variant_residue}' is not present "
                "in the codon table."
            )

        if reference_residue == variant_residue:

            original_codon = reference_dna[
                position * 3 : (position + 1) * 3
            ]

            modified_dna.append(
                original_codon
            )

        else:

            modified_dna.append(
                codon_table[variant_residue]
            )

    return "".join(modified_dna)


# ============================================================================
# MAIN
# ============================================================================

def main():

    input_file = INPUT_FILE

    df = pd.read_excel(
        input_file
    )

    if SEQUENCE_COLUMN not in df.columns:
        raise ValueError(
            f"The input file must contain a column named "
            f"'{SEQUENCE_COLUMN}'."
        )

    variant_sequences = (
        df[SEQUENCE_COLUMN]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )

    if not variant_sequences:
        raise ValueError(
            "No amino-acid variant sequences were found."
        )

    dna_sequences = []

    for variant_sequence in variant_sequences:

        try:
            dna_sequence = generate_modified_dna(
                reference_aa=REFERENCE_AA_SEQUENCE,
                reference_dna=REFERENCE_DNA_SEQUENCE,
                variant_aa=variant_sequence,
                codon_table=CODON_TABLE,
            )

            dna_sequences.append(
                dna_sequence
            )

        except ValueError as error:

            dna_sequences.append(
                f"ERROR: {error}"
            )

    output_df = pd.DataFrame({
        "aa_sequence": variant_sequences,
        "dna_sequence": dna_sequences,
    })

    output_df.to_excel(
        OUTPUT_FILE,
        index=False,
    )

    print(
        f"Generated DNA sequences for "
        f"{len(variant_sequences)} protein variants."
    )

    print(
        f"Output saved to: {OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
