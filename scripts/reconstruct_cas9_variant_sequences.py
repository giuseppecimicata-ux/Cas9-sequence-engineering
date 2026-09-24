from pathlib import Path

import pandas as pd


# ============================================================================
# USER INPUT
# ============================================================================

FILE_PATH = "/path/to/variant_sequences.xlsx"
OUTPUT_PATH = "reconstructed_cas9_variants.xlsx"

# Column containing the internal protein variant sequences
SEQUENCE_COLUMN = 0

# Flanking residues required to reconstruct the complete engineered region
N_TERMINAL_FLANK = "ELDKAGFIKRQ"
C_TERMINAL_FLANK = "EVQTGGF"


# ============================================================================
# LOAD AND PROCESS SEQUENCES
# ============================================================================

def load_sequences(file_path, sequence_column=0):
    """
    Load amino-acid sequences from an Excel file.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the input Excel file.

    sequence_column : int
        Zero-based index of the column containing the amino-acid sequences.

    Returns
    -------
    list of str
        Input amino-acid sequences.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {file_path}"
        )

    df = pd.read_excel(
        file_path,
        header=None,
    )

    if df.shape[1] <= sequence_column:
        raise ValueError(
            f"The input file does not contain column "
            f"{sequence_column}."
        )

    sequences = []

    for value in df.iloc[:, sequence_column]:

        if pd.isna(value):
            continue

        sequence = str(value).strip()

        if not sequence:
            continue

        sequences.append(sequence)

    if not sequences:
        raise ValueError(
            "No amino-acid sequences were found in the input file."
        )

    return sequences


def reconstruct_sequences(
    sequences,
    n_terminal_flank,
    c_terminal_flank,
):
    """
    Add N- and C-terminal flanking residues to each variant sequence.

    Returns
    -------
    list of str
        Reconstructed protein sequences.
    """
    return [
        f"{n_terminal_flank}{sequence}{c_terminal_flank}"
        for sequence in sequences
    ]


def save_sequences(
    output_path,
    sequences,
):
    """
    Save reconstructed protein sequences to an Excel file.
    """
    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_df = pd.DataFrame(
        {"protein_sequence": sequences}
    )

    output_df.to_excel(
        output_path,
        index=False,
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    sequences = load_sequences(
        FILE_PATH,
        sequence_column=SEQUENCE_COLUMN,
    )

    reconstructed_sequences = reconstruct_sequences(
        sequences,
        n_terminal_flank=N_TERMINAL_FLANK,
        c_terminal_flank=C_TERMINAL_FLANK,
    )

    save_sequences(
        OUTPUT_PATH,
        reconstructed_sequences,
    )

    print(
        f"Reconstructed {len(reconstructed_sequences)} "
        f"protein variant sequences."
    )

    print(
        f"Output saved to: {OUTPUT_PATH}"
    )


if __name__ == "__main__":
    main()
