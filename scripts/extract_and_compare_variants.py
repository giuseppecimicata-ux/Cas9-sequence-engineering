from pathlib import Path

# ============================================================================
# USER INPUT
# ============================================================================

IDS_FILE_PATH = "/path/to/selected_variant_ids.txt"
FASTA_FILE_PATH = "/path/to/input_sequences.fasta"
OUTPUT_FASTA_PATH = "selected_variants.fasta"


# ============================================================================
# INPUT PARSING
# ============================================================================

def read_selected_ids(file_path):
    """
    Read selected variant identifiers from a text file.

    Expected format:
        VARIANT_GROUP: VARIANT_GROUP_1, VARIANT_GROUP_4, ...

    The corresponding reference sequence (<VARIANT_GROUP>_0) is
    automatically included for comparison.

    Returns
    -------
    group_name : str
        Name of the variant group.

    selected_ids : set
        Set of sequence identifiers to extract.

    reference_id : str
        Identifier of the reference sequence.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Selected-ID file not found: {file_path}"
        )

    with open(
        file_path,
        "r",
        encoding="utf-8",
    ) as file:

        lines = [
            line.strip()
            for line in file
            if line.strip()
        ]

    if not lines:
        raise ValueError(
            "The selected-ID file is empty."
        )

    # The current input format contains one group definition.
    header, ids_string = lines[0].split(":", 1)

    group_name = header.strip()

    selected_ids = {
        identifier.strip()
        for identifier in ids_string.split(",")
        if identifier.strip()
    }

    reference_id = f"{group_name}_0"
    selected_ids.add(reference_id)

    return group_name, selected_ids, reference_id


def read_fasta(file_path):
    """
    Read sequences from a FASTA file.

    Returns
    -------
    dict
        Mapping from sequence identifiers to sequences.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"FASTA file not found: {file_path}"
        )

    sequences = {}

    current_id = None
    current_sequence = []

    with open(
        file_path,
        "r",
        encoding="utf-8",
    ) as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:
                    sequences[current_id] = "".join(
                        current_sequence
                    )

                current_id = line[1:].strip()
                current_sequence = []

            else:
                if current_id is None:
                    raise ValueError(
                        "Sequence data found before the first FASTA header."
                    )

                current_sequence.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(
            current_sequence
        )

    if not sequences:
        raise ValueError(
            "No sequences were found in the FASTA file."
        )

    return sequences


# ============================================================================
# SEQUENCE COMPARISON
# ============================================================================

def count_sequence_differences(
    sequence_1,
    sequence_2,
):
    """
    Count residue differences between two aligned sequences.

    Differences are counted position-by-position. The sequences are
    expected to have equal length.
    """
    if len(sequence_1) != len(sequence_2):
        raise ValueError(
            "Sequences must have the same length for positional comparison."
        )

    return sum(
        residue_1 != residue_2
        for residue_1, residue_2 in zip(
            sequence_1,
            sequence_2,
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    (
        group_name,
        selected_ids,
        reference_id,
    ) = read_selected_ids(
        IDS_FILE_PATH
    )

    sequences = read_fasta(
        FASTA_FILE_PATH
    )

    if reference_id not in sequences:
        raise ValueError(
            f"Reference sequence '{reference_id}' "
            "was not found in the FASTA file."
        )

    reference_sequence = sequences[
        reference_id
    ]

    output_path = Path(
        OUTPUT_FASTA_PATH
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    extracted_count = 0

    with open(
        output_path,
        "w",
        encoding="utf-8",
    ) as output_file:

        for sequence_id in sorted(
            selected_ids,
            key=lambda identifier: (
                int(identifier.rsplit("_", 1)[-1])
                if identifier.rsplit("_", 1)[-1].isdigit()
                else float("inf")
            ),
        ):

            if sequence_id not in sequences:
                print(
                    f"Warning: sequence ID not found in FASTA: "
                    f"{sequence_id}"
                )
                continue

            sequence = sequences[
                sequence_id
            ]

            differences = count_sequence_differences(
                sequence,
                reference_sequence,
            )

            output_file.write(
                f">{sequence_id} | differences from "
                f"{reference_id}: {differences}\n"
            )

            output_file.write(
                sequence + "\n"
            )

            extracted_count += 1

    print(
        f"Variant group: {group_name}"
    )

    print(
        f"Reference sequence: {reference_id}"
    )

    print(
        f"Sequences written: {extracted_count}"
    )

    print(
        f"Output FASTA: {output_path}"
    )


if __name__ == "__main__":
    main()
