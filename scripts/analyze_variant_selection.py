from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap


# ============================================================================
# USER INPUT
# ============================================================================

FILE_PATH = "/path/to/your/variant_selection_matrix.xlsx"

OUTPUT_IDS_PATH = "selected_variant_ids.txt"

# Plot settings
COLORMAP = "binary"
ZERO_COLOR = "white"
VMIN = 0.0
VMAX = 1.0

FIGURE_SIZE = (6, 5)


# ============================================================================
# LOAD SELECTION MATRIX
# ============================================================================

def load_selection_matrix(file_path):
    """
    Load a variant-selection matrix from an Excel file.

    Columns correspond to design groups or regions.
    Rows correspond to individual variant indices.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {file_path}"
        )

    df = pd.read_excel(file_path)

    if df.empty:
        raise ValueError("The input Excel file is empty.")

    return df


# ============================================================================
# IDENTIFY SELECTED VARIANTS
# ============================================================================

def extract_selected_variant_ids(df):
    """
    Identify non-zero entries in the selection matrix.

    Each non-zero value is converted into a variant identifier using
    the format:

        <column>_<row index>

    Row indices are one-based.

    Returns
    -------
    selected_ids : dict
        Dictionary mapping each column/group to its selected variant IDs.

    total_selected : int
        Total number of non-zero entries.
    """
    selected_ids = {}
    total_selected = 0

    for column in df.columns:

        ids = []

        for row_index, value in enumerate(df[column], start=1):

            if pd.notna(value) and value != 0:
                ids.append(
                    f"{column}_{row_index}"
                )

        if ids:
            selected_ids[column] = ids
            total_selected += len(ids)

    return selected_ids, total_selected


# ============================================================================
# SAVE SELECTED VARIANT IDS
# ============================================================================

def save_selected_variant_ids(
    output_path,
    selected_ids,
    total_selected,
):
    """
    Save selected variant IDs to a text file.

    Each column/group is written as one line containing its selected IDs.
    """
    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(
        output_path,
        "w",
        encoding="utf-8",
    ) as output_file:

        for group, ids in selected_ids.items():
            output_file.write(
                f"{group}: {', '.join(ids)}\n"
            )

        output_file.write(
            f"\nTotal non-zero entries: {total_selected}\n"
        )

    return output_path


# ============================================================================
# CREATE SELECTION HEATMAP
# ============================================================================

def plot_selection_matrix(
    df,
    colormap=COLORMAP,
    zero_color=ZERO_COLOR,
    vmin=VMIN,
    vmax=VMAX,
    figure_size=FIGURE_SIZE,
):
    """
    Visualize the variant-selection matrix as a heatmap.

    Zero values are displayed as white.
    """
    data = df.to_numpy(dtype=float)

    base_cmap = plt.get_cmap(
        colormap,
        256,
    )

    colors = base_cmap(
        np.linspace(0, 1, 256)
    )

    colors[0, :] = np.array(
        [1, 1, 1, 1]
    )

    custom_cmap = ListedColormap(colors)

    fig, ax = plt.subplots(
        figsize=figure_size
    )

    image = ax.imshow(
        data,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
    )

    colorbar = fig.colorbar(
        image,
        ax=ax,
    )

    colorbar.set_label(
        "Selection score"
    )

    ax.set_xlabel(
        "Variant group"
    )

    ax.set_ylabel(
        "Variant index"
    )

    ax.set_xticks(
        np.arange(len(df.columns))
    )

    ax.set_xticklabels(
        df.columns,
        rotation=90,
    )

    ax.grid(False)

    fig.tight_layout()

    plt.show()


# ============================================================================
# MAIN
# ============================================================================

def main():

    file_path = Path(FILE_PATH)

    df = load_selection_matrix(
        file_path
    )

    selected_ids, total_selected = (
        extract_selected_variant_ids(df)
    )

    print(
        f"Total non-zero entries: {total_selected}"
    )

    for group, ids in selected_ids.items():
        print(
            f"{group}: {', '.join(ids)}"
        )

    output_path = save_selected_variant_ids(
        OUTPUT_IDS_PATH,
        selected_ids,
        total_selected,
    )

    print(
        f"\nSelected variant IDs saved to: "
        f"{output_path}"
    )

    plot_selection_matrix(
        df
    )


if __name__ == "__main__":
    main()
