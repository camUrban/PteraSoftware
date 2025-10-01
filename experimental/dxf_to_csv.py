"""DXF to CSV converter.

This script extracts outline coordinates from a DXF file and saves them to a CSV
file. It then plots the outline using matplotlib.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate as sp_interp


def _parse_dxf(filepath):
    """Parse a DXF file and extract geometric entities.

    :param filepath: Path to the DXF file.

    :return: List of entities, where each entity is a dict with 'type' and 'coords'.
    """
    these_entities = []

    with open(filepath, "r") as f:
        lines = [line.strip() for line in f.readlines()]

    i = 0
    while i < len(lines):
        # Look for entity markers.
        if lines[i] == "0":
            i += 1
            if i >= len(lines):
                break
            entity_type = lines[i]

            if entity_type == "LWPOLYLINE":
                this_entity = _parse_lwpolyline(lines, i)
                if this_entity:
                    these_entities.append(this_entity)
                    i = (
                        this_entity["end_index"] - 1
                    )  # Subtract 1 because we'll increment at the end.
            elif entity_type == "LINE":
                this_entity = _parse_line(lines, i)
                if this_entity:
                    these_entities.append(this_entity)
                    i = (
                        this_entity["end_index"] - 1
                    )  # Subtract 1 because we'll increment at the end.
        i += 1

    return these_entities


def _parse_lwpolyline(lines, start_index):
    """Parse a LWPOLYLINE entity from DXF file lines.

    :param lines: List of all lines in the DXF file.

    :param start_index: Starting index for this entity.

    :return: Dictionary with entity type and coordinates.
    """
    these_coords = []
    i = start_index + 1  # Skip the entity type line.
    current_x = None

    while i < len(lines) - 1:
        code = lines[i]
        value = lines[i + 1]

        # Check for end of entity (group code 0).
        if code == "0":
            break

        # Code 10 is x-coordinate.
        if code == "10":
            current_x = float(value)
        # Code 20 is y-coordinate.
        elif code == "20" and current_x is not None:
            these_coords.append([current_x, float(value)])
            current_x = None

        i += 2  # Move to next code-value pair.

    if these_coords:
        return {"type": "LWPOLYLINE", "coords": these_coords, "end_index": i}
    return None


def _parse_line(lines, start_index):
    """Parse a LINE entity from DXF file lines.

    :param lines: List of all lines in the DXF file.
    :param start_index: Starting index for this entity.
    :return: Dictionary with entity type and coordinates.
    """
    i = start_index + 1  # Skip the entity type line.
    x1 = y1 = x2 = y2 = None

    while i < len(lines) - 1:
        code = lines[i]
        value = lines[i + 1]

        # Check for end of entity (group code 0).
        if code == "0":
            break

        # Parse coordinates.
        if code == "10":
            x1 = float(value)
        elif code == "20":
            y1 = float(value)
        elif code == "11":
            x2 = float(value)
        elif code == "21":
            y2 = float(value)

        i += 2  # Move to next code-value pair.

    if x1 is not None and y1 is not None and x2 is not None and y2 is not None:
        these_coords = [[x1, y1], [x2, y2]]
        return {"type": "LINE", "coords": these_coords, "end_index": i}
    return None


def _extract_all_coordinates(these_entities):
    """Extract all coordinates from parsed entities in order.

    :param these_entities: List of entity dictionaries.
    :return: List of [x, y] coordinate pairs.
    """
    all_coords = []

    for this_entity in these_entities:
        all_coords.extend(this_entity["coords"])

    return all_coords


def _transform_coordinates(these_coords):
    """Transform coordinates by swapping X and Y, then negating new X.

    :param these_coords: List of [x, y] coordinate pairs.

    :return: List of [new_x, new_y, new_z] coordinate triples where new_x = -old_y,
    new_y = old_x, and new_z = 0.
    """
    these_transformed_coords = []

    for x, y in these_coords:
        new_x = -y
        new_y = x
        new_z = 0.0
        these_transformed_coords.append([new_x, new_y, new_z])

    return these_transformed_coords


def _classify_edges(these_entities):
    """Classify entities into wing edges based on their geometric properties.

    :param these_entities: List of entity dictionaries with transformed coordinates.

    :return: Dictionary with keys 'root', 'leading', 'tip', 'trailing',
    each containing list of [x, y, z] coordinates.
    """
    these_edges = {"root": [], "leading": [], "tip": [], "trailing": []}

    for this_entity in these_entities:
        these_coords = this_entity["coords"]
        if len(these_coords) < 2:
            continue

        this_coords_array = np.array(these_coords, dtype=float)
        x_vals = this_coords_array[:, 0]
        y_vals = this_coords_array[:, 1]

        this_x_min, this_x_max = float(np.min(x_vals)), float(np.max(x_vals))
        this_y_min, this_y_max = float(np.min(y_vals)), float(np.max(y_vals))
        this_x_range = this_x_max - this_x_min
        this_y_range = this_y_max - this_y_min

        # Classify based on primary direction of variation.
        if this_y_range > this_x_range:
            # Varies primarily in Y (spanwise), so it's a chordwise edge.
            if np.mean(x_vals) < 0:
                these_edges["leading"].extend(these_coords)
            else:
                these_edges["trailing"].extend(these_coords)
        else:
            # Varies primarily in X (chordwise), so it's a spanwise edge.
            if np.mean(y_vals) < 0.005:
                these_edges["root"].extend(these_coords)
            else:
                these_edges["tip"].extend(these_coords)

    return these_edges


def _sort_edge_root_to_tip(these_edge_coords):
    """Sort edge coordinates from root to tip (ascending Y).

    :param these_edge_coords: List of [x, y, z] coordinates.

    :return: Sorted list of [x, y, z] coordinates.
    """
    if len(these_edge_coords) == 0:
        return these_edge_coords

    this_coords_array = np.array(these_edge_coords, dtype=float)
    these_sorted_indices = np.argsort(this_coords_array[:, 1])
    these_sorted_coords = this_coords_array[these_sorted_indices].tolist()

    return these_sorted_coords


def _resample_edge(these_edge_coords, y_values):
    """Resample edge coordinates at specified Y values using PCHIP interpolation.

    :param these_edge_coords: List of [x, y, z] coordinates sorted by Y.

    :param y_values: (N,) ndarray of floats, Y values at which to sample.

    :return: List of [x, y, z] resampled coordinates.
    """
    this_coords_array = np.array(these_edge_coords, dtype=float)
    y_original = this_coords_array[:, 1]
    x_original = this_coords_array[:, 0]

    # Create PCHIP interpolator for x as a function of y.
    x_interpolator = sp_interp.PchipInterpolator(y_original, x_original)

    # Evaluate at the new y values.
    x_resampled = x_interpolator(y_values)

    # Z values are all zero (flat planform).
    z_resampled = np.zeros_like(y_values, dtype=float)

    # Combine into [x, y, z] coordinates.
    resampled_coords = np.column_stack([x_resampled, y_values, z_resampled])

    return resampled_coords.tolist()


def _save_to_csv(these_coords, output_filepath):
    """Save coordinates to a CSV file.

    :param these_coords: List of [x, y, z] coordinate triples.

    :param output_filepath: Path to output CSV file.
    """
    with open(output_filepath, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["X", "Y", "Z"])
        writer.writerows(these_coords)

    print(f"Saved {len(these_coords)} coordinates to {output_filepath}")


def process_dxf_to_wing_section_data(this_dxf_filepath, this_num_sections):
    """Process a DXF file to extract wing section data.

    :param this_dxf_filepath: str, path to the DXF file.

    :param this_num_sections: int, number of spanwise sections.

    :return: (N, 4) ndarray of floats, where N = num_sections + 1. Each row contains
    [displacement_x, displacement_y, displacement_z, chord] for a wing section.
    """
    # Parse the DXF file.
    these_entities = _parse_dxf(this_dxf_filepath)

    # Transform entity coordinates for classification.
    these_transformed_entities = []
    for this_entity in these_entities:
        these_transformed_entity_coords = _transform_coordinates(this_entity["coords"])
        these_transformed_entities.append(
            {
                "type": this_entity["type"],
                "coords": these_transformed_entity_coords,
                "end_index": this_entity["end_index"],
            }
        )

    # Classify edges.
    these_edges = _classify_edges(these_transformed_entities)

    # Sort leading and trailing edges from root to tip.
    these_edges["leading"] = _sort_edge_root_to_tip(these_edges["leading"])
    these_edges["trailing"] = _sort_edge_root_to_tip(these_edges["trailing"])

    # Get max y value.
    these_leading_coords = np.array(these_edges["leading"], dtype=float)
    these_trailing_coords = np.array(these_edges["trailing"], dtype=float)
    this_y_max = max(these_leading_coords[-1, 1], these_trailing_coords[-1, 1])

    # Resample leading and trailing edges.
    this_y_uniform = np.linspace(0, this_y_max, this_num_sections + 1, dtype=float)
    these_edges["leading"] = _resample_edge(these_edges["leading"], this_y_uniform)
    these_edges["trailing"] = _resample_edge(these_edges["trailing"], this_y_uniform)

    # Calculate chord lengths.
    this_leading_resampled = np.array(these_edges["leading"], dtype=float)
    this_trailing_resampled = np.array(these_edges["trailing"], dtype=float)
    these_chord_lengths = this_trailing_resampled[:, 0] - this_leading_resampled[:, 0]

    # Calculate displacements.
    result = np.zeros((this_num_sections + 1, 4), dtype=float)

    for i in range(len(this_leading_resampled)):
        if i == 0:
            # First point: displacement from origin (0, 0, 0).
            displacement_x = this_leading_resampled[i, 0] - 0.0
            displacement_y = this_leading_resampled[i, 1] - 0.0
            displacement_z = this_leading_resampled[i, 2] - 0.0
        else:
            # Subsequent points: displacement from previous point.
            displacement_x = (
                this_leading_resampled[i, 0] - this_leading_resampled[i - 1, 0]
            )
            displacement_y = (
                this_leading_resampled[i, 1] - this_leading_resampled[i - 1, 1]
            )
            displacement_z = (
                this_leading_resampled[i, 2] - this_leading_resampled[i - 1, 2]
            )

        result[i, 0] = displacement_x
        result[i, 1] = displacement_y
        result[i, 2] = displacement_z
        result[i, 3] = these_chord_lengths[i]

    return result


def save_leading_edge_with_chords(
    these_leading_coords, these_chord_lengths, output_filepath
):
    """Save leading edge displacements with chord lengths to a CSV file.

    :param these_leading_coords: List of [x, y, z] coordinate triples for leading edge.

    :param these_chord_lengths: (N,) ndarray of floats, chord lengths at each point.

    :param output_filepath: Path to output CSV file.
    """
    this_coords_array = np.array(these_leading_coords, dtype=float)

    with open(output_filepath, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Displacement_X", "Displacement_Y", "Displacement_Z", "Chord"])

        for i in range(len(this_coords_array)):
            if i == 0:
                # First point: displacement from origin (0, 0, 0).
                displacement_x = this_coords_array[i, 0] - 0.0
                displacement_y = this_coords_array[i, 1] - 0.0
                displacement_z = this_coords_array[i, 2] - 0.0
            else:
                # Subsequent points: displacement from previous point.
                displacement_x = this_coords_array[i, 0] - this_coords_array[i - 1, 0]
                displacement_y = this_coords_array[i, 1] - this_coords_array[i - 1, 1]
                displacement_z = this_coords_array[i, 2] - this_coords_array[i - 1, 2]

            row = [
                displacement_x,
                displacement_y,
                displacement_z,
                these_chord_lengths[i],
            ]
            writer.writerow(row)

    n_disps = len(these_leading_coords)

    print(
        f"Saved {n_disps} leading edge displacements with chords to {output_filepath}"
    )


def plot_outline(these_coords):
    """Plot the outline using matplotlib.

    Plots with +y pointing right (horizontal axis) and +x pointing down (vertical axis).

    :param these_coords: List of [x, y] coordinate pairs.
    """
    this_coords_array = np.array(these_coords, dtype=float)

    plt.figure(figsize=(10, 8))
    plt.plot(
        this_coords_array[:, 1],
        this_coords_array[:, 0],
        "b-",
        linewidth=1.5,
        label="Outline",
    )
    plt.plot(
        this_coords_array[:, 1],
        this_coords_array[:, 0],
        "ro",
        markersize=2,
        label="Points",
    )
    plt.xlabel("Y")
    plt.ylabel("X", labelpad=0)
    plt.title("DXF Outline")
    plt.gca().invert_yaxis()
    plt.grid(True, alpha=0.3)
    plt.axis("equal")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_edges(edges_dict):
    """Plot wing edges in different colors.

    Plots with +y pointing right (horizontal axis) and +x pointing down (vertical axis).

    :param edges_dict: Dictionary with keys 'root', 'leading', 'tip', 'trailing',
    each containing coordinates.
    """
    colors = {
        "root": "red",
        "leading": "blue",
        "tip": "green",
        "trailing": "orange",
    }

    plt.figure(figsize=(10, 8))

    for this_edge_name, these_coords in edges_dict.items():
        if len(these_coords) == 0:
            continue

        this_coords_array = np.array(these_coords, dtype=float)
        plt.plot(
            this_coords_array[:, 1],
            this_coords_array[:, 0],
            "-",
            color=colors[this_edge_name],
            linewidth=2,
            label=f"{this_edge_name.capitalize()} edge",
        )
        plt.plot(
            this_coords_array[:, 1],
            this_coords_array[:, 0],
            "o",
            color=colors[this_edge_name],
            markersize=2,
        )

    plt.xlabel("Y")
    plt.ylabel("X", labelpad=0)
    plt.title("Wing Planform Edges")
    plt.gca().invert_yaxis()
    plt.grid(True, alpha=0.3)
    plt.axis("equal")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    dxf_filepath = "gammabot_approximate_wing.dxf"
    csv_filepath = "gammabot_approximate_wing_outline.csv"
    num_sections = 20

    # Parse the DXF file.
    print(f"Parsing {dxf_filepath}...")
    entities = _parse_dxf(dxf_filepath)
    print(f"Found {len(entities)} entities")

    # Extract the coordinates.
    coords = _extract_all_coordinates(entities)
    print(f"Extracted {len(coords)} coordinate points")

    # Transform the coordinates.
    transformed_coords = _transform_coordinates(coords)
    print(f"Transformed coordinates")

    # Transform entity coordinates for classification.
    transformed_entities = []
    for entity in entities:
        transformed_entity_coords = _transform_coordinates(entity["coords"])
        transformed_entities.append(
            {
                "type": entity["type"],
                "coords": transformed_entity_coords,
                "end_index": entity["end_index"],
            }
        )

    # Classify edges.
    edges = _classify_edges(transformed_entities)
    print(f"Classified edges:")
    for edge_name, edge_coords in edges.items():
        print(f"  {edge_name}: {len(edge_coords)} points")

    # Sort leading and trailing edges from root to tip.
    edges["leading"] = _sort_edge_root_to_tip(edges["leading"])
    edges["trailing"] = _sort_edge_root_to_tip(edges["trailing"])
    print(f"Sorted leading and trailing edges from root to tip")

    # Validate edge properties.
    leading_coords = np.array(edges["leading"], dtype=float)
    trailing_coords = np.array(edges["trailing"], dtype=float)

    print(f"\nValidation checks:")
    print(f"  First leading edge point: {leading_coords[0]}")
    print(f"  Is at origin (0, 0, 0): {np.allclose(leading_coords[0], [0, 0, 0])}")

    le_y_start = leading_coords[0, 1]
    le_y_end = leading_coords[-1, 1]
    te_y_start = trailing_coords[0, 1]
    te_y_end = trailing_coords[-1, 1]

    print(f"  Leading edge Y range: {le_y_start:.6f} to {le_y_end:.6f}")
    print(f"  Trailing edge Y range: {te_y_start:.6f} to {te_y_end:.6f}")
    print(f"  Initial Y coordinates match: {np.isclose(le_y_start, te_y_start)}")
    print(f"  Final Y coordinates match: {np.isclose(le_y_end, te_y_end)}")

    # Resample leading and trailing edges.
    y_max = max(le_y_end, te_y_end)
    y_uniform = np.linspace(0, y_max, num_sections + 1, dtype=float)

    edges["leading"] = _resample_edge(edges["leading"], y_uniform)
    edges["trailing"] = _resample_edge(edges["trailing"], y_uniform)
    print(f"\nResampled leading and trailing edges to {num_sections + 1} points")
    print(f"  Y values: {y_uniform[0]:.6f} to {y_uniform[-1]:.6f}")

    # Validate that leading edge is forward of trailing edge.
    leading_resampled_check = np.array(edges["leading"], dtype=float)
    trailing_resampled_check = np.array(edges["trailing"], dtype=float)
    x_leading = leading_resampled_check[:, 0]
    x_trailing = trailing_resampled_check[:, 0]

    leading_ahead = np.all(x_leading < x_trailing)
    print(f"  Leading edge X < trailing edge X for all points: {leading_ahead}")
    if not leading_ahead:
        print(
            f"    WARNING: Leading edge is not ahead of trailing edge at some points!"
        )

    # Calculate chord length at each section.
    chord_lengths = x_trailing - x_leading
    print(f"\nChord lengths at each section:")
    print(f"  Root chord: {chord_lengths[0]:.6f}")
    print(f"  Tip chord: {chord_lengths[-1]:.6f}")
    print(f"  Min chord: {np.min(chord_lengths):.6f}")
    print(f"  Max chord: {np.max(chord_lengths):.6f}")

    # Save leading edge points with chord lengths to CSV.
    save_leading_edge_with_chords(edges["leading"], chord_lengths, csv_filepath)

    # Create four subplots for visualization.
    plt.rcParams.update({"font.size": 7})
    fig, axes = plt.subplots(1, 4, figsize=(32, 8))

    # Subplot 1: Original outline.
    ax1 = axes[0]
    coords_array = np.array(transformed_coords, dtype=float)
    ax1.plot(coords_array[:, 1], coords_array[:, 0], "b-", linewidth=1, label="Outline")
    ax1.plot(coords_array[:, 1], coords_array[:, 0], "ro", markersize=2, label="Points")
    ax1.set_xlabel("Y")
    ax1.set_ylabel("X", labelpad=0)
    ax1.set_title("Original Outline")
    ax1.invert_yaxis()
    ax1.grid(True, alpha=0.3)
    ax1.axis("equal")
    ax1.legend()
    ax1.tick_params(axis="y", pad=1)

    # Subplot 2: Classified edges (before resampling).
    ax2 = axes[1]
    edge_colors = {
        "root": "g",
        "leading": "c",
        "tip": "m",
        "trailing": "y",
    }
    # Use the sorted but not resampled edges.
    edges_before_resample = {
        "root": edges["root"],
        "leading": _sort_edge_root_to_tip(
            _classify_edges(transformed_entities)["leading"]
        ),
        "tip": edges["tip"],
        "trailing": _sort_edge_root_to_tip(
            _classify_edges(transformed_entities)["trailing"]
        ),
    }
    for edge_name, edge_coords in edges_before_resample.items():
        if len(edge_coords) == 0:
            continue
        edge_array = np.array(edge_coords, dtype=float)
        ax2.plot(
            edge_array[:, 1],
            edge_array[:, 0],
            "-",
            color=edge_colors[edge_name],
            linewidth=1,
            label=f"{edge_name.capitalize()} edge",
        )
        ax2.plot(
            edge_array[:, 1],
            edge_array[:, 0],
            "o",
            color=edge_colors[edge_name],
            markersize=2,
        )
    ax2.set_xlabel("Y")
    ax2.set_ylabel("X", labelpad=0)
    ax2.set_title("Classified Edges")
    ax2.invert_yaxis()
    ax2.grid(True, alpha=0.3)
    ax2.axis("equal")
    ax2.legend()
    ax2.tick_params(axis="y", pad=1)

    # Subplot 3: Resampled leading and trailing edges.
    ax3 = axes[2]
    leading_resampled = np.array(edges["leading"], dtype=float)
    trailing_resampled = np.array(edges["trailing"], dtype=float)
    ax3.plot(
        leading_resampled[:, 1],
        leading_resampled[:, 0],
        "c-",
        linewidth=1,
        label="Leading edge",
    )
    ax3.plot(
        leading_resampled[:, 1],
        leading_resampled[:, 0],
        "co",
        markersize=2,
    )
    ax3.plot(
        trailing_resampled[:, 1],
        trailing_resampled[:, 0],
        "y",
        linewidth=1,
        label="Trailing edge",
    )
    ax3.plot(
        trailing_resampled[:, 1],
        trailing_resampled[:, 0],
        "yo",
        markersize=2,
    )
    ax3.set_xlabel("Y")
    ax3.set_ylabel("X", labelpad=0)
    ax3.set_title(f"Resampled Edges ({num_sections + 1} points)")
    ax3.invert_yaxis()
    ax3.grid(True, alpha=0.3)
    ax3.axis("equal")
    ax3.legend()
    ax3.tick_params(axis="y", pad=1)

    # Subplot 4: Chord length distribution.
    ax4 = axes[3]
    ax4.plot(y_uniform, chord_lengths, "k-", linewidth=1, label="Chord")
    ax4.plot(y_uniform, chord_lengths, "ko", markersize=2)
    ax4.set_xlabel("Y (Spanwise Position)")
    ax4.set_ylabel("Chord Length", labelpad=0)
    ax4.set_title("Chord Distribution")
    ax4.grid(True, alpha=0.3)
    ax4.axis("equal")
    ax4.legend()
    ax4.tick_params(axis="y", pad=1)

    plt.subplots_adjust(
        left=0.05, right=0.98, top=0.92, bottom=0.15, wspace=0.3, hspace=0.3
    )
    plt.show()
