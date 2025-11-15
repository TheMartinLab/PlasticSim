#!/usr/bin/env python3
"""
Convert individual .xray files to PNG images.

This script reads all .xray files from a specified directory and converts
each one to a PNG image file.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def read_xray_file(file_path):
    """
    Read an .xray file and return pixel coordinates and intensities.

    Args:
        file_path: Path to the .xray file

    Returns:
        tuple: (x_pixels, y_pixels, intensities) as numpy arrays
    """
    data_lines = []
    skip_next_blank = False

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            # Skip comment lines
            if line.startswith("#"):
                continue
            # Skip the "xyI_column" line
            if line == "xyI_column":
                skip_next_blank = True
                continue
            # Skip the blank line after "xyI_column"
            if skip_next_blank and line == "":
                skip_next_blank = False
                continue
            # Skip blank lines
            if line == "":
                continue
            # Collect data lines
            data_lines.append(line)

    # Parse the tab-separated data
    x_pixels = []
    y_pixels = []
    intensities = []

    for line in data_lines:
        parts = line.split("\t")
        if len(parts) >= 3:
            x_pixels.append(float(parts[0]))
            y_pixels.append(float(parts[1]))
            intensities.append(float(parts[2]))

    return np.array(x_pixels), np.array(y_pixels), np.array(intensities)


def plot_intensity_map(
    x_pixels, y_pixels, intensities, output_path, title=None, log_scale=False
):
    """
    Create and save a 2D intensity map plot.

    Args:
        x_pixels: numpy array of X pixel coordinates
        y_pixels: numpy array of Y pixel coordinates
        intensities: numpy array of intensity values
        output_path: Path to save the PNG file
        title: Optional custom title
    """
    # Create 2D grid for visualization
    unique_x = np.unique(x_pixels)
    unique_y = np.unique(y_pixels)

    # Create mapping from pixel values to indices
    x_to_idx = {val: idx for idx, val in enumerate(unique_x)}
    y_to_idx = {val: idx for idx, val in enumerate(unique_y)}

    # Create intensity grid
    intensity_grid = np.zeros((len(unique_y), len(unique_x)))
    for x, y, intensity in zip(x_pixels, y_pixels, intensities):
        x_idx = x_to_idx[x]
        y_idx = y_to_idx[y]
        intensity_grid[y_idx, x_idx] = np.log(intensity) if log_scale else intensity

    # Create the plot
    plt.figure(figsize=(12, 10))
    plt.imshow(
        intensity_grid,
        extent=[x_pixels.min(), x_pixels.max(), y_pixels.min(), y_pixels.max()],
        origin="lower",
        aspect="auto",
        cmap="viridis",
        interpolation="nearest",
    )
    plt.colorbar(label="Intensity" + (" (log scale)" if log_scale else ""))
    plt.xlabel("X Pixel")
    plt.ylabel("Y Pixel")

    if title is None:
        plt.title("X-ray Diffraction Intensity")
    else:
        plt.title(title)

    plt.tight_layout()

    # Save the plot
    output_path = Path(output_path)
    if log_scale:
        output_path = output_path.with_suffix(".log.png")
    else:
        output_path = output_path.with_suffix(".png")
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {output_path}")
    plt.close()


def convert_xray_files(folder_path, output_folder=None):
    """
    Convert all .xray files in a folder to PNG images.

    Args:
        folder_path: Path to folder containing .xray files
        output_folder: Path to folder for output PNG files (default: same as input folder)
        log_scale: If True, also create log-scale versions
    """
    folder = Path(folder_path)

    if not folder.exists():
        raise ValueError(f"Folder does not exist: {folder_path}")

    # Find all .xray files
    xray_files = sorted(folder.glob("*.xray"))

    if len(xray_files) == 0:
        raise ValueError(f"No .xray files found in {folder_path}")

    print(f"Found {len(xray_files)} .xray files")

    # Determine output folder
    if output_folder is None:
        output_folder = folder
    else:
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

    # Process each file
    for i, xray_file in enumerate(xray_files, 1):
        print(f"Processing {xray_file.name} ({i}/{len(xray_files)})...")

        # Read the file
        x_pixels, y_pixels, intensities = read_xray_file(xray_file)

        # Determine output filename
        output_name = xray_file.stem + ".png"
        output_path = output_folder / output_name

        # Create title from filename
        title = f"X-ray Diffraction Intensity: {xray_file.name}"

        # Plot and save
        plot_intensity_map(
            x_pixels, y_pixels, intensities, output_path, title, log_scale=False
        )

        # Also create log-scale version if requested
        plot_intensity_map(
            x_pixels,
            y_pixels,
            intensities,
            output_path,
            title + " (log scale)",
            log_scale=True,
        )

    print(f"\nConverted {len(xray_files)} file(s) to PNG")


def main():
    parser = argparse.ArgumentParser(description="Convert .xray files to PNG images")
    parser.add_argument(
        "folder", type=str, help="Path to folder containing .xray files"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output folder for PNG files (default: same as input folder)",
    )
    # parser.add_argument(
    #     "--log-scale",
    #     action="store_true",
    #     help="Also create log-scale versions of the images",
    # )

    args = parser.parse_args()

    try:
        convert_xray_files(args.folder, args.output)
    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
