#!/usr/bin/env python3
"""
Sum multiple .xray files and save the result as a PNG image.

This script reads all .xray files from a specified directory, sums their
intensities pixel-by-pixel, and saves the result as a PNG image.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict


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


def sum_xray_files(folder_path):
    """
    Sum all .xray files in a folder.

    Args:
        folder_path: Path to folder containing .xray files

    Returns:
        tuple: (x_pixels, y_pixels, summed_intensities, num_files)
            - x_pixels: numpy array of X pixel coordinates
            - y_pixels: numpy array of Y pixel coordinates
            - summed_intensities: numpy array of summed intensity values
            - num_files: number of files processed
    """
    folder = Path(folder_path)

    if not folder.exists():
        raise ValueError(f"Folder does not exist: {folder_path}")

    # Find all .xray files
    xray_files = sorted(folder.glob("*.xray"))

    if len(xray_files) == 0:
        raise ValueError(f"No .xray files found in {folder_path}")

    print(f"Found {len(xray_files)} .xray files")

    # Dictionary to accumulate intensities by (x, y) pixel coordinates
    intensity_sum = defaultdict(float)

    # Process each file
    for i, xray_file in enumerate(xray_files, 1):
        print(f"Processing {xray_file.name} ({i}/{len(xray_files)})...")
        x_pixels, y_pixels, intensities = read_xray_file(xray_file)

        # Sum intensities for each pixel
        for x, y, intensity in zip(x_pixels, y_pixels, intensities):
            intensity_sum[(x, y)] += intensity

    # Convert to arrays
    x_pixels = []
    y_pixels = []
    summed_intensities = []

    for (x, y), intensity in intensity_sum.items():
        x_pixels.append(x)
        y_pixels.append(y)
        summed_intensities.append(intensity)

    x_pixels = np.array(x_pixels)
    y_pixels = np.array(y_pixels)
    summed_intensities = np.array(summed_intensities)

    print(f"Total pixels: {len(x_pixels)}")
    print(
        f"Intensity range: {summed_intensities.min():.2e} to {summed_intensities.max():.2e}"
    )

    return x_pixels, y_pixels, summed_intensities, len(xray_files)


def plot_intensity_map(
    x_pixels, y_pixels, intensities, num_files, output_path, title=None, log_scale=False
):
    """
    Create and save a 2D intensity map plot.

    Args:
        x_pixels: numpy array of X pixel coordinates
        y_pixels: numpy array of Y pixel coordinates
        intensities: numpy array of intensity values
        num_files: number of files that were summed (for title)
        output_path: Path to save the PNG file
        title: Optional custom title (default: uses num_files)
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
    plt.colorbar(label="Summed Intensity")
    plt.xlabel("X Pixel")
    plt.ylabel("Y Pixel")

    if title is None:
        plt.title(f"Summed X-ray Diffraction Intensity ({num_files} files)")
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
    print(f"Saved summed image to: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Sum multiple .xray files and save as PNG image"
    )
    parser.add_argument(
        "folder", type=str, help="Path to folder containing .xray files"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output PNG file path (default: folder/summed_xray.png)",
    )

    args = parser.parse_args()

    try:
        # Sum the files
        x_pixels, y_pixels, summed_intensities, num_files = sum_xray_files(args.folder)

        # Determine output path
        if args.output is None:
            folder = Path(args.folder)
            output_path = folder / "summed_xray.png"
        else:
            output_path = args.output

        # Plot and save
        plot_intensity_map(
            x_pixels, y_pixels, summed_intensities, num_files, output_path
        )
        # Plot and save
        plot_intensity_map(
            x_pixels,
            y_pixels,
            summed_intensities,
            num_files,
            output_path,
            log_scale=True,
        )
    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
