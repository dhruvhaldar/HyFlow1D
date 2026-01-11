import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
import sys
import argparse
import pathlib

def is_safe_path(path_str):
    """
    Security: Validate path to prevent directory traversal.
    Returns False if path contains '..' components.
    """
    if not path_str:
        return True
    try:
        path = pathlib.Path(path_str)
        if ".." in path.parts:
            return False
        return True
    except Exception:
        return False

def plot_all():
    parser = argparse.ArgumentParser(
        description="Visualize HyFlow1D simulation results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_dir",
        nargs="?",
        default=None,
        help="Directory containing solution_*.csv files. Defaults to 'output/' or current directory."
    )
    parser.add_argument(
        "-o", "--output",
        default="advection_plot.png",
        help="Output filename for the plot image."
    )

    args = parser.parse_args()

    # Security: Validate inputs
    if not is_safe_path(args.output):
        print(f"❌ Error: Invalid output path '{args.output}'. Path traversal ('..') is not allowed.", file=sys.stderr)
        sys.exit(1)

    if args.input_dir and not is_safe_path(args.input_dir):
        print(f"❌ Error: Invalid input directory '{args.input_dir}'. Path traversal ('..') is not allowed.", file=sys.stderr)
        sys.exit(1)

    # Determine search path
    if args.input_dir:
        search_path = os.path.join(args.input_dir, "solution_*.csv")
        display_dir = args.input_dir
    else:
        # Auto-detect: check 'output/' first, then current directory
        if glob.glob("output/solution_*.csv"):
            search_path = "output/solution_*.csv"
            display_dir = "output/"
        else:
            search_path = "solution_*.csv"
            display_dir = "./"

    files = sorted(glob.glob(search_path), key=lambda f: int(''.join(filter(str.isdigit, os.path.basename(f))) or 999999))
    
    if not files:
        print(f"❌ No solution files found in: {display_dir}")
        print("   (Looking for 'solution_*.csv')")
        return

    print(f"📊 Found {len(files)} solution files in '{display_dir}'")

    # Plot initial, middle, and final
    # Limit to a few frames
    if len(files) > 5:
        indices = [0, len(files)//4, len(files)//2, 3*len(files)//4, len(files)-1]
        files_to_plot = [files[i] for i in indices]
    else:
        files_to_plot = files

    plt.figure(figsize=(10, 6))
    
    for f in files_to_plot:
        try:
            data = pd.read_csv(f)
            # Extract step number for label
            step_num = ''.join(filter(str.isdigit, os.path.basename(f)))
            label = f"Step {step_num}" if step_num else os.path.basename(f)
            plt.plot(data['x'], data['u'], label=label)
        except Exception as e:
            print(f"⚠️  Warning: Could not read {f}: {e}")
        
    plt.axvline(x=0.5, color='k', linestyle='--', label='Interface (FV | DG)')
    plt.xlabel('Position (x)')
    plt.ylabel('Value (u)')
    plt.title('Hybrid FV-DG Linear Advection Simulation')
    plt.legend()
    plt.grid(True, alpha=0.3)

    try:
        plt.savefig(args.output, dpi=150)
        print(f"✅ Plot saved to: {os.path.abspath(args.output)}")
    except Exception as e:
        print(f"❌ Error saving plot: {e}")

if __name__ == "__main__":
    plot_all()
