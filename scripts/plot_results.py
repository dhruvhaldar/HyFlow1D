import sys
import glob
import os
import argparse
from pathlib import Path
from itertools import cycle

# Palette: Setup colors for better UX (respects NO_COLOR standard)
class Colors:
    _enabled = sys.stdout.isatty() and os.getenv("NO_COLOR") is None

    RESET   = "\033[0m" if _enabled else ""
    BOLD    = "\033[1m" if _enabled else ""
    RED     = "\033[31m" if _enabled else ""
    GREEN   = "\033[32m" if _enabled else ""
    YELLOW  = "\033[33m" if _enabled else ""
    BLUE    = "\033[34m" if _enabled else ""

    BOLD_RED     = "\033[1;31m" if _enabled else ""
    BOLD_GREEN   = "\033[1;32m" if _enabled else ""

try:
    import matplotlib.pyplot as plt
    import pandas as pd
except ImportError as e:
    print(f"\n{Colors.BOLD_RED}❌ Error: Missing required dependencies for visualization.{Colors.RESET}")
    print(f"   Reason: {e}")
    print(f"\n{Colors.YELLOW}💡 Please install them with:{Colors.RESET}")
    print(f"   pip install matplotlib pandas\n")
    sys.exit(1)

def is_safe_path(path_str):
    if not path_str:
        return False
    path = Path(path_str)
    for part in path.parts:
        if part == "..":
            return False
    return True

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

    # Security: Validate paths
    if args.input_dir and not is_safe_path(args.input_dir):
        print("❌ Error: Invalid input directory. Path traversal ('..') is not allowed.")
        sys.exit(1)

    if not is_safe_path(args.output):
        print("❌ Error: Invalid output filename. Path traversal ('..') is not allowed.")
        sys.exit(1)

    # Determine search path
    if args.input_dir:
        search_path = os.path.join(args.input_dir, "solution_*.csv")
        display_dir = args.input_dir
    else:
        # Auto-detect: check 'output/' first, then 'build/output/', then current directory
        if glob.glob("output/solution_*.csv"):
            search_path = "output/solution_*.csv"
            display_dir = "output/"
        elif glob.glob("build/output/solution_*.csv"):
            search_path = "build/output/solution_*.csv"
            display_dir = "build/output/"
        else:
            search_path = "solution_*.csv"
            display_dir = "./"

    files = sorted(glob.glob(search_path), key=lambda f: int(''.join(filter(str.isdigit, os.path.basename(f))) or 999999))
    
    if not files:
        print(f"{Colors.BOLD_RED}❌ No solution files found in: {display_dir}{Colors.RESET}")
        print(f"   (Looking for 'solution_*.csv')")
        return

    print(f"{Colors.BOLD}📊 Found {Colors.BLUE}{len(files)}{Colors.RESET}{Colors.BOLD} solution files in '{display_dir}'{Colors.RESET}")

    # Plot initial, middle, and final
    # Limit to a few frames
    if len(files) > 5:
        indices = [0, len(files)//4, len(files)//2, 3*len(files)//4, len(files)-1]
        files_to_plot = [files[i] for i in indices]
    else:
        files_to_plot = files

    plt.figure(figsize=(10, 6))

    # Cycle through line styles for better accessibility (colorblind friendly)
    line_styles = cycle(['-', '--', '-.', ':'])
    
    for f in files_to_plot:
        try:
            data = pd.read_csv(f)
            # Extract step number for label
            step_num = ''.join(filter(str.isdigit, os.path.basename(f)))
            label = f"Step {step_num}" if step_num else os.path.basename(f)
            plt.plot(data['x'], data['u'], label=label, linestyle=next(line_styles), linewidth=2)
        except Exception as e:
            print(f"{Colors.YELLOW}⚠️  Warning: Could not read {f}: {e}{Colors.RESET}")
        
    plt.axvline(x=0.5, color='gray', linestyle='--', alpha=0.7, label='Interface (FV | DG)')
    plt.xlabel('Position (x)')
    plt.ylabel('Value (u)')
    plt.title('Hybrid FV-DG Linear Advection Simulation')
    plt.legend()
    plt.grid(True, alpha=0.3)

    try:
        plt.savefig(args.output, dpi=150)
        print(f"{Colors.BOLD_GREEN}✅ Plot saved to: {Colors.RESET}{Colors.BOLD}{os.path.abspath(args.output)}{Colors.RESET}")
    except Exception as e:
        print(f"{Colors.BOLD_RED}❌ Error saving plot: {e}{Colors.RESET}")

if __name__ == "__main__":
    plot_all()
