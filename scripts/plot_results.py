#!/usr/bin/env python3
import sys
import glob
import os
import argparse
import difflib
import subprocess
import platform
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

# Palette: Smart Argument Parser for "Did you mean?" suggestions
class SmartArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"{Colors.BOLD_RED}❌ Error: {message}{Colors.RESET}\n")

        if "unrecognized arguments" in message:
            # Get all valid flags from the parser actions
            valid_flags = []
            for action in self._actions:
                valid_flags.extend(action.option_strings)

            # Check arguments provided in sys.argv for potential typos
            # Note: This is a heuristic scan of sys.argv.
            for arg in sys.argv[1:]:
                if arg.startswith('-') and arg not in valid_flags:
                    # It's an unknown flag, try to find a match
                    matches = difflib.get_close_matches(arg, valid_flags, n=1, cutoff=0.6)
                    if matches:
                        sys.stderr.write(f"       Did you mean '{Colors.YELLOW}{matches[0]}{Colors.RESET}'?\n")

        print(f"\n{Colors.BOLD}Usage:{Colors.RESET}")
        self.print_usage(sys.stderr)
        sys.exit(2)

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

def get_time_from_file(filepath):
    """Extracts time from the first line comment like '# t=0.123'."""
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if first_line.startswith("# t="):
                return float(first_line.strip().split('=')[1])
    except Exception:
        return None
    return None

def get_open_command():
    system = platform.system()
    if system == 'Darwin':
        return 'open'
    elif system == 'Windows':
        return 'start'
    elif system == 'Linux':
        return 'xdg-open'
    return None

def open_file(filepath):
    if platform.system() == 'Windows':
        try:
            if hasattr(os, 'startfile'):
                os.startfile(filepath)
            else:
                # Fallback: empty title arg required when path is quoted
                subprocess.run(['cmd', '/c', 'start', '', filepath], check=True)
            print(f"{Colors.BLUE}👀 Opening preview...{Colors.RESET}")
        except Exception as e:
            print(f"{Colors.YELLOW}⚠️  Warning: Failed to open preview: {e}{Colors.RESET}")
        return

    cmd = get_open_command()
    if not cmd:
        print(f"{Colors.YELLOW}⚠️  Warning: Could not detect how to open files on this OS.{Colors.RESET}")
        return

    try:
        subprocess.run([cmd, filepath], check=True, stderr=subprocess.DEVNULL)
        print(f"{Colors.BLUE}👀 Opening preview...{Colors.RESET}")
    except Exception as e:
        print(f"{Colors.YELLOW}⚠️  Warning: Failed to open preview: {e}{Colors.RESET}")

def plot_all():
    parser = SmartArgumentParser(
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
    parser.add_argument(
        "-p", "--preview",
        action="store_true",
        help="Automatically open the generated plot file."
    )

    args = parser.parse_args()

    # Security: Validate paths
    if args.input_dir and not is_safe_path(args.input_dir):
        print("❌ Error: Invalid input directory. Path traversal ('..') is not allowed.")
        sys.exit(1)

    if not is_safe_path(args.output):
        print("❌ Error: Invalid output filename. Path traversal ('..') is not allowed.")
        sys.exit(1)

    # Security: Prevent overwriting symbolic links to prevent Symlink Attacks
    if Path(args.output).is_symlink():
        print(f"{Colors.BOLD_RED}❌ Error: Output file '{args.output}' is a symbolic link.{Colors.RESET}")
        print(f"   Refusing to overwrite symbolic links to prevent security risks.")
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
        print(f"\n{Colors.YELLOW}💡 Tip: Run the simulation first to generate results:{Colors.RESET}")
        print(f"   ./hyflow1d  (or build/hyflow1d)")
        return

    print(f"{Colors.BOLD}📊 Found {Colors.BLUE}{len(files)}{Colors.RESET}{Colors.BOLD} solution files in '{display_dir}'{Colors.RESET}")

    # Plot initial, middle, and final
    # Limit to a few frames
    if len(files) > 5:
        indices = [0, len(files)//4, len(files)//2, 3*len(files)//4, len(files)-1]
        files_to_plot = [files[i] for i in indices]
    else:
        files_to_plot = files

    # UX: List the specific snapshots being plotted so the user knows what they are looking at
    print(f"{Colors.BOLD}🎨 Plotting {len(files_to_plot)} snapshots:{Colors.RESET}")

    loaded_data = []
    for f in files_to_plot:
        t_val = get_time_from_file(f)
        fname = os.path.basename(f)
        t_str = f"t={t_val:.4f}s" if t_val is not None else "t=?"

        try:
            df = pd.read_csv(f, comment='#')
            min_u, max_u = df['u'].min(), df['u'].max()
            stats = f"[min: {min_u:+.2f}, max: {max_u:+.2f}]"
            loaded_data.append((f, df, t_val))
            print(f"   • {Colors.BLUE}{fname:<20}{Colors.RESET} ({t_str}) {Colors.YELLOW}{stats}{Colors.RESET}")
        except Exception as e:
            print(f"   • {Colors.BLUE}{fname:<20}{Colors.RESET} ({t_str}) {Colors.RED}[read error]{Colors.RESET}")
            loaded_data.append((f, None, t_val))

    plt.figure(figsize=(10, 6))

    # Cycle through line styles for better accessibility (colorblind friendly)
    line_styles = cycle(['-', '--', '-.', ':'])
    
    for f, data, time_val in loaded_data:
        if data is None:
            continue

        try:
            # Extract step number for fallback/auxiliary label
            step_num = ''.join(filter(str.isdigit, os.path.basename(f)))

            if time_val is not None:
                if step_num:
                    label = f"t = {time_val:.2f} (Step {step_num})"
                else:
                    label = f"t = {time_val:.2f}"
            else:
                label = f"Step {step_num}" if step_num else os.path.basename(f)

            plt.plot(data['x'], data['u'], label=label, linestyle=next(line_styles), linewidth=2)
        except Exception as e:
            print(f"{Colors.YELLOW}⚠️  Warning: Could not plot {f}: {e}{Colors.RESET}")
        
    plt.axvline(x=0.5, color='gray', linestyle='--', alpha=0.7, label='Interface (FV | DG)')
    plt.xlabel('Position (x)')
    plt.ylabel('Value (u)')
    plt.title('Hybrid FV-DG Linear Advection Simulation')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Security: Ensure sensitive data visualization is protected
    # Set umask to 0o177 to ensure the file is created with 0600 permissions (rw-------)
    # We save the old umask to restore it later.
    # 0o177 masks group/other rwx and user x.
    old_umask = os.umask(0o177)
    try:
        plt.savefig(args.output, dpi=150)

        # Explicitly enforce 0600 permissions as a second layer of defense (Defense in Depth)
        # This ensures that even if umask was somehow ineffective, we restrict access.
        try:
            os.chmod(args.output, 0o600)
        except OSError as e:
             # If we don't own the file (but could write to it), chmod might fail.
             print(f"{Colors.YELLOW}⚠️  Warning: Could not set secure permissions (0600) on '{args.output}': {e}{Colors.RESET}")

        print(f"{Colors.BOLD_GREEN}✅ Plot saved to: {Colors.RESET}{Colors.BOLD}{os.path.abspath(args.output)}{Colors.RESET}")

        if args.preview:
            open_file(os.path.abspath(args.output))
        else:
            cmd = get_open_command()
            if cmd:
                print(f"{Colors.YELLOW}💡 Tip: View with: {Colors.RESET}{cmd} {args.output}")
                print(f"       Or run with {Colors.BOLD}--preview{Colors.RESET} next time.")

    except Exception as e:
        print(f"{Colors.BOLD_RED}❌ Error saving plot: {e}{Colors.RESET}")
    finally:
        # Always restore the original umask
        os.umask(old_umask)

if __name__ == "__main__":
    plot_all()
