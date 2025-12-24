import matplotlib.pyplot as plt
import pandas as pd
import glob
import os

def plot_all():
    # Check if output dir exists, else use current
    search_path = "output/solution_*.csv"
    if not glob.glob(search_path):
        search_path = "solution_*.csv"

    files = sorted(glob.glob(search_path), key=lambda f: int(''.join(filter(str.isdigit, f)) or 999999))
    
    if not files:
        print("No solution files found.")
        return

    # Plot initial, middle, and final
    # Limit to a few frames
    if len(files) > 5:
        indices = [0, len(files)//4, len(files)//2, 3*len(files)//4, len(files)-1]
        files_to_plot = [files[i] for i in indices]
    else:
        files_to_plot = files

    plt.figure(figsize=(10, 6))
    
    for f in files_to_plot:
        data = pd.read_csv(f)
        plt.plot(data['x'], data['u'], label=f"File {f}")
        
    plt.axvline(x=0.5, color='k', linestyle='--', label='Interface (FV | DG)')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('Hybrid FV-DG Linear Advection')
    plt.legend()
    plt.grid(True)
    plt.savefig('advection_plot.png')
    print("Plot saved to advection_plot.png")

if __name__ == "__main__":
    plot_all()
