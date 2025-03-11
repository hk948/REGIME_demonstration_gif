import sys

sys.stderr = open('error.log', 'w')

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import random
import matplotlib

matplotlib.use("TkAgg")  # Use Tkinter backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

# Global variables (to store simulation results)
matrix_C_stock = None
mean_series = None
example_series = None

# Default parameter values (simulation_iteration is fixed)
default_params = {
    "dstb_interval_1": "40",
    "dstb_severity_1": "0.4",
    "dstb_interval_2": "130",
    "dstb_severity_2": "1.0",
    "initial_C_stock": "0",
    "C_influx": "2",
    "residence_time": "50",
    "simulation_length": "1000",
    "simulation_iteration": "100"  # Fixed value
}

# Descriptive labels for each parameter
param_labels = {
    "dstb_interval_1": "Return interval of disturbance type A (years)",
    "dstb_severity_1": "Severity of disturbance type A (fraction)",
    "dstb_interval_2": "Return interval of disturbance type B (years)",
    "dstb_severity_2": "Severity of disturbance type B (fraction)",
    "initial_C_stock": "Initial carbon storage (Mg C ha⁻¹)",
    "C_influx": "Carbon input (Mg C ha⁻¹ year⁻¹)",
    "residence_time": "Carbon turnover time (years)",
    "simulation_length": "Simulation length (years)"
}

# Create main window (parameter settings)
root = tk.Tk()
root.title("Parameters for REGIME")

# Create parameter input frame
param_frame = ttk.Frame(root, padding="10")
param_frame.grid(row=0, column=0, sticky="nsew")

entries = {}


# Create BooleanVars for disturbance types A and B
consider_A = tk.BooleanVar(value=True)
consider_B = tk.BooleanVar(value=True)

# Place checkbuttons at the top of the parameter input frame (rows 0 and 1)
check_A = ttk.Checkbutton(param_frame, text="Consider disturbance type A", variable=consider_A,
                           command=lambda: toggle_entry_state(consider_A, entries["dstb_interval_1"], entries["dstb_severity_1"]))
check_A.grid(row=0, column=0, columnspan=2, sticky="w", pady=2)

check_B = ttk.Checkbutton(param_frame, text="Consider disturbance type B", variable=consider_B,
                           command=lambda: toggle_entry_state(consider_B, entries["dstb_interval_2"], entries["dstb_severity_2"]))
check_B.grid(row=1, column=0, columnspan=2, sticky="w", pady=2)

# Then, when creating the numeric entry fields, start the row index from 2:
row = 2
for key, default_val in default_params.items():
    if key in ["dstb_interval_1", "dstb_severity_1", "dstb_interval_2", "dstb_severity_2",
               "initial_C_stock", "C_influx", "residence_time", "simulation_length"]:
        label_text = param_labels.get(key, key)
        ttk.Label(param_frame, text=label_text).grid(row=row, column=0, sticky="w", pady=2)
        entry = ttk.Entry(param_frame, width=20)
        entry.insert(0, default_val)

        entry.grid(row=row, column=1, pady=2)
        entries[key] = entry
        row += 1

def toggle_entry_state(var, interval_entry, severity_entry):
    # If not checked, disable the corresponding entries; else, enable them.
    if var.get():
        interval_entry.state(["!disabled"])
        severity_entry.state(["!disabled"])
    else:
        interval_entry.state(["disabled"])
        severity_entry.state(["disabled"])

# "Run Simulation" button
run_button = ttk.Button(param_frame, text="Run Simulation")
run_button.grid(row=row, column=0, columnspan=2, pady=10)

# Window for simulation results and animation (Toplevel window)
sim_window = None  # For animation window

def run_simulation():
    global matrix_C_stock, mean_series, example_series, sim_window
    try:
        dstb_interval_1 = int(entries["dstb_interval_1"].get()) if consider_A.get() else None
        dstb_severity_1 = float(entries["dstb_severity_1"].get()) if consider_A.get() else None
        dstb_interval_2 = int(entries["dstb_interval_2"].get()) if consider_B.get() else None
        dstb_severity_2 = float(entries["dstb_severity_2"].get()) if consider_B.get() else None
        initial_C_stock = float(entries["initial_C_stock"].get())
        C_influx = float(entries["C_influx"].get())
        residence_time = float(entries["residence_time"].get())
        simulation_length = int(entries["simulation_length"].get())
        simulation_iteration = 100
    except Exception as e:
        messagebox.showerror("Input Error", f"Error in input values:\n{e}")
        return

    if consider_A.get() and dstb_interval_1 == 0:
        messagebox.showerror("Input Error",
                             "Disturbance type A return interval cannot be 0.")
        return

    if consider_B.get() and dstb_interval_2 == 0:
        messagebox.showerror("Input Error",
                             "Disturbance type B return interval cannot be 0.")
        return

    if consider_A.get() and (dstb_severity_1 > 1 or dstb_severity_1 < 0):
        messagebox.showerror("Input Error",
                             "Disturbance type A severity should be between 0 and 1.")
        return

    if consider_B.get() and (dstb_severity_2 > 1 or dstb_severity_2 < 0):
        messagebox.showerror("Input Error",
                             "Disturbance type B severity should be between 0 and 1.")
        return

    # Check: if a type is not considered, its fields are not used.
    # If both are not considered, then set real_C = max_C.
    if not consider_A.get() and not consider_B.get():
        compound_interval = None
        compound_severity = None
        real_C = C_influx * residence_time  # max_C
    elif consider_A.get() and not consider_B.get():
        compound_interval = dstb_interval_1
        compound_severity = dstb_severity_1
        real_C = (C_influx * residence_time * compound_interval) / (compound_interval + residence_time * compound_severity)
        # Here real_C becomes less meaningful; you might choose to just use the chosen type's values.
    elif not consider_A.get() and consider_B.get():
        compound_interval = dstb_interval_2
        compound_severity = dstb_severity_2
        real_C = (C_influx * residence_time * compound_interval) / (compound_interval + residence_time * compound_severity)
    else:
        compound_interval = 1 / (1/dstb_interval_1 + 1/dstb_interval_2)
        compound_severity = (dstb_severity_1/dstb_interval_1 + dstb_severity_2/dstb_interval_2) / (1/dstb_interval_1 + 1/dstb_interval_2)
        real_C = (C_influx * residence_time * compound_interval) / (compound_interval + residence_time * compound_severity)
        theta_1 = dstb_severity_1 / dstb_interval_1
        theta_2 = dstb_severity_2 / dstb_interval_2
        theta = theta_1 + theta_2
        I_1_p = round(theta_1/theta * 100, 1)
        I_2_p = round(theta_2/theta * 100, 1)

    max_C = C_influx * residence_time

    matrix_C_stock = np.zeros((simulation_iteration, simulation_length))
    matrix_C_stock[:, 0] = initial_C_stock

    for year in range(1, simulation_length):
        for itr in range(simulation_iteration):
            # Use considered types only. If a type is not considered, its chance is 0.
            chance_A = (1.0 / dstb_interval_1) if consider_A.get() else 0
            chance_B = (1.0 / dstb_interval_2) if consider_B.get() else 0
            r = random.random()
            if r < chance_A:
                matrix_C_stock[itr, year] = matrix_C_stock[itr, year - 1] * (1 - dstb_severity_1)
            elif r < chance_A + chance_B:
                matrix_C_stock[itr, year] = matrix_C_stock[itr, year - 1] * (1 - dstb_severity_2)
            else:
                matrix_C_stock[itr, year] = C_influx + matrix_C_stock[itr, year - 1] * (1 - 1.0 / residence_time)

    mean_series = matrix_C_stock.mean(axis=0)
    # Use the specified example iteration
    example_series = matrix_C_stock[99, :]

    global sim_window
    if sim_window is None or not tk.Toplevel.winfo_exists(sim_window):
        sim_window = tk.Toplevel(root)
        sim_window.title("Simulated impacts of disturbance regimes on carbon storage dynamics with REGIME-defined potentials")
    else:
        sim_window.lift()

    # Build metadata string from user inputs
    metadata_str = "\n"
    if consider_A.get():
        metadata_str += f"Disturbance A: Interval = {dstb_interval_1} years; Severity = {dstb_severity_1}\n"
    else:
        metadata_str += "Disturbance A: Not considered\n"
    if consider_B.get():
        metadata_str += f"Disturbance B: Interval = {dstb_interval_2} years; Severity = {dstb_severity_2}\n"
    else:
        metadata_str += "Disturbance B: Not considered\n"

    metadata_str += f"Initial carbon storage = {initial_C_stock} Mg C ha⁻¹\n"
    metadata_str += f"Carbon input = {C_influx} Mg C ha⁻¹ year⁻¹\n"
    metadata_str += f"Carbon turnover time = {residence_time} years\n"

    if consider_A.get() and consider_B.get():
        metadata_str += f"Relative contribution to potential reduction:\n A = {I_1_p} %; B = {I_2_p} %"
    else:
        metadata_str += f"Relative contribution to potential reduction:\nNot considered"


    run_animation(sim_window, max_C, real_C, simulation_length, metadata_str)

def run_animation(window, max_C, real_C, simulation_length, metadata_str):
    # Create Figure with a custom gridspec so that the left plot is wider.
    fig, (ax_time, ax_grid) = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [3, 2]})

    # Increase margins to prevent clipping on screen
    fig.subplots_adjust(left=0.08, right=0.92, top=0.92, bottom=0.08)
    # Configure time-series plot (legend fixed at upper right)
    line_mean, = ax_time.plot([], [], 'k-', linewidth=2, label="Mean of 100 cells")
    line_iter, = ax_time.plot([], [], 'r-', linewidth=1, label="Example cell (Red bordered)")
    ax_time.axhline(y=max_C, color='blue', linestyle='dashed', linewidth=2,
                    label=f"Maximal potential ({max_C:.1f} Mg C ha⁻¹)")
    ax_time.axhline(y=real_C, color='black', linestyle='dashed', linewidth=1.5,
                    label=f"Realizable potential ({real_C:.1f} Mg C ha⁻¹)")
    ax_time.set_xlim(0, simulation_length)
    ax_time.set_ylim(0, max_C * 1.1)
    ax_time.set_xlabel("Year")
    ax_time.set_ylabel("AGBC (Mg C ha⁻¹)")
    ax_time.set_title("Time Series")
    ax_time.legend(loc='upper right')

    # Configure grid image: use viridis_r colormap, and set extent so each cell is 1x1
    grid_img = ax_grid.imshow(np.zeros((10, 10)), cmap="viridis_r", vmin=0, vmax=max_C,
                              extent=[-0.5, 9.5, 9.5, -0.5])
    ax_grid.set_title("Spatial Pattern (10x10 grid)")
    ax_grid.axis('off')

    # Add a colorbar to the right grid image
    cbar = fig.colorbar(grid_img, ax=ax_grid, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label("AGBC (Mg C ha⁻¹)")

    # Since we want the down-right cell
    red_rect = patches.Rectangle((9 - 0.5, 9.5 - (0 + 1)), 0.96, 0.96, linewidth=3, edgecolor='red',
                                 facecolor='none')
    ax_grid.add_patch(red_rect)

    def init():
        line_iter.set_data([], [])
        line_mean.set_data([], [])
        grid_img.set_data(np.zeros((10, 10)))
        return line_iter, line_mean, grid_img

    def update(t):
        years = np.arange(t + 1)
        line_iter.set_data(years, example_series[:t + 1])
        line_mean.set_data(years, mean_series[:t + 1])
        grid_data = matrix_C_stock[:, t].reshape(10, 10)
        grid_img.set_data(grid_data)
        ax_grid.set_title(f"Spatial Pattern (Year: {t})")
        ax_time.set_title(f"Time Series (up to Year: {t})")
        return line_iter, line_mean, grid_img

    anim = FuncAnimation(fig, update, frames=simulation_length, init_func=init, interval=50, blit=False)

    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)

    # Add metadata text below the left plot using fig.text
    fig.text(0.62, 0.02, metadata_str, fontsize=10, ha='left', va='bottom', wrap=True)

    def save_gif():
        file_path = filedialog.asksaveasfilename(defaultextension=".gif",
                                                 filetypes=[("GIF files", "*.gif"), ("All Files", "*.*")])
        if file_path:
            try:
                fig.subplots_adjust(left=0.08, right=0.92, top=0.92, bottom=0.08)
                # To reduce file size, skip frames: use every 2nd frame.
                skip = 2
                save_frames = range(0, simulation_length, skip)
                # Create a new animation with skipped frames for saving
                anim2 = FuncAnimation(fig, update, frames=save_frames, init_func=init, interval=50, blit=False)

                anim2.save(file_path, writer="pillow", fps=20/skip)
                messagebox.showinfo("Save GIF",
                                    f"Animation saved as GIF:\n{file_path}")
            except Exception as e:
                messagebox.showerror("Save Error", f"Error saving GIF:\n{e}")

    save_button = ttk.Button(window, text="Save GIF", command=save_gif)
    save_button.pack(pady=5)


run_button.config(command=run_simulation)

root.mainloop()

