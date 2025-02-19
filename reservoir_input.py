
import os
import pandas as pd
import numpy as np
from openpyxl import load_workbook
from openpyxl.worksheet.datavalidation import DataValidation


# Declare global variables
num_x_blocks = int(input("Enter number of grid blocks in X-direction: "))
num_y_blocks = int(input("Enter number of grid blocks in Y-direction: "))
num_z_blocks = int(input("Enter number of grid blocks in Z-direction: "))


def add_dropdown(worksheet, cell_range, options):
    """
    Adds a dropdown validation to a specified range in a worksheet.

    :param worksheet: The target worksheet to add the dropdown.
    :param cell_range: The Excel range (e.g., "B2:B100") to apply the dropdown.
    :param options: A comma-separated string of dropdown options.
    """
    dv = DataValidation(type="list", formula1=f'"{options}"', allow_blank=True)
    worksheet.add_data_validation(dv)
    dv.add(cell_range)

def get_grid_block_sizes(total_length, num_blocks, custom_sizes=None):
    """
    Computes grid block sizes based on user input.

    :param total_length: Total length of the reservoir in one dimension.
    :param num_blocks: Number of grid blocks.
    :param custom_sizes: List of custom sizes for each block (optional).
    :return: List of grid block sizes.
    """
    if custom_sizes:
        if len(custom_sizes) != num_blocks:
            raise ValueError(f"Custom sizes must match the number of grid blocks ({num_blocks}).")
        if sum(custom_sizes) != total_length:
            raise ValueError(f"Custom sizes must sum up to the total reservoir length ({total_length} ft).")
        return custom_sizes
    else:
        # Uniform grid block size
        return [total_length / num_blocks] * num_blocks

def create_reservoir_input_file(output_file):
    """
    Creates an Excel workbook for reservoir input with PVT, fluid density, gridblock properties,
    and dropdowns for Well Model.

    :param output_file: Path to save the Excel workbook.
    """
    global num_x_blocks, num_y_blocks, num_z_blocks
    # Check if the file exists
    if os.path.exists(output_file):
        print(f"The file '{output_file}' already exists.")
        user_choice = input("Would you like to overwrite it? (yes/no/rename): ").strip().lower()
        if user_choice == "no":
            print("Operation aborted. The file will not be overwritten.")
            return
        elif user_choice == "rename":
            new_name = input("Enter a new file name (with .xlsx extension): ").strip()
            output_file = new_name

    # Create a workbook
    with pd.ExcelWriter(output_file) as writer:
        # PVT Worksheet
        pvt_data = {
            "Pressure (psi)": [],
            "Oil FVF (rb/stb)": [],
            "Gas FVF (rb/scf)": [],
            "Water FVF (rb/stb)": [],
            "Oil Viscosity (cp)": [],
            "Water Viscosity (cp)": [],
            "Gas Viscosity (cp)": [],
            "Oil Compressibility (1/psi)": [],
            "Water Compressibility (1/psi)": [],
            "Gas Compressibility (1/psi)": [],
            'Solution Gas-Oil Ratio (scf/stb)': [],
        }
        pd.DataFrame(pvt_data).to_excel(writer, sheet_name="PVT", index=False)

        # Fluid Density Worksheet
        density_data = {
            "Fluid": ["Oil", "Water", "Gas"],
            "Density (lbm/ft3)": ["", "", ""]
        }
        pd.DataFrame(density_data).to_excel(writer, sheet_name="Fluid Density", index=False)
        
        # Reference Depth and Pressure 
        Reference_depth = {'Reference Depth_OWC (ft)':[], 'Reference Pressure_OWC (psi)':[], 
                           'Reference Depth_GOC (ft)':[], 'Reference Pressure_GOC (psi)':[], 
                           'Rock_Compressibility':[]}
        pd.DataFrame(Reference_depth).to_excel(writer, sheet_name="Reference Values", index = False)
        
        # Numerical Worksheet
        numerical_data = {
            "Simulation Time (days)": [],
            "Timestep Size (days)": []
        }
        pd.DataFrame(numerical_data).to_excel(writer, sheet_name="Numerical", index=False)

        # Petrophysical Properties Worksheet
        petrophysical_data = {'Oil-Gas IFT':[], 'Oil-Water IFT':[], 'Gas-Water IFT': [], 
                              'Pe_ow':[], 'Swr':[], 'Sgr':[], 'Krow_max':[], 'Krw_max':[], 
                              'Krg_max':[], 'nw':[], 'no':[], 'ng':[], 'now':[], 'nog':[], 
                              'Krog_max':[], 'Sorw':[], 'Sorg':[]}
        pd.DataFrame(petrophysical_data).to_excel(writer, sheet_name="Petrophysical Data", index = False)

        # Gridblock Properties Worksheet
        reservoir_length = float(input("Enter reservoir length in X-direction (ft): "))
        reservoir_width = float(input("Enter reservoir width in Y-direction (ft): "))
        reservoir_height = float(input("Enter reservoir height in Z-direction (ft): "))
        num_x_blocks = int(input("Enter number of grid blocks in X-direction: "))
        num_y_blocks = int(input("Enter number of grid blocks in Y-direction: "))
        num_z_blocks = int(input("Enter number of grid blocks in Z-direction: "))

        # Custom sizes for X-direction
        custom_x_sizes = input("Enter custom grid sizes in X-direction (comma-separated, or leave blank for uniform): ")
        custom_x_sizes = (
            [float(size) for size in custom_x_sizes.split(",")] if custom_x_sizes else None
        )
        grid_size_x = get_grid_block_sizes(reservoir_length, num_x_blocks, custom_x_sizes)

        # Custom sizes for Y-direction
        custom_y_sizes = input("Enter custom grid sizes in Y-direction (comma-separated, or leave blank for uniform): ")
        custom_y_sizes = (
            [float(size) for size in custom_y_sizes.split(",")] if custom_y_sizes else None
        )
        grid_size_y = get_grid_block_sizes(reservoir_width, num_y_blocks, custom_y_sizes)

        # Custom sizes for Z-direction
        custom_z_sizes = input("Enter custom grid sizes in Z-direction (comma-separated, or leave blank for uniform): ")
        custom_z_sizes = (
            [float(size) for size in custom_z_sizes.split(",")] if custom_z_sizes else None
        )
        grid_size_z = get_grid_block_sizes(reservoir_height, num_z_blocks, custom_z_sizes)

        grid_data = {
            "Grid Block": range(1, (num_x_blocks * num_y_blocks * num_z_blocks) + 1),
            "Porosity": np.random.uniform(0.2, 0.3, num_x_blocks * num_y_blocks * num_z_blocks),
            "Permeability X (mD)": np.random.uniform(50, 200, num_x_blocks * num_y_blocks * num_z_blocks),
            "Permeability Y (mD)": np.random.uniform(50, 200, num_x_blocks * num_y_blocks * num_z_blocks),
            "Permeability Z (mD)": np.random.uniform(50, 200, num_x_blocks * num_y_blocks * num_z_blocks),
            "Grid Size X (ft)": np.tile(grid_size_x, num_y_blocks * num_z_blocks),
            "Grid Size Y (ft)": np.repeat(np.tile(grid_size_y, num_z_blocks), num_x_blocks),
            "Grid Size Z (ft)": np.repeat(grid_size_z, num_x_blocks * num_y_blocks),
            "Grid Depth (ft)": np.random.uniform(1000, 2000, num_x_blocks * num_y_blocks * num_z_blocks) 
        }
        pd.DataFrame(grid_data).to_excel(writer, sheet_name="Gridblock Properties", index=False)

        # Well Model Worksheet
        well_data = {
            "Well ID": [],
            "Well Type (Injector/Producer)": [],
            "Well Deviation (Vertical/Horizontal)": [],
            "Wellbore Radius (ft)": [],
            "Skin": [],
            "Control Mode (Rate/BHP)": [],
            "Rate (STB/day)": [],
            "BHP (psi)": [],
            "Perforated Gridblocks": [],
        }
        pd.DataFrame(well_data).to_excel(writer, sheet_name="Well Model", index=False)

        # Boundary Conditions Worksheet
        boundary_conditions_data = {
            "Boundary Location": ["(0, L)", "(L, 0)", "(0, W)", "(W, 0)", "(0, H)", "(H, 0)"],
            "Boundary Condition": ["", "", "", "", "", ""],
            "Value": ["", "", "", "", "", ""]
        }
        pd.DataFrame(boundary_conditions_data).to_excel(writer, sheet_name="Boundary Conditions", index=False)

    # Open the created workbook to add dropdowns
    workbook = load_workbook(output_file)
    boundary_ws = workbook["Boundary Conditions"]
    well_ws = workbook["Well Model"]

    # Add dropdowns for Boundary Conditions
    add_dropdown(boundary_ws, "B2:B7", "Dirichlet,Neumann,No Flow")
    
    # Add dropdowns for the Well Model
    add_dropdown(boundary_ws, "B2:B7", "Dirichlet,Neumann,No Flow")
    add_dropdown(well_ws, "B2:B100", "Injector,Producer")
    add_dropdown(well_ws, "C2:C100", "Vertical,Horizontal")
    add_dropdown(well_ws, "F2:F100", "Rate,BHP")

    # Add dropdowns for Boundary Conditions
    add_dropdown(boundary_ws, "B2:B7", "Dirichlet,Neumann,No Flow")

    # Save the workbook with dropdowns
    workbook.save(output_file)
    print(f"Reservoir input file created with dropdowns and grid customization: {output_file}")

create_reservoir_input_file("Res_Input.xlsx")
