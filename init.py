
import pandas as pd
import numpy as np
from reservoir_input import create_reservoir_input_file


def initialize_reservoir_with_phase_saturations(input_file, output_file):
    """
    Initializes reservoir by calculating pressures, capillary pressures, phase saturations,
    and relative permeabilities using the Stone I model.

    :param input_file: Path to the Excel file containing reservoir input data.
    :param output_file: Path to save the reservoir initialization file.
    """
   
    # Constants
    GRAVITY = 32.174  # ft/s^2

    # Load data from the preprocessed file
    ref_values_df = pd.read_excel(input_file, sheet_name="Reference Values")
    fluid_density_df = pd.read_excel(input_file, sheet_name="Fluid Density")
    gridblock_df = pd.read_excel(input_file, sheet_name="Gridblock Properties")
    petrophysics_df = pd.read_excel(input_file, sheet_name="Petrophysical Data")

    # Extract reference values
    ref_depth_OWC = ref_values_df["Reference Depth_OWC (ft)"].iloc[0]
    ref_pressure_OWC = ref_values_df["Reference Pressure_OWC (psi)"].iloc[0]
    ref_depth_GOC = ref_values_df["Reference Depth_GOC (ft)"].iloc[0]
    ref_pressure_GOC = ref_values_df["Reference Pressure_GOC (psi)"].iloc[0]

    # Extract fluid densities
    densities = fluid_density_df.set_index("Fluid")["Density (lbm/ft3)"]
    oil_density = densities["Oil"]
    water_density = densities["Water"]
    gas_density = densities["Gas"]

    # Extract petrophysical values
    Pe_ow = petrophysics_df["Pe_ow"].values[0]  # Entry pressure for oil-water
    Oil_Water_IFT = petrophysics_df["Oil-Water IFT"].values[0]
    Oil_Gas_IFT = petrophysics_df["Oil-Gas IFT"].values[0]
    Gas_Water_IFT = petrophysics_df["Gas-Water IFT"].values[0]
    S_wr = petrophysics_df["Swr"].iloc[0]
    S_orw = petrophysics_df["Sorw"].iloc[0]  # Residual oil saturation to water
    S_org = petrophysics_df["Sorg"].iloc[0]  # Residual oil saturation to gas
    S_gr = petrophysics_df["Sgr"].iloc[0]  # Residual gas saturation
    k_rw_max = petrophysics_df["Krw_max"].iloc[0]
    k_rg_max = petrophysics_df["Krg_max"].iloc[0]
    k_row_max = petrophysics_df["Krow_max"].iloc[0]
    k_rog_max = petrophysics_df["Krog_max"].iloc[0]
    n_w = petrophysics_df["nw"].iloc[0]
    n_o = petrophysics_df["no"].iloc[0]
    n_g = petrophysics_df["ng"].iloc[0]

    # Calculate Pe_og
    Pe_og = Pe_ow * (Oil_Gas_IFT / Oil_Water_IFT)
    Pe_gw = Pe_ow * (Gas_Water_IFT / Oil_Water_IFT)
    

    # Calculate pressures for each gridblock
    gridblock_df["Pressure Oil (psi)"] = (
        (ref_pressure_OWC + Pe_ow) + (oil_density * (gridblock_df["Grid Depth (ft)"] - ref_depth_OWC) / 144)
    )
    gridblock_df["Pressure Water (psi)"] = (
        ref_pressure_OWC + (water_density * (gridblock_df["Grid Depth (ft)"] - ref_depth_OWC) / 144)
    )
    gridblock_df["Pressure Gas (psi)"] = (
        (ref_pressure_GOC + Pe_og) + (gas_density * (gridblock_df["Grid Depth (ft)"] - ref_depth_GOC) / 144)
    )

    # Calculate capillary pressures
    gridblock_df["Capillary Pressure OW (psi)"] = abs(
        gridblock_df["Pressure Oil (psi)"] - gridblock_df["Pressure Water (psi)"]
    )
    gridblock_df["Capillary Pressure OG (psi)"] = abs(
        gridblock_df["Pressure Gas (psi)"] - gridblock_df["Pressure Oil (psi)"]
    )
    gridblock_df["Capillary Pressure GW (psi)"] = abs(
        gridblock_df["Pressure Gas (psi)"] - gridblock_df["Pressure Water (psi)"]
    )

    # Function to calculate saturations based on depth
    def calculate_saturations(row):
        depth = row["Grid Depth (ft)"]
        Pc_og = row["Capillary Pressure OG (psi)"]
        Pc_ow = row["Capillary Pressure OW (psi)"]
        Pc_gw = row["Capillary Pressure GW (psi)"]

        if depth < ref_depth_GOC:
            S_w = S_wr + (1 - S_wr) * ((Pc_gw / Pe_gw) ** -2)
            S_o = S_org
            S_g = 1 - S_w - S_o
            
            
        elif ref_depth_GOC <= depth < ref_depth_OWC:
            S_w = S_wr + (1 - S_wr) * ((Pc_ow / Pe_ow) ** -2)
            S_g = S_gr
            S_o = 1 - S_w - S_gr
        else:
            S_w = 1.0
            S_g = 0.0
            S_o = 0.0

        

        return pd.Series([S_w, S_o, S_g])

    # Apply the saturation calculation function
    gridblock_df[["Saturation Water", "Saturation Oil", "Saturation Gas"]] = gridblock_df.apply(
        calculate_saturations, axis=1
    )

    def calculate_relative_permeability_stone1(row):
            
        """
        Calculate relative permeability for a grid block using the Stone I model.

        :param row: A row of the DataFrame containing saturations and grid properties.
        :return: Series of krw (water), kro (oil), krg (gas).
        """
        # Extract saturations
        S_w = row["Saturation Water"]
        S_o = row["Saturation Oil"]
        S_g = row["Saturation Gas"]

        # Compute normalized saturations
        S_wD = (S_w - S_wr) / (1 - S_wr - S_orw)
        S_gD = (S_g - S_gr) / (1 - S_wr - S_org - S_gr)
        S_oD = (S_o - S_orw) / (1 - S_wr - S_orw - S_gr)


        # Two-phase relative permeabilities
        k_rw = k_rw_max * (S_wD ** n_w)
        if S_gD >= 1:
            k_rg = k_rg_max 
        elif S_gD <= 0:
             k_rg = 0
        else:
             k_rg = k_rg_max * ((S_gD) ** n_g)
                
         
        k_row = k_row_max * ((1 - S_wD) ** n_w) 
        k_rog = k_rog_max * ((1 - S_gD) ** n_g) 
        # Compute "a" and "S_om" for oil saturation
        a = S_g / (1 - S_wr - S_org) 
        S_om = (1 - a) * S_orw + a * S_org

        # Updated normalized saturations
        S_o_star = (S_o - S_om) / (1 - S_wr - S_om - S_gr) 
        S_w_star = (S_w - S_wr) / (1 - S_wr - S_om - S_gr) 
        S_g_star = (S_g - S_gr) / (1 - S_wr - S_om - S_gr) 

        

        # Three-phase oil relative permeability
        if S_oD <= 0:
            k_ro = 0.0
        elif S_oD >= 1:
            k_ro = k_row_max
        else:
            k_ro = (S_o_star * k_row * k_rog) / (
                    k_row_max * (1 - S_w_star) * (1 - S_g_star))
        
        # Return the relative permeabilities
        return pd.Series([k_rw, k_ro, k_rg])
    
    def ref_pressure(row):
        """ Calculate reservoir pressure for each gridblock by weighting phase pressure
        with phase saturation
        """
        # Extract saturations
        S_w = row["Saturation Water"]
        S_o = row["Saturation Oil"]
        S_g = row["Saturation Gas"]
        
        #Extract Phase_Pressure
        P_w = row['Pressure Water (psi)']
        P_o = row['Pressure Oil (psi)']
        P_g = row['Pressure Gas (psi)']
        
        Reservoir_Pressure = P_w*S_w + P_o*S_o + P_g*S_g
        
        # Return the reservoir pressure
        return pd.Series([Reservoir_Pressure])
    

        

    gridblock_df[['Reservoir_Pressure (psia)']] = gridblock_df.apply(ref_pressure, axis=1)

    # Apply the Stone I relative permeability calculation
    gridblock_df[["krw", "kro", "krg"]] = gridblock_df.apply(calculate_relative_permeability_stone1, axis=1)

    # Prepare output DataFrame
    output_df = gridblock_df[[
        "Grid Block", "Grid Depth (ft)", "Pressure Oil (psi)", "Pressure Water (psi)", "Pressure Gas (psi)",
        "Capillary Pressure OW (psi)", "Capillary Pressure GW (psi)","Reservoir_Pressure (psia)","Capillary Pressure OG (psi)",
        "Saturation Water", "Saturation Oil", "Saturation Gas",
        "krw", "kro", "krg"
    ]]

    # Save to a new Excel file
    output_df.to_excel(output_file, index=False)
    print(f"Reservoir initialization with Stone I relative permeabilities completed. Results saved to: {output_file}")

# Example Usage
if __name__ == "__main__":
    input_file = "Res_Input.xlsx"
    output_file = "Res_Init_With_Phase_Saturations.xlsx"

    # Run the reservoir initialization with phase saturations and relative permeabilities
    initialize_reservoir_with_phase_saturations(input_file, output_file)





