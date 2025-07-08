import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d

# Introduction 
print("\n//////////////////////////////////////////////////////////////////////////////////")
print("INSTITUTO TECNOLÓGICO DE COSTA RICA")
print("Escuela de Ingeniería en Construcción")
print("https://www.tec.ac.cr")
print("Session: FLUJO NO-UNIFORME/ FLUJO CRÍTICO")

print("\nM.Sc. Eng. Maikel Méndez M")
print("Water Resources + GIS + DataScience")
print("Instituto Tecnológico de Costa Rica")
print("https://www.tec.ac.cr")
print("https://orcid.org/0000-0003-1919-141X")
print("https://www.scopus.com/authid/detail.uri?authorId=51665581300")
print("https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en")
print("https://www.youtube.com/c/maikelmendez")
print("https://github.com/maikelonu")
print("//////////////////////////////////////////////////////////////////////////////////")

print("\n////////////////////////////////////////////////////////")
print("BLOCK: Declarations")
print("////////////////////////////////////////////////////////")
base_m = 0.086  # Hydraulic flume base (m)
q = 4.00  # Water flow (m3/h)
q /= 3600  #Water flow (m3/s)

# Triangular Weir dimensions and vectors are defined
block_ID = ["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9"] # Control-points IDs
block_X = np.array([15.050, 15.505, 15.550, 15.602, 15.670, 15.740, 15.810, 15.870, 16.120]) # Control-points X distance (m)
block_DZ = np.array([0.45, 0.45, 0.45, 0.45, 0.50, 0.55, 0.55, 0.55, 0.45]) # Hydraulic flume DeltaZ to bottom (cm)
block_ELEV = np.array([0.0, 0.0, 2.7, 6.2, 4.6, 3.0, 1.5, 0.0, 0.0]) # Relative elevation of the Triangular Weir (cm)

# Experimental Data !!!!!!!!!!!!!!
block_Y = np.array([10.40, 10.35, 10.25, 9.65, 7.10, 5.10, 3.50, 2.15, 1.90]) # experimental water depth (cm)

# Vectors are transformed
block_X = block_X - 15.05  # Normalized to start in 0 (m)
block_DZ /= 100  # Transformed to (m)
block_ELEV /= 100  # Transformed to (m)
block_Y /= 100  # Transformed to (m)
block_YEFF = block_Y - block_DZ  # Effective water depth is calculated for plotting (m)
block_YEFF2 = block_YEFF - block_ELEV  # effective water depth is calculated for numerica analysis (m)

# A dataFrame object is created   
df_weir = pd.DataFrame({
    "block_ID": block_ID,
    "block_X": block_X,
    "block_DZ": block_DZ,
    "block_ELEV": block_ELEV,
    "block_Y": block_Y,
    "block_YEFF": block_YEFF,
    "block_YEFF2": block_YEFF2
})

# Static energy is calculated (m)
df_weir["block_energ_EST"] = df_weir["block_YEFF2"]

# Dynamic energy is calculated (m)
df_weir["block_energ_DYM"] = (q**2) / (2 * 9.81 * (base_m**2) * (df_weir["block_energ_EST"]**2))

# Total energy is calculated (m)
df_weir["block_energ_TOTAL"] = df_weir["block_energ_EST"] + df_weir["block_energ_DYM"]

# Theoretical total-energy is calculated (m) 
df_weir["block_energ_TOTAL_TEO"] = df_weir["block_energ_TOTAL"].iloc[0]

# Energy losses are calculated (m)
df_weir["block_energ_LOSS"] = df_weir["block_energ_TOTAL_TEO"] - df_weir["block_energ_TOTAL"]

# Energy losses are calculated for plotting(m) 
df_weir["block_energ_PLOT"] = df_weir["block_energ_TOTAL"] + df_weir["block_ELEV"]

# Hydraulic area is calculated (m2)
df_weir["area"] = df_weir["block_YEFF2"] * base_m

# Hydraulic perimeter is calculated (m)
df_weir["perimeter"] = (df_weir["block_YEFF2"] * 2) + base_m

# Water velocity is calculated (m/s)
df_weir["vel"] = q / df_weir["area"]

# Froude number is calculated
df_weir["Froude"] = df_weir["vel"] / np.sqrt(df_weir["area"] * 9.81 / base_m)

# If-statement & For-loop for Froude number
df_weir["rule"] = np.where(df_weir["Froude"] > 1, "SUPER", "SUB")

# Approxfun {stats} is used to interpolate yc (m)
tfun01 = interp1d(df_weir["Froude"], df_weir["block_energ_EST"], kind='linear', fill_value="extrapolate") # Froude number = 1

# Approxfun {stats} is used to interpolate yc position (m)
tfun02 = interp1d(df_weir["Froude"], df_weir["block_X"], kind='linear', fill_value="extrapolate") # Froude number = 1

# Approxfun {stats} is used to interpolate Ec (m)
tfun03 = interp1d(df_weir["Froude"], df_weir["block_energ_TOTAL"], kind='linear', fill_value="extrapolate") # Froude number = 1

# Yc is estimated (m)
yc_inter = np.round(tfun01(1), 4)

# Yc distance is estimated (m)
yc_inter_dist = np.round(tfun02(1), 4)

# Ec is estimated (m)
Ec = np.round(tfun03(1), 4)

# Yc is calculated based on ec.03 (m)
yc_ec03 = round(((q**2) / ((base_m**2) * 9.81)) ** (1/3), 4)

# Yc is calculated based on ec.04 (m)
yc_ec04 = (2/3) * Ec

# Colors for plotting
cols = {
    "pto_control": "firebrick",
    "Energ_Est": "deepskyblue",
    "base_canal": "black",
    "vertedor": "darkgreen",
    "Energ_Total_Teo": "magenta",
    "Energ_Total_Exp": "lightsalmon",
    "yc": "#009999",
    "y_effec": "#00cc00"
}

# Plot 1: Energy Profile with 'p1' to 'p9' labels
plt.figure(figsize=(10, 6))
plt.plot(df_weir["block_X"], df_weir["block_ELEV"], color=cols["vertedor"], linewidth=1.5, label="Vertedor")
plt.axhline(0.0, color=cols["base_canal"], linewidth=1.5, label="Base del Canal")
plt.axhline(df_weir["block_energ_TOTAL_TEO"].iloc[0], color=cols["Energ_Total_Teo"], linewidth=0.75, linestyle="--", label="Energía Total Teórica")
plt.plot(df_weir["block_X"], df_weir["block_YEFF"], color=cols["Energ_Est"], linewidth=1.5, label="Energía Estática")
plt.plot(df_weir["block_X"], df_weir["block_energ_PLOT"], color=cols["Energ_Total_Exp"], linewidth=1.0, label="Energía Total Experimental")

# Adding labels for 'p1' to 'p9'
for i, txt in enumerate(df_weir['block_ID']):
    plt.text(df_weir["block_X"].iloc[i], df_weir["block_ELEV"].iloc[i] + 0.005, txt, fontsize=10, color='gray', ha='center')

# Plot settings
plt.axvline(yc_inter_dist, color=cols["yc"], linewidth=1.25, linestyle="--", label="yc (Froude=1)")
plt.xlabel("Distancia (m)")
plt.ylabel("Energía (m)")
plt.title("Perfil Energético. Vertedor Triangular")
plt.legend(loc="upper right")
plt.grid(True)

# Show the plot
plt.show()

# Plot 2: Specific Energies with 'p1' to 'p9' labels
plt.figure(figsize=(10, 6))
plt.scatter(df_weir["block_energ_TOTAL"], df_weir["block_energ_EST"], color='#0000cc', s=50, label="Static Energy (Exp)")
plt.plot(df_weir["block_energ_TOTAL"], df_weir["block_energ_EST"], color='#0000cc', linewidth=0.75)

# Adding labels for 'p1' to 'p9'
for i, txt in enumerate(df_weir['block_ID']):
    plt.text(df_weir["block_energ_TOTAL"].iloc[i], df_weir["block_energ_EST"].iloc[i] + 0.005, txt, fontsize=10, color='gray', ha='center')

# Plot settings
plt.axhline(yc_ec03, color='#009999', linewidth=0.75, linestyle="--", label="yc")
plt.xlabel("Total Energy (m)")
plt.ylabel("Piezometric Energy (m)")
plt.title("Specific Energies. Triangular Weir")
plt.legend(loc="upper right")
plt.grid(True)

# Show the plot
plt.show()

print("\n/////////////////////////////////////////////////////////////")
print("END OF SCRIPT")
print("/////////////////////////////////////////////////////////////")