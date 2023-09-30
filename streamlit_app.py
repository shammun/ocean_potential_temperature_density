import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
import math
import seaborn as sns

# (Your functions: density0, density, adiabatic_lapse_rate, etc... remain unchanged)

def density0(s, t):
    A = 1.001685e-04 + t * (-1.120083e-06 + t * 6.536332e-09)
    A = 999.842594 + t * (6.793952e-02 + t * (-9.095290e-03 + t * A))
    B = 7.6438e-05 + t * (-8.2467e-07 + t * 5.3875e-09)
    B = 0.824493 + t * (-4.0899e-03 + t * B)
    C = -5.72466e-03 + t * (1.0227e-04 - t * 1.6546e-06)
    D = 4.8314e-04
    dens0 = A + s * (B + C * math.sqrt(s) + D * s)
    return dens0

def density(s, t, p):
    t2 = t * t
    t3 = t2 * t
    t4 = t3 * t

    d0 = density0(s, t)
    E = 19652.21 + 148.4206 * t - 2.327105 * t2 + 1.360477e-2 * t3 - 5.155288e-5 * t4
    F = 54.6746 - 0.603459 * t + 1.09987e-2 * t2 - 6.1670e-5 * t3
    G = 7.944e-2 + 1.6483e-2 * t - 5.3009e-4 * t2
    H = 3.239908 + 1.43713e-3 * t + 1.16092e-4 * t2 - 5.77905e-7 * t3
    I = 2.2838e-3 - 1.0981e-5 * t - 1.6078e-6 * t2
    J = 1.91075e-4
    M = 8.50935e-5 - 6.12293e-6 * t + 5.2787e-8 * t2
    N = -9.9348e-7 + 2.0816e-8 * t + 9.1697e-10 * t2

    s1p5 = s * math.sqrt(s)
    pb = p / 10
    K = (E + F * s + G * s1p5) + (H + I * s + J * s1p5) * pb + (M + N * s) * pb * pb
    d = d0 / (1 - pb/K)
    return d

def adiabatic_lapse_rate(s, t, p):
    ds = s - 35.0
    atg = ((-2.1687e-16 * t + 1.8676e-14) * t - 4.6206e-13) * p * p
    atg += (2.7759e-12 * t - 1.1351e-10) * ds * p
    atg += (((-5.4481e-14 * t + 8.7330e-12) * t - 6.7795e-10) * t + 1.8741e-8) * p
    atg += (-4.2393e-8 * t + 1.8932e-6) * ds
    atg += ((6.6228e-10 * t - 6.8360e-8) * t + 8.5258e-6) * t + 3.5803e-5
    return atg

def potential_temperature(s, t0, p0):
    p = p0
    t = t0
    h = 0 - p
    xk = h * adiabatic_lapse_rate(s, t, p)
    t += 0.5 * xk
    q = xk
    p += 0.5 * h
    xk = h * adiabatic_lapse_rate(s, t, p)
    t += 0.29289322 * (xk - q)
    q = 0.58578644 * xk + 0.121320344 * q
    xk = h * adiabatic_lapse_rate(s, t, p)
    t += 1.707106781 * (xk - q)
    q = 3.414213562 * xk - 4.121320344 * q
    p += 0.5 * h
    xk = h * adiabatic_lapse_rate(s, t, p)
    theta = t + (xk - 2.0 * q) / 6.0
    return theta

def density_anomaly(S, T): # here T means potential Temperature 
    row_st0 = (999.842594) + (6.793952E-2 *T) - (9.095290E-3 * (T**2)) - (1.001685E-4 * (T**3)) - (1.120083E-6 * (T**4)) + (6.536332E-9 * (T**5))+(S * (0.824493 - (4.0899E-3 * T) + (7.6438E-5 * (T**2)) - (8.2467E-7 * (T**3)) + (5.3875E-9 * (T**4)))) +((S**(3/2)) * (-5.72466E-3 + (1.0227E-4 * T) - (1.6546E-6 * (T**2)))) + ( 4.8314E-4 * (S**2))
    potential_density = row_st0 -1000
    return potential_density

def potential_density(S, T): # here T means potential Temperature 
    row_st0 = (999.842594) + (6.793952E-2 *T) - (9.095290E-3 * (T**2)) - (1.001685E-4 * (T**3)) - (1.120083E-6 * (T**4)) + (6.536332E-9 * (T**5))+(S * (0.824493 - (4.0899E-3 * T) + (7.6438E-5 * (T**2)) - (8.2467E-7 * (T**3)) + (5.3875E-9 * (T**4)))) +((S**(3/2)) * (-5.72466E-3 + (1.0227E-4 * T) - (1.6546E-6 * (T**2)))) + ( 4.8314E-4 * (S**2))
    potential_density = row_st0 -1000
    return potential_density


def compute_density_from_df(df):
    # Check if the required columns 'salinity', 'temperature', 'pressure' exist in the dataframe
    if all(col in df.columns for col in ['salinity', 'temperature', 'pressure']):
        # Compute the density values using the provided functions and append to a new column
        df['density'] = df.apply(lambda row: density(row['salinity'], row['temperature'], row['pressure']), axis=1)
    else:
        raise ValueError("The required columns 'salinity', 'temperature', 'pressure' must be present in the CSV file.")
    
    return df




def compute_potential_temperature_from_df(df):
    # Check if the required columns 'salinity', 'temperature', 'pressure' exist in the dataframe
    if all(col in df.columns for col in ['salinity', 'temperature', 'pressure']):
        # Compute the density values using the provided functions and append to a new column
        df['potential_temperature'] = df.apply(lambda row: potential_temperature(row['salinity'], row['temperature'], row['pressure']), axis=1)
    else:
        raise ValueError("The required columns 'salinity', 'temperature', 'pressure' must be present in the CSV file.")
    
    return df



def compute_density_anomaly_from_df(df): # here T means potential Temperature 
    df['density_anomaly'] = df.apply(lambda row: density_anomaly(row['salinity'], row['temperature']), axis=1)
    return df



def compute_potential_density_and_E_from_df(df): # here T means potential Temperature 
    df['potential_density'] = df.apply(lambda row: density_anomaly(row['salinity'], row['potential_temperature']), axis=1)
    
    E = []
    potential_density = df['potential_density']
    depth = df['pressure']
    density = df['density']

    for i in range(len(depth)-1):
        delta_rho = potential_density[i] - potential_density[i+1]
        delta_z = (-1*depth[i]) - (-1*depth[i+1])
        avg_density = 0.5*(density[i] + density[i+1])
        e_val = - (1 / avg_density) * (delta_rho / delta_z)
        E.append(e_val)
    E.append(np.nan)
    df['E'] = E
    
    return df



def intermediate_plot(df):
    # Set Seaborn style for more appealing visuals
    sns.set_style("whitegrid", {'grid.linestyle': '--'})

    # Increase the size of all labels and titles
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    
    # Create the plot with increased height
    fig, ax1 = plt.subplots(figsize=(10, 11))

    # Plot Temperature on x-axis with enhanced visuals
    color = sns.color_palette("colorblind")[0]
    ax1.set_xlabel('Temperature (°C)', color=color, weight='bold')
    ax1.plot(df['temperature'], df['pressure'], color=color, linewidth=3)
    ax1.tick_params(axis='x', labelcolor=color)
    ax1.invert_yaxis()  # Invert y axis so depth increases downwards
    
    # Enhanced Potential Temperature plot
    color = sns.color_palette("colorblind")[1]
    ax2 = ax1.twiny()
    ax2.xaxis.tick_bottom()
    ax2.xaxis.set_label_position('bottom') 
    ax2.spines['bottom'].set_position(('outward', 60))
    ax2.set_xlabel('Potential Temperature (°C)', color=color, weight='bold')
    ax2.plot(df['potential_temperature'], df['pressure'], color=color, linewidth=3)
    ax2.tick_params(axis='x', labelcolor=color, direction='out', which='both')
    
    # Enhanced Salinity plot
    color = sns.color_palette("colorblind")[2]
    ax3 = ax1.twiny()
    ax3.spines['top'].set_position(('outward', 0))
    ax3.set_xlabel('Salinity (psu)', color=color, weight='bold')
    ax3.plot(df['salinity'], df['pressure'], color=color, linewidth=3)
    ax3.tick_params(axis='x', labelcolor=color)

    # Enhanced Density plot
    color = sns.color_palette("colorblind")[3]
    ax4 = ax1.twiny()
    ax4.spines['top'].set_position(('outward', 60))
    ax4.set_xlabel('Sigma-t (kg/m³)', color=color, weight='bold')
    ax4.plot(df['density_anomaly'], df['pressure'], color=color, linestyle='--', linewidth=3.5)
    ax4.tick_params(axis='x', labelcolor=color)
    
    # Enhanced Potential Density plot
    color = sns.color_palette("colorblind")[7]
    ax5 = ax1.twiny()
    ax5.spines['top'].set_position(('outward', 120))
    ax5.set_xlabel('Potential density (kg/m³)', color=color, weight='bold')
    ax5.plot(df['potential_density'], df['pressure'], color=color, linewidth=3)
    ax5.tick_params(axis='x', labelcolor=color)

    # Grid and title setup
    ax1.grid(True, which='both', linestyle='--', linewidth=1)
    plt.title('Profiles of Temperature, Potential Temperature, Salinity, Sigma-t and Potential Density against Depth', y=1.4, fontsize=18, fontweight='bold')
    
    # Display the enhanced plot
    plt.tight_layout()
    plt.show()

    
def final_plot(df):
    depth = df['pressure']
    T = df['temperature']
    S = df['salinity']
    Potential_temperature = df['potential_temperature']
    Sigma_tee = df['density_anomaly']
    potential_density = df['potential_density']
    density = df['density']
    E = df["E"]

    # Set global font size and style for all plots
    plt.rcParams['axes.labelsize'] = 16  # Increase label font size
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 22
    plt.rcParams['legend.fontsize'] = 18

    # Plotting: 2 rows and 4 columns
    fig, axs = plt.subplots(2, 4, figsize=(20, 15))

    # Plot Temperature
    # For the first plot
    axs[0, 0].plot(T, -depth, label='Temperature', color='b', linewidth=2.5)
    axs[0, 0].set_xlabel('Temperature')
    axs[0, 0].set_ylabel('Depth (m)')
    axs[0, 0].legend()
    
    # For the first plot
    #handles, labels = axs[0, 0].get_legend_handles_labels()
    #labels[0] = r'$\mathbf{Temperature}$'  # Using a bold LaTeX font here
    #axs[0, 0].legend(handles, labels, fontsize=18)  # Set the fontsize to your desired size


    # Plot Potential Temperature
    axs[0, 1].plot(Potential_temperature, -depth, label='Potential temperature', color='r', linewidth=2.5)
    axs[0, 1].set_xlabel('Potential Temperature')
    axs[0, 1].legend()

    # Plot Salinity
    axs[0, 2].plot(S, -depth, label='Salinity', color='g', linewidth=2.5)
    axs[0, 2].set_xlabel('Salinity')
    axs[0, 2].legend()

    # Plot Density
    axs[0, 3].plot(density, -depth, label='Density', color='purple', linewidth=2.5)
    axs[0, 3].set_xlabel('Density')
    axs[0, 3].legend()

    # Plot Density Anomaly
    axs[1, 0].plot(Sigma_tee, -depth, label='Sigma-t', color='cyan', linewidth=2.5)
    axs[1, 0].set_xlabel('Density Anomaly')
    axs[1, 0].set_ylabel('Depth (m)')
    axs[1, 0].legend()

    # Plot Potential Density
    axs[1, 1].plot(potential_density, -depth, label='Potential Density', color='orange', linewidth=2.5)
    axs[1, 1].set_xlabel('Potential Density')
    axs[1, 1].legend()
    
    # Plot static stability
    axs[1, 2].plot(E, -depth, label='E', color='black', linewidth=2.5)
    axs[1, 2].set_xlabel('Static Stability')
    axs[1, 2].legend()
    
    # Turn off the axis for the empty subplot
    axs[1, 3].axis('off')

    # Add title
    fig.suptitle("Profiles of Temperature, Potential Temperature, Salinity, Density, Sigma-t, Potential Density and E against Depth", fontsize=24, fontweight='bold', y=1.08)
    
    plt.tight_layout()
    plt.show()


def all_values_in_one_table(df):
    new_data = compute_density_from_df(df)
    new_data2 = compute_potential_temperature_from_df(new_data)
    new_data3 = compute_density_anomaly_from_df(new_data2)
    new_data4 = compute_potential_density_and_E_from_df(new_data3)
    return new_data4


def app():

    st.title("Seawater Potential Temperature, Density, Potential Density, Density Anomaly, and Static Stability Calculator")

    # Display your name
    st.write("Shammunul Islam")

    # Display your LinkedIn profile as a clickable link
    st.markdown("[LinkedIn Profile](https://www.linkedin.com/in/shammunul/)")
    
    # Display your picture
    st.image('ocean_wave_shorter.gif', use_column_width=True)

    # Embed audio with autoplay and loop using raw HTML
    st.markdown(f'<audio src="{audio_file}" autoplay loop></audio>', unsafe_allow_html=True)
    
    # Inject custom CSS to style file uploader
    st.markdown("""
    <style>
        .stFileUploader > div > label {
            font-size: 20px;
            color: red;
        }
    </style>
    """, unsafe_allow_html=True)

    uploaded_file = st.file_uploader("Choose a CSV file with at least three column names: 'temperature', 'salinity', and 'pressure'")

    # uploaded_file = st.file_uploader("Choose a CSV file with at least three column names: 'temperature', 'salinity', and 'pressure'", type="csv")
    
    if uploaded_file is not None:
        # Read the CSV data
        df = pd.read_csv(uploaded_file)

        # Ensure the DataFrame is not empty
        if df.empty:
            st.error("The uploaded CSV file does not contain any data. Please upload a valid CSV file.")
            return

        # Process the DataFrame using the provided functions
        try:
            result_df = all_values_in_one_table(df)
            
            st.markdown("<span style='color: blue; font-size: 20px'>Please find a new table below with the newly calculated columns of potential temperature, density, sigma-t, potential density and E</span>", unsafe_allow_html=True)
            # st.write("Please find a new table below with the newly calculated columns of potential temperature, density, sigma-t, potential Density and E")
            # Display DataFrames using Streamlit
            st.write(result_df)

            # Download button
            csv = result_df.to_csv(index=False)
            st.download_button(
                label="Download the data as a CSV file",
                data=csv,
                file_name="processed_data.csv",
                mime="text/csv",
            )
            
            # Plotting intermediate and final plots
            st.markdown("<span style='color: blue; font-size: 20px'>Plot of temperature, potential temperature, salinity, sigma-t and potential density against depth in one plot</span>", unsafe_allow_html=True)
            #st.write("Plot of Temperature, Potential Temperature, Salinity, Sigma-t and Potential Density against Depth in one plot")
            intermediate_plot(result_df)
            st.pyplot(plt.gcf())  # Display the current figure

            st.markdown("<span style='color: blue; font-size: 20px'>Plot of temperature, potential temperature, salinity, density, sigma-t, potential density and E against depth in separate plots</span>", unsafe_allow_html=True)
            # st.write("Plot of Temperature, Potential Temperature, Salinity, Density, Sigma-t, Potential Density and E against Depth in separate plots")
            final_plot(result_df)
            st.pyplot(plt.gcf())  # Display the current figure

        except Exception as e:
            st.error(f"An error occurred: {e}")

if __name__ == "__main__":
    app()
