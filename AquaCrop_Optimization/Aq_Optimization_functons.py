# add soil data if there's internet, add flowering gdd

from aquacrop import AquaCropModel, Soil, Crop, InitialWaterContent,IrrigationManagement
from aquacrop.utils import prepare_weather, get_filepath
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
import math
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.stats import pearsonr
plt.rcParams["figure.dpi"] = 300



def run_aquacrop_model( planting_date, harvest_date, smt, soil, wdf, sim_start ='2006/01/01', sim_end = '2022/12/30',
                       Emergence =80,  Maturity = 1700, MaxIrr =12, MaxIrrSeason =  600,MaxRooting =1409,HIstart =880,
                       Senescence = 1400, CCx = .96, WP = 33.7, Kcb =1.05, HI0 = 0.48, a_HI = 7.0, Zmax = 2.3):
    

             
    irr_mngt = IrrigationManagement(irrigation_method=1, SMT=smt, MaxIrrSeason=MaxIrrSeason, MaxIrr= MaxIrr)
        
    
    # Create Crop, InitialWaterContent, and IrrigationManagement instances
#     if crop is None:
      
    crop = Crop(c_name='MaizeGDD', 
                Name='MaizeGDD', planting_date=planting_date, harvest_date=harvest_date,GDDmethod = 3,
               Emergence =Emergence,  Maturity = Maturity, MaxRooting =MaxRooting,
                CCx =CCx, WP = WP, Kcb = Kcb, HI0 = HI0, a_HI = a_HI, Zmax = Zmax)
#     print('cc0', crop.CC0)
    crop.CGC = (
            np.log(
                (((0.98 * crop.CCx) - crop.CCx) * crop.CC0)
                / (-0.25 * (crop.CCx**2))
            )
        ) / (-(705 - crop.Emergence))
    
    
    
    tCD = crop.MaturityCD - crop.SenescenceCD
    if tCD <= 0:
        tCD = 1

    CCi = crop.CCx * (1 - 0.05 * (np.exp(((3.33 * crop.CDC_CD) / (crop.CCx + 2.29)) * tCD) - 1))
    if CCi < 0:
        CCi = 0

    tGDD = crop.Maturity - crop.Senescence
    if tGDD <= 0:
        tGDD = 5

    crop.CDC = ((crop.CCx + 2.29) * np.log((((CCi/crop.CCx) - 1) / -0.05) + 1)) / (3.33 * tGDD)
#     print("maturity", crop.Maturity)
#     print(WP)
#     print("SMT", smt)
#     print("CGC:", crop.CGC) 
#     print("CDC:",crop.CDC)
    initWC = InitialWaterContent(value=['FC'])
    
# MaxIrr= 6.5, MaxIrrSeason=600
    # Create AquaCropModel instance
    model = AquaCropModel(sim_start, sim_end, wdf, soil, crop,
                          initial_water_content=initWC, irrigation_management=irr_mngt, off_season=True)
    # Run the model
    model.run_model(till_termination=True)
    
    # Save results
    field_df = model._outputs.water_flux
    date_range = pd.date_range(start=sim_start, end=sim_end).strftime("%Y-%m-%d").tolist()
    field_df["Date"] = date_range
    field_df["ET_aqua"] = field_df.Es + field_df.Tr
    field_df["Date"] = pd.to_datetime(field_df["Date"])
    field_df.index = field_df["Date"]

    field_yld = model._outputs.final_stats
    gdd_df =model._outputs.crop_growth
    x= pd.Series([crop.Maturity], name='Maturity')
    
    return gdd_df, field_df, field_yld, x




# ---------------------------for All fields -------------------------


def run_aquacrop_model_for_fid(df, fid, sim_start='01/01', sim_end='12/30', smt=[55, 75, 45, 30], Emergence =80,
                               Maturity = 1700, MaxIrr =12,
                               CCx = .96, WP = 33.7, Kcb =1.05, HI0 = 0.48, a_HI = 7.0, Zmax = 2.3 ):

    # Subset the DataFrame for the current fid
    fid_df = df[df['fid'] == fid].copy()

    # Create WaterDataFrame (wdf) based on your data
    wdf = pd.DataFrame(fid_df[['MinTemp', 'MaxTemp', 'Precipitation', 'ReferenceET', 'Date']])
        # Default values for simulation parameters if not provided
   
 # maximum irrigation df

    sam_df =pd.read_csv("Data/FieldData_AllFieldsCompiled-Annual.csv")
    sam_df= sam_df[sam_df['cropType'].str.lower() == "corn"]
    maxirrdf = sam_df.groupby("FieldID")["irrigation_mm"].max().reset_index()   

    if fid in maxirrdf['FieldID'].values:
        # Get the corresponding max irrigation season
        max_irr_season = maxirrdf[maxirrdf['FieldID'] == fid]['irrigation_mm'].values[0]
    else:
        # Set default max irrigation season to 550
        max_irr_season = 600
    
    irr_mngt = IrrigationManagement(irrigation_method=1, SMT=smt, MaxIrrSeason=max_irr_season, MaxIrr= MaxIrr)
    

#     add gdd df
    gdd_df =pd.read_csv("Data/all_fields_GDD.csv") 

    default_maturity = 1700
    default_max_rooting = 1409
    default_hI_start = 880
    default_senescence = 1400
    soil = Soil(soil_type='SiltLoam')
    
    # Define planting and harvest dates based on FieldID
    if fid.startswith("SW"):
        
        Emergence = gdd_df.loc[gdd_df["FieldID"]=="SW", "Emergence" ].item()
        Maturity =  gdd_df.loc[gdd_df["FieldID"]=="SW", "Maturity" ].item()
        planting_date = '05/10'
        harvest_date = '10/07'
        soil = Soil(soil_type = 'SiltClayLoam')
        
    elif fid.startswith("NB") or fid.startswith("NW"):
        
        Emergence = gdd_df.loc[gdd_df["FieldID"]=="NW", "Emergence" ].item()
        Maturity =  gdd_df.loc[gdd_df["FieldID"]=="NW", "Maturity" ].item()
        
        planting_date = '05/18'
        harvest_date = '10/17'
        
    elif fid.startswith("NC"):
        
        Emergence = gdd_df.loc[gdd_df["FieldID"]=="NC", "Emergence" ].item()
        Maturity =  gdd_df.loc[gdd_df["FieldID"]=="NC", "Maturity" ].item()
        
        planting_date = '05/07'
        harvest_date = '10/12'
        
    elif fid.startswith("WC"):
        Emergence = gdd_df.loc[gdd_df["FieldID"]=="WC", "Emergence" ].item()
        Maturity =  gdd_df.loc[gdd_df["FieldID"]=="WC", "Maturity" ].item()
        
        planting_date = '05/16'
        harvest_date = '10/11'
        
    else:
        # Default planting and harvest dates if FieldID does not match any condition
        planting_date = '05/08'
        harvest_date = '10/07'
        
    percentage_change = (Maturity  - default_maturity) / default_maturity * 100
    MaxRooting = default_max_rooting + (percentage_change / 100) * default_max_rooting
#     print('maxrooting:',  MaxRooting)
    HIstart = default_hI_start + (percentage_change / 100) * default_hI_start
    Senescence = default_senescence + (percentage_change / 100) * default_senescence

    # Run AquaCrop model for the current fid
    gdd_df, field_df, field_yld, x = run_aquacrop_model(planting_date, harvest_date,smt =smt, soil =soil, wdf=wdf,                                                                     sim_start = sim_start, sim_end = sim_end, Emergence = Emergence,
                                                       Maturity = Maturity, MaxIrrSeason = max_irr_season,
                                                       MaxRooting = MaxRooting, HIstart = HIstart, Senescence = Senescence,
                                                       CCx =CCx, WP = WP, Kcb = Kcb, HI0 = HI0, a_HI = a_HI, Zmax = Zmax 
                                                       )
    # Add fid column to identify the data
    field_df['FieldID'] = fid
    field_yld['FieldID'] = fid
    gdd_df['FieldID'] = fid
    x['FieldID'] = fid
    
    return gdd_df, field_df, field_yld, x



                          # ********************************************
def for_objf (smt =[50, 65, 40, 30], CCx = 0.96, WP = 33.7, Kcb =1.05, HI0 = 0.48, 
              a_HI = 7.0, Zmax = 2.3,df_type ="train", no_et = True ):

    # run all fields aquacrop model simulation
    df = pd.read_csv('Data/Updated_all_fieldsclimate.csv')
    df["Date"] = df.Date.str[:8]
    df["Date"] = pd.to_datetime(df["Date"])
    df['Year'] = df['Date'].dt.year
#     df =df[~df['fid'].isin(['NW6', 'NW7'])]
#     df = df[~df["fid"].isin(["NW6", "NW7"]) & ~df["fid"].str.startswith("WC")]


    unique_fids = df['fid'].unique()
    # DataFrames to store results for all fids
    all_ET_df = pd.DataFrame()
    all_yld_df = pd.DataFrame()

    # Iterate over each fid and run AquaCrop model
    for fid in unique_fids:
        model, field_df, field_yld, x = run_aquacrop_model_for_fid(df, fid, 
                                                                   sim_start='2006/01/01', 
                                                                   sim_end='2023/12/30', 
                                                                   smt = smt , CCx = CCx, WP = WP, Kcb = Kcb,
                                                                   HI0 = HI0, a_HI = a_HI, Zmax = Zmax
                                                                   )

        # Append results to the aggregated DataFrames
        all_ET_df = pd.concat([all_ET_df, field_df], ignore_index=True)
        all_yld_df = pd.concat([all_yld_df, field_yld], ignore_index=True)

    all_yld_df["Year"] = all_yld_df["Harvest Date (YYYY/MM/DD)"].dt.year
    all_ET_df["Year"] = all_ET_df.Date.dt.year
    

    
#     

    #for train run
    
    if df_type == "test":
        ydff  = pd.read_csv("/Users/m089r172/Library/CloudStorage/OneDrive-UniversityofKansas/old_codes_And_new_mac_codes/Python_docs/Aquacrop/Data_for_Crop_ML_model/test_data.csv")

    elif df_type == "train":
        ydff  = pd.read_csv("/Users/m089r172/Library/CloudStorage/OneDrive-UniversityofKansas/old_codes_And_new_mac_codes/Python_docs/Aquacrop/Data_for_Crop_ML_model/train_data.csv")
    else:
        ydf= pd.read_csv ("/Users/m089r172/Library/CloudStorage/OneDrive-UniversityofKansas/AquaCrop documents/All_reported_data.csv")

        bad_fields = pd.read_csv("../Data_for_Crop_ML_model/bad_fields.csv")
        ydf = ydf.merge(bad_fields, on=["Year", "FieldID"], how="left", indicator=True)
        ydff = ydf[ydf["_merge"] == "left_only"].drop(columns=["_merge"])# full data filtered should be here
              
        
#     ydff = ydff[~((ydff['FieldID'].isin(['NW3', 'NW5'])) & (ydff['Year'] == 2022))]
    print(len(ydff))
        
    if no_et:
        simul_reported = pd.merge(all_yld_df, ydff, on=['FieldID', 'Year'], how='inner')
        print("without_ET:",len(simul_reported ))
        test_df = "Not available as we are using full data"

            
    else:
        #     OpenET data
        GEEdf =pd.read_csv("AquaCrop_Optimization/Data/Full_year_ET2.csv")
        #         Yearly sum
        sim_et= all_ET_df.groupby(["FieldID" ,"Year"])["ET_aqua"].sum().reset_index()
        all_ET_df =pd.merge(GEEdf,  sim_et, on=['FieldID', 'Year'], how='inner')
        simul_reported =  pd.merge(ydff,all_yld_df, on=['FieldID', 'Year'], how='inner') 
        simul_reported  = pd.merge ( all_ET_df ,  simul_reported, on=['FieldID', 'Year'], how='inner')
        print("100 prcnt_with_ET:",len(simul_reported ))



    
    return all_yld_df, simul_reported, all_ET_df 





def obj_func (param, no_et =True):
    
    """
    you can remove any parameter from "for_objf" function
    to avaoid calibration of that parameter. for that lb and up adjustment is necessary.
    """

    smt = [param[0], param[1], param[2], param[3]]
    CCx = param[4]
    WP = param[5] 
    Kcb = param[6]
    HI0 = param[7]
    a_HI = param[8]
    Zmax = param[7]

# , a_HI = a_HI, Zmax = Zmax
    print("smt:", smt)
    print("CCx:", CCx)
    print("WP:",  WP)
    print("Kcb:", Kcb)
    print("HI0:", HI0)
    print("a_HI:", a_HI)
    print("Zmax:",Zmax)
    
    all_yld_df,  df = for_objf(smt =smt, CCx =CCx, WP = WP, Kcb = Kcb,
                                                     HI0 = HI0, a_HI = a_HI, Zmax = Zmax,
                                                     no_et = no_et)
    
    df.rename(columns={"Seasonal irrigation (mm)": "Simulated_Irrigation",
                                        "Yield (tonne/ha)": "Simulated_Yield"},inplace =True)
    # Step 1: Calculate residuals
    df["Irrigation_Residual"] = (df["Reported_Irrigation"] - df["Simulated_Irrigation"]).abs()
    df["Yield_Residual"] = (df["Reported_Yield"] - df["Simulated_Yield"]).abs()



    # Step 2: Normalize residuals using Min-Max scaling
    def min_max_scaling(series):
        return (series - series.min()) / (series.max() - series.min())

    df["Normalized_Irrigation_Residual"] = min_max_scaling(df["Irrigation_Residual"])
    df["Normalized_Yield_Residual"] = min_max_scaling(df["Yield_Residual"])


    # Step 3: Compute mean of normalized residuals
    mean_normalized_irrigation_residual = df["Normalized_Irrigation_Residual"].mean()
    mean_normalized_yield_residual = df["Normalized_Yield_Residual"].mean()

    if no_et:
        weight_yield = 6
        weight_irrigation =4


        # Objective function
        fitness = (
            weight_irrigation * mean_normalized_irrigation_residual +
            weight_yield * mean_normalized_yield_residual
        )

        # Output results
#         print("Mean Normalized Irrigation Residual:", mean_normalized_irrigation_residual)
#         print("Mean Normalized Yield Residual:", mean_normalized_yield_residual)
        print("Irrigation contribution:",weight_irrigation * mean_normalized_irrigation_residual)
        print("Yield contribution:",  weight_yield *mean_normalized_yield_residual)
#         print("Objective Function Value:", fitness)

    else:
        df["ET_Residual"] =  (df['Ensemble_ET'] - df['ET_aqua']).abs()
        df["Normalized_ET_Residual"] = min_max_scaling(df["ET_Residual"])
        mean_normalized_ET_residual = df["Normalized_ET_Residual"].mean()
        
        weight_yield = 5
        weight_irrigation =3
        weight_et  = 2

        # Objective function
        fitness = (
            weight_irrigation * mean_normalized_irrigation_residual +
            weight_yield * mean_normalized_yield_residual +
            weight_et * mean_normalized_ET_residual
        )

        # Output results
#         print("Mean Normalized Irrigation Residual:", mean_normalized_irrigation_residual)
#         print("Mean Normalized Yield Residual:", mean_normalized_yield_residual)
#         print("ET_resid:", mean_normalized_ET_residual)
        print("Irrigation contribution:",weight_irrigation * mean_normalized_irrigation_residual)
        print("Yield contribution:",  weight_yield *mean_normalized_yield_residual)
        print("ET contribution:",  weight_et * mean_normalized_ET_residual)
#         print("Objective Function Value:", fitness)

    print("loss:", fitness)
    print("------------------------")
    
    return fitness





# ------------------ comparison plot ///////////////////////////


def pbias(predicted, observed):
    """Percent bias"""
    return 100.0 * np.sum(observed - predicted) / np.sum(observed)

def nmae(predicted, observed):
    """Normalized Mean Absolute Error"""
    return 100* mean_absolute_error(observed, predicted) / np.mean(observed)


def plot_yield_comparison(
    df, 
    x_col: str = 'Reported_Yield', 
    y_col: str = 'Simulated_Yield', 
    df_type: str = "Test", 
    hue: str = None,
    style: str = None,
    metric: bool = None,
    x_axis: str = None,
    y_axis: str = None,
    xlim: tuple = (-200, 1100), 
    ylim: tuple = (-200, 1100),
    bound: tuple = None 
    # New parameter for bounding region
) -> tuple:
    plt.rcParams["figure.dpi"] = 300
    fig, axes = plt.subplots(figsize=(8, 5))

    # 1:1 line
    axes.axline((0, 0), slope=1, label='1:1 line', color='gray', linestyle='--')
    
    # Scatter plot
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue, ax=axes, palette='tab10', style=style, s=50, edgecolor='k')
    
    # Regression line
    sns.regplot(data=df, x=x_col, y=y_col, scatter=False, ci=None, ax=axes, color='red', label='Regression line', line_kws={'label': None})

    # Add shaded region if `bound` is specified
    if bound:
        # Create a continuous range of x values spanning the x-axis limits
        x_range = np.linspace(df[x_col].min(), df[x_col].max(), len(df))
        
        # Fill the area between the bounds for the entire x-range
        axes.fill_between(x_range, bound[0], bound[1], color='gray', alpha=0.3, label=f'Y range: {bound[0]} to {bound[1]}')

    # Set axis labels if provided
    if x_axis is not None:
        axes.set_xlabel(x_axis, fontsize=12)
    else:
        axes.set_xlabel(x_col, fontsize=12)  # Default to x_col if x_axis is not provided

    if y_axis is not None:
        axes.set_ylabel(y_axis, fontsize=12)
    else:
        axes.set_ylabel(y_col, fontsize=12)  # Default to y_col if y_axis is not provided
        
    # Set axis limits if provided
    if xlim:
        axes.set_xlim(xlim)
    if ylim:
        axes.set_ylim(ylim)

    # Set limits
    if x_col == "Reported_Yield":
        axes.set_ylim(0, 19)
        axes.set_xlim(0, 19)
    elif x_col == "Ensemble_ET":
        axes.set_ylim(600, 922)
        axes.set_xlim(600, 922)
    else:
        axes.set_ylim(140, 756)
        axes.set_xlim(140, 756)

    # Calculate metrics
    observed = df[x_col].values
    predicted = df[y_col].values

    corr, _ = pearsonr(observed, predicted)
    mae = mean_absolute_error(observed, predicted)
    rmse = np.sqrt(mean_squared_error(observed, predicted))
    r2 = r2_score(observed, predicted)
    pbias_value = pbias(predicted, observed)
    nmae_value = nmae(predicted, observed)

    # Annotate metrics
#     axes.annotate(f"MAE: {mae:.2f}", xy=(0.05, 0.95), xycoords='axes fraction', fontsize=11)
    axes.annotate(f"RMSE: {rmse:.2f}", xy=(0.05, 0.90), xycoords='axes fraction', fontsize=11)
    axes.annotate(f"PBIAS: {pbias_value:.2f}%", xy=(0.05, 0.85), xycoords='axes fraction', fontsize=11)
    axes.annotate(f"NMAE: {nmae_value:.2f}%", xy=(0.05, 0.80), xycoords='axes fraction', fontsize=11)
    axes.annotate(f"$R^2$: {corr**2:.2f}", xy=(0.05, 0.75), xycoords='axes fraction', fontsize=11)
#     axes.annotate(f"$r$: {corr:.2f}", xy=(0.05, 0.70), xycoords='axes fraction', fontsize=11)

    # Legend & title
    axes.legend( loc='upper right')
    axes.set_title(f"{df_type}")

    return fig, axes



# --------------------------------------------------
def calculate_single_metric(df, observed_col, simulated_col):
    """
    Calculate MAE, RMSE, and R^2 for a single pair of observed and simulated columns,
    while handling NaN values.
    
    Parameters:
    - df: DataFrame containing the data.
    - observed_col: Column name for observed values.
    - simulated_col: Column name for simulated values.
    
    Returns:
    - Dictionary containing MAE, RMSE, and R^2.
    """
    # Drop NaN values for the current pair of columns
    filtered_df = df[[observed_col, simulated_col]].dropna()

    mae = mean_absolute_error(filtered_df[observed_col], filtered_df[simulated_col])
    rmse = np.sqrt(mean_squared_error(filtered_df[observed_col], filtered_df[simulated_col]))
    r2 = r2_score(filtered_df[observed_col], filtered_df[simulated_col])

    return {"MAE": mae, "RMSE": rmse, "R^2": r2}

def calculate_metrics_to_dataframe(df,d_set ="train",prm ="prms", iteration =1,no_et =True):
    
    df.rename(columns={"Seasonal irrigation (mm)": "Simulated_Irrigation",
                                        "Yield (tonne/ha)": "Simulated_Yield"},inplace =True)
    # Initialize a list to hold metric rows
    metric_rows = []

    # Define pairs of observed and simulated columns
    pairs = [
        ("Reported_Irrigation", "Simulated_Irrigation"),
        ("Reported_Yield", "Simulated_Yield"),
        ("Ensemble_ET", "ET_aqua")
    ]
    if no_et:
            pairs = [
        ("Reported_Irrigation", "Simulated_Irrigation"),
        ("Reported_Yield", "Simulated_Yield")
    ]
            
    for obs, sim in pairs:
        # Drop NaN values specific to the pair
        filtered_df = df[[obs, sim]].dropna()
        
        # Calculate metrics
        mae = mean_absolute_error(filtered_df[obs], filtered_df[sim])
        rmse = np.sqrt(mean_squared_error(filtered_df[obs], filtered_df[sim]))
        r2 = r2_score(filtered_df[obs], filtered_df[sim] )
        
        # Calculate Pearson correlation coefficient
        pearson_corr, _ = pearsonr(filtered_df[obs], filtered_df[sim])
        
        # Append the results as a dictionary
        metric_rows.append({
            "iteration": iteration,
            "obj_variable": obs.split("_")[1],
            "MAE": mae,
            "RMSE": rmse,
            "R^2": r2,
            "Pearson r": pearson_corr,
             "Data_set": d_set,
             "prms":prm
        })
    
    # Convert the list of dictionaries to a DataFrame
    ms = pd.DataFrame(metric_rows)
    return ms




def plot_comparison(df, reported_col="Reported_Yield", simulated_col="Simulated_Yield", ylim=(0, 19),
                    cols=2, loc='upper center'):
    """
    Plots a comparison of simulated and reported values for each unique field,
    and displays RMSE, PBIAS, NMAE, and R² inside each subplot.

    Parameters:
        df (pd.DataFrame): Dataset containing fields and comparison data.
        reported_col (str): Column name for reported values.
        simulated_col (str): Column name for simulated values.
        ylim (tuple): Y-axis limits. Default is (0, 19).
        cols (int): Number of columns for subplot grid. Default is 2.
    """
    df = df[df.Year < 2023]
    plt.rcParams["figure.dpi"] = 300
    unique_fields = df['FieldID'].nunique()
    rows = math.ceil(unique_fields / cols)
    fig = plt.figure(figsize=(14, rows * 3))

    for i, (field_id, group) in enumerate(df.sort_values('Year').groupby('FieldID')):
        group_no_nan = group.dropna(subset=[reported_col, simulated_col])

        observed = group_no_nan[reported_col].values
        predicted = group_no_nan[simulated_col].values

        # Metrics
        rmse = mean_squared_error(observed, predicted, squared=False)
        pbias_val = pbias(predicted, observed)
        nmae_val = nmae(predicted, observed)
#         r = np.corrcoef(observed, predicted)[0, 1]
#         r2 = r ** 2

        ax = fig.add_subplot(rows, cols, i + 1)

        # Plot simulated and reported values
        ax.plot(group['Year'], group[simulated_col], label='Simulated', marker='o')
        ax.plot(group['Year'], group[reported_col], label='Reported', marker='o')

        # Title
        ax.set_title(f'FieldID {field_id}', fontsize=11)

        # Add metrics inside plot
        ax.annotate(
            f"RMSE: {rmse:.2f}\n"
            f"PBIAS: {pbias_val:.2f}%\n"
            f"NMAE: {nmae_val:.2f}%\n",
#             f"$R^2$: {r2:.2f}",
            xy=(0.05, 0.95), xycoords='axes fraction',
            fontsize=9, verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white")
        )

        if reported_col == "Reported_Yield":
            ax.set_xlabel('Year', fontsize=12)
            ax.set_ylabel('Yield (tonne/ha)', fontsize=12)
            ax.set_ylim(*ylim)
        else:
            ax.set_xlabel('Year', fontsize=12)
            ax.set_ylabel('Irrigation (mm)', fontsize=12)
            ax.set_ylim(90, 610)

        ax.legend(loc=loc)

    plt.tight_layout()
    plt.show()




# def plot_comparison(df, reported_col="Reported_Yield", simulated_col="Simulated_Yield", ylim=(0, 19), cols=2):
#     """
#     Plots a comparison of simulated and reported values for each unique field, 
#     including correlation and RMSE in the title.

#     Parameters:
#         df (pd.DataFrame): The dataset containing the fields and comparison data.
#         reported_col (str): The column name for reported values.
#         simulated_col (str): The column name for simulated values.
#         ylim (tuple): Limits for the y-axis. Default is (0, 17).
#         cols (int): The number of columns for the subplot grid. Default is 2.
#     """
#     df = df[df.Year<2023]
#     plt.rcParams["figure.dpi"] = 300
#     unique_fields = df['FieldID'].nunique()  # Count unique fields
#     rows = math.ceil(unique_fields / cols)  # Calculate rows dynamically
# #     df["Year"] = pd.to_datetime(df["Year"])
#     fig = plt.figure(figsize=(14, rows * 3))  # Adjust figure size

#     # Iterate over each FieldID group
#     for i, (field_id, group) in enumerate(df.sort_values('Year').groupby('FieldID')):
#         # Calculate correlation
#         group_no_nan = group.dropna(subset =['Reported_Yield'])
#         correlation = group_no_nan[simulated_col].corr(group_no_nan[reported_col])
#         # Calculate RMSE
# #         rmse = mean_squared_error(group[reported_col], group[simulated_col], squared=False)
#         ax = fig.add_subplot(rows, cols, i + 1)  # Dynamically determine position
#         # Plot simulated and reported values
#         ax.plot(group['Year'], group[simulated_col], label='Simulated', marker='o')
#         ax.plot(group['Year'], group[reported_col], label='Reported', marker='o')

#         # Add title with correlation and RMSE
#         ax.set_title(f'FieldID {field_id}\nCorrelation: {correlation:.2f}')
# #         ax.set_title(f'FieldID {field_id}\nCorrelation: {correlation:.2f}, RMSE: {rmse:.2f}')
        
#         if reported_col == "Reported_Yield":
#             ax.set_xlabel('Year', fontsize=14, fontweight='bold')
#             ax.set_ylabel('Yield (tonne/ha)', fontsize=14, fontweight='bold')
#             ax.set_ylim(*ylim)
#             ax.legend()
#         else:
#             ax.set_xlabel('Year', fontsize=14, fontweight='bold')
#             ax.set_ylabel('Irrigation (mm)', fontsize=14, fontweight='bold')
#             ax.set_ylim(90,610)
#             ax.legend()

#     plt.tight_layout()
#     plt.show()

