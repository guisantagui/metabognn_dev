################################################################################
# MetaboPert: align averaged delta dataset to Schuhknecht metabolomics dataset.#
################################################################################

import pandas as pd
import numpy as np
import os
import re
import copy as cp
import argparse

# Functions
################################################################################

COL_RE = re.compile(r"^(.*)_([\d.]+) uM_(\d+) h$")

def parse_metchange_columns(metchange_df):
    parsed = {}
    for col in metchange_df.columns:
        m = COL_RE.match(col)
        if m:
            compound, dose, time = m.groups()
            parsed[col] = (compound, float(dose), float(time))
    return (parsed)

def available_timepoints(parsed_cols):
    return (
        sorted(
            set(t for _, _, t in parsed_cols.values())
        )
    )
    
def select_best_column(
    compound,
    parsed_cols,
    target_dose=10.0,
    target_time=48.0,
    max_log_dose_dist=np.log10(3)
):
    candidates = [(col, d, t) for col, (c, d, t) in parsed_cols.items() if c == compound]
    if not candidates:
        return None, None, None, "no_data_for_compound"
    times_here = sorted(set(t for _, _, t in candidates))
    chosen_time = target_time if target_time in times_here else \
        min(times_here, key=lambda t: abs(t - target_time))
    time_dist = abs(chosen_time - target_time)
    same_time = [(col, d) for col, d, t in candidates if t == chosen_time]
    best_col, best_dose, best_log_dist = None, None, np.inf
    for col, d in same_time:
        ld = abs(np.log10(d) - np.log10(target_dose))
        if ld < best_log_dist:
            best_col, best_dose, best_log_dist = col, d, ld
    if best_log_dist > max_log_dose_dist:
        return None, best_log_dist, time_dist, "dose_too_far"
    status = "exact" if (best_log_dist == 0 and time_dist == 0) else \
             "nearest_dose" if time_dist == 0 else \
             "nearest_time" if best_log_dist == 0 else \
             "nearest_time_and_dose"
    return best_col, best_log_dist, time_dist, status

def build_match_table(
    metchange_df,
    compounds,
    target_dose=10.0,
    target_time=48.0,
    max_log_dose_dist=np.log10(3)
):
    parsed = parse_metchange_columns(metchange_df)
    print("Timepoints available in MetChange data:", available_timepoints(parsed))
    rows = []
    for c in compounds:
        col, dose_dist, time_dist, status = select_best_column(
            c, parsed, target_dose, target_time, max_log_dose_dist
        )
        rows.append(
            {
                "compound": c,
                "matched_column": col,
                "log_dose_distance": dose_dist,
                "time_distance_h": time_dist,
                "match_status": status
            }
        )
        match_df = pd.DataFrame(rows)
    print(match_df["match_status"].value_counts())
    return (match_df)

# Terminal argument parser
################################################################################
parser = argparse.ArgumentParser(
    description='Align averaged delta dataset to Schuhknecht metabolomics dataset.'
)

parser.add_argument(
        '--mc_delta_avg',
        type=str,
        help='MetChange average deltas file.'
)

parser.add_argument(
        '--mc_mdat_avg',
        type=str,
        help='Metadata file of the MetChange average delta file.'
)

parser.add_argument(
        '--drug_map',
        type=str,
        help='Drug equivalence table, obtained with SMILES.'
)

parser.add_argument(
        '--cpd_map',
        type=str,
        help='Compound map, obtained with GEM table.'
)

parser.add_argument(
        '--sch',
        type=str,
        help='Excel file from Schuhknecht paper.'
)

parser.add_argument(
        '--compartments',
        type=str,
        help='Compartments whose metabolites should be retained. Either one, or several, separated by commas.'
)

parser.add_argument(
        '--out_dir', 
	type=str,
        help='Directory where files will be stored.', 
        default=os.getcwd()
)

args = parser.parse_args()


# Directory stuff
################################################################################
#mc_delta_avg_f = "/home/users/gsantamaria/projects/metabo_pert/results/mc_deltas/A549/mixPnlty/zFPKM_penalty_gem_filt_metChange_rpmi_medium_concat_delta_avg.csv"
#mc_mdat_avg_f = "/home/users/gsantamaria/projects/metabo_pert/results/mc_deltas/A549/mixPnlty/mdat_penalties_delta_avg.csv"
#drug_map_f = "/home/users/gsantamaria/projects/metabo_pert/results/lincs_to_schuhknecht_mapping/common_drugs_mapping.csv"
#cpd_map_f = "/home/users/gsantamaria/projects/metabo_pert/results/gem_cpd_mapping/cpd_info_wSchuhknecht.csv"
#sch_f = "/home/users/gsantamaria/projects/metabo_pert/data/schuhknecht/41587_2024_2524_MOESM2_ESM.xlsx"
#compartment = "c"
#out_dir = "/home/users/gsantamaria/projects/metabo_pert/results/mc_schuhk_assoc/mixPenlty/"

mc_delta_avg_f = args.mc_delta_avg
mc_mdat_avg_f = args.mc_mdat_avg
drug_map_f = args.drug_map
cpd_map_f = args.cpd_map
sch_f = args.sch
compartment = args.compartments
out_dir = args.out_dir

os.makedirs(out_dir, exist_ok=True)

compartment = compartment.split(',')
# Load data
################################################################################

# Schuhknecht
sch_xlsx = pd.ExcelFile(sch_f)

sch_meanlog2FC = sch_xlsx.parse("mean log2-FC")
sch_pvals = sch_xlsx.parse("p-values")
sch_qvals = sch_xlsx.parse("q-values")
sch_z = sch_xlsx.parse("z-scores")

# Metchange
mc_delta_avg = pd.read_csv(mc_delta_avg_f, index_col = 0)
mc_mdat_avg = pd.read_csv(mc_mdat_avg_f)

mc_delta_avg.shape
mc_mdat_avg.shape

# Drug map
drug_map = pd.read_csv(drug_map_f, index_col = 0)

# Compound (metabolite) map
cpd_map = pd.read_csv(cpd_map_f, index_col = 0)

# Drug matching
################################################################################

# Filter mdat and delta dataframes to keep only perturbations that appear in
# Schuhknecht dataset.
keep_perts = [x in drug_map["cmap_name"].to_list() for x in mc_mdat_avg["pertname"]]
mc_mdat_avg = mc_mdat_avg.loc[keep_perts, :]
mc_delta_avg = mc_delta_avg.loc[:, keep_perts]

# Change index names and remove unnecessary columns
sch_meanlog2FC.index = sch_meanlog2FC["metabolites/drug drug m/z"].to_list()
sch_meanlog2FC = sch_meanlog2FC.iloc[:, 2:sch_meanlog2FC.shape[1]]

sch_pvals.index = sch_pvals["metabolites/drug drug m/z"].to_list()
sch_pvals = sch_pvals.iloc[:, 2:sch_pvals.shape[1]]

sch_qvals.index = sch_qvals["metabolites/drug drug m/z"].to_list()
sch_qvals = sch_qvals.iloc[:, 2:sch_qvals.shape[1]]

# Get ordered Schuhknecht names
mc_mdat_sch_names = mc_mdat_avg["pertname"].map(
    drug_map.set_index("cmap_name")["schuhknecht_name"]
).to_list()

keep_drugs_schuhknecht = [x in mc_mdat_sch_names for x in sch_meanlog2FC.columns]

sch_meanlog2FC = sch_meanlog2FC.loc[:, keep_drugs_schuhknecht]
sch_pvals = sch_pvals.loc[:, keep_drugs_schuhknecht]
sch_qvals = sch_qvals.loc[:, keep_drugs_schuhknecht]


# Metabolite matching
################################################################################

# Filter compartment
cpd_map = cpd_map.loc[[x in compartment for x in cpd_map["compartment"].to_list()], :]

# Keep only compounds that only have one match to Schuhknecht
cpd_map = cpd_map.loc[cpd_map["dup_after_shuhknecht_map"] == False, :]

# Keep only compounds that are in metchange file (some compounds are not there
# due to having no production in LP1).
cpd_map = cpd_map.loc[
    [x in mc_delta_avg.index for x in cpd_map["id"].to_list()]
]

mc_delta_avg = mc_delta_avg.loc[cpd_map["id"], :]

# Matching to best timepoint and best dose
################################################################################
mc_mdat_avg["dose_float"] = [float(x.split()[0]) for x in mc_mdat_avg["dose"]]
mc_mdat_avg["timpnt_int"] = [
    int(x.split()[0]) for x in mc_mdat_avg["timepoint"]
]

sch_meanlog2FC = sch_meanlog2FC.loc[cpd_map["Schuhknecht_mapping"], :]
sch_pvals = sch_pvals.loc[cpd_map["Schuhknecht_mapping"], :]
sch_qvals = sch_qvals.loc[cpd_map["Schuhknecht_mapping"], :]

# Assign variables for Schuhknecht dosage and timepoint (from the paper)

# Dosage in uM
sch_dosage = 10
# Timepoint in h
sch_timpnt = 48

cpds = drug_map.set_index("schuhknecht_name").loc[
    sch_meanlog2FC.columns
]["cmap_name"].to_list()

drug_match_table = build_match_table(mc_delta_avg, cpds)

drug_match_table["schuhknecht_name"] = sch_meanlog2FC.columns

# Keep only exact matches and nearest time.
zip_drugmatch = zip(
    drug_match_table["matched_column"].to_list(),
    drug_match_table["match_status"].to_list()
)

match_stat_keep = ["exact", "nearest_time"]

keep_cols = [
    col != None and stat in match_stat_keep for col, stat in zip_drugmatch
]

drug_match_table = drug_match_table.loc[keep_cols, :]

sch_meanlog2FC = sch_meanlog2FC.loc[
    :, drug_match_table["schuhknecht_name"].to_list()
]

sch_pvals = sch_pvals.loc[
    :, drug_match_table["schuhknecht_name"].to_list()
]

sch_qvals = sch_qvals.loc[
    :, drug_match_table["schuhknecht_name"].to_list()
]

mc_delta_avg = mc_delta_avg.loc[
    :, drug_match_table["matched_column"].to_list()
]

# Homogenize drug names
mc_delta_avg.columns = sch_meanlog2FC.columns

# Sanity checks
print(
    "Correct metabolite name in Schuhknect: %s"
    % (sch_meanlog2FC.index.to_list() == cpd_map["Schuhknecht_mapping"].to_list())
)

print(
    "Correct metabolite name in MetChange: %s"
    % (mc_delta_avg.index.to_list() == cpd_map["id"].to_list())
)

# Homogenize metabilite names
mc_delta_avg.index = cpd_map["name"].to_list()
sch_meanlog2FC.index = cpd_map["name"].to_list()
sch_pvals.index = cpd_map["name"].to_list()
sch_qvals.index = cpd_map["name"].to_list()

# Add colum of whether the drugs in the match table are validated metabolic
# drugs

validated_metabolic = [
    # --------------------------------------------------------------------------
    # 1. ANTIDIABETICS / GLUCOSE-INSULIN AXIS & INCRETINS
    # Direct modulators of AMPK, PPARγ, K_ATP channels, DPP-4, α-glucosidase, etc.
    # --------------------------------------------------------------------------
    "Metformin hydrochloride",
    "Phenformin hydrochloride",
    "Tolbutamide",
    "Chlorpropamide",
    "Tolazamide",
    "Glipizide",
    "Glibenclamide",
    "Gliclazide",
    "Glimepiride",
    "Gliquidone",
    "Repaglinide",
    "Nateglinide",
    "Acarbose",
    "Miglitol",
    "Voglibose",
    "Pioglitazone",
    "Troglitazone",
    "Rosiglitazone hydrochloride",
    "Linagliptin",

    # --------------------------------------------------------------------------
    # 2. LIPID & LIPOPROTEIN METABOLISM (HMG-CoA Reductase, PPARα, NPC1L1, MTP)
    # --------------------------------------------------------------------------
    "Atorvastatin",
    "Rosuvastatin",
    "Pitavastatin calcium",
    "Gemfibrozil",
    "Fenofibrate",
    "Bezafibrate",
    "Clofibrate",
    "Clofibric acid",
    "Ciprofibrate",
    "Ezetimibe",
    "Lomitapide",
    "Niacin",
    "Acipimox",
    "Probucol",

    # --------------------------------------------------------------------------
    # 3. FATTY-ACID OXIDATION & CARDIAC ENERGY METABOLISM MODULATORS
    # CPT-1 inhibitors, 3-KAT inhibitors, or metabolic shifts from FA to glucose
    # --------------------------------------------------------------------------
    "Trimetazidine dihydrochloride",
    "Ranolazine",
    "Perhexiline maleate",

    # --------------------------------------------------------------------------
    # 4. PURINE / PYRIMIDINE / URIC ACID METABOLISM
    # Xanthine oxidase inhibitors, uricosurics (URAT1 inhibitors), xanthines
    # --------------------------------------------------------------------------
    "Allopurinol",
    "Probenecid",
    "Sulfinpyrazone",
    "Benzbromarone",
    "Caffeine",
    "Theophylline monohydrate",
    "Theobromine",

    # --------------------------------------------------------------------------
    # 5. THYROID AXIS & STEROID/HORMONAL METABOLISM
    # Thyroid hormones, peroxidase inhibitors, aromatase & 5α-reductase inhibitors
    # --------------------------------------------------------------------------
    "(L) Thyroxine",
    "Propylthiouracil",
    "Carbimazole",
    "Anastrozole",
    "Letrozole",
    "Bicalutamide",
    "Flutamide",
    "Nilutamide",
    "Tamoxifen citrate",
    "Toremifene",
    "Raloxifene hydrochloride",

    # --------------------------------------------------------------------------
    # 6. APPETITE / WEIGHT / MONOAMINERGIC METABOLIC REGULATION
    # CB1 antagonists, reuptake inhibitors, sympathomimetics affecting metabolic rate
    # --------------------------------------------------------------------------
    "Sibutramine hydrochloride",
    "Dexfenfluramine hydrochloride",
    "Phentermine hydrochloride",
    "Rimonabant",

    # --------------------------------------------------------------------------
    # 7. FOLATE / ONE-CARBON & ENZYMATIC AMINO ACID METABOLISM
    # Dihydrofolate reductase (DHFR) and thymidylate synthase inhibitors
    # --------------------------------------------------------------------------
    "Amethopterin",  # Methotrexate
    "Trimethoprim",
    "Pyrimethamine",
    "Folic acid",
    "5-fluorouracil",
    "Floxuridine",
    "Capecitabine",
    "Tegafur",
    "Carmofur",

    # --------------------------------------------------------------------------
    # 8. ALCOHOL / ALDEHYDE / ENZYME MODULATORS & MISCELLANEOUS METABOLIC TARGETS
    # ALDH inhibitors, ADH inhibitors, mitochondrial respirators, HDACs
    # --------------------------------------------------------------------------
    "Disulfiram",
    "Fomepizole",
    "Vorinostat",     # HDAC inhibitor -> global histone/non-histone acetylation metabolic shifts
    "Idebenone",      # CoQ10 analog affecting mitochondrial electron transport chain
    "Valproic acid",  # Carnitine depletion, mitochondrial fatty acid beta-oxidation inhibition
]

drug_match_table[
    "validated_metabolic"
] = [
    x in validated_metabolic for x in drug_match_table["schuhknecht_name"].to_list()
]

# Save results
################################################################################

drug_match_table.to_csv(os.path.join(out_dir, "drug_match_table.csv"))
print("drug_match_table.csv saved in %s." % out_dir)
mc_delta_avg.to_csv(os.path.join(out_dir, "mc_delta_avg_matched.csv"))
print("mc_delta_avg_matched.csv saved in %s." % out_dir)
sch_meanlog2FC.to_csv(os.path.join(out_dir, "sch_log2FC_matched.csv"))
print("sch_log2FC_matched.csv saved in %s." % out_dir)
sch_pvals.to_csv(os.path.join(out_dir, "sch_pvals_matched.csv"))
print("sch_pvals_matched.csv saved in %s." % out_dir)
sch_qvals.to_csv(os.path.join(out_dir, "sch_qvals_matched.csv"))
print("sch_qvals_matched.csv saved in %s." % out_dir)