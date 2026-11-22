# All the knobs live here.
from pathlib import Path

DATA = Path("data")
DATA.mkdir(exist_ok=True)
SEED = 42

# --- Who counts as "tumor" and "healthy"? Adjust after running step 1. ---
TUMOR_FILTER = ("tissue_general == 'breast' and disease == 'breast cancer' "
                "and cell_type == 'malignant cell' and is_primary_data == True")
HEALTHY_FILTER = "disease == 'normal' and is_primary_data == True"

# Healthy tissues we must NOT hurt. Names that don't exist in the Census are skipped (and printed).
HEALTHY_TISSUES = [
    "breast", "lung", "heart", "brain", "liver", "kidney", "colon", "small intestine",
    "intestine", "pancreas", "stomach", "skin of body", "muscle organ", "blood",
    "bone marrow", "spleen", "lymph node", "thymus", "urinary bladder", "prostate gland",
    "uterus", "ovary", "eye", "adipose tissue", "esophagus", "thyroid gland",
    "adrenal gland", "tongue", "vasculature", "nose", "salivary gland", "testis",
]

# --- Size limits (keep Colab T4 / ~12 GB RAM happy) ---
MAX_TUMOR_CELLS = 20000
HEALTHY_PER_TISSUE = 2000

# --- Definitions ---
MIN_COUNT = 1             # a gene is "positive" in a cell if raw count >= this
MIN_TUMOR_PREVALENCE = 0.05   # keep genes positive in >=5% of TRAIN tumor cells
OFFTARGET_CAP = 0.05      # max fraction of cells in ANY healthy tissue we tolerate hitting
TEST_FRACTION = 0.3       # fraction of DONORS (patients) held out

SURFY_URL = "http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx"