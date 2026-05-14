import sys
from pyprojroot.here import here
import pandas as pd

# Set basedir with here(), based on presence of .git file
here()
# Add script dir to path to import gcep.py and gcep_config.py
sys.path.insert(0, here('scripts'))
from clingenAPI import clingenAPI
import clingenAPI_config
from datetime import datetime
import os

print(here())

# Load api information from gcep_config
api_dict = {'pird':clingenAPI_config.api_key_pird, 'scid':clingenAPI_config.api_key_scid}
affiliation_dict = {'pird':clingenAPI_config.affiliation_pird, 'scid':clingenAPI_config.affiliation_scid}
   
# Set gcep of interest (scid or pird)    
# active_gcep = 'scid'
active_gcep = 'pird'

# # combine = "BMA"
# combine = "funSimAvg"

# start_date="2020-11-01"
# end_date="2021-12-31"
start_date = "2024-12-01"
end_date = "2025-01-31"

# Query GCEP for HPO data using gcep class in gcep.py
clingen_query = clingenAPI(
    api_key = api_dict[active_gcep], 
    clingen_url = clingenAPI_config.clingen_url,
    status = "approved",
    affiliation=affiliation_dict[active_gcep],
    start = start_date, 
    end = end_date
)

# Generate HPO table from query using hpoTable() method
df_probands = clingen_query.hpoTable()

# Save data
## Create timestamped directory to save data
dt_str = datetime.now().strftime('%Y%m%d_%H%M%S')
save_dir = here(f'data/clingen/api/{dt_str}')
os.makedirs(save_dir, exist_ok=True)
## Create save path with gcep and query date range
date_range = f's{start_date.replace("-", "")}_e{end_date.replace("-", "")}'
save_path = os.path.join(save_dir, f'{active_gcep}_hpo_{date_range}.csv.gz')
print(f'Saving to {save_path}')

df_probands.to_csv(here(save_path))