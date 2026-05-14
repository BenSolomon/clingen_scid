import sys
from pyprojroot.here import here
import argparse
import pandas as pd
import gzip

here()
sys.path.insert(0, here('scripts'))
from scripts.clingenScrape import clingen_scrape


def main():
    parser = argparse.ArgumentParser(description="Scrape ClinGen probands")
    parser.add_argument('HGNC_ID', type=str, nargs='?', default="HGNC:12731", help='The gene HGNC ID to scrape')
    parser.add_argument('SAVE_DIR', type=str, nargs='?', default = here(), help='Directory to save output files')
    args = parser.parse_args()
    
    print(f'### STARTING {args.HGNC_ID}', file=sys.stderr)
    
    
    clingen_query = clingen_scrape(args.HGNC_ID)
    print(f'{args.HGNC_ID}\t{clingen_query.valid_entry}')

        
    if clingen_query.table is not None:
        with gzip.open(f'{args.SAVE_DIR}/{args.HGNC_ID.replace(":", "_")}.pkl.gz', 'wb') as f:
            clingen_query.table.to_pickle(f)
        df_hpo=clingen_query.hpoTable()
        if not df_hpo.empty:
            df_hpo.to_csv(f'{args.SAVE_DIR}/{args.HGNC_ID.replace(':', '_')}_hpo.csv' , index=False)


if __name__ ==  "__main__":
    main()
