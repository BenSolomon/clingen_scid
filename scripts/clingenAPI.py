import requests
import clingenAPI_config
import pandas as pd
import re

# API documentation: https://vci-gci-docs.clinicalgenome.org/vci-gci-docs/gci-help/gci-api#affiliations-list

class clingenAPI:
    """
    Class that takes parameters query to the ClinGen API and parses response data
    Includes methods to return data in a variety of formats
    """    
    def __init__(self, api_key, clingen_url, status, affiliation, start, end):  
        self.api_key = api_key
        self.clingen_url = clingen_url
        self.params = {
            "target": "gci",
            "status": status,
            "affiliation": affiliation,
            "start": start,
            "end": end
        }
        self.response_proband = self._apiGet()
        self.response_summary = self._apiGet(summary = True)
        self.table = self._makeDataFrame()
        self.genes = self.table['Gene'].unique()
        self.n_genes = len(self.genes)
    
    def _apiGet(self, summary = False):
        """
        Internal method used by __init__ to make a get request to the ClinGen API
        
        Args:
            summary (bool): False if proband data is requested, True if summary data is requested. 
                If True, the name of the parameter is set to "summary" to get summary data.
        Returns:
            _type_: _description_
        """
        if summary:
            self.params['name'] = "summary"        
                
        response = requests.get(
            f"{self.clingen_url}/snapshots",
            headers={"x-api-key": self.api_key},
            params=self.params)
        return response
    
    def _makeDataFrame(self):
        """
        Internal method used by __init__ to convert API response to a pandas dataframe
        Joins the proband response and the summary table response into a single dataframe

        Returns:
            
            pd.DataFrame: DataFrame with columns:
                `Gene`: Gene name
                `Disease`: Disease name
                `Mode of Inheritance`: Mode of inheritance
                `Status`: Whether gene entry is published or approved
                `Approval Date`: When the entry was approved
                `Approval Review Date`: When the approval was reviewed
                `Published date`: When the entry was published
                `probands`: JSON tree of proband data
                `GCEP Affiliation`: Affiliation of the GCEP
                `Final Classification`: Final classification
                `Final Classification Date`: Date of final classification
                `Genetic Total Points`: Total points for the genetic data
                `Experimental Total Points`: Total points for the experimental data
                `Total Points`: Total points for the genetic and experimental data
                `Proband Count`: Number of probands
                `Scored Proband Count`: Number of scored probands
                `Earliest PMID Year`: Year of the earliest PMID
                `Most Recent PMID Year`: Year of the most recent PMID
                `GDM UUID`: Unique identifier for the GDM
        """        
        json_proband = self.response_proband.json()
        json_summary = self.response_summary.json()
        df_proband = pd.DataFrame(json_proband)
        df_proband = df_proband.sort_values(by='Approval Date', ascending=False).groupby(['Gene', 'Disease', 'Mode of Inheritance']).first().reset_index()
        df_summary = pd.DataFrame(json_summary)
        df_summary = df_summary.sort_values(by='Final Classification Date', ascending=False).groupby(['Gene', 'Disease', 'Mode of Inheritance']).first().reset_index()
        df_merged = pd.merge(
            df_proband, 
            df_summary, 
            on=['Gene', 'Disease', 'Mode of Inheritance'], 
            how='outer'
        )
        print(f'df_proband shape: {df_proband.shape}')
        print(f'df_summary shape: {df_summary.shape}')
        print(f'df_merged shape: {df_merged.shape}')
        return df_merged            
    
    def _singleGeneProbandTable(self, row):
        """
        Internal method used as pandas apply function to format proband data for a single gene
        Parses the proband sub-json format and joins it with the original Gene and Disease data

        Args:
            row (_type_): _description_

        Returns:
            _type_: _description_
        """        
        df = pd.DataFrame(row['probands'])
        gene = row['Gene']
        disease = row['Disease']
        df['Gene'] = gene
        df['Disease'] = disease
        df = df[['Gene', 'Disease'] + [col for col in df.columns if col not in ['Gene', 'Disease']]]
        return df
    
    
    def _formatHpoString(self, hpo_str):
        """
        Internal method used to reformat how HPO terms are stored in ClinGen json data
        Intake format from ClinGen JSON follows pattern "HPO term (HPO_ID)"
        Output format is a dictionary with keys "HPO_ID" and "HPO_term"
        This allows the HPO_ID to be used with HPO3 module for further analysis
        
        Args:
            hpo_str (_type_): _description_

        Returns:
            _type_: _description_
        """        
        try: 
            match = re.match(r"(.+?)\s*\((HP:\d+)\)", hpo_str)
            if match:
                hpo_dict = {
                    "HPO_ID": match.group(2),
                    "HPO_term": match.group(1)
                }
                return hpo_dict
            else:
                return hpo_str, None
        except:
            return hpo_str, None
            
    def probandTable(self):
        """
        Takes data from self.table and returns a new table with each proband in the dataset.
        Each entry in self.table['probands'] is a json tree of proband data that is converted to a dataframe
        The original Gene and Disease from the original row is added to this new dataframe
        The resulting proband dataframes are concatenated into a single dataframe
        """        
        df = self.table.apply(self._singleGeneProbandTable, axis=1)
        df = pd.concat(df.to_list(), axis = 0, ignore_index=True)
        return df
    
    def hpoTable(self):
        """
        Takes data from self.table and returns a new table for each HPO term for every 
        Gene-Disease-proband combination
        Includes application of _formatHpoString to reformat HPO terms

        Returns:
            _type_: _description_
        """        
        df = self.probandTable()
        df = df[['Gene', 'Disease', 'label','HPO terms']].explode("HPO terms", ignore_index=True)
        df['HPO terms'] = df['HPO terms'].apply(self._formatHpoString) 
        df = df.join(pd.json_normalize(df.pop('HPO terms')))
        df = df.dropna(subset=['HPO_ID'])
        return df
    
    #TODO: Add a method that extracts all proband variant information
    
                
if __name__ == "__main__":
    # If run as script to test, create a clingenAPI object and print some data 
    clingen_query = clingenAPI(
        api_key = clingenAPI_config.api_key_pird, 
        clingen_url = clingenAPI_config.clingen_url,
        status = "approved",
        affiliation=clingenAPI_config.affiliation_pird,
        start = "2024-12-01", 
        end = "2025-01-31"
    )
    
    print(clingen_query.n_genes)
    print(clingen_query.genes)
    print(clingen_query.table)
    print(pd.DataFrame(clingen_query.table['probands']))
    # print(clingen_query.probandTable())
    # print(clingen_query.hpoTable())
