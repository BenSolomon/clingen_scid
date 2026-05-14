import requests
from bs4 import BeautifulSoup
import pandas as pd
import re
from io import StringIO

class clingenScrape:
    """
    Class that takes parameters scrape the ClinGen HTML site and parses response data
    Includes methods to return data in a variety of formats
    """    
    def __init__(self, hgnc_id):
        self.hgnc_id = hgnc_id
        self.clingen_base_url = 'https://search.clinicalgenome.org'
        self.clingen_gene_url = f"{self.clingen_base_url}/kb/genes/{self.hgnc_id}"
        self.gene_response = requests.get(self.clingen_gene_url)
        self.valid_gene = self._validGene() 
        self.disease_entries = self._getClingenDiseaseEntries()
        self.valid_entry = self.disease_entries is not None    
        
        if self.valid_entry:
            self.disease_responses = [requests.get(x) for x in self.disease_entries]
            table = [self._get_table(x) for x in self.disease_responses]
            table = [x for x in table if x is not None]
            self.table = pd.concat([self._get_table(x) for x in self.disease_responses], ignore_index=True) if len(table) > 0 else None
        else:
            self.disease_responses = None
            self.table = None
            
    def _validGene(self):
        """
        Internal method to determine if the HGNC ID is valid

        Returns:
            _type_: _description_
        """
        if self.gene_response.status_code != 200:
            return False
        else: 
            soup = BeautifulSoup(self.gene_response.text, 'html.parser')
            error_header = soup.find('h1', string='Error retrieving Gene details')
            error_div = soup.find('div', {'class': 'alert-danger'})
            
            return not (error_header or error_div)

    def _getClingenDiseaseEntries(self):
        """If an HGNC ID is valid, determine if there are any gene-disease curations
        in ClinGen and return their URLs

        Returns:
            _type_: _description_
        """        
        if not self.valid_gene:
            return None
        else:
            soup = BeautifulSoup(self.gene_response.text, 'html.parser')
            links = soup.find_all('a', class_='btn btn-xs btn-success btn-block btn-report')
            output = [link.get('href') for link in links if link.get('href').startswith('/kb/gene-validity/CGGV:assertion')]
            output = [f'{self.clingen_base_url}{link}' for link in output]
            if len(output) == 0:
                return None
            return output
        
    def _parseHpoFreeText(self, input_string):
        """
        Takes a string of format "HPO terms(s): HPO term (HPO_ID), HPO term (HPO_ID) Free text: Free text"
        or "Free text: Free text HPO terms(s): HPO term (HPO_ID), HPO term (HPO_ID)"
        and returns a dictionary with keys "phenotype_hpo" and "phenotype_text"

        Args:
            input_string (_type_): _description_

        Returns:
            _type_: _description_
        """
        if type(input_string) != str:
            return None
                
        hpo_terms_start = "HPO terms(s):"
        free_text_start = "Free text:"

        hpo_terms_index = input_string.find(hpo_terms_start)
        free_text_index = input_string.find(free_text_start)

        hpo_terms = None
        free_text = None

        if hpo_terms_index != -1:
            if free_text_index != -1:
                hpo_terms = input_string[hpo_terms_index + len(hpo_terms_start):free_text_index].strip()
            else:
                hpo_terms = input_string[hpo_terms_index + len(hpo_terms_start):].strip()

        if free_text_index != -1:
            free_text = input_string[free_text_index + len(free_text_start):].strip()

        return {"phenotype_hpo": hpo_terms, "phenotype_text": free_text}
    
    def _formatHpoList(self, hpo_string):
        """
        Converts a string of HPO terms into a list of HPO terms.

        Args:
            hpo_string (str): A string containing HPO terms.

        Returns:
            list: A list of HPO terms.
        """
        if hpo_string is None:
            return None
        else:
            pattern = r"([^()]+?\(HP:\d+\))"
            matches = re.findall(pattern, hpo_string)
            return matches
    
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
                hpo_dict = {
                    "HPO_ID": None,
                    "HPO_term": None
                }
                return hpo_dict
        except:
            hpo_dict = {
                    "HPO_ID": None,
                    "HPO_term": None
                }
            return hpo_dict
        
    def _getDiseaseFeatures(self, disease_response):
        """
        Get additional data for gene-disease relationship that occurs
        outside of proband table. This inlcudes the disease name, the
        disease MONDO ID, the mode of inheritance, and the expert panel 
        the gene was associated with, and the classificaiton status of
        the gene-disease relationship

        Args:
            disease_response (_type_): _description_

        Returns:
            _type_: _description_
        """        
        
        soup = BeautifulSoup(disease_response.text, 'html.parser')
        
        # Get MONDO disease ID
        mondo_tag = soup.find('div', string=re.compile(r"MONDO:\d+"))
        if mondo_tag:
            mondo =  mondo_tag.text
        else:
            mondo = None
        
        # Get gene-disease calssification status
        target_href = "https://www.clinicalgenome.org/docs/gene-disease-validity-classification-information/"
        classification_tag = soup.find("a", href=target_href)
        if classification_tag:
            classification = classification_tag.get_text(strip=True)
        else:
            classification = None
        
        # Get additional stats in <dt><dd> structure    
        dt_output = []
        for dt in soup.find_all('dt'):
            entry = ' '.join(dt.find_all(string=True, recursive=False)).strip()
            entry = entry.split(':')[0]
            entry = entry.replace('\n', '')
            entry = entry.rstrip()
            dt_output.append(entry)
            
        dd_output = []
        for dd in soup.find_all('dd'):
            entry = ' '.join(dd.find_all(string=True, recursive=False)).strip()
            entry = entry.split(':')[0]
            entry = entry.replace('\n', '')
            entry = entry.rstrip()
            dd_output.append(entry)
            
        # Compile data into dictionary   
        output_dict = dict(zip(dt_output, dd_output))
        output_dict = {k:output_dict.get(k) for k in output_dict.keys() if k in ['Gene', 'Disease', 'Mode of Inheritance', 'Expert Panel']}
        output_dict.update({"MONDO": mondo, "classification": classification})
        return output_dict
    
    def _get_table(self, disease_response):
        disease_features = self._getDiseaseFeatures(disease_response)
                
        soup = BeautifulSoup(disease_response.text, 'html.parser')
        table = soup.find('table', {'id': 'geclv'})
        if table is None:
            return None
        
        df_proband = pd.read_html(StringIO(str(table)))[0]
        
        # Split phenotype string in to HPO terms and free text
        df_hpo = df_proband['Proband Phenotypes'].apply(self._parseHpoFreeText)        
        df_proband['HPO terms'] = df_hpo.apply(lambda x: x['phenotype_hpo'] if x else None)
        df_proband['HPO terms'] = df_proband['HPO terms'].apply(self._formatHpoList) 
        df_proband['HPO free text'] = df_hpo.apply(lambda x: x['phenotype_text'] if x else None)
        df_proband.drop(columns=['Proband Phenotypes'], inplace=True)
        
        df_proband.insert(0, 'Gene', disease_features['Gene'])
        df_proband.insert(1, 'HGNC', self.hgnc_id)
        df_proband.insert(2, 'Disease', disease_features['Disease'])
        df_proband.insert(3, 'MONDO', disease_features['MONDO'])
        df_proband.insert(4, 'Inheritance', disease_features['Mode of Inheritance'])
        df_proband.insert(5, 'GCEP', disease_features['Expert Panel'])
        df_proband.insert(6, 'Classification', disease_features['classification'])
        
        df_proband.rename(columns={'Proband Label': 'label'}, inplace=True)
        
        if df_proband.empty:
            return None
        
        return(df_proband)
    
    def hpoTable(self):
        """
        Takes data from self.table and returns a new table for each HPO term for every 
        Gene-Disease-proband combination
        Includes application of _formatHpoString to reformat HPO terms

        Returns:
            _type_: _description_
        """
        if self.table is None:
            return None
        else:                 
            df = self.table
            df = df[['Gene', 'Disease', "MONDO", 'label','HPO terms']].explode("HPO terms", ignore_index=True)
            df['HPO terms'] = df['HPO terms'].apply(self._formatHpoString) 
            df = df.join(pd.json_normalize(df.pop('HPO terms')))
            df = df.dropna(subset=['HPO_ID'])
            return df
    
def clingen_diagnostics(HGNC_ID):
    clingen_query = clingenScrape(HGNC_ID)
    print(clingen_query.hgnc_id)
    print(clingen_query.gene_response)
    print(clingen_query.valid_gene)
    print(clingen_query.disease_entries)
    print(clingen_query.table)
    print(clingen_query.hpoTable())
    print("\n")

                
if __name__ == "__main__":

    # If run as script to test, create a clingen object and print some data     
    print("Query for FASL")
    clingen_diagnostics("HGNC:11936")
    
    print("Query for STAT3 - Two curation entries")
    clingen_diagnostics("HGNC:11364")
    
    print("Query for gene with no curation")
    clingen_diagnostics("HGNC:5")
    
    print("Query for invalid HGNC")
    clingen_diagnostics("HGNC:123412341234")
    
    print("Query for A2ML1 - Currated but no HPO")
    clingen_diagnostics("HGNC:23336")
    
    print("Query for AARS1 - ")
    clingen_diagnostics("HGNC:20")
    
    print("Query for ABCC8 - Has a ClinGen page, but no proband table")
    clingen_diagnostics("HGNC:59")
    
    print("Query for ACAT2 - Has proband table, but table is empty ")
    clingen_diagnostics("HGNC:94")
    
