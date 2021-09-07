
from bioservices import EUtils
e = EUtils(email="silas.kieser@unige.ch")

import pandas as pd

from tqdm import tqdm
import os
import warnings


def parse_sample(sample):

    if 'SAMPLE_ATTRIBUTES' in sample:
        parsed_sample_atributes=sample['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']

        if type(parsed_sample_atributes)== list:

            sample_atributes= pd.DataFrame(parsed_sample_atributes)
            sample_atributes.index=sample_atributes.TAG
            sample_atributes=sample_atributes.VALUE
        else:
            sample_atributes=pd.Series(parsed_sample_atributes['VALUE'],
                                       index=[parsed_sample_atributes['TAG']]
                                      )
    else:
        
        sample_atributes = pd.Series()

    sample_data=pd.concat((
                pd.Series(sample['IDENTIFIERS']),
                pd.Series(sample['SAMPLE_NAME']),
                sample_atributes
                ))
    try:
        sample_data['Sample Title']=sample['TITLE']
    except KeyError:
        pass



    return sample_data


def parse_set_of_runs(parsed_xml):


    have_multiple_runs=False
    D=pd.DataFrame()

    for i,r in enumerate(parsed_xml['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']):
        #sample data
        run_data= parse_sample(r['SAMPLE'])
        

        if not run_data.index.is_unique:

            #print(f'Duplicated values for run {run_data.name}, I wil drop the first values')

            #print(run_data.loc[run_data.index.duplicated(keep=False)])

            run_data=run_data.loc[run_data.index.duplicated(keep='last')]

        # Project
        try:
            run_data['BioProject']=r['STUDY']['IDENTIFIERS']['EXTERNAL_ID']

            run_data['BioProject Title']=r['STUDY']['DESCRIPTOR']['STUDY_TITLE']
        except KeyError:
            warnings.warn('Run has no associated bioproject id')

        # Run
        run_set=r['RUN_SET']['RUN']
        if type(run_set) ==list:
            have_multiple_runs=True
        else:
            #else print(f"I have multiple runs for sample ")

            run_set = [run_set]



        for run_info in run_set:


            run_data['RunID']= run_info['IDENTIFIERS']['PRIMARY_ID']
            run_data.name= run_data['RunID']

            D=D.append(run_data)


            #assert ~pd.isnull(run_data.BioProject)
            
        if have_multiple_runs:
            warnings.warn("Some samples have multiple runs")

    return D








def download_sra_metadata(biosamples,xml_folder,database='sra'):
    """download metadata for biosamples in chuncks of 200"""

    n_max=200

    if os.path.exists(xml_folder):
        print("XML folder existsts already, will continue downloading missing files")
    else:
        os.makedirs(xml_folder)

    if type(biosamples) != str:
        biosamples= list(biosamples)


    for i,ret_Start in tqdm(list(enumerate(range(0,len(biosamples),n_max)))):

        xml_file=f'{xml_folder}/metadata_{i}.xml'

        if not os.path.exists(xml_file):

            handle=e.EFetch(database,id=biosamples[ret_Start:ret_Start+n_max],retmode='xml')

            if type(handle)==int:
                print(f'got error {handle.decode()} for ids:',biosamples[ret_Start:ret_Start+n_max])

            with open(xml_file,'wb') as outf:
                outf.write(handle)

def parse_sra_metadata(xml_folder):


    parsed=[]

    for xml_file in tqdm(os.listdir(xml_folder)):
        
        tsv_file= os.path.join(xml_folder,xml_file.replace('.xml','.tsv'))
        
        if os.path.exists(tsv_file):
            table = pd.read_table(tsv_file,index_col=0)
        else:
            parsed_xml =e.parse_xml(open(os.path.join(xml_folder,xml_file),'rb').read())

            table =parse_set_of_runs(parsed_xml)
            table.to_csv(tsv_file, sep='\t')
            
        parsed.append(table)

    df= pd.concat(parsed)
    del parsed

    return df

from bioservices.eutils import EUtilsParser

def get_attribute_from_xml(parsed_xml,key,expect_list=True):
    """SRA data structures often contain nested drectoreis with list as output
        but sometimes only 1 level of directory is present with a string.

        This function assures that wlways a list is returned.
    """


    first_level=parsed_xml[key]

    dtype= type(first_level)

    if type(first_level) == EUtilsParser :
        keys=first_level.keys()
        if len(keys)>1:
            TypeError(f"I didn't expect multiple keys on first level for {key}. Got {keys}")
        else:
            values= list(first_level.values())[0]

    elif dtype == type(None):
        raise KeyError(f"Returned NoneType for key {key}")

    elif dtype==str:
        values=str
    else:
        raise TypeError(f"Did't expected type {dtype} for key {key}")

    if expect_list:


        dtype=type(values)

        if  dtype==list:
            return values
        else :
            return [values]
    else:
        return values

def search_indicies(db,accessions):

    queries=[]
    mapped_ids=[]

    for acc in tqdm(accessions):
        res=e.ESearch(db,acc)

        idlist= res['idlist']

        if len(idlist)> 1:
            print(f'Got multiple Ids for search {acc}')

        for id in idlist:
            queries.append(acc)
            mapped_ids.append(id)

    return pd.Series(queries,index=mapped_ids,name='Accession')
