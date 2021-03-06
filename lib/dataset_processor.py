#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import datetime
import os
import json
import re
import time
import json
import requests
from ftplib import FTP

from response import Response


class DatasetProcessor:

    # Class variables


    ###############################################################################################
    # Constructor
    def __init__(self):
        self.status = 'OK'
        self.verbose = None
        self.response = None
        self.tasks_todo = []
        self.base_dir = "/proteomics/peptideatlas2/archive/Arabidopsis"
        self.state = { 'processing_state': 'Unknown', 'todo': 'assess' }
        self.datasets = { 'identifiers': {} }
        self.compressed_extension = 'zip'

        response = Response()
        self.response = response


    ###############################################################################################
    def add_dataset(self, dataset_id):
        """Add a dataset to the dict of tracked datasets along with basic default attributes
        so that it may get processed
        """

        # Set the self state and begin examining what we have
        response = self.response
        response.info(f"Add dataset_id {dataset_id} to the list of tracked datasets")

        dataset = { 'status': 'QUEUED', 'state': { 'processing_state': 'Queued', 'message': 'Ready to begin setup' },
            'dataset_id': dataset_id, 'metadata': { 'location': f"{self.base_dir}/{dataset_id}" } }
        self.datasets['identifiers'][dataset_id] = dataset

        return response


    ###############################################################################################
    def process(self):
        """Top level method to loop over all datasets in the system to process

        """

        # Set the self state and begin examining what we have
        response = self.response
        #response.info(f"Begin processing all datasets")

        for dataset_id in self.datasets['identifiers']:
            self.process_dataset(dataset_id)

        return response


    ###############################################################################################
    def process_dataset(self, dataset_id):
        """Finish the current step or take the next step in processing one dataset
        """

        #### Set the self state and begin examining what we have
        response = self.response

        # Get the dataset handle and set status
        dataset = self.datasets['identifiers'][dataset_id]

        #response.info(f"Processing dataset {dataset_id}. Currently in state {dataset['state']['processing_state']}")
        dataset['status'] = 'PROCESSING'

        #### Apply needed processing tasks to the dataset
        loop_counter = 0
        while 1:

            #### If in an error state, return
            if dataset['status'] == 'ERROR':
                return response

            # Consult the processing_state to figure out where to go
            if dataset['state']['processing_state'] == 'Queued':
                self.assess_setup(dataset_id)

            elif dataset['state']['processing_state'] == 'Set up':
                self.assess_setup(dataset_id)

            elif dataset['state']['processing_state'] == 'Ready to download':
                self.assess_download(dataset_id)
                #if dataset['state']['processing_state'] == 'Ready to download':
                #    self.download_dataset(dataset_id)

            elif dataset['state']['processing_state'] == 'Downloading':
                self.assess_download(dataset_id)

            elif dataset['state']['processing_state'] == 'Ready to convert':
                dataset['state']['processing_state'] = 'Wait'
                #self.assess_conversion(dataset_id)
                #if dataset['state']['processing_state'] == 'Ready to convert':
                #       self.convert_dataset(dataset_id)

            elif dataset['state']['processing_state'] == 'Converting':
                self.assess_conversion(dataset_id)

            elif dataset['state']['processing_state'] == 'Ready to compress':
                self.assess_conversion(dataset_id)
                if dataset['state']['processing_state'] == 'Ready to compress':
                    self.compress_dataset(dataset_id)

            elif dataset['state']['processing_state'] == 'Compressing':
                self.assess_conversion(dataset_id)

            elif dataset['state']['processing_state'] == 'Wait':
                #response.info(f"Done with what we are able to do so far")
                dataset['status'] = 'READY'
                return response

            else:
                strange_state = dataset['state']['processing_state']
                dataset['status'] = 'ERROR'
                #dataset['state']['processing_state'] = 'UnrecognizedState'
                #response.error(f"Unrecognized processing state {strange_state}", error_code=dataset['state']['processing_state'])
                #response.warning(f"Unrecognized processing state {strange_state}")
                return response

            # Don't need a loop here after all??
            break

            #### Provide a mechanism to escape an infinite loop
            loop_counter += 1
            if loop_counter > 100:
                response.error(f"process_dataset loop_counter reached {loop_counter}. Likely runaway infinite loop", error_code='InfiniteLoopExit')
                return response


    ###############################################################################################
    def assess_setup(self, dataset_id):
        """Assess the current state of the setup of the dataset
        """

        # Set the self state and begin examining what we have
        response = self.response
        response.info(f"Assessing setup of dataset {dataset_id}")

        #### Ensure that we have a dataset_id and its in our list
        if dataset_id is None:
            self.status = 'ERROR'
            self.state['processing_state'] = 'MissingDatasetId'
            response.error(f"No dataset_id has been specified. Set dataset_id", error_code=self.state['processing_state'])
            return response
        if dataset_id not in self.datasets['identifiers']:
            self.status = 'ERROR'
            self.state['processing_state'] = 'DatasetIdNotInList'
            response.error(f"The dataset_id {dataset_id} is not in the processing list", error_code=self.state['processing_state'])
            return response

        # Get the dataset handle and set status
        dataset = self.datasets['identifiers'][dataset_id]
        dataset['status'] = 'PROCESSING'

        # Ensure that we have a location
        response.info(f"Verify storage location")
        if dataset['metadata']['location'] is None:
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'MissingLocation'
            response.error(f"No location has been specified. Set location", error_code=dataset['state']['processing_state'])
            return response

        # Ensure that the location ends with the dataset_id
        match = re.match(r'.+?/([A-Z\d]+)$', dataset['metadata']['location'])
        if match:
            if match.group(1) != dataset_id:
                dataset['status'] = 'ERROR'
                dataset['state']['processing_state'] = 'InvalidLocation'
                response.error(f"Location does not end with the dataset_id. It must", error_code=dataset['state']['processing_state'])
                return response
        else:
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'InvalidLocation'
            response.error(f"Dataset location cannot be parsed", error_code=dataset['state']['processing_state'])
            return response

        # Set the mode to either assess (as a check on work that may have been done previously)
        # or verify (to verify that work that was just done has completed)
        mode = 'assess'
        if dataset['state']['processing_state'] == 'Set up':
            mode = 'verify'

        # See if the location exists
        target_path = dataset['metadata']['location']
        response.info(f"Check dataset location {dataset['metadata']['location']}")
        if not os.path.exists(target_path):
            if mode == 'assess':
                self.create_destination(dataset_id)
                if dataset['status'] == 'ERROR':
                    return
            if mode == 'verify':
                dataset['status'] = 'ERROR'
                dataset['state']['processing_state'] = 'LocationMissing'
                response.error(f"Tried to create location but it is not found", error_code=dataset['state']['processing_state'])
                return response

        # See if the ProteomeXchange record exists
        response.info(f"Check for ProteomeXchange metadata file")
        target_path = f"{dataset['metadata']['location']}/data/ProteomeXchange.json"
        if not os.path.exists(target_path):
            if mode == 'assess':
                self.fetch_px_record(dataset_id)
                if dataset['status'] == 'ERROR':
                    return
            if mode == 'verify':
                dataset['status'] = 'ERROR'
                dataset['state']['processing_state'] = 'PXRecordMissing'
                response.error(f"Tried to verify PX record but it is not found", error_code=dataset['state']['processing_state'])
                return response

        # See if the PX data is parsable that that we have it in memory
        response.info(f"Check parsing of the ProteomeXchange metadata file")
        with open(target_path,'r') as infile:
            try:
                px_data = json.load(infile)
            except Exception as error: 
                dataset['status'] = 'ERROR'
                dataset['state']['processing_state'] = 'CannotParsePXJSON'
                response.error(f"Unable to parse the ProteomeXchange JSON file {target_path}: {error}", error_code=dataset['state']['processing_state'])
                return response

        # See if the information is already stored
        if 'px_data' in dataset['metadata'] is not None:
            response.info(f"ProteomeXchange information is already stored. Will not overwrite")
        else:
            response.info(f"Importing ProteomeXchange record information")
            dataset['metadata']['px_data'] = px_data


        # Extract some key information from the px_data

        # Extract the FTP location from the PX record
        ftp_location = None
        dataset['metadata']['ftp_location'] = None
        if 'fullDatasetLinks' in px_data:
            for term in px_data['fullDatasetLinks']:
                if term['name'] == 'Dataset FTP location':
                    ftp_location = term['value']
        if ftp_location is not None:
            dataset['metadata']['ftp_location'] = ftp_location
        else:
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'CannotFindfullDatasetLinks'
            response.error(f"Unable to find fullDatasetLinks in PX record for {dataset_id}", error_code=dataset['state']['processing_state'])


        # If we got this far, then we're ready to download
        response.info(f"Ready to move to download")
        dataset['state']['processing_state'] = 'Ready to download'
        return response


    ###############################################################################################
    # Create the dataset destination
    def create_destination(self, dataset_id):
        """Create the destination area for the dataset
        """

        #### Set the self state and begin examining what we have
        response = self.response
        response.info(f"Creating dataset destination")

        # Get the dataset handle
        dataset = self.datasets['identifiers'][dataset_id]

        #### See if the location exists. If not, create it
        target_path = dataset['metadata']['location']
        if not os.path.exists(target_path):
            response.info(f"Create dataset directory '{target_path}'")
            try: 
                os.mkdir(target_path)
            except Exception as error:
                dataset['status'] = 'ERROR'
                dataset['state']['processing_state'] = 'CannotCreateLocation'
                response.error(f"Unable to create dataset location {target_path}: {error}", error_code=dataset['state']['processing_state'])
                return response

        #### Double check if the location exists
        if not os.path.exists(target_path):
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'LocationNotCreated'
            response.error(f"Dataset location was not created in {target_path}", error_code=dataset['state']['processing_state'])
            return response

        #### See if the data subdir exists. If not, create it
        target_path = f"{dataset['metadata']['location']}/data"
        response.debug(f"Check data subdirectory '{target_path}'")
        if not os.path.exists(target_path):
            response.info(f"Create data subdirectory '{target_path}'")
            try:
                os.mkdir(target_path)
            except Exception as error:
                dataset['status'] = 'ERROR'
                dataset['state']['processing_state'] = 'CannotCreateLocation'
                response.error(f"Unable to create dataset location {target_path}: {error}", error_code=dataset['state']['processing_state'])
                return response

        return response


    ###############################################################################################
    def fetch_px_record(self, dataset_id):
        """Fetch the PX record and store it in the specified location
        """

        #### Create the URL for ProteomeCentral
        response = self.response
        response.info(f"Fetching dataset record from ProteomeXchange")

        # Get the dataset handle and set status
        dataset = self.datasets['identifiers'][dataset_id]

        url = f"http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={dataset_id}&outputMode=json"
        response_content = requests.get(url, headers={'accept': 'application/json'})

        #### Examine response
        status_code = response_content.status_code
        if status_code != 200:
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'CannotFetchPXRecord'
            print(response_content)
            response.error(f"Unable to fetch ProteomeCentral record for dataset '{dataset_id}', status_code={status_code}. Not a publicly released dataset?", error_code=dataset['state']['processing_state'])
            return response

        with open(f"{dataset['metadata']['location']}/data/ProteomeXchange.json", "w", encoding="utf-8") as outfile:
            outfile.write(str(response_content.text))


    ###############################################################################################
    def assess_download(self, dataset_id):
        """Assess to see if the desired files are there and downloaded and what needs to be done

        """

        verify_by_curl_continue = False

        response = self.response
        response.info(f"Assess the dataset directory to see what has been done and is to do")

        # Get the dataset handle and set status
        dataset = self.datasets['identifiers'][dataset_id]
        px_data = dataset['metadata']['px_data']

        # Set the mode to either assess (as a check on work that may have been done previously)
        # or verify (to verify that work that was just done has completed)
        mode = 'assess'


        # This code never gets used! It is never in this state!!!!!!!!!!!!!!!!!!!! FIXME
        if dataset['state']['processing_state'] == 'Ready for conversion':
            mode = 'verify'

        # Get the location of the remote FTP folder
        ftp_dir = dataset['metadata']['ftp_location']
        if ftp_dir is None:
            ftp_dir = '????'

        # Check on the dataset manifest
        response.info(f"Checking for PRIDE manifest")

        # If there's no manifest entry in metadata, then create it
        if 'manifest' not in dataset['metadata']:
            dataset['metadata']['manifest'] = { 'status': 'UNKNOWN' }

        # If we're already in a ready state, then there's nothing to do
        print(dataset['metadata']['manifest'])
        if dataset['metadata']['manifest']['status'] == 'READY':
            response.info(f"PRIDE manifest is READY")

        #### If UNAVAILABLE, then wait some more
        elif dataset['metadata']['manifest']['status'] == 'UNAVAILABLE':
            response.info(f"Manifest is is not available at the source")
            pass

        # Otherwise, figure out what to do
        else:

            # If the file is already there, then record it
            destination_filepath = f"{dataset['metadata']['location']}/data/README.txt"
            if os.path.exists(destination_filepath):
                response.info(f"PRIDE manifest (README.txt) file is READY")
                dataset['metadata']['manifest']['status'] = 'READY'
                dataset['metadata']['manifest']['file'] = { 'status': 'READY', 'fileroot': 'README',
                    'filename': 'README.txt', 'full_path': f"{dataset['metadata']['location']}/data/README.txt",
                    'location': f"{dataset['metadata']['location']}/data",
                    'expected_size': None, 'current_size': None, 'uri': f"{ftp_dir}/README.txt",
                    'local_age': None, 'is_complete': True, 'filetype': 'txt' }

            # If the file is not there
            else:

                # If the status is UNKNOWN, then this is the first we'vbe considered it and need to download it
                if dataset['metadata']['manifest']['status'] == 'UNKNOWN':
                    dataset['metadata']['manifest']['file'] = { 'status': 'TODO', 'fileroot': 'README',
                        'filename': 'README.txt', 'full_path': f"{dataset['metadata']['location']}/data/README.txt",
                        'location': f"{dataset['metadata']['location']}/data",
                        'expected_size': None, 'current_size': None, 'uri': f"{ftp_dir}/README.txt",
                        'local_age': None, 'is_complete': False, 'filetype': 'txt' }
                    self.tasks_todo.append( { 'command': 'download_file', 'file_metadata': dataset['metadata']['manifest']['file'] } )
                    dataset['metadata']['manifest']['status'] = 'DOWNLOADING'

                #### If DOWNLOADING, then wait some more
                elif dataset['metadata']['manifest']['status'] == 'DOWNLOADING':
                    if dataset['metadata']['manifest']['file']['status'] == 'UNAVAILABLE':
                        dataset['metadata']['manifest']['status'] = 'UNAVAILABLE'
                    else:
                        response.info(f"Manifest is still not available. Keep waiting.")
                        pass


        # Check on individual runs
        ms_runs = []
        previous_msruns = {}
        missing_msruns = {}
        if 'ms_runs' in dataset['metadata'] and dataset['metadata']['ms_runs'] is not None:
            response.debug(f"There are already previous MS runs imported")
            for msrun_name in dataset['metadata']['ms_runs']:
                have_previous_msruns = True
                previous_msruns[msrun_name] = 1
        else:
            response.debug(f"Creating container for MS runs")
            dataset['metadata']['ms_runs'] = {}


        #### If there's no information on the runs yet, then create it
        if len(previous_msruns) == 0 and len(missing_msruns) == 0:

            #### If there is available information in the PX record datasetFiles
            if 'datasetFiles' in px_data:
                response.info(f"Found datasetFiles in the PX record")
                #### Loop over all the datasetFiles to find the MS Runs
                for dataset_file in px_data['datasetFiles']:
                    if dataset_file['name'] == 'Associated raw file URI':
                        uri = dataset_file['value']
                        match = re.match(r'(.+)/(.+)?$',uri)
                        if match:
                            filename = match.group(2)
                            ms_runs.append({ 'filename': filename, 'uri': uri })
                        else:
                            dataset['status'] = 'ERROR'
                            dataset['state']['processing_state'] = 'FailedFilenameMatch'
                            response.error(f"Unable to get the file file name in uri {uri}", error_code=dataset['state']['processing_state'])
                            return response

            #### If the files are not listed in the PX record, try getting a listing at the source via FTP
            else:
                response.warning(f"Unable to find the datasetFiles in the PX record. Trying to get an FTP listing at the source.")
                ms_runs = self.get_ftp_dir_listing(dataset_id, previous_msruns, missing_msruns, verify_by_curl_continue)
                if ms_runs is None or len(ms_runs) == 0:
                    dataset['status'] = 'ERROR'
                    dataset['state']['processing_state'] = 'DatasetFilesFailedFTPDirListing'
                    return response

                response.warning(f"Obtained a file listing at the source. This is great, but if there are subfolders, this won't work")

            #### Loop over the runs we have and process
            for ms_run in ms_runs:
                filename = ms_run['filename']
                uri = ms_run['uri']
                match = re.match(r'(.+)\.(.+)?$',filename)
                if match:
                    fileroot = match.group(1)
                    if fileroot in dataset['metadata']['ms_runs']:
                        if dataset['metadata']['ms_runs'][fileroot]['raw_file']['status'] == 'READY':
                            response.info(f"MS Run {fileroot} is READY")
                            del previous_msruns[fileroot]
                        else:
                            #response.info(f"MS Run {fileroot} is still downloading")
                            pass
                    else:
                        destination_filepath = f"{dataset['metadata']['location']}/data/{filename}"
                        file_info = { 'status': 'TODO',
                            'fileroot': fileroot,
                            'filename': filename,
                            'full_path': f"{dataset['metadata']['location']}/data/{filename}",
                            'location': f"{dataset['metadata']['location']}/data",
                            'expected_size': None,
                            'current_size': None,
                            'uri': uri,
                            'local_age': None,
                            'is_complete': False,
                            'filetype': match.group(2) }

                        if os.path.exists(destination_filepath):
                            response.info(f"Found MS Run raw file {filename} untracked but already present")
                            file_info['status'] = 'READY'
                            file_info['is_complete'] = True

                        if verify_by_curl_continue:
                            response.info(f"But verify_by_curl_continue is set, so perform a curl continue anyway")
                            file_info['status'] = 'TODO'
                            file_info['is_complete'] = False

                        dataset['metadata']['ms_runs'][fileroot] = {}
                        dataset['metadata']['ms_runs'][fileroot]['raw_file'] = file_info

                        if file_info['status'] == 'TODO':
                            self.tasks_todo.append( { 'command': 'download_file', 'file_metadata': dataset['metadata']['ms_runs'][fileroot]['raw_file'] } )
                            missing_msruns[fileroot] = 1
                            dataset['status'] = 'PROCESSING'
                            dataset['state']['processing_state'] = 'Downloading'

                else:
                    dataset['status'] = 'ERROR'
                    dataset['state']['processing_state'] = 'FailedFileRootMatch'
                    response.error(f"Unable to get the file root for {filename}", error_code=dataset['state']['processing_state'])
                    return response



        if len(previous_msruns) != 0 or len(missing_msruns) != 0:
            if mode == 'assess':
                dataset['state']['processing_state'] == 'Ready to download'
                response.info(f"MS Run files are not present. Need to download")
                return response
            else:
                dataset['status'] = 'ERROR'
                dataset['state']['processing_state'] = 'MissingMSRunFiles'
                #response.error(f"Unable to download some MS Runs", error_code=dataset['state']['processing_state'])
                response.warning(f"Unable to download some MS Runs")
                return response

        # If we got this far, then everything is good
        response.info(f"MS Runs seem all accounted for and in order. Ready for conversion")
        dataset['state']['processing_state'] = 'Ready to convert'
        return response


    ###############################################################################################
    def download_dataset(self, dataset_id):
        """Public method that reads the PX dataset information and checks/updates extracted metadata

        """

        response = self.response
        response.info(f"Preparing to download manifest and MS Run raw files")

        # Get the dataset handle and set status
        dataset = self.datasets['identifiers'][dataset_id]
        dataset['status'] = 'PROCESSING'
        px_data = dataset['metadata']['px_data']

        # If we don't already have the manifest, queue a fetch for that
        #if 'manifest' not in dataset['metadata']:
        #    ftp_dir = None
        #    dataset['metadata']['manifest'] = {}
        #    if 'fullDatasetLinks' in px_data:
        #        for term in px_data['fullDatasetLinks']:
        #            if term['name'] == 'Dataset FTP location':
        #                ftp_dir = term['value']
        #    else:
        #        dataset['status'] = 'ERROR'
        #        dataset['state']['processing_state'] = 'CannotFindfullDatasetLinks'
        #        response.error(f"Unable to find fullDatasetLinks in PX record for {dataset_id}", error_code=dataset['state']['processing_state'])
        #    if ftp_dir is not None:
        #        response.info(f"Queueing manifest (README.txt) for download")
        #        dataset['metadata']['manifest']['file'] = { 'status': 'TODO', 'fileroot': 'README',
        #            'filename': 'README.txt', 'full_path': f"{dataset['metadata']['location']}/data/README.txt",
        #            'location': f"{dataset['metadata']['location']}/data",
        #            'expected_size': None, 'current_size': None, 'uri': f"{ftp_dir}/README.txt",
        #            'local_age': None, 'is_complete': False, 'filetype': 'txt' }
        #        self.tasks_todo.append( { 'command': 'download_file', 'file_metadata': dataset['metadata']['manifest']['file'] } )

        # Check the previous list of ms runs
        have_previous_msruns = False
        previous_msruns = {}
        if 'ms_runs' in dataset['metadata'] and dataset['metadata']['ms_runs'] is not None:
            response.debug(f"Collecting previous MS runs")
            for msrun_name in dataset['metadata']['ms_runs']:
                have_previous_msruns = True
                previous_msruns[msrun_name] = 1
        else:
            response.debug(f"Creating container for MS runs")
            dataset['metadata']['ms_runs'] = {}

        # Check that there are datasetFiles specified
        if 'datasetFiles' in px_data:
            #### Loop over all the datasetFiles to find the MS Runs
            for dataset_file in px_data['datasetFiles']:
                if dataset_file['name'] == 'Associated raw file URI':
                    uri = dataset_file['value']
                    match = re.match(r'(.+)/(.+)?$',uri)
                    if match:
                        filename = match.group(2)
                        match = re.match(r'(.+)\.(.+)?$',filename)
                        if match:
                            fileroot = match.group(1)
                            if fileroot in dataset['metadata']['ms_runs']:
                                del previous_msruns[fileroot]
                            else:
                                response.info(f"Queueing MS Run raw file {filename} for download")
                                dataset['metadata']['ms_runs'][fileroot] = {}
                                dataset['metadata']['ms_runs'][fileroot]['raw_file'] = { 'status': 'TODO', 'fileroot': fileroot,
                                    'filename': filename, 'full_path': f"{dataset['metadata']['location']}/data/{filename}",
                                    'location': f"{dataset['metadata']['location']}/data",
                                    'expected_size': None, 'current_size': None, 'uri': uri,
                                    'local_age': None, 'is_complete': False, 'filetype': match.group(2) }
                                self.tasks_todo.append( { 'command': 'download_file', 'file_metadata': dataset['metadata']['ms_runs'][fileroot]['raw_file'] } )
                                if have_previous_msruns:
                                    response.warning(f"Previous catalog of MS run did not have {fileroot}")

            # Complain if there are leftovers
            #for fileroot in self.metadata['ms_runs']

        dataset['status'] = 'PROCESSING'
        dataset['state']['processing_state'] = 'Downloading'


    ###############################################################################################
    def get_ftp_dir_listing(self, dataset_id, previous_msruns, missing_msruns, verify_by_curl_continue):
        """Public method that reads the contents of the remote FTP directory

        """

        response = self.response
        response.info(f"Checking the contents of the remote FTP directory")

        # Get the dataset handle and set status
        dataset = self.datasets['identifiers'][dataset_id]
        dataset['status'] = 'PROCESSING'
        px_data = dataset['metadata']['px_data']
        ms_runs = []

        # Get the FTP location
        ftp_url = None
        if 'fullDatasetLinks' in px_data:
            for term in px_data['fullDatasetLinks']:
                if term['name'] == 'Dataset FTP location':
                    ftp_url = term['value']
        else:
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'CannotFindfullDatasetLinks'
            response.error(f"Unable to find fullDatasetLinks in PX record for {dataset_id}", error_code=dataset['state']['processing_state'])
            return

        # Try to decompose the FTP location
        response.info(f"Remote FTP location is {ftp_url}")
        match = re.match(r'ftp://(.+?)/(.+)$',ftp_url)
        if match:
            ftp_host = match.group(1)
            ftp_dir = match.group(2)

        else:
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'CannotReadRemoteFTPDir'
            response.error(f"Unable to read the remote FTP directory {ftp_url}", error_code=dataset['state']['processing_state'])
            return

        try:
            response.info(f"Connecting to {ftp_host}")
            ftp_session = FTP(ftp_host)
            response.info(f"Login as anonymous")
            ftp_session.login()
            response.info(f"Set CWD to {ftp_dir}")
            ftp_session.cwd(ftp_dir)
            response.info(f"Executing nlst")
            files = ftp_session.nlst()
        except:
            dataset['status'] = 'ERROR'
            dataset['state']['processing_state'] = 'FailedFTPDirListing'
            response.error(f"Unable to get the dir listing at {ftp_url}", error_code=dataset['state']['processing_state'])
            return

        #### Loop over all the files to guess the MS Runs
        for filename in files:
            if filename.endswith('.raw') or filename.endswith('.RAW'):
                ms_runs.append( { 'filename': filename, 'uri': f"{ftp_url}/{filename}" } )

        return ms_runs


    ###############################################################################################
    def assess_conversion(self, dataset_id):
        """Assess if the RAW files have been converted to mzML

        """

        response = self.response
        dataset = self.datasets['identifiers'][dataset_id]
        response.info(f"Check if mzML files for {dataset_id} have already been created")

        # Set the mode to either assess (as a check on work that may have been done previously)
        # or verify (to verify that work that was just done has completed)
        mode = 'assess'
        if dataset['state']['processing_state'] == 'Ready to compress':
            mode = 'verify'
        if dataset['state']['processing_state'] == 'Ready to assess':
            mode = 'verify'

        # Loop over the MS runs and make a dict of things we might need to do
        ms_runs_to_convert = {}
        ms_runs_to_compress = {}
        for fileroot in dataset['metadata']['ms_runs']:
            ms_runs_to_convert['fileroot'] = 1
            ms_runs_to_compress['fileroot'] = 1

        # Loop over the MS runs and queue a conversion job
        for fileroot in dataset['metadata']['ms_runs']:
            have_mzML_record = False
            have_mzML_gz_record = False
            if 'mzML_file' in dataset['metadata']['ms_runs'][fileroot]:
                have_mzML_record = True
            if 'mzML_gz_file' in dataset['metadata']['ms_runs'][fileroot]:
                have_mzML_gz_record = True

            have_mzML_file = False
            have_mzML_gz_file = False
            if os.path.exists(f"{dataset['metadata']['location']}/data/{filename}.mzML"):
                have_mzML_file = True
            if os.path.exists(f"{dataset['metadata']['location']}/data/{filename}.mzML.{self.compressed_extension}"):
                have_mzML_gz_file = True

            if not have_mzML_record and have_mzML_file:
                filename = f"{fileroot}.mzML"
                destination_filepath = f"{dataset['metadata']['location']}/data/{filename}"
                response.info(f"Found MS Run raw file {filename}.mzML untracked but already present")
                dataset['metadata']['ms_runs'][fileroot]['mzML_file'] = { 'status': 'READY', 'fileroot': fileroot,
                    'filename': filename, 'full_path': f"{dataset['metadata']['location']}/data/{filename}",
                    'location': f"{dataset['metadata']['location']}/data",
                    'expected_size': None, 'current_size': None,
                    'local_age': None, 'is_complete': True, 'filetype': 'mzML' }
                have_mzML_record = True

            if not have_mzML_gz_record and have_mzML_gz_file:
                filename = f"{fileroot}.mzML.{self.compressed_extension}"
                destination_filepath = f"{dataset['metadata']['location']}/data/{filename}"
                response.info(f"Found MS Run raw file {filename}.mzML.{self.compressed_extension} untracked but already present")
                dataset['metadata']['ms_runs'][fileroot]['mzML_gz_file'] = { 'status': 'READY', 'fileroot': fileroot,
                    'filename': filename, 'full_path': f"{dataset['metadata']['location']}/data/{filename}",
                    'location': f"{dataset['metadata']['location']}/data",
                    'expected_size': None, 'current_size': None,
                    'local_age': None, 'is_complete': True, 'filetype': 'mzML' }
                have_mzML_gz_record = True

            if have_mzML_gz_record:
                if dataset['metadata']['ms_runs'][fileroot]['mzML_gz_file']['status'] == 'READY':
                    del ms_runs_to_convert['fileroot']
                    del ms_runs_to_compress['fileroot']
                else:
                    pass
            elif have_mzML_record:
                if dataset['metadata']['ms_runs'][fileroot]['mzML_file']['status'] == 'READY':
                    del ms_runs_to_compress['fileroot']
                else:
                    pass
            else:
                pass


        if dataset['state']['processing_state'] == 'Ready to convert':
            if mode == 'assess':
                if len(ms_runs_to_convert) == 0:
                    dataset['state']['processing_state'] == 'Ready to compress'
                    return
                else:
                    # stay ready to convert and the conversion will start
                    return
            else:
                response.warning(f"I should not be here in Ready to convert and not assess")
                return
        elif dataset['state']['processing_state'] == 'Converting':
            if mode == 'assess':
                if len(ms_runs_to_convert) == 0:
                    dataset['state']['processing_state'] == 'Ready to compress'
                    return
                else:
                    response.warning(f"Still waiting for conversion to complete")
                    return
            else:
                response.warning(f"I should not be here in Ready to convert and not assess")
                return
        elif dataset['state']['processing_state'] == 'Ready to compress':
            if mode == 'assess':
                if len(ms_runs_to_convert) == 0:
                    dataset['state']['processing_state'] == 'Ready to compress'
                    return
                else:
                    response.warning(f"Still waiting for conversion to complete")
                    return
            else:
                response.warning(f"I should not be here in Ready to convert and not assess")
                return


                # Continue here with assess conversion! FIXME


    ###############################################################################################
    def assess_mzML(self, dataset_id):
        """Check if there is an mzML file and if not queue a conversion

        """

        response = self.response
        response.info(f"Convert to mzML")

        # Get the dataset handle and set status
        dataset = self.datasets['identifiers'][dataset_id]
        dataset['status'] = 'PROCESSING'

        # Loop over the MS runs and queue a conversion job
        for fileroot in dataset['metadata']['ms_runs']:
            have_mzML = False
            have_mzML_gz = False
            if 'mzML_file' in dataset['metadata']['ms_runs'][fileroot]: have_mzML = True
            if 'mzML_gz_file' in dataset['metadata']['ms_runs'][fileroot]: have_mzML_gz = True

            # If we don't have either uncompressed or compressed, queue up the conversion
            if not have_mzML and not have_mzML_gz:
                response.info(f"Queuing conversion to mzML for MS run {fileroot}")
                filename = f"{fileroot}.mzML"
                dataset['metadata']['ms_runs'][fileroot]['mzML_file'] = { 'status': 'TODO', 'fileroot': fileroot,
                    'filename': filename, 'full_path': f"{dataset['metadata']['location']}/data/{filename}",
                    'location': f"{dataset['metadata']['location']}/data",
                    'expected_size': None, 'current_size': None,
                    'local_age': None, 'is_complete': False, 'filetype': 'mzML' }
                self.tasks_todo.append( { 'command': 'convert_to_mzML', 'file_metadata': dataset['metadata']['ms_runs'][fileroot]['mzML_file'] } )

            # Or if we have the mzMLs but not mzML.gz, then queue the READY ones


    ###############################################################################################
    def show(self, level='high'):
        """Public method that shows the current status of the processor

        """

        response = self.response
        #response.debug(f"Show DatasetProcessor status")

        # Get a listing of available datasets and show their status
        dataset_ids = self.datasets['identifiers']
        buffer = ''
        buffer += f"DatasetProcessor: status: {self.status}, state: {self.state}\n"
        if len(dataset_ids) > 0:
            buffer += '  - Datasets:\n'
            for dataset_id in dataset_ids:
                status = self.datasets['identifiers'][dataset_id]['status']
                n_files = 0
                if 'ms_runs' in self.datasets['identifiers'][dataset_id]['metadata']:
                    n_files = len(self.datasets['identifiers'][dataset_id]['metadata']['ms_runs'])
                buffer += f"      {dataset_id} - {status} - {n_files} files\n"
        else:
            buffer += '  - No datasets'

        # Show some other key things
        buffer += f"  - Tasks to do: {len(self.tasks_todo)}"


        return buffer




##########################################################################################
def main():

    # Parse command line options
    import argparse
    argparser = argparse.ArgumentParser(description='Command line interface to the Response class. Runs tests when run without the --example=N flag')
    argparser.add_argument('--verbose', action='count', help='If set, print out messages to STDERR as they are generated' )
    argparser.add_argument('dataset_id', type=str, help='PXDnnnnnn identifier to process')
    params = argparser.parse_args()

    # Create our DatasetProcessor
    processor = DatasetProcessor()
    response = Response()

    # Set verbosity
    if params.verbose is not None:
        Response.output = 'STDERR'
        processor.verbose = 1

    # Set the dataset_id
    processor.add_dataset(params.dataset_id)

    # Process the dataset
    result = processor.process()
    if result.status != 'OK':
        response.merge(result)
        print(response.show(level=Response.DEBUG))
        return response

    response.merge(result)
    print(response.show(level=Response.DEBUG))

    print(json.dumps(processor.datasets,sort_keys=True,indent=2))
    print(json.dumps(processor.tasks_todo,sort_keys=True,indent=2))

    return response


if __name__ == "__main__": main()
