#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import datetime
import os
import json
import re
import time
import subprocess

from response import Response
from dataset_processor import DatasetProcessor


class AutomationAgent:

    # Class variables


    # Constructor
    def __init__(self):
        self.status = 'OK'
        self.verbose = None
        self.response = None
        self.config = None
        self.state = { 'status': 'Starting', 'command_pointer': 0 }
        self.tasks_state = { 'previous_show_buffer': '' }
        self.jobs = { }
        self.job_control = { 'job_index': 1, 'n_running_jobs': 0, 'n_jobs': 0, 'n_running_jobs_by_type': {} }
        self.dataset_processor = DatasetProcessor()
        self.start_directory = os.getcwd()

    # Destructor
    def __del__(self):
        #### Remove our PID file
        pid_file = self.start_directory+"/PID"
        if os.path.exists(pid_file):
            os.remove(pid_file)

        # Remove our STOP file
        stop_file = self.start_directory+"/STOP"
        if os.path.exists(stop_file):
            os.remove(stop_file)


    # Start the agent
    def start(self, startstop=None):
        """Public method that starts the agent

        :param startstop: If not None, then start everything, show status, and exit.
        :type startstop: int
        """

        # Create a response object
        response = Response()
        self.response = response

        # Read the config file and prepare the current state
        self.configure()
        self.prepare_state()
        if self.response.status != 'OK':
            print(self.show(level='full'))
            return

        # If startstop is set, then show our status and exit
        if startstop is not None:
            self.stop()
            print(self.show(level='full'))
            return

        # Otherwise run
        self.run()
        self.stop()
        print(self.show(level='full'))


    # Prepare the agent configuration parameters
    def configure(self):
        """Public method that prepares a default configuration and reads the agent's config file

        """

        self.response.debug(f"Setting up agent configuration")
        self.state['status'] = 'Configuring'

        # Set up a default configuration
        self.config = {
            'sleep_interval': 3,
            'heartbeat_interval': 60,
            'data_path': "/proteomics/peptideatlas2/archive/Arabidopsis",
            'max_running_jobs': 2,
            'max_running_jobs_by_type': { 'download': 2 },
        }

        # Try to find a config file and read it
        config_file = "./agent_config.json"
        if os.path.exists(config_file):
            try:
                with open(config_file) as infile:
                    input_config = json.load(infile)
            except Exception as error:
                self.response.error(f"Error reading config file {config_file} - {error}", error_code='ConfigFileReadError')
                return
        else:
            self.response.warning(f"Did not find a local config file {config_file}. Creating one with current defaults")
            with open(config_file,'w') as outfile:
                outfile.write(json.dumps(self.config, indent=4, sort_keys=True)+'\n')
            return

        self.response.info(f"Read agent config file {config_file}")
        for key,value in input_config.items():
            if key in self.config:
                self.config[key] = value
            else:
                self.response.error(f"Local config file has unrecognized key {key} that is not supported. Check spelling", error_code='ConfigFileKeyError')
                return

        # Set the data processors data path, too
        self.dataset_processor.base_dir = self.config['data_path']


    # Prepare the agent state
    def prepare_state(self):
        """Public method that prepares basic low-level state of agent

        """

        self.response.debug(f"Preparing the agent state")
        self.state['status'] = 'Preparing'

        # See if there is a PID file
        pid_file = self.start_directory+"/PID"
        if os.path.exists(pid_file):
            self.response.error(f"There is already a PID file {pid_file}. Another instance is running.", error_code='PIDFileAlreadyExists')
            return

        # See if there is a STOP file
        stop_file = self.start_directory+"/STOP"
        if os.path.exists(stop_file):
            try:
                os.remove(stop_file)
            except Exception as error:
                self.response.error(f"Error trying to delete STOP signal file {stop_file} - {error}", error_code='CannotDeleteSTOPFile')
                return

        # Create the PID file of the current process
        self.state['pid'] = os.getpid()
        try:
            with open(pid_file,'w+') as outfile:
                outfile.write(f"{self.state['pid']}\n")
        except Exception as error:
            self.response.error(f"Unable to write PID file {pid_file}. Write permission error? - {error}", error_code='CannotWritePIDFile')
            return

        # Read the command pointer file
        self.read_command_pointer_file()

        # Verify the data path
        self.verify_data_path()


    # Read or prepare the command pointer file
    def verify_data_path(self):
        """Public method that verifies the configured data path

        """

        self.response.debug(f"Verifying the data_path")

        # Try to find a pointer file and read it
        data_path = self.config['data_path']
        if not os.path.exists(data_path):
            self.response.error(f"Specified data_path {data_path} is missing", error_code='MissingDataPath')
            return


    # Read or prepare the command pointer file
    def read_command_pointer_file(self):
        """Public method that reads the agent's command pointer file

        """

        self.response.debug(f"Reading command pointer file")
        self.state['status'] = 'Configuring'

        # Try to find a pointer file and read it
        pointer_file = os.path.dirname(os.path.abspath(__file__))+"/agent_commands.pointer"
        if os.path.exists(pointer_file):
            try:
                with open(pointer_file) as infile:
                    for line in infile:
                        line = line.rstrip()
                        self.state['command_pointer'] = int(line)
            except Exception as error:
                self.response.error(f"Error reading command pointer file {pointer_file} - {error}", error_code='PointerFileReadError')
                return
        else:
            self.update_command_pointer_file()
            return


    # Update the command pointer file
    def update_command_pointer_file(self):
        """Public method that updates the agent's command pointer file

        """

        self.response.debug(f"Updating command pointer file")

        # Try to find a pointer file and read it
        pointer_file = os.path.dirname(os.path.abspath(__file__))+"/agent_commands.pointer"
        try:
            with open(pointer_file,'w') as outfile:
                outfile.write(f"{self.state['command_pointer']}\n")
        except Exception as error:
            self.response.error(f"Error writing command pointer file {pointer_file} - {error}", error_code='PointerFileWriteError')
            return


    # Read the next command from the command file
    def read_command(self):
        """Public method that reads the next command from the agent's command file

        """

        #self.response.debug(f"Reading command file")

        # Try to find a pointer file and read it
        command_file = os.path.dirname(os.path.abspath(__file__))+"/agent_commands.txt"
        if os.path.exists(command_file):
            try:
                with open(command_file) as infile:
                    pointer = 0
                    for line in infile:
                        pointer += len(line)
                        line = line.rstrip()
                        if pointer > self.state['command_pointer']:
                            self.state['command_pointer'] = pointer
                            self.update_command_pointer_file()
                            return line
            except Exception as error:
                self.response.error(f"Error reading command file {command_file} - {error}", error_code='CommandFileReadError')
                return
        else:
            try:
                with open(command_file,'w'):
                    pass
            except Exception as error:
                self.response.error(f"Error writing command file {command_file}", error_code='CommandFileCreateError')
                return



    # Return a text summary of the current state of the Response
    def show(self, level='compact'):
        """Public method that returns a string buffer of a nice plain text rendering of the state of the agent.

        :param level: Level of detail to provide: 'full' or 'compact'.
        :type level: str
        :return: A string buffer (with newlines) suitable for plain-text printing
        :rtype: str
        """
        if level == 'compact':
            buffer = f"Status: {self.status}, State: {self.state['status']}\n"
        else:
            buffer = f"Automation agent: status: {self.status}, state: {self.state['status']}\n"
            buffer += self.response.show(level=self.response.DEBUG)
        return buffer


    # Run the agent in a loop
    def run(self):
        """Public method that runs the agent in an endless loop until a STOP file is seen

        """

        self.response.info(f"Agent entering run mode")
        self.state['status'] = 'Running'

        #### Define the STOP file to watch for
        stop_file = self.start_directory+"/STOP"
        heartbeat_time = 0

        while 1:

            # Run the main task of the agent
            self.main_task()
            if self.response.status != 'OK':
                print(self.response.show(level=Response.DEBUG))
                return self.response

            # Sleep for the prescribed interval
            #self.response.debug(f"Sleeping {self.config['sleep_interval']} seconds")
            time.sleep(self.config['sleep_interval'])

            # If the heartbeat time is reached, then send a message
            heartbeat_time += self.config['sleep_interval']
            if heartbeat_time > self.config['heartbeat_interval']:
                self.response.info(f"Agent is alive")
                heartbeat_time = 0

            # Check on the STOP file
            if os.path.exists(stop_file):
                self.response.info(f"Detected agent STOP file. Shutting down.")
                break


    # Stop the agent
    def stop(self):
        """Public method that performs cleanup tasks for the agent

        """

        self.response.info(f"Stopping agent")
        self.state['status'] = 'Stopping'

        #### Remove our PID file
        pid_file = self.start_directory+"/PID"
        if os.path.exists(pid_file):
            os.remove(pid_file)

        # Remove our STOP file
        stop_file = self.start_directory+"/STOP"
        if os.path.exists(stop_file):
            os.remove(stop_file)


    # Execute a presumed quick command that will block the current process
    def execute_blocking(self, command, location):
        """Public method that executes a provided command

        """

        self.response.debug(f"Executing command '{command}'")
        cwd = os.getcwd()
        os.chdir(location)
        result = subprocess.run(command, capture_output=True)
        os.chdir(cwd)
        eprint(result)


    # Add a new job to the queue
    def add_job(self, job):
        """Public method that adds a job to the queue

        """

        job_index = self.job_control['job_index']
        self.response.debug(f"Adding new job {job_index} to the queue")
        self.jobs[job_index] = job
        self.job_control['n_jobs'] += 1
        self.job_control['job_index'] += 1


    # Execute a presumed quick command that will block the current process
    def show_jobs(self):
        """Public method that summarizes the jobs in the queue

        """

        eprint(f"  n_jobs={self.job_control['n_jobs']}, n_running_jobs={self.job_control['n_running_jobs']}")
        for job_id,job in self.jobs.items():
            output_status = 'none'
            if 'expected_output_file' in job:
                if os.path.exists(job['expected_output_file']):
                    mtime = os.path.getmtime(job['expected_output_file'])
                    now = time.time()
                    age = age = int(now-mtime)
                    output_status = f"file age: {age} s"
            if job['status'] == 'run':
                eprint(f"    - {job_id}: status={job['status']}, type={job['status']}, cmd={' '.join(job['args'])}, output: {output_status}")


    # See if there are any queued jobs that can be launched
    def launch_jobs(self):
        """Public method that determines which waiting jobs to launch, if any

        """

        # If there is nothing to do, return immediately
        if self.job_control['n_jobs'] == 0:
            return
        if self.job_control['n_running_jobs'] >= self.config['max_running_jobs']:
            return

        # Loop through the jobs to determine priority
        urgent_jobs_list = []
        other_jobs_list = []
        for job_id,job in self.jobs.items():
            if job['status'] == 'run':
                continue
            elif job['status'] == 'redo':
                urgent_jobs_list.append(job_id)
            else:
                other_jobs_list.append(job_id)

        jobs_list = urgent_jobs_list
        for job_id in other_jobs_list:
            jobs_list.append(job_id)

        # Loop through the jobs and see if any can be started
        for job_id in jobs_list:
            job = self.jobs[job_id]
            job_type = job['type']
            if job_type not in self.job_control['n_running_jobs_by_type']:
                self.job_control['n_running_jobs_by_type'][job_type] = 0
            if self.job_control['n_running_jobs_by_type'][job_type] >= self.config['max_running_jobs_by_type'][job_type]:
                continue

            self.response.info(f"Launching job '{job_id}'")
            cwd = os.getcwd()
            os.chdir(job['location'])
            proc = subprocess.Popen(job['args'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            job['handle'] = proc
            os.chdir(cwd)
            job['pid'] = proc.pid
            job['status'] = 'run'
            job['launch_timestamp'] = time.time()
            self.job_control['n_running_jobs'] += 1
            self.job_control['n_running_jobs_by_type'][job_type] += 1


    # See if there are any running jobs that are done
    def poll_jobs(self):
        """Public method that reviews the list of running jobs for anything to complete

        """

        # If there is nothing to do, return immediately
        if self.job_control['n_running_jobs'] == 0:
            return

        # Loop through the jobs and poll any running ones
        job_ids_to_delete = []
        jobs_to_restart = []
        for job_id,job in self.jobs.items():
            if job['status'] != 'run':
                continue

            proc = job['handle']
            return_code = proc.poll()

            #### If still running, try to determine if still productive
            if return_code is None:
                if 'retry_staleness' in job and 'expected_output_file' in job and job['retry_staleness'] > 0:
                    if os.path.exists(job['expected_output_file']):
                        mtime = os.path.getmtime(job['expected_output_file'])
                        now = time.time()
                        file_age = int(now - mtime)
                        job_age = int(now - job['launch_timestamp'])
                        if file_age > job['retry_staleness'] and job_age > job['retry_staleness']:
                            self.response.warning(f"Maximum staleness {job['retry_staleness']} reached for file {job['expected_output_file']}. Kill and restart.")
                            if 'n_retries' not in job: job['n_retries'] = 0
                            if 'max_retries' not in job: job['max_retries'] = 10
                            job['handle'].terminate()
                            self.job_control['n_jobs'] -= 1
                            self.job_control['n_running_jobs'] -= 1
                            self.job_control['n_running_jobs_by_type'][job['type']] -= 1
                            time.sleep(5)
                            job_ids_to_delete.append(job_id)
                            #os.rename(job['expected_output_file'],f"{job['expected_output_file']}-{job['n_retries']}")

                            #### See if we should restart it
                            job['n_retries'] += 1
                            if job['n_retries'] > job['max_retries']:
                                self.response.error(f"Max retries {job['max_retries']} reached for file {job['expected_output_file']}", error_code='MaxRetriesReached')
                            else:
                                jobs_to_restart.append(job)

                # Back to the top if it was still running
                continue

            # If the job has exited, close it off
            self.response.info(f"Job {job_id} is complete with return code {return_code}")
            job_ids_to_delete.append(job_id)
            self.job_control['n_jobs'] -= 1
            self.job_control['n_running_jobs'] -= 1
            self.job_control['n_running_jobs_by_type'][job['type']] -= 1

            # If this was a file download job, update some key information
            if job['type'] == 'download':
                if os.path.exists(job['file_handle']['full_path']):
                    job['file_handle']['status'] = 'READY'
                    job['file_handle']['is_complete'] = True
                    job['file_handle']['current_size'] = os.path.getsize(job['file_handle']['full_path'])
                else:
                    self.response.warning(f"File download was supposedly complete, the but file isn't there! Requeue it.")
                    jobs_to_restart.append(job)

        # Delete any finished jobs from the jobs queue. Must be done here at the end
        # since it is not permissable to delete inside the above loop
        for job_id in job_ids_to_delete:
            del self.jobs[job_id]

        # If there are any jobs to restart, put them back in the queue
        for new_job in jobs_to_restart:
            new_job['pid'] = None
            new_job['handle'] = None
            new_job['status'] = 'redo'
            self.add_job(new_job)


    # Queue the tasks from the worker
    def queue_tasks(self):
        """Loop through the tasks called out by the worker and queue them to execution

        """

        # Loop over the requested tasks and add them to the scheduler
        for task in self.dataset_processor.tasks_todo:

            # Process command download_file
            if task['command'] == 'download_file':
                uri = task['file_metadata']['uri']
                location = task['file_metadata']['location']
                expected_output_file = task['file_metadata']['full_path']
                new_job = { 'pid': None, 'type': 'download', 'args': [ "curl", "-R", "-O", "-C", "-", uri ], 'retry_staleness': 30,
                    'n_retries': 0, 'max_retries': 10, 'file_handle': task['file_metadata'],
                    'location': location, 'status': 'qw', 'handle': None, 'expected_output_file': expected_output_file }
                self.add_job(new_job)

            # Process command convert_to_mzML
            elif task['command'] == 'convert_to_mzML':
                location = task['file_metadata']['location']
                full_path = task['file_metadata']['full_path']
                expected_output_file = re.sub(r'\.raw$','.mzML',full_path, flags=re.I)
                new_job = { 'pid': None, 'type': 'download',
                    'args': [ "C:/Users/ericd/Documents/Software/Thermo/ThermoRawFileParser/ThermoRawFileParser", "-m", "0", "-f", "2", "-i", full_path ],
                    'retry_staleness': 30,
                    'n_retries': 0, 'max_retries': 10, 'file_handle': task['file_metadata'],
                    'location': location, 'status': 'qw', 'handle': None, 'expected_output_file': expected_output_file }
                self.add_job(new_job)

            # If not handled yet, then this is unrecognized
            else:
                self.response.error(f"Unrecognized task from worker: {task['command']}", error_code='UnrecognizedTaskName')

        # Clear out the tasks
        self.dataset_processor.tasks_todo = []


    ##########################################################################################

    # Main task of the agent. Is called every sleep_interval seconds
    def main_task(self):
        """Public method that runs the main task of agent. Override with the true work

        """

        # First look for new work to put in the queue
        command = self.read_command()
        if command is not None:
            self.response.info(f"Received command '{command}'")
            got_match = False
            match = re.match(r'get\s+(.+)$',command)
            if match:
                new_job = { 'pid': None, 'type': 'download', 'args': [ "curl", "-R", "-O", match.group(1) ],
                    'location': self.config['data_path'], 'status': 'qw', 'handle': None }
                self.add_job(new_job)
                got_match = True

            match = re.match(r'add_dataset\s+(.+)$',command)
            if got_match is False and match:
                self.dataset_processor.add_dataset(match.group(1))
                got_match = True

            else:
                self.response.warning(f"Unable to interpret received command '{command}'")

        # Run the DatasetProcessor for a cycle
        result = self.dataset_processor.process()
        if result.status != 'OK':
            self.response.merge(result)
            print(self.response.show(level=Response.DEBUG))
            return self.response

        # Show the DatasetProcessor status
        status_buffer = self.dataset_processor.show(level='high')
        if status_buffer != self.tasks_state['previous_show_buffer']:
            print(status_buffer)
            self.tasks_state['previous_show_buffer'] = status_buffer

        # If there are tasks to perform from the DatasetProcessor, queue them
        if len(self.dataset_processor.tasks_todo) > 0:
            self.queue_tasks()

        # Show jobs
        if self.job_control['n_jobs'] > 0:
            self.show_jobs()

        # Then see if there is something to launch
        self.launch_jobs()

        # Check in on running jobs
        self.poll_jobs()




##########################################################################################
def main():

    # Parse command line options
    import argparse
    argparser = argparse.ArgumentParser(description='Command line interface to the Response class. Runs tests when run without the --example=N flag')
    argparser.add_argument('--verbose', action='count', help='If set, print out messages to STDERR as they are generated' )
    argparser.add_argument('--startstop', action='count', help='If set, perform all startup work and describe the current state and exit' )
    params = argparser.parse_args()

    # Create our AutomationAgent
    agent = AutomationAgent()

    # Set verbosity
    if params.verbose is not None:
        Response.output = 'STDERR'
        agent.verbose = 1

    # Run the agent
    agent.start(startstop=params.startstop)


if __name__ == "__main__": main()
