import json
import os
from src.invertedV_fitting.user_toggles import FileToggles

class ExecutableClasses:

    def generate_run_directories(self,dict_executable):

        # Check if the mission name directory exists, if not create it
        check_path = f'{FileToggles.PATH_TO_FOLDER}{FileToggles.mission_name}'
        if not os.path.exists(check_path):
            os.makedirs(FileToggles.PATH_TO_FOLDER)

        # Check if the instrument, date-specific folder exists, if not create it
        if not os.path.exists(FileToggles.RUN_PATH):
            os.makedirs(FileToggles.RUN_PATH)

        # Check if all the subroutine folders exist
        folders_paths = [
            'fit_data',
            'primary_beam_fits',
            'primary_beam_fits/plots',
            'backscatter'
        ]

        for path in folders_paths:
            check_path = f'{FileToggles.RUN_PATH}/{path}/'
            if not os.path.exists(check_path):
                os.makedirs(check_path)

    # ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    def retrieve_ace_data(self,year, month, day, spacecraft, level, local_dir=None):

        from bs4 import BeautifulSoup
        import requests
        import sys

        """
        This routine retrieves ACE data from the UIowa server.

        Required inputs are year, month, and day of the file you wish to retreive, as well as the spacecraft
        name ('1' or '2') and data product level (l2 or l3).

        This routine will download the most recent file for that date to a local subdirectory, which
        is indicated as the "local_dir" variable. This is defaulted to ./data/TS1(2)/ACE/LL/YYYY/MM.
        """

        level = level.lower()
        base_url = f'https://tracers-portal.physics.uiowa.edu/teams/flight/ACE/ts{spacecraft}/{level}'
        date_url = f'{base_url}/{year}/{month}'

        page = requests.get(date_url, auth=('tracers-sot', 'SciOpsTeamFlight!'))
        data = page.text
        soup = BeautifulSoup(data, "html.parser")
        ds = f'{year}{month}{day}'
        all_strings = soup.find_all('a',href=True)
        idx = []
        for i in range(len(all_strings)):
            string_name = all_strings[i].get('href')
            if ds in string_name:
                idx.append(i)
        if len(idx) > 0:
            day_file = all_strings[idx[-1]].get('href')
            sys.stdout.write('\nDownloading ' + f'{day_file}' + '\n')
            file_url_path = date_url + '/' + day_file
            local_file_path = local_dir + '/' + day_file
            r = requests.get(file_url_path, auth=('tracers-sot', 'SciOpsTeamFlight!'))
            with open(local_file_path, 'wb') as df:
                df.write(r.content)
        else:
            ymd = f'{year}-{month}-{day}'
            raise Exception(f"No ACE {spacecraft.upper()} files for {ymd}")

        return None


    def load_fit_data(self):

        ################################################################################################
        # check if the fitting data has been downloaded, if not download it
        from glob import glob
        data_files = glob(f'{FileToggles.RUN_PATH}/fit_data/*.cdf*')

        # FILE I/O
        print('Finding local fit data...', end='')

        if len(data_files) == 0:

            if FileToggles.mission_name == 'TRACERS':
                print('Data Files not found.')

                # Download the ACE data and store it
                self.retrieve_ace_data(year=FileToggles.data_year,
                                       month=FileToggles.data_month,
                                       day=FileToggles.data_day,
                                       spacecraft=FileToggles.payload_designator,
                                       level=FileToggles.data_level,
                                       local_dir=f'{FileToggles.RUN_PATH}/fit_data/')
            elif FileToggles.mission_name == 'ACESII':
                raise Exception('Need to Implement ACESII Code!')
        else:
            print(f'Data Found! {data_files[0]}')

    def generate_run_JSON(self, filename='run_config.json'):
        """
        Write every variable from every toggle class in user_toggles.py to a JSON file.

        The JSON is grouped by class name, and each key matches the variable name
        in its class. Classes appear in the order they are defined in the file,
        and variables in the order they are defined in each class.

        Parameters
        ----------
        filename : str
            Name of the output file.

        Returns
        -------
        str
            Full path to the written JSON file.
        """
        import datetime as dt
        import inspect
        import json
        import os
        import numpy as np
        import src.invertedV_fitting.user_toggles as user_toggles

        # Handle types the json module can't serialize on its own
        def json_default(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.generic):  # np.int64, np.float32, np.bool_, etc.
                return obj.item()
            if isinstance(obj, (dt.datetime, dt.date, dt.time)):
                return obj.isoformat()
            if isinstance(obj, dt.timedelta):
                return obj.total_seconds()
            return str(obj)  # last resort so the dump never fails

        output_dir = f'{user_toggles.FileToggles().RUN_PATH}/'
        os.makedirs(output_dir, exist_ok=True)

        config_dict = {}
        for class_name, cls in vars(user_toggles).items():
            # Only classes defined in user_toggles.py, not ones imported into it
            if not inspect.isclass(cls) or cls.__module__ != user_toggles.__name__:
                continue

            config_dict[class_name] = {
                var_name: value
                for var_name, value in vars(cls).items()
                if not var_name.startswith('__')
                   and not callable(value)
                   and not isinstance(value, (staticmethod, classmethod, property))
            }

        outpath = os.path.join(output_dir, filename)
        with open(outpath, 'w') as outfile:
            json.dump(config_dict, outfile, indent=3, default=json_default)

        return outpath
