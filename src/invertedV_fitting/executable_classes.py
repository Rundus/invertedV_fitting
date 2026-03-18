import json
import os
from src.invertedV_fitting.user_toggles.user_toggles import UserToggles

class ExecutableClasses:

    def generate_run_JSON(self):

        # Check if this specific run has been made
        if not os.path.exists(UserToggles.run_folder_path):
            os.makedirs(UserToggles.run_folder_path)

        # Get all the values from UserInputs
        members = [attr for attr in dir(UserToggles) if not callable(getattr(UserToggles, attr)) and not attr.startswith("__")]
        values = [getattr(UserToggles, attr) for attr in dir(UserToggles) if not callable(getattr(UserToggles, attr)) and not attr.startswith("__")]
        config_dict = {
            f'{mVal}':f'{vVal}' for mVal,vVal in zip(members,values)
        }

        # JSON I/O
        outpath = f'{UserToggles.run_folder_path}/run_config.json'
        with open(outpath, 'w') as outfile:
            json.dump(config_dict, outfile, indent=3)

        return
