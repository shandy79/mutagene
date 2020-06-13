import datetime
import yaml
from pathlib import Path


def inputs_to_yml(input_dict):
    # Construct file names here, set datetime
    input_dict_cp = input_dict.copy()

    if isinstance(input_dict_cp['input_files'], list):
        for i in range(len(input_dict_cp['input_files'])):
            if not isinstance(input_dict_cp['input_files'][i], (str, bytes, Path)):
                file_name = input_dict_cp['input_files'][i].name
                input_dict_cp['input_files'][i] = file_name
    else:
        if not isinstance(input_dict_cp['input_files'], (str, bytes, Path)):
            file_name = input_dict_cp['input_files'].name
            input_dict_cp['input_files'] = file_name

    if not isinstance(input_dict_cp['output_file'], (str, bytes, Path)):
        file_name = input_dict_cp['output_file'].name
        input_dict_cp['output_file'] = file_name

    now = datetime.datetime.now()
    input_dict_cp['datetime'] = now
    datetime_str = now.strftime('%Y%m%d_%H%M%S')

    with open(f"mutagene-{input_dict_cp['command']}-{datetime_str}.yml", 'w') as f:
        yaml.dump(input_dict_cp, f)
