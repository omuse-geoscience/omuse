from omuse import units
def partition_forcing_file(forcing_file):
    pass


def add_forcing_to_ext_file(partitioned_file_paths, ext_file_path):
    pass


def ext_force_file_creator(ext_file_path, forcing_file_list, number_of_workers):
    for forcing_file in forcing_file_list:
        partitioned_file_paths = partition_forcing_file(forcing_file)
        add_forcing_to_ext_file(partitioned_file_paths, ext_file_path)
