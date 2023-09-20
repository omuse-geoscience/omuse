import xarray as xr
import numpy as np
from pathlib import Path

def partition_forcing_file(forcing_file_path, num_parts, output_dir):
    """
    Partitions a NetCFD `forcing_file` into `num_parts` equal parts.
    """
    forcing_file_path = Path(forcing_file_path)
    ds = xr.open_dataset(forcing_file_path)
    longitudes = ds.coords['longitude'].values
    lon_part = np.floor(len(longitudes)/num_parts)
    for i in range(num_parts):
        lon_start = int(i*lon_part)
        lon_end = int((i+1)*lon_part)
        lon_range = longitudes[lon_start:lon_end]
        if i == range(num_parts)[-1]:
            lon_range = longitudes[lon_start:]

        output_path = Path(output_dir) / f"{forcing_file_path.stem}_part{i}.nc"
        print(ds.sel(longitude=lon_range))
        print(output_path)
        ds.to_netcdf(output_path)


def add_forcing_to_ext_file(partitioned_file_paths, operators, ext_file_path):
    """
    Adds `file_paths` as external forcings to the Delft3D .ext file.
    """
    pass


def partition_and_add_forcing(ext_file_path, forcing_file_list, num_parts):
    for forcing_file_path in forcing_file_list:
        partitioned_file_paths = partition_forcing_file(forcing_file_path, num_parts)
        add_forcing_to_ext_file(partitioned_file_paths, ext_file_path)
