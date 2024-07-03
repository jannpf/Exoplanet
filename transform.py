import pandas as pd
from glob import glob
import os

folder = "./data/EXOPLANET_ARCHIVE"
output_folder = f"{folder}_OUT"
files = sorted(glob(folder + "/**.tbl"))


if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for f in files:
    ll = open(f, 'rt').readlines()

    # parse parameter lines
    params = {k.strip(): v.strip().strip('\"\'')
              for k, v in (l[1:].split("=")[:2] for l in ll if l.startswith("\\") and "=" in l)}
    
    star_id = str.replace(params["STAR_ID"], " ", "_")
    print(f"Opened {f} with star_id {star_id}")

    # parse header section
    col_names, types, units = [
        list(map(str.strip, l[1:-2].split("|")))[:3] for l in ll if l.startswith("|")][:3]

    # treat every line starting with a number as data
    data = [list(map(str.strip, l.strip().split()))[:3]
            for l in ll if l.strip()[0].isdigit()]

    # read dataframe
    df = pd.DataFrame(data, columns=col_names)

    # set types according to header
    for c, t in list(zip(col_names, types)):
        df[c] = df[c].astype(t)
    df = df.set_index(col_names[0])

    # file names tend to be more descriptive than the id with MARVELS
    output_file_name = os.path.basename(f) if star_id.startswith("2MASS") else star_id

    # dont merge, as later datasets sometimes contain newer information //TODO check for others
    # # merge existing
    # if os.path.isfile(f"{output_folder}/{output_file_name}.txt"):
    #     df_old = pd.read_table(f"{output_folder}/{output_file_name}.txt", header=None, sep=" ", index_col=0)
    #     if len(df_old.columns) == len(df.columns):
    #         df_old.columns = df.columns
    #         df = df.combine_first(df_old).sort_index()
    #     else:
    #         print(f"Shape mismatch at {output_file_name}, old: {len(df_old.columns)}, new: {len(df.columns)}. Skipping...")
    #         continue

    # write out in std format
    df.to_csv(
        f"{output_folder}/{output_file_name}.txt",
        header=False,
        sep=" ",
        index=True,
        float_format='%.18e',
        mode='w'
    )
