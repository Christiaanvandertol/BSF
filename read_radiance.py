import pandas as pd
from glob import glob

files = glob('*.csv') # gather all *.csv files in the current working folder

for file in files: # loop over all *.csv files
    f = open(file)
    lines = f.readlines() # fully read the file in memory
    f.close() # release the file resource; it is no longer needed
    lines = [line.strip() for line in lines] # remove the newline character from each line
    header = lines[0]
    if len(header.split(',')) > len(header.split(';')): # detect the csv delimiter: ; or ,
        delimiter = ','
    else:
        delimiter = ';'
    header = header.split(delimiter)
    # Extract times and bands to get an indication of the content
    times = [time.strip('"').lstrip('X').replace('_',':') for time in header[1:]]
    bands = [line.split(delimiter)[0] for line in lines[1:]]
    # Print some metadata
    print(file)
    print('Nr times:', len(times))
    print('First and last time:', times[0], times[-1])
    print('Nr bands:', len(bands))
    print('First and last band:', bands[0], bands[-1])
    # Transform the data to a pandas DataFrame
    cols = [header[0].strip('"')] + times
    data = [line.split(delimiter) for line in lines[1:]]
    df = pd.DataFrame(data, columns = cols)
    df = df.apply(pd.to_numeric, errors='coerce') # convert all data from string to float; all unreadable fields (e.g. #N/D) are converted to pythonic NaN
    df = df.set_index(df.columns[0]) # mark the first column (the spectra) as the index
    # All data is now in 'df', which is a pandas DataFrame object
    # pandas is a python library that allows Excel-like calculations on 2D tables
    # The 'df' object now has all data as float.
    # The df columns are named, by the corresponding time (the name is still a string)
    # The df index is the spectral bands
    # This way, df can easily be transposed if needed
    # Do something with df ....
    print()
