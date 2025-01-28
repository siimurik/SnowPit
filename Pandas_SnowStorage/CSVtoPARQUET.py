import pandas as pd
import chardet

input_filename = 'Surface_Melt_Rate_Data.csv'
#with open('C:\\Users\\sipuga\\Documents\\SnowStorageSolvers\\Python_based\\Snow_Storage_Data.csv', 'rb') as f:
with open(input_filename, 'rb') as f:
    result = chardet.detect(f.read())
    encoding = result['encoding']

print(encoding)

# Drop the first row of the CSV file
df = pd.read_csv(input_filename, skiprows=1, encoding=encoding)

# Convert the DataFrame to a Parquet file
output_filename = 'SurfaceMeltRateDATA.parquet'
df.to_parquet(output_filename)    # PyArrow and fastparquet are necessary

print(f"CSV data has been successfully converted to Parquet format and saved as '{output_filename}'.")