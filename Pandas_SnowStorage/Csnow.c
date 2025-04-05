#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>

typedef struct {
    char*** data;
    int rows;
    int cols;
    int* cols_to_keep;
    int cols_kept;
} CSVData;

// Structure to hold converted data with error information
typedef struct {
    double* values;    // Array of converted values
    bool* is_valid;    // Array indicating valid conversions
    int length;        // Number of elements
} NumericArray;


typedef struct {
    char** strings;  // Array of strings
    int length;      // Number of strings
} StringArray;

// Function prototypes
//const char* detect_encoding(const char* file_path);
//int count_columns(const char* line);
//char** parse_csv_line(const char* line, int max_columns);
//void free_csv_row(char** row, int cols);
//void free_csv_data(CSVData* data);
//CSVData* read_csv_with_encoding(const char* file_path, int* columns_to_keep, int num_cols_to_keep);

const char* detect_encoding(const char* file_path) {
    (void)file_path;  // Explicitly mark as unused
    return "ISO-8859-1";
}

int count_columns(const char* line) {
    int count = 0;
    bool in_quotes = false;
    const char* p = line;

    while (*p) {
        if (*p == '"') {
            in_quotes = !in_quotes;
        } else if (*p == ',' && !in_quotes) {
            count++;
        }
        p++;
    }
    return count + 1; // Columns = separators + 1
}

void free_csv_row(char** row, int cols) {
    if (!row) return;
    for (int i = 0; i < cols; i++) {
        free(row[i]);
    }
    free(row);
}

// Helper function to trim all whitespace including newlines
char* trim_whitespace(char* str) {
    if (!str) return NULL;
    
    // Trim leading space
    while(isspace(*str)) str++;
    
    if (*str == 0) return str; // All spaces
    
    // Trim trailing space
    char* end = str + strlen(str) - 1;
    while(end > str && isspace(*end)) end--;
    
    // Write new null terminator
    *(end+1) = '\0';
    
    return str;
}

char** parse_csv_line(const char* line, int max_columns) {
    char** columns = calloc(max_columns, sizeof(char*));
    if (!columns) return NULL;

    char* buffer = strdup(line);
    if (!buffer) {
        free(columns);
        return NULL;
    }

    // First remove any trailing newline
    char* nl = strchr(buffer, '\n');
    if (nl) *nl = '\0';

    int col = 0;
    bool in_quotes = false;
    char* start = buffer;
    char* ptr = buffer;

    while (*ptr && col < max_columns) {
        if (*ptr == '"') {
            in_quotes = !in_quotes;
            ptr++;
        } else if (*ptr == ',' && !in_quotes) {
            *ptr = '\0';
            
            // Clean the field
            char* field = trim_whitespace(start);
            columns[col++] = strdup(field[0] ? field : " ");
            
            start = ptr + 1;
            ptr++;
        } else {
            ptr++;
        }
    }

    // Handle last column
    if (col < max_columns) {
        char* field = trim_whitespace(start);
        columns[col++] = strdup(field[0] ? field : " ");
    }

    // Fill remaining columns with empty strings if needed
    while (col < max_columns) {
        columns[col++] = strdup(" ");
    }

    free(buffer);
    return columns;
}

CSVData* read_csv_with_encoding(const char* file_path, int* columns_to_keep, int num_cols_to_keep) {
    const char* encoding = detect_encoding(file_path);
    printf("\nAssuming encoding: %s\n", encoding);

    FILE* file = fopen(file_path, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    // First pass to count rows and columns
    int rows = 0;
    int max_columns = 0;
    char line[100000];

    while (fgets(line, sizeof(line), file)) {
        rows++;
        int cols_in_line = count_columns(line);
        if (cols_in_line > max_columns) max_columns = cols_in_line;
    }

    // Handle column selection
    int* cols_kept = NULL;
    int cols_kept_count = 0;

    if (columns_to_keep == NULL || num_cols_to_keep == 0) {
        cols_kept_count = max_columns;
        cols_kept = malloc(cols_kept_count * sizeof(int));
        for (int i = 0; i < cols_kept_count; i++) {
            cols_kept[i] = i;
        }
    } else {
        cols_kept_count = num_cols_to_keep;
        cols_kept = malloc(cols_kept_count * sizeof(int));
        memcpy(cols_kept, columns_to_keep, cols_kept_count * sizeof(int));
    }

    // Allocate main structure
    CSVData* csv_data = calloc(1, sizeof(CSVData));
    if (!csv_data) {
        fclose(file);
        free(cols_kept);
        return NULL;
    }

    csv_data->rows = rows;
    csv_data->cols = max_columns;
    csv_data->cols_to_keep = cols_kept;
    csv_data->cols_kept = cols_kept_count;
    csv_data->data = calloc(rows, sizeof(char**));

    // Second pass to read data
    rewind(file);
    for (int i = 0; i < rows; i++) {
        if (!fgets(line, sizeof(line), file)) break;

        char** all_columns = parse_csv_line(line, max_columns);
        if (!all_columns) {
            // Cleanup
            for (int j = 0; j < i; j++) free_csv_row(csv_data->data[j], csv_data->cols_kept);
            free(csv_data->data);
            free(cols_kept);
            free(csv_data);
            fclose(file);
            return NULL;
        }

        csv_data->data[i] = calloc(cols_kept_count, sizeof(char*));
        for (int j = 0; j < cols_kept_count; j++) {
            int col_idx = cols_kept[j];
            if (col_idx < max_columns && all_columns[col_idx]) {
                csv_data->data[i][j] = strdup(all_columns[col_idx]);
            } else {
                csv_data->data[i][j] = strdup(" ");
            }
        }

        free_csv_row(all_columns, max_columns);
    }

    fclose(file);
    return csv_data;
}

void free_csv_data(CSVData* data) {
    if (!data) return;
    for (int i = 0; i < data->rows; i++) {
        if (data->data[i]) {
            free_csv_row(data->data[i], data->cols_kept);
        }
    }
    free(data->data);
    free(data->cols_to_keep);
    free(data);
}

// Function to print a specific column from the first n rows
void print_column(CSVData* data, int column_index, int max_rows) {
    if (!data || column_index < 0 || column_index >= data->cols_kept) {
        printf("Invalid column index\n");
        return;
    }
    
    int rows_to_print = (max_rows < data->rows) ? max_rows : data->rows;
    for (int i = 0; i < rows_to_print; i++) {
        printf("%s\n", data->data[i][column_index]);
    }
}

// Function to find first empty value in a column
int find_first_empty(CSVData* data, int column_index) {
    if (!data || column_index < 0 || column_index >= data->cols_kept) {
        return -1;
    }
    
    for (int i = 0; i < data->rows; i++) {
        if (strlen(data->data[i][column_index]) == 0 || 
            strcmp(data->data[i][column_index], " ") == 0) {
            return i;
        }
    }
    return -1;
}

// Function to find index of a specific value in a column
int find_index(CSVData* data, int column_index, const char* value) {
    if (!data || column_index < 0 || column_index >= data->cols_kept) {
        return -1;
    }
    
    for (int i = 0; i < data->rows; i++) {
        if (strcmp(data->data[i][column_index], value) == 0) {
            return i;
        }
    }
    return -1;
}

// Function to extract a column as an array of strings
char** extract_column(CSVData* data, int column_index) {
    if (!data || column_index < 0 || column_index >= data->cols_kept) {
        return NULL;
    }
    
    char** column = malloc(data->rows * sizeof(char*));
    if (!column) return NULL;
    
    for (int i = 0; i < data->rows; i++) {
        column[i] = strdup(data->data[i][column_index]);
    }
    
    return column;
}

// Similar to printVec
void print_column_stats(CSVData* data, int column_index, const char* column_name) {
    if (!data || column_index < 0 || column_index >= data->cols_kept) {
        printf("Invalid column\n");
        return;
    }
    
    //printf("\n");
    int length = data->rows;
    
    if (length > 10) {
        for (int i = 0; i < 5; i++) {
            printf("%d\t%s\n", i, data->data[i][column_index]);
        }
        printf("...\n");
        for (int i = length-5; i < length; i++) {
            printf("%d\t%s\n", i, data->data[i][column_index]);
        }
    } else {
        for (int i = 0; i < length; i++) {
            printf("%d\t%s\n", i, data->data[i][column_index]);
        }
    }
    
    printf("Name: %s, Length: %d\n", column_name, length);
}

// Similar to printAny but for a range of rows
void print_data_range(CSVData* data, int start, int end, const char* column_name) {
    if (!data || start < 0 || end >= data->rows || start > end) {
        printf("Invalid range\n");
        return;
    }
    
    printf("\t");
    int length = end - start + 1;
    
    if (length > 10) {
        for (int i = start; i < start+5; i++) {
            for (int j = 0; j < data->cols_kept; j++) {
                printf("%s\t", data->data[i][j]);
            }
            //printf("\n");
        }
        printf("...\n\t");
        for (int i = end-4; i <= end; i++) {
            for (int j = 0; j < data->cols_kept; j++) {
                printf("%s\t", data->data[i][j]);
            }
            //printf("\n");
        }
    } else {
        for (int i = start; i <= end; i++) {
            for (int j = 0; j < data->cols_kept; j++) {
                printf("%s\t", data->data[i][j]);
            }
            //printf("\n");
        }
    }
    
    printf("Name: %s, Length: %d\n", column_name, length);
}



// Function to convert string array to double array
NumericArray convert_string_array_to_double(StringArray* string_arr) {
    NumericArray result;
    result.values = NULL;
    result.is_valid = NULL;
    result.length = 0;
    
    if (!string_arr || string_arr->length == 0) {
        return result;
    }
    
    result.values = malloc(string_arr->length * sizeof(double));
    result.is_valid = malloc(string_arr->length * sizeof(bool));
    if (!result.values || !result.is_valid) {
        free(result.values);
        free(result.is_valid);
        result.length = 0;
        return result;
    }
    
    result.length = string_arr->length;
    
    for (int i = 0; i < string_arr->length; i++) {
        char* cleaned = trim_whitespace(string_arr->strings[i]);
        char* endptr;
        result.values[i] = strtod(cleaned, &endptr);
        
        // Check if conversion was successful
        result.is_valid[i] = (*cleaned != '\0' && *endptr == '\0');
        if (!result.is_valid[i]) {
            fprintf(stderr, "Warning: Could not convert '%s' to double\n", string_arr->strings[i]);
            result.values[i] = 0.0; // Default value for invalid
        }
    }
    
    return result;
}

// Function to convert string array to int array
NumericArray convert_to_int(char** string_array, int length) {
    NumericArray result;
    result.values = malloc(length * sizeof(double)); // Using double to match the struct
    result.is_valid = malloc(length * sizeof(bool));
    result.length = length;
    
    if (!result.values || !result.is_valid) {
        // Handle memory allocation failure
        if (result.values) free(result.values);
        if (result.is_valid) free(result.is_valid);
        result.length = 0;
        return result;
    }
    
    for (int i = 0; i < length; i++) {
        char* endptr;
        long value = strtol(string_array[i], &endptr, 10);
        
        // Check if conversion was successful
        if (endptr == string_array[i] || *endptr != '\0') {
            result.values[i] = 0;  // Default value for invalid conversions
            result.is_valid[i] = false;
            fprintf(stderr, "Warning: Could not convert '%s' to int\n", string_array[i]);
        } else {
            result.values[i] = (double)value; // Store as double to match the struct
            result.is_valid[i] = true;
        }
    }
    
    return result;
}

// Function to free a NumericArray
void free_numeric_array(NumericArray* arr) {
    if (arr) {
        free(arr->values);
        free(arr->is_valid);
        arr->values = NULL;
        arr->is_valid = NULL;
        arr->length = 0;
    }
}

// Initialize an empty StringArray
StringArray create_empty_string_array() {
    StringArray arr;
    arr.strings = NULL;
    arr.length = 0;
    return arr;
}

// Free all memory used by a StringArray
void free_string_array(StringArray* arr) {
    if (arr) {
        for (int i = 0; i < arr->length; i++) {
            free(arr->strings[i]);
        }
        free(arr->strings);
        arr->length = 0;
    }
}

// Function to extract a column from a CSVData range
StringArray extract_column_range(CSVData* data, int column_index, int start, int end) {
    StringArray result = create_empty_string_array();
    
    if (!data || column_index < 0 || column_index >= data->cols_kept || 
        start < 0 || end >= data->rows || start > end) {
        return result;
    }
    
    result.length = end - start + 1;
    result.strings = malloc(result.length * sizeof(char*));
    if (!result.strings) {
        result.length = 0;
        return result;
    }
    
    for (int i = 0; i < result.length; i++) {
        result.strings[i] = strdup(data->data[start + i][column_index]);
        if (!result.strings[i]) {
            // If allocation fails, clean up what we've allocated so far
            for (int j = 0; j < i; j++) {
                free(result.strings[j]);
            }
            free(result.strings);
            result.length = 0;
            return result;
        }
    }
    
    return result;
}

// Function to print a NumericArray (similar to printVec)
void print_numeric_array(NumericArray* array, const char* column_name) {
    if (!array || array->length == 0) {
        printf("Empty or invalid array\n");
        return;
    }
    
    printf("\n");
    int length = array->length;
    
    if (length > 10) {
        for (int i = 0; i < 5; i++) {
            if (array->is_valid[i]) {
                printf("%d\t%.6f\n", i, array->values[i]);
            } else {
                printf("%d\tINVALID\n", i);
            }
        }
        printf("...\n");
        for (int i = length-5; i < length; i++) {
            if (array->is_valid[i]) {
                printf("%d\t%.6f\n", i, array->values[i]);
            } else {
                printf("%d\tINVALID\n", i);
            }
        }
    } else {
        for (int i = 0; i < length; i++) {
            if (array->is_valid[i]) {
                printf("%d\t%.6f\n", i, array->values[i]);
            } else {
                printf("%d\tINVALID\n", i);
            }
        }
    }
    
    printf("Name: %s, Length: %d\n", column_name, length);
}



void print_all_headers(CSVData* data) {
    if (!data || !data->data || data->rows == 0) return;
    
    printf("\nColumn headers:\n");
    for (int i = 0; i < data->cols_kept; i++) {
        printf("%d: %s\n", i, data->data[0][i]);
    }
}

void debug_string_array(StringArray* arr, const char* name) {
    printf("\nDebugging %s (%d items):\n", name, arr->length);
    for (int i = 0; i < (arr->length > 5 ? 5 : arr->length); i++) {
        printf("[%d] Raw: '%s' (len:%zu)\n", i, arr->strings[i], strlen(arr->strings[i]));
    }
    if (arr->length > 5) printf("...\n");
    for (int i = arr->length-5; i < arr->length; i++) {
        printf("[%d] Raw: '%s' (len:%zu)\n", i, arr->strings[i], strlen(arr->strings[i]));
    }
}

// Fast calculation with separate allocations
NumericArray calculate_solar_air_temp_fast(const NumericArray* glob_solir_vec,
    const NumericArray* air_temp_vec,
    double alpha,
    double h,
    double T_cor_fact) 
{
    NumericArray result = {NULL, NULL, 0};

    // Input validation
    if (!glob_solir_vec || !air_temp_vec || 
        !glob_solir_vec->values || !air_temp_vec->values ||
        !glob_solir_vec->is_valid || !air_temp_vec->is_valid ||
        glob_solir_vec->length != air_temp_vec->length) {
        return result;
    }

    const int length = glob_solir_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    result.is_valid = malloc(length * sizeof(bool));

    if (!result.values || !result.is_valid) {
        free(result.values);
        free(result.is_valid);
        return result;
    }

    result.length = length;

    // Performance-critical loop
    const double* restrict glob_solir = glob_solir_vec->values;
    const double* restrict air_temp = air_temp_vec->values;
    const bool* restrict valid1 = glob_solir_vec->is_valid;
    const bool* restrict valid2 = air_temp_vec->is_valid;

    double* restrict out_values = result.values;
    bool* restrict out_valid = result.is_valid;

    for (int i = 0; i < length; i++) {
        const bool valid = valid1[i] && valid2[i];
        out_valid[i] = valid;
        out_values[i] = valid ? (alpha * glob_solir[i] / h + air_temp[i] - T_cor_fact) : 0.0;
    }

    return result;
}

NumericArray calculate_solar_air_temp_simple(const NumericArray* glob_solir_vec,
    const NumericArray* air_temp_vec,
    double alpha,
    double h,
    double T_cor_fact) 
{
    // Initialize empty result
    NumericArray result = {NULL, NULL, 0};

    // Basic validation
    if (!glob_solir_vec || !air_temp_vec ||             // NULL pointer check
        glob_solir_vec->length != air_temp_vec->length) // Vector length check
    {  
        return result;
    }

    const int length = glob_solir_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    result.is_valid = malloc(length * sizeof(bool));
    if (!result.values || !result.is_valid) {
        free(result.values);
        free(result.is_valid);
        return result;
    }
    result.length = length;

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = alpha * glob_solir_vec->values[i] / h 
        + air_temp_vec->values[i] 
        - T_cor_fact;
        result.is_valid[i] = true; // Assuming all inputs are valid
    }

    return result;
}

NumericArray calculate_surface_power(const NumericArray* T_sol_air_vec,
    double A_surf,
    double lam_i,
    double d_ins)     
{
    // Initialize empty result
    NumericArray result = {NULL, NULL, 0};

    // Basic validation
    if (!T_sol_air_vec || !T_sol_air_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = T_sol_air_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    result.is_valid = malloc(length * sizeof(bool));
    if (!result.values || !result.is_valid) {
        free(result.values);
        free(result.is_valid);
    return result;
    }
    result.length = length;

    // Pre-calculate the constant factor
    const double heat_transfer_coeff = A_surf * lam_i / d_ins;

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = heat_transfer_coeff * T_sol_air_vec->values[i];
        result.is_valid[i] = T_sol_air_vec->is_valid[i]; // Inherit validity
    }

    return result;
}

NumericArray calculate_surf_meltrate(const NumericArray* Q_surf_vec,
    double L_f,
    double rho_snow)     
{
    // Initialize empty result
    NumericArray result = {NULL, NULL, 0};

    // Basic validation
    if (!Q_surf_vec || !Q_surf_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = Q_surf_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    result.is_valid = malloc(length * sizeof(bool));
    if (!result.values || !result.is_valid) {
        free(result.values);
        free(result.is_valid);
    return result;
    }
    result.length = length;

    // Pre-calculate the constant factor
    const double melt_rate_coeff = 1.0/(L_f*rho_snow);

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = melt_rate_coeff * Q_surf_vec->values[i];
        result.is_valid[i] = Q_surf_vec->is_valid[i]; // Inherit validity
    }

    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    // Keep all columns (0-17)
    //int columns_to_keep[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
    int columns_to_keep[] = {5, 6, 7, 12, 15, 17};
    int num_cols = sizeof(columns_to_keep)/sizeof(columns_to_keep[0]);
    
    CSVData* data = read_csv_with_encoding("Snow_Storage_Data.csv", columns_to_keep, num_cols);
    if (!data) {
        printf("Failed to read CSV\n");
        return 1;
    }

    printf("Total columns detected: %d\n", data->cols);

    // 1. Print first 5 rows of first YEAR column (index 0)
    printf("\nFirst 5 rows of YEAR column:\n");
    print_column(data, 0, 5);

    // 2. Find first empty value in YEAR column
    int first_empty = find_first_empty(data, 0);
    printf("\nFirst empty YEAR at index: %d\n", first_empty);

    // 3. Get time data and find period - using first time column (index 6)
    // 5,  6, 7, 12, 15, 17
    // 0, [1], 2,  3, 4,  5
    int time_col_index = 1;
    int period_start = find_index(data, time_col_index, "2023-04-01T00:00");
    int period_end = find_index(data, time_col_index, "2023-08-31T23:00");
    
    printf("\nPeriod from %d to %d\n", period_start, period_end);

    // 4. Print time column stats
    printf("\nTime column stats:\n");
    print_column_stats(data, time_col_index, "time");

    // 5. Print data for the period
    printf("\nData for the period:\n");
    print_data_range(data, period_start, period_end, "Period Data");

    // Get air temperature data (index 7)
    // 5, 6,  7, 12, 15, 17
    // 0, 1, [2],  3, 4,  5
    int air_temp_col = 2;
    StringArray air_temp_raw = extract_column_range(data, air_temp_col, period_start, period_end);
    
    // Convert to double
    NumericArray air_temp_vec = convert_string_array_to_double(&air_temp_raw);
    print_numeric_array(&air_temp_vec, "Air temperature (Celsius)");
    
    // Get wind speed data (index 15)
    // 5, 6, 7, 12, 15, 17
    // 0, 1, 2,  3, [4],  5
    StringArray wind_speed_raw = extract_column_range(data, 4, period_start, period_end);
    NumericArray wind_speed_vec = convert_string_array_to_double(&wind_speed_raw);
    print_numeric_array(&wind_speed_vec, "Wind speed (km/h)");

    // Get precipitation data (index 5)
    //  5, 6, 7, 12, 15, 17
    // [0], 1, 2,  3, 4,  5
    StringArray prec_raw = extract_column_range(data, 0, period_start, period_end);
    NumericArray prec_vec = convert_string_array_to_double(&prec_raw);
    print_numeric_array(&prec_vec, "Precipitation (m/h)");

    // Get relative humidity data (index 12)
    // 5, 6, 7, 12, 15, 17
    // 0, 1, 2, [3], 4, 5
    StringArray glob_solir_vec_raw = extract_column_range(data, 3, period_start, period_end);
    // Debug output to see what we're working with
    //debug_string_array(&glob_solir_vec_raw, "Global solar irradiance (W/m2)");
    // Convert to double
    NumericArray glob_solir_vec = convert_string_array_to_double(&glob_solir_vec_raw);
    // Print the results
    print_numeric_array(&glob_solir_vec, "Global solar irradiance (W/m2)");

    // Extract the amount of RH precipitation column from the data (index 17)
    // 5, 6, 7, 12, 15, 17
    // 0, 1, 2,  3, 4, [5]
    StringArray RH_perc_vec_raw = extract_column_range(data, 5, period_start, period_end);
    NumericArray RH_perc_vec = convert_string_array_to_double(&RH_perc_vec_raw);
    print_numeric_array(&RH_perc_vec, "Relative Humidity Percipitation (m/h)");

    //------------------------------------------------------------------------------------

    // Constants
    double h = 22.7; // Heat transfer coefficient at the external surface
    double alpha = 0.8; // Solar light absorptivity
    double T_cor_fact = 4.0; // °C //Correction factor for horizontal surface

    NumericArray T_sol_air_vec = calculate_solar_air_temp_simple(
        &glob_solir_vec, 
        &air_temp_vec,
        alpha, h, T_cor_fact
    );

    print_numeric_array(&T_sol_air_vec, "Solar-Air (Celsius)");

    //------------------------------------------------------------------------------------
   
    double d_ins = 0.1; // m        // Insulation layer thickness
    double lam_i = 0.32; // W/(mK)  // Thermal conductivity for the insulating material
    double A_surf = 210.0; // m^2   // The surface area (m2) of the pile of snow

    NumericArray Q_surf_vec = calculate_surface_power(&T_sol_air_vec, A_surf, lam_i, d_ins);
    print_numeric_array(&Q_surf_vec, "Surface power (W)");

    //------------------------------------------------------------------------------------

    double L_f = 333.4E03;   // J/kg; latent heat of fusion
    double rho_snow = 411.0; // kg/m^3; density of snow

    NumericArray f_srf_melt_vec = calculate_surf_meltrate(&Q_surf_vec, L_f, rho_snow);
    print_numeric_array(&f_srf_melt_vec, "Surface melt rate (m^3/s)");

    //------------------------------------------------------------------------------------


    // Cleanup
    free_string_array(&air_temp_raw);
    free_numeric_array(&air_temp_vec);
    free_string_array(&wind_speed_raw);
    free_numeric_array(&wind_speed_vec);
    free_string_array(&prec_raw);
    free_numeric_array(&prec_vec);
    free_string_array(&glob_solir_vec_raw);
    free_numeric_array(&glob_solir_vec);
    free_string_array(&RH_perc_vec_raw);
    free_numeric_array(&RH_perc_vec);

    free_numeric_array(&T_sol_air_vec);
    free_numeric_array(&Q_surf_vec);
    free_numeric_array(&f_srf_melt_vec);

    free_csv_data(data);

    return 0;
}