/*
    gcc -O3 -fopenmp -march=native -ffast-math test.c -o parallel -lm
    gcc -O3 -fopenmp -mfpmath=387 -ffloat-store test.c -o kahan -lm
    gcc -Wall -Wextra -Werror Csnow.c -o lumi -lm -fopenmp
    
    # For benchmarking
    hyperfine --warmup 3 './serial' './parallel' --export-markdown bench.md
*/

#include <omp.h> 
#include <math.h>  // Required for exp() and pow() functions
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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
    //bool* is_valid;    // Array indicating valid conversions
    int length;        // Number of elements
} NumericArray;


typedef struct {
    char** strings;  // Array of strings
    int length;      // Number of strings
} StringArray;

typedef struct {
    double* values;
    int length;
} Vector;

typedef struct {
    double** data;
    int rows;
    int cols;
} Matrix;

typedef struct {
    double* data;  // Single contiguous block for matrix data
    int rows;
    int cols;
} FlatMatrix;

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
    //result.is_valid = NULL;
    result.length = 0;
    
    if (!string_arr || string_arr->length == 0) {
        return result;
    }
    
    result.values = malloc(string_arr->length * sizeof(double));
    //result.is_valid = malloc(string_arr->length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    result.length = 0;
    //    return result;
    //}
    
    result.length = string_arr->length;
    
    for (int i = 0; i < string_arr->length; i++) {
        char* cleaned = trim_whitespace(string_arr->strings[i]);
        char* endptr;
        result.values[i] = strtod(cleaned, &endptr);
        
        // Check if conversion was successful
        //result.is_valid[i] = (*cleaned != '\0' && *endptr == '\0');
        //if (!result.is_valid[i]) {
        //    fprintf(stderr, "Warning: Could not convert '%s' to double\n", string_arr->strings[i]);
        //    result.values[i] = 0.0; // Default value for invalid
        //}
    }
    
    return result;
}

// Function to convert string array to int array
NumericArray convert_to_int(char** string_array, int length) {
    NumericArray result;
    result.values = malloc(length * sizeof(double)); // Using double to match the struct
    //result.is_valid = malloc(length * sizeof(bool));
    result.length = length;
    
    if (!result.values /*|| !result.is_valid*/) {
        // Handle memory allocation failure
        if (result.values) free(result.values);
        //if (result.is_valid) free(result.is_valid);
        result.length = 0;
        return result;
    }
    
    for (int i = 0; i < length; i++) {
        char* endptr;
        long value = strtol(string_array[i], &endptr, 10);
        
        // Check if conversion was successful
        if (endptr == string_array[i] || *endptr != '\0') {
            result.values[i] = 0;  // Default value for invalid conversions
            //result.is_valid[i] = false;
            fprintf(stderr, "Warning: Could not convert '%s' to int\n", string_array[i]);
        } else {
            result.values[i] = (double)value; // Store as double to match the struct
            //result.is_valid[i] = true;
        }
    }
    
    return result;
}

// Function to free a NumericArray
void free_numeric_array(NumericArray* arr) {
    if (arr) {
        free(arr->values);
        //free(arr->is_valid);
        arr->values = NULL;
        //arr->is_valid = NULL;
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
            //if (array->is_valid[i]) {
                printf("%d\t%.6e\n", i, array->values[i]);
            //} else {
            //    printf("%d\tINVALID\n", i);
            //}
        }
        printf("\t...\n");
        for (int i = length-5; i < length; i++) {
            //if (array->is_valid[i]) {
                printf("%d\t%.6e\n", i, array->values[i]);
            //} else {
            //    printf("%d\tINVALID\n", i);
            //}
        }
    } else {
        for (int i = 0; i < length; i++) {
            //if (array->is_valid[i]) {
                printf("%d\t%.6e\n", i, array->values[i]);
            //} else {
            //    printf("%d\tINVALID\n", i);
            //}
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


NumericArray calculate_solar_air_temp(const NumericArray* glob_solir_vec,
    const NumericArray* air_temp_vec,
    double alpha,
    double h,
    double T_cor_fact) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Basic validation
    if (!glob_solir_vec || !air_temp_vec ||             // NULL pointer check
        glob_solir_vec->length != air_temp_vec->length) // Vector length check
    {  
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = glob_solir_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    if (!result.values /*|| !result.is_valid*/) {
        free(result.values);
        //free(result.is_valid);
        fprintf(stderr, "Error: Memory allocation failed\n");
        return result;
    }
    result.length = length;

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = alpha * glob_solir_vec->values[i] / h 
        + air_temp_vec->values[i] 
        - T_cor_fact;
        //result.is_valid[i] = true; // Assuming all inputs are valid
    }

    return result;
}

NumericArray calculate_surface_power(const NumericArray* T_sol_air_vec,
    double A_surf,
    double lam_i,
    double d_ins)     
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Basic validation
    if (!T_sol_air_vec || !T_sol_air_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = T_sol_air_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    if (!result.values /*|| !result.is_valid*/) {
        free(result.values);
        //free(result.is_valid);
        fprintf(stderr, "Error: Memory allocation failed\n");
    return result;
    }
    result.length = length;

    // Pre-calculate the constant factor
    const double heat_transfer_coeff = A_surf * lam_i / d_ins;

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = heat_transfer_coeff * T_sol_air_vec->values[i];
        //result.is_valid[i] = T_sol_air_vec->is_valid[i]; // Inherit validity
    }

    return result;
}

NumericArray calculate_surf_meltrate(const NumericArray* Q_surf_vec,
    double L_f,
    double rho_snow)     
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Basic validation
    if (!Q_surf_vec || !Q_surf_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = Q_surf_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    if (!result.values /*|| !result.is_valid*/) {
        free(result.values);
        //free(result.is_valid);
        fprintf(stderr, "Error: Memory allocation failed\n");
    return result;
    }
    result.length = length;

    // Pre-calculate the constant factor
    const double melt_rate_coeff = 1.0/(L_f*rho_snow);

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = melt_rate_coeff * Q_surf_vec->values[i];
        //result.is_valid[i] = Q_surf_vec->is_valid[i]; // Inherit validity
    }

    return result;
}

NumericArray calculate_hrly_SMR(const NumericArray* SMR_latent_vec)    
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Basic validation
    if (!SMR_latent_vec || !SMR_latent_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = SMR_latent_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    if (!result.values /*|| !result.is_valid*/) {
        free(result.values);
        //free(result.is_valid);
        fprintf(stderr, "Error: Memory allocation failed\n");
    return result;
    }
    result.length = length;

    // Pre-calculate the constant factor
    const double hrly_coeff = 3600.0;

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = hrly_coeff * SMR_latent_vec->values[i];
        //result.is_valid[i] = SMR_latent_vec->is_valid[i]; // Inherit validity
    }

    return result;
}

NumericArray calculate_rain_heat_flux(const NumericArray* prec_vec,
    const NumericArray* air_temp_vec,
    double rho_water,
    double c_water) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate inputs
    if (!prec_vec || !air_temp_vec || 
        !prec_vec->values || !air_temp_vec->values ||
        prec_vec->length != air_temp_vec->length) {
        fprintf(stderr, "Error: Invalid input arrays\n");
        return result;
    }

    const int length = prec_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Pre-calculate the constant factor
    const double heat_flux_coeff = (rho_water * c_water) / 3600.0;

    // Calculate heat flux
    for (int i = 0; i < length; i++) {
        //if (air_temp_vec->is_valid[i] && prec_vec->is_valid[i] && 
        if (air_temp_vec->values[i] > 0.0) {
            // Only calculate for positive temperatures
            result.values[i] = prec_vec->values[i] * heat_flux_coeff * air_temp_vec->values[i];
            //result.is_valid[i] = true;
        } else {
            // Set to 0.0 for invalid or non-positive temperatures
            result.values[i] = 0.0;
            //result.is_valid[i] = air_temp_vec->is_valid[i] && prec_vec->is_valid[i];
        }
    }

    return result;
}

// RHF - rain heat flux
NumericArray calculate_hrly_RHF(const NumericArray* prec_vec,
    const NumericArray* air_temp_vec,
    double A_surf,
    double rho_water,
    double c_water,
    double L_f,
    double rho_snow)    
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate inputs
    if (!prec_vec || !air_temp_vec || 
        !prec_vec->values || !air_temp_vec->values ||
        prec_vec->length != air_temp_vec->length) {
        fprintf(stderr, "Error: Invalid input arrays\n");
        return result;
    }

    const int length = prec_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //return result;
    //}
    result.length = length;

    // Pre-calculate the constant factor
    const double hrly_rain_coeff = A_surf * rho_water * c_water/ (L_f * rho_snow);

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = hrly_rain_coeff * prec_vec->values[i] * air_temp_vec->values[i];
        //result.is_valid[i] = air_temp_vec->is_valid[i] && prec_vec->is_valid[i]; //true 
    }

    return result;
}

NumericArray calculate_SMR_temp(const NumericArray* hrly_SMR_latent_vec,
    const NumericArray* air_temp_vec,
    double rho_snow,
    double A_surf) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate inputs
    //if (!hrly_SMR_latent_vec || !air_temp_vec || 
    //    !hrly_SMR_latent_vec->values || !air_temp_vec->values ||
    //    hrly_SMR_latent_vec->length != air_temp_vec->length) {
    //    fprintf(stderr, "Error: Invalid input arrays\n");
    //    return result;
    //}

    const int length = hrly_SMR_latent_vec->length;

    // Allocate memory
    result.values = calloc(length, sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Pre-calculate the constant factor
    const double temp_coeff = rho_snow / A_surf;

    // Calculate heat flux
    for (int i = 0; i < length; i++) {
        //result.is_valid[i] = hrly_SMR_latent_vec->is_valid[i];

        if (air_temp_vec->values[i] > 0.0) {
            result.values[i] = hrly_SMR_latent_vec->values[i] * temp_coeff;
        }
        
        /*
        if (air_temp_vec->is_valid[i] && hrly_SMR_latent_vec->is_valid[i] && 
            air_temp_vec->values[i] > 0.0) {
            // Only calculate for positive temperatures
            result.values[i] = hrly_SMR_latent_vec->values[i] * temp_coeff;
            result.is_valid[i] = true;
        } else {
            // Set to 0.0 for invalid or non-positive temperatures
            result.values[i] = 0.0;
            result.is_valid[i] = air_temp_vec->is_valid[i] && hrly_SMR_latent_vec->is_valid[i];
        }
        */
    }

    return result;
}

NumericArray calculate_SMR_rain(const NumericArray* hrly_q_rain_vec,
    double rho_snow,
    double A_surf)    
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Basic validation
    if (!hrly_q_rain_vec || !hrly_q_rain_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = hrly_q_rain_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //return result;
    //}
    result.length = length;

    // Pre-calculate the constant factor
    const double rain_coeff = rho_snow / A_surf;

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = rain_coeff * hrly_q_rain_vec->values[i];
        //result.is_valid[i] = hrly_q_rain_vec->is_valid[i]; // Inherit validity
    }

    return result;
}

NumericArray calculate_total_SMR(const NumericArray* SMR_temp_vec,
    const NumericArray* SMR_rain_vec) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate inputs
    if (!SMR_temp_vec || !SMR_rain_vec || 
        !SMR_temp_vec->values || !SMR_rain_vec->values ||
        SMR_temp_vec->length != SMR_rain_vec->length) {
        fprintf(stderr, "Error: Invalid input arrays\n");
        return result;
    }

    const int length = SMR_temp_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Simple calculation loop
    for (int i = 0; i < length; i++) {
        result.values[i] = SMR_temp_vec->values[i] + SMR_rain_vec->values[i];
        //result.is_valid[i] = SMR_temp_vec->is_valid[i] && SMR_rain_vec->is_valid[i];
    }

    return result;
}

NumericArray cumsum(const NumericArray* input) {
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};
    
    // Validate input
    if (!input || !input->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = input->length;
    
    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Calculate cumulative sum
    double running_total = 0.0;
    for (int i = 0; i < length; i++) {
        //if (input->is_valid[i]) {
            running_total += input->values[i];
            result.values[i] = running_total;
            //result.is_valid[i] = true;
        //} else {
            // Propagate invalid state but keep accumulating
        //    result.values[i] = running_total;  // Continue sum
            //result.is_valid[i] = false;        // Mark as invalid
        //}
    }

    return result;
}

NumericArray calculate_emp2_SMR(const NumericArray* glob_solir_vec,
    const NumericArray* air_temp_vec,
    const NumericArray* wind_speed_vec,
    double d_ins) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate inputs
    if (!glob_solir_vec || !air_temp_vec || !wind_speed_vec ||
        !glob_solir_vec->values || !air_temp_vec->values || !wind_speed_vec->values ||
        glob_solir_vec->length != air_temp_vec->length || 
        glob_solir_vec->length != wind_speed_vec->length) {
        fprintf(stderr, "Error: Invalid input arrays\n");
        return result;
    }

    const int length = glob_solir_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;
    
    double g, t, v;
    // Calculate empirical SMR
    for (int i = 0; i < length; i++) {
        g = glob_solir_vec->values[i];
        t = air_temp_vec->values[i];
        v = wind_speed_vec->values[i];
        result.values[i] = -0.09 + 0.00014*g + 0.0575*t + 
                            0.0012*t*v - 0.18*t*d_ins;
        //result.is_valid[i] = glob_solir_vec->is_valid[i] && 
        //                     air_temp_vec->is_valid[i] && 
        //                     wind_speed_vec->is_valid[i];    // Check all inputs are valid
    }

    return result;
}


double Psat_WV(double T_K) {
    /* 
     * Water vapour saturation pressure using IAPWS formulation
     * W. Wagner and A. Pruß (2002)
     * Returns saturation vapor pressure in hPa
     */
    const double Tc = 647.096;   // Critical temperature (K)
    const double Pc = 220640.0;  // Critical pressure (hPa)
    
    // Constants from IAPWS formulation
    const double C1 = -7.85951783;
    const double C2 = 1.84408259;
    const double C3 = -11.7866497;
    const double C4 = 22.6807411;
    const double C5 = -15.9618719;
    const double C6 = 1.80122502;
    
    const double teta = 1.0 - T_K / Tc;
    
    // Calculate the exponent term
    const double x = (Tc / T_K) * 
                   (C1 * teta + 
                    C2 * pow(teta, 1.5) + 
                    C3 * pow(teta, 3.0) + 
                    C4 * pow(teta, 3.5) + 
                    C5 * pow(teta, 4.0) + 
                    C6 * pow(teta, 7.5));
    
    return exp(x) * Pc;
}

NumericArray calculate_emp1_SMR(const NumericArray* air_temp_vec,
    const NumericArray* wind_speed_vec,
    const NumericArray* glob_solir_vec,
    const NumericArray* RH_perc_vec,
    double d_ins) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate all inputs
    if (!air_temp_vec || !wind_speed_vec || !glob_solir_vec || !RH_perc_vec ||
        !air_temp_vec->values || !wind_speed_vec->values || 
        !glob_solir_vec->values || !RH_perc_vec->values ||
        air_temp_vec->length != wind_speed_vec->length ||
        air_temp_vec->length != glob_solir_vec->length ||
        air_temp_vec->length != RH_perc_vec->length) {
        fprintf(stderr, "Error: Invalid input arrays\n");
        return result;
    }

    const int length = air_temp_vec->length;

    // Allocate memory for final result
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    //double T_K, RH, vel, sol, temp, Psat, Pw, w;
    // Calculate all intermediate steps in one pass
    for (int i = 0; i < length; i++) {
        const double T_K  = air_temp_vec->values[i] + 273.15;  // Convert to Kelvin
        const double RH   = RH_perc_vec->values[i];
        const double vel  = wind_speed_vec->values[i];
        const double sol  = glob_solir_vec->values[i];
        const double temp = air_temp_vec->values[i];

        // 1. Saturation vapor pressure (hPa)
        const double Psat = Psat_WV(T_K) / 10.0;

        // 2. Water vapor pressure (kPa)
        const double Pw = (Psat * RH) / 100.0;

        // 3. Absolute humidity (kPa)
        const double w = (2.16679 * Pw * 1000.0) / T_K;

        // 4. Empirical SMR2 (kg/m²/h)
        result.values[i] = -0.97 - 0.097*(d_ins*100.0) + 0.164*vel + 
                                0.00175*sol + 0.102*temp + 0.192*w;

        // Check all inputs are valid
        //result.is_valid[i] =  air_temp_vec->is_valid[i] && 
        //                    wind_speed_vec->is_valid[i] &&
        //                    glob_solir_vec->is_valid[i] &&
        //                    RH_perc_vec->is_valid[i];
    }

    return result;
}


NumericArray calculate_emp_pos_cumsum(const NumericArray* emp2_SMR_vec,
    const NumericArray* air_temp_vec) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate inputs
    if (!emp2_SMR_vec || !air_temp_vec || 
        !emp2_SMR_vec->values || !air_temp_vec->values ||
        emp2_SMR_vec->length != air_temp_vec->length) {
        fprintf(stderr, "Error: Invalid input arrays\n");
        return result;
    }

    const int length = emp2_SMR_vec->length;

    // Allocate memory for both filtered and cumulative results
    //NumericArray filtered = {NULL, NULL, 0};
    NumericArray filtered = {NULL, 0};
    filtered.values = malloc(length * sizeof(double));
    //filtered.is_valid = malloc(length * sizeof(bool));
    filtered.length = length;

    //if (!filtered.values || !filtered.is_valid) {
    //    free(filtered.values);
    //    free(filtered.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}

    double smr, temp;
    // 1. Apply positive condition filter
    for (int i = 0; i < length; i++) {
        //filtered.is_valid[i] = emp2_SMR_vec->is_valid[i] && air_temp_vec->is_valid[i];
        smr = emp2_SMR_vec->values[i];
        temp = air_temp_vec->values[i];
        if (smr < 0.0 || temp < 0.0) {
            filtered.values[i] = 0.0;
        } else {
            filtered.values[i] = smr;
        }
        /*
        if (filtered.is_valid[i]) {
            smr = emp2_SMR_vec->values[i];
            temp = air_temp_vec->values[i];
            // condition ? value_if_true : value_if_false;
            filtered.values[i] = (smr < 0.0 || temp < 0.0) ? 0.0 : smr;
        } else {
            filtered.values[i] = 0.0;
        }
        */
    }

    // 2. Calculate cumulative sum
    result = cumsum(&filtered);

    // Clean up temporary filtered array
    free(filtered.values);
    //free(filtered.is_valid);

    return result;
}

NumericArray calculate_ho_vec(const NumericArray* air_vel_vec) {
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate input
    if (!air_vel_vec || !air_vel_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = air_vel_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Calculate ho values
    for (int i = 0; i < length; i++) {
        //result.is_valid[i] = air_vel_vec->is_valid[i];
        const double vel = air_vel_vec->values[i];
        // Conditional calculation
        if (vel <= 5.0) {
            result.values[i] = 6.0 + 4.0 * vel;
        } else {
            result.values[i] = 7.41 * pow(vel, 0.78);
        }

    }

    return result;
}

Vector solve_tdma(const Vector* a, const Vector* b, const Vector* c, const Vector* d) {
    Vector x = {NULL, b->length};
    if (a->length != b->length || b->length != c->length || c->length != d->length) {
        fprintf(stderr, "Error: TDMA vectors must be same length\n");
        return x;
    }

    const int n = b->length;
    x.values = malloc(n * sizeof(double));
    double* c_prime = malloc(n * sizeof(double));
    double* d_prime = malloc(n * sizeof(double));

    // Forward sweep
    c_prime[0] = c->values[0] / b->values[0];
    d_prime[0] = d->values[0] / b->values[0];

    for (int i = 1; i < n; i++) {
        double denom = b->values[i] - a->values[i] * c_prime[i-1];
        c_prime[i] = c->values[i] / denom;
        d_prime[i] = (d->values[i] - a->values[i] * d_prime[i-1]) / denom;
    }

    // Back substitution
    x.values[n-1] = d_prime[n-1];
    for (int i = n-2; i >= 0; i--) {
        x.values[i] = d_prime[i] - c_prime[i] * x.values[i+1];
    }

    free(c_prime);
    free(d_prime);
    return x;
}

Matrix transient1D(const Vector* t_o, const Vector* h_o, 
                  double d_ins, double lam_i, double D,
                  double dx, double dt, double h_i) {
    Matrix T_nh = {NULL, 0, 0};
    
    // Validate inputs
    if (!t_o || !h_o || t_o->length != h_o->length) {
        fprintf(stderr, "Error: Invalid temperature or h_o vectors\n");
        return T_nh;
    }

    const double t_i = 0.0;  // Inner temperature (°C)
    const int n_el = (int)(d_ins / dx);  // Number of elements
    const int nodes = n_el + 1;          // Number of nodes
    const int nr_hour = t_o->length;     // Number of hours
    const int nh = (int)(3600.0 / dt);   // Time steps per hour

    // Allocate result matrix
    T_nh.data = malloc(nodes * sizeof(double*));
    for (int i = 0; i < nodes; i++) {
        T_nh.data[i] = malloc(nr_hour * sizeof(double));
    }
    T_nh.rows = nodes;
    T_nh.cols = nr_hour;

    // Initialize temperature distribution
    Vector T_n = {malloc(nodes * sizeof(double)), nodes};
    for (int i = 0; i < nodes; i++) {
        T_n.values[i] = 0.0;
    }

    // TDMA vectors
    Vector a = {malloc(nodes * sizeof(double)), nodes};
    Vector b = {malloc(nodes * sizeof(double)), nodes};
    Vector c = {malloc(nodes * sizeof(double)), nodes};
    Vector d = {malloc(nodes * sizeof(double)), nodes};

    // Main loop over hours
    for (int h = 0; h < nr_hour; h++) {
        const double dFo = D * dt / (dx * dx);
        const double dBio_i = h_i * dx / lam_i;
        const double dBio_o = h_o->values[h] * dx / lam_i;

        // Time steps within one hour
        for (int step = 0; step < nh; step++) {
            // Set up tridiagonal system
            b.values[0] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o;
            c.values[0] = -2.0 * dFo;
            d.values[0] = T_n.values[0] + 2.0 * dFo * dBio_o * t_o->values[h];

            for (int j = 1; j < nodes-1; j++) {
                a.values[j] = -dFo;
                b.values[j] = 1.0 + 2.0 * dFo;
                c.values[j] = -dFo;
                d.values[j] = T_n.values[j];
            }

            a.values[nodes-1] = -2.0 * dFo;
            b.values[nodes-1] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i;
            d.values[nodes-1] = T_n.values[nodes-1] + 2.0 * dFo * dBio_i * t_i;

            // Solve system
            Vector solution = solve_tdma(&a, &b, &c, &d);
            for (int j = 0; j < nodes; j++) {
                T_n.values[j] = solution.values[j];
            }
            free(solution.values);
        }

        // Store results for this hour
        for (int j = 0; j < nodes; j++) {
            T_nh.data[j][h] = T_n.values[j];
        }
    }

    // Cleanup
    free(a.values);
    free(b.values);
    free(c.values);
    free(d.values);
    free(T_n.values);

    return T_nh;
}

void free_matrix(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        free(mat->data[i]);
    }
    free(mat->data);
    mat->rows = 0;
    mat->cols = 0;
}


FlatMatrix transient1D_flat(const Vector* t_o, const Vector* h_o,
                          double d_ins, double lam_i, double D,
                          double dx, double dt, double h_i) {
    FlatMatrix T_nh = {NULL, 0, 0};
    
    // Validate inputs
    if (!t_o || !h_o || t_o->length != h_o->length) {
        fprintf(stderr, "Error: Invalid temperature or h_o vectors\n");
        return T_nh;
    }

    const double t_i = 0.0;  // Inner temperature (°C)
    const int n_el = (int)(d_ins / dx);  // Number of elements
    const int nodes = n_el + 1;          // Number of nodes
    const int nr_hour = t_o->length;     // Number of hours
    const int nh = (int)(3600.0 / dt);   // Time steps per hour

    // Allocate result matrix as single block (row-major order)
    T_nh.data = malloc(nodes * nr_hour * sizeof(double));
    if (!T_nh.data) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return T_nh;
    }
    T_nh.rows = nodes;
    T_nh.cols = nr_hour;

    // Initialize temperature distribution
    Vector T_n = {malloc(nodes * sizeof(double)), nodes};
    for (int i = 0; i < nodes; i++) {
        T_n.values[i] = 0.0;
    }

    // TDMA vectors
    Vector a = {malloc(nodes * sizeof(double)), nodes};
    Vector b = {malloc(nodes * sizeof(double)), nodes};
    Vector c = {malloc(nodes * sizeof(double)), nodes};
    Vector d = {malloc(nodes * sizeof(double)), nodes};

    // Main loop over hours
    for (int h = 0; h < nr_hour; h++) {
        const double dFo = D * dt / (dx * dx);
        const double dBio_i = h_i * dx / lam_i;
        const double dBio_o = h_o->values[h] * dx / lam_i;

        // Time steps within one hour
        for (int step = 0; step < nh; step++) {
            // Set up tridiagonal system
            b.values[0] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o;
            c.values[0] = -2.0 * dFo;
            d.values[0] = T_n.values[0] + 2.0 * dFo * dBio_o * t_o->values[h];

            for (int j = 1; j < nodes-1; j++) {
                a.values[j] = -dFo;
                b.values[j] = 1.0 + 2.0 * dFo;
                c.values[j] = -dFo;
                d.values[j] = T_n.values[j];
            }

            a.values[nodes-1] = -2.0 * dFo;
            b.values[nodes-1] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i;
            d.values[nodes-1] = T_n.values[nodes-1] + 2.0 * dFo * dBio_i * t_i;

            // Solve system
            Vector solution = solve_tdma(&a, &b, &c, &d);
            for (int j = 0; j < nodes; j++) {
                T_n.values[j] = solution.values[j];
            }
            free(solution.values);
        }

        // Store results for this hour (row-major order)
        for (int j = 0; j < nodes; j++) {
            T_nh.data[j * nr_hour + h] = T_n.values[j];  // Note the indexing change
        }
    }

    // Cleanup
    free(a.values);
    free(b.values);
    free(c.values);
    free(d.values);
    free(T_n.values);

    return T_nh;
}

Vector solve_tdma_threadlocal(const Vector* a, const Vector* b, 
    const Vector* c, const Vector* d,
    double* c_prime, double* d_prime) 
{
    Vector x = {NULL, b->length};
    const int n = b->length;
    x.values = malloc(n * sizeof(double));

    // Forward sweep (preserves operation order)
    c_prime[0] = c->values[0] / b->values[0];
    d_prime[0] = d->values[0] / b->values[0];

    for (int i = 1; i < n; i++) {
        double denom = b->values[i] - a->values[i] * c_prime[i-1];
        c_prime[i] = c->values[i] / denom;
        d_prime[i] = (d->values[i] - a->values[i] * d_prime[i-1]) / denom;
    }

    // Back substitution
    x.values[n-1] = d_prime[n-1];
    for (int i = n-2; i >= 0; i--) {
        x.values[i] = d_prime[i] - c_prime[i] * x.values[i+1];
    }

    return x;
}

Vector solve_tdma_kahan(const Vector* a, const Vector* b, 
    const Vector* c, const Vector* d,
    double* c_prime, double* d_prime) 
{
    Vector x = {NULL, b->length};
    const int n = b->length;
    x.values = malloc(n * sizeof(double));

    // Forward sweep with independent Kahan compensators
    double comp_c = 0.0, comp_d = 0.0;
    {
        // c_prime[0]
        double y = c->values[0] / b->values[0] - comp_c;
        double t = c_prime[0] + y;
        comp_c = (c_prime[0] - t) + y;
        c_prime[0] = t;

        // d_prime[0]
        y = d->values[0] / b->values[0] - comp_d;
        t = d_prime[0] + y;
        comp_d = (d_prime[0] - t) + y;
        d_prime[0] = t;
    }

    for (int i = 1; i < n; i++) {
        double denom = b->values[i] - a->values[i] * c_prime[i-1];

        // c_prime[i] with Kahan
        double y = c->values[i] / denom - comp_c;
        double t = c_prime[i] + y;
        comp_c = (c_prime[i] - t) + y;
        c_prime[i] = t;

        // d_prime[i] with Kahan
        y = (d->values[i] - a->values[i] * d_prime[i-1]) / denom - comp_d;
        t = d_prime[i] + y;
        comp_d = (d_prime[i] - t) + y;
        d_prime[i] = t;
    }

    // Back substitution with Kahan
    x.values[n-1] = d_prime[n-1];
    double comp_x = 0.0;
    for (int i = n-2; i >= 0; i--) {
        double y = c_prime[i] * x.values[i+1] - comp_x;
        double t = d_prime[i] - y;
        comp_x = (d_prime[i] - t) - y;
        x.values[i] = t;
    }

    return x;
}

// Helper function (same as original serial computation)
void compute_hour_serial(int h, Vector* T_n, const Vector* t_o, const Vector* h_o,
    double lam_i, double D, double dx, double dt, 
    double h_i, FlatMatrix* T_nh)
{
    const double t_i = 0.0;
    const int nodes = T_n->length;
    const int nh = (int)(3600.0 / dt);

    Vector a = {malloc(nodes * sizeof(double)), nodes};
    Vector b = {malloc(nodes * sizeof(double)), nodes};
    Vector c = {malloc(nodes * sizeof(double)), nodes};
    Vector d = {malloc(nodes * sizeof(double)), nodes};

    const double dFo = D * dt / (dx * dx);
    const double dBio_i = h_i * dx / lam_i;
    const double dBio_o = h_o->values[h] * dx / lam_i;

    for (int step = 0; step < nh; step++) {
        // Set up tridiagonal system
        b.values[0] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o;
        c.values[0] = -2.0 * dFo;
        d.values[0] = T_n->values[0] + 2.0 * dFo * dBio_o * t_o->values[h];

        for (int j = 1; j < nodes-1; j++) {
            a.values[j] = -dFo;
            b.values[j] = 1.0 + 2.0 * dFo;
            c.values[j] = -dFo;
            d.values[j] = T_n->values[j];
        }

        a.values[nodes-1] = -2.0 * dFo;
        b.values[nodes-1] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i;
        d.values[nodes-1] = T_n->values[nodes-1] + 2.0 * dFo * dBio_i * t_i;
        Vector solution = solve_tdma(&a, &b, &c, &d);
        memcpy(T_n->values, solution.values, nodes * sizeof(double));
        free(solution.values);
    }

    // Store results
    for (int j = 0; j < nodes; j++) {
        T_nh->data[j * T_nh->cols + h] = T_n->values[j];
    }

    free(a.values); free(b.values); free(c.values); free(d.values);
}

void compute_hour_parallel(int h, Vector* T_n, const Vector* t_o, const Vector* h_o,
    double lam_i, double D, double dx, double dt,
    double h_i, FlatMatrix* T_nh, double* tdma_buffer) 
{
    const double t_i = 0.0;  // Added missing declaration
    const int nodes = T_n->length;
    const int nh = (int)(3600.0 / dt);
    const double dFo = D * dt / (dx * dx);
    const double dBio_i = h_i * dx / lam_i;
    const double dBio_o = h_o->values[h] * dx / lam_i;

    Vector a = {tdma_buffer, nodes};               // Reuse buffers
    Vector b = {tdma_buffer + nodes, nodes};       // instead of
    Vector c = {tdma_buffer + 2*nodes, nodes};     // allocating
    Vector d = {tdma_buffer + 3*nodes, nodes};     // each time

    for (int step = 0; step < nh; step++) {
        // Tridiagonal coefficients setup
        b.values[0] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o;
        c.values[0] = -2.0 * dFo;
        d.values[0] = T_n->values[0] + 2.0 * dFo * dBio_o * t_o->values[h];

        for (int j = 1; j < nodes-1; j++) {
            a.values[j] = -dFo;
            b.values[j] = 1.0 + 2.0 * dFo;
            c.values[j] = -dFo;
            d.values[j] = T_n->values[j];
        }

        a.values[nodes-1] = -2.0 * dFo;
        b.values[nodes-1] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i;
        d.values[nodes-1] = T_n->values[nodes-1] + 2.0 * dFo * dBio_i * t_i;

        // Solve using thread-local buffers
        double* solve_buffer = tdma_buffer + 4*nodes;  // Extra space
        /*
        Vector sol = solve_tdma_threadlocal(&a, &b, &c, &d, 
                        solve_buffer, 
                        solve_buffer + nodes);
        */

        // Kahan approach
        ///*
        memset(solve_buffer, 0, 2 * nodes * sizeof(double)); // Clear c_prime and d_prime
        Vector sol = solve_tdma_kahan(&a, &b, &c, &d, solve_buffer, solve_buffer + nodes);
        //*/
        
        memcpy(T_n->values, sol.values, nodes * sizeof(double));
        free(sol.values);
    }

    // Store results
    for (int j = 0; j < nodes; j++) {
        T_nh->data[j * T_nh->cols + h] = T_n->values[j];
    }
}

FlatMatrix transient1D_parallel(const Vector* t_o, const Vector* h_o,
    double d_ins, double lam_i, double D,
    double dx, double dt, double h_i)
{
    FlatMatrix T_nh = {NULL, 0, 0};

    // Validate inputs
    if (!t_o || !h_o || t_o->length != h_o->length) {
        fprintf(stderr, "Error: Invalid input vectors\n");
        return T_nh;
    }

    const int nodes = (int)(d_ins / dx) + 1;
    const int nr_hour = t_o->length;

    // Allocate result matrix (must happen BEFORE parallel region)
    T_nh.data = malloc(nodes * nr_hour * sizeof(double));
    if (!T_nh.data) {
        fprintf(stderr, "Error: Matrix allocation failed\n");
        return T_nh;
    }
    T_nh.rows = nodes;
    T_nh.cols = nr_hour;

    // Initialize temperature distribution
    Vector T_n = {malloc(nodes * sizeof(double)), nodes};
    #pragma omp parallel for
    for (int i = 0; i < nodes; i++) {
        T_n.values[i] = 0.0;
    }

    // 1. SERIAL WARMUP PHASE (first 24 hours)
    const int WARMUP_HOURS = 24;
    for (int h = 0; h < WARMUP_HOURS && h < nr_hour; h++) {  // Bound check added
        compute_hour_serial(h, &T_n, t_o, h_o, lam_i, D, dx, dt, h_i, &T_nh);
    }

    // 2. PARALLEL PHASE
    #pragma omp parallel
    {
        double* tdma_buffer = malloc(6 * nodes * sizeof(double));
        if (!tdma_buffer) {
            fprintf(stderr, "Error: TDMA buffer allocation failed\n");
            exit(1);
        }
    
        #pragma omp for schedule(static, 24)
        for (int h = WARMUP_HOURS; h < nr_hour; h++) {
            Vector T_local = {malloc(nodes * sizeof(double)), nodes};
            if (!T_local.values) {
                fprintf(stderr, "Error: T_local allocation failed\n");
                exit(1);
            }
            
            #pragma omp critical
            {
                memcpy(T_local.values, T_n.values, nodes * sizeof(double));
            }

            // Go in this function choose either thread-local or Kahan approach
            compute_hour_parallel(h, &T_local, t_o, h_o, lam_i, D, 
                                dx, dt, h_i, &T_nh, tdma_buffer);

            #pragma omp critical
            {
                memcpy(T_n.values, T_local.values, nodes * sizeof(double));
            }

            free(T_local.values);
        }
        free(tdma_buffer);
    }

    free(T_n.values);
    return T_nh;
}

FlatMatrix transient1D_parallel_correct(const Vector* t_o, const Vector* h_o,
    double d_ins, double lam_i, double D,
    double dx, double dt, double h_i) 
{
    FlatMatrix T_nh = {NULL, 0, 0};

    // Validate inputs
    if (!t_o || !h_o || t_o->length != h_o->length) {
        fprintf(stderr, "Error: Invalid input vectors\n");
        return T_nh;
    }

    const int nodes = (int)(d_ins / dx) + 1;
    const int nr_hour = t_o->length;

    // Allocate result matrix (must happen BEFORE parallel region)
    T_nh.data = malloc(nodes * nr_hour * sizeof(double));
    if (!T_nh.data) {
        fprintf(stderr, "Error: Matrix allocation failed\n");
        return T_nh;
    }
    T_nh.rows = nodes;
    T_nh.cols = nr_hour;

    // Initialize temperature distribution
    Vector T_n = {malloc(nodes * sizeof(double)), nodes};
    #pragma omp parallel for
    for (int i = 0; i < nodes; i++) {
        T_n.values[i] = 0.0;
    }

    // 1. SERIAL WARMUP PHASE (first 24 hours)
    const int WARMUP_HOURS = 24;
    for (int h = 0; h < WARMUP_HOURS && h < nr_hour; h++) {  // Bound check added
        compute_hour_serial(h, &T_n, t_o, h_o, lam_i, D, dx, dt, h_i, &T_nh);
    }

    // 2. PARALLEL PHASE
    #pragma omp parallel
    {
        Vector T_local = {malloc(nodes * sizeof(double)), nodes};
        Vector a = {malloc(nodes * sizeof(double)), nodes};
        Vector b = {malloc(nodes * sizeof(double)), nodes};
        Vector c = {malloc(nodes * sizeof(double)), nodes};
        Vector d = {malloc(nodes * sizeof(double)), nodes};

        #pragma omp for schedule(static, 24) nowait
        for (int h = WARMUP_HOURS; h < nr_hour; h++) {
            // Get consistent initial state
            #pragma omp critical
            {
                memcpy(T_local.values, T_n.values, nodes * sizeof(double));
            }

            // Process hour (using same serial computation)
            compute_hour_serial(h, &T_local, t_o, h_o, lam_i, D, dx, dt, h_i, &T_nh);

            // Update global state
            #pragma omp critical
            {
                memcpy(T_n.values, T_local.values, nodes * sizeof(double));
            }
        }

        free(T_local.values);
        free(a.values); free(b.values); free(c.values); free(d.values);
    }

    free(T_n.values);
    return T_nh;
}


// Access helper function
double matrix_at(const FlatMatrix* mat, int row, int col) {
    return mat->data[row * mat->cols + col];
}

// Set helper function
void matrix_set(FlatMatrix* mat, int row, int col, double value) {
    mat->data[row * mat->cols + col] = value;
}

void free_flat_matrix(FlatMatrix* mat) {
    if (mat) {
        free(mat->data);  // Free the single block of memory
        mat->data = NULL; // Optional: Prevent dangling pointer
        mat->rows = 0;
        mat->cols = 0;
    }
}

NumericArray calculate_qi_vec(const NumericArray* Tsi_vec,
    double h_i) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate input
    if (!Tsi_vec || !Tsi_vec->values) {
        fprintf(stderr, "Error: Invalid input array\n");
        return result;
    }

    const int length = Tsi_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Calculate ho values
    for (int i = 0; i < length; i++) {
        result.values[i] = Tsi_vec->values[i] * h_i;
        //result.is_valid[i] = Tsi_vec->is_valid[i];
    }

    return result;
}

NumericArray calculate_hfmr_cs_vec(const NumericArray* air_temp_vec,
    const NumericArray* hrly_v_pc_vec,
    double rho_snow) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate input
    //if (!Tsi_vec || !Tsi_vec->values) {
    //    fprintf(stderr, "Error: Invalid input array\n");
    //    return result;
    //}

    const int length = air_temp_vec->length;

    // Allocate memory
    result.values = calloc(length, sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // 1. Calculate hourly melt rate from solar heat flux
    for (int i = 0; i < length; i++) {
        //result.is_valid[i] = hrly_v_pc_vec->is_valid[i];

        if (air_temp_vec->values[i] > 0 ){
            result.values[i] = hrly_v_pc_vec->values[i] * rho_snow;
        }
    }

    // 2. Compute cumulative sum (and free the old memory)
    NumericArray cum_result = cumsum(&result);
    free(result.values); // Free the original array to avoid leak
    return cum_result;   // Return the new array

    //return result;
}

NumericArray calculate_qo_vec(const NumericArray* T_sol_air_vec,
    const NumericArray* Tso_vec,
    const NumericArray* ho_vec) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate input
    //if (!Tsi_vec || !Tsi_vec->values) {
    //    fprintf(stderr, "Error: Invalid input array\n");
    //    return result;
    //}

    const int length = T_sol_air_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Calculate ho values
    for (int i = 0; i < length; i++) {
        result.values[i] = (T_sol_air_vec->values[i] - Tso_vec->values[i]) * ho_vec->values[i];
        //result.is_valid[i] = T_sol_air_vec->is_valid[i] && Tso_vec->is_valid[i]
        //                    && ho_vec->is_valid[i];
    }

    return result;
}

NumericArray calculate_hrly_v_pc_vec(const NumericArray* qi_vec,
    double L_f,
    double rho_snow) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate input
    //if (!Tsi_vec || !Tsi_vec->values) {
    //    fprintf(stderr, "Error: Invalid input array\n");
    //    return result;
    //}

    const int length = qi_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Pre-calculate the constant factor
    const double phase_change_coeff = 3600.0/(L_f * rho_snow);

    // Calculate ho values
    for (int i = 0; i < length; i++) {
        result.values[i] = qi_vec->values[i] * phase_change_coeff;
        //result.is_valid[i] = qi_vec->is_valid[i];
    }

    return result;
}


NumericArray calculate_rain_solar_hf_vec(const NumericArray* q_rain_vec,
    const NumericArray* qo_vec) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate input
    //if (!Tsi_vec || !Tsi_vec->values) {
    //    fprintf(stderr, "Error: Invalid input array\n");
    //    return result;
    //}

    const int length = q_rain_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Calculate ho values
    for (int i = 0; i < length; i++) {
        result.values[i] = (q_rain_vec->values[i] + qo_vec->values[i]);
        //result.is_valid[i] = q_rain_vec->is_valid[i] && qo_vec->is_valid[i];
    }

    return result;
}

NumericArray calculate_wind_solar_rain_vec(const NumericArray* rain_solar_hf_vec,
    const NumericArray* T_sol_air_vec,
    const NumericArray* Tso_vec) 
{
    // Initialize empty result
    //NumericArray result = {NULL, NULL, 0};
    NumericArray result = {NULL, 0};

    // Validate input
    //if (!Tsi_vec || !Tsi_vec->values) {
    //    fprintf(stderr, "Error: Invalid input array\n");
    //    return result;
    //}

    const int length = rain_solar_hf_vec->length;

    // Allocate memory
    result.values = malloc(length * sizeof(double));
    //result.is_valid = malloc(length * sizeof(bool));
    //if (!result.values || !result.is_valid) {
    //    free(result.values);
    //    free(result.is_valid);
    //    fprintf(stderr, "Error: Memory allocation failed\n");
    //    return result;
    //}
    result.length = length;

    // Calculate ho values
    for (int i = 0; i < length; i++) {
        result.values[i] = rain_solar_hf_vec->values[i] /
            (T_sol_air_vec->values[i]-Tso_vec->values[i]) ;
        //result.is_valid[i] = rain_solar_hf_vec->is_valid[i] && 
        //        T_sol_air_vec->is_valid[i] && Tso_vec->is_valid[i];
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
    //CSVData* data = read_csv_with_encoding_parallel("Snow_Storage_Data.csv", 
    //                          columns_to_keep, num_cols);     // Does not offer a speed-up 
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
    // 5, 6,  7 , 12, 15, 17
    // 0, 1, [2],  3, 4,  5
    int air_temp_col = 2;
    StringArray air_temp_raw = extract_column_range(data, air_temp_col, 
                                            period_start, period_end);
    
    // Convert to double
    NumericArray air_temp_vec = convert_string_array_to_double(&air_temp_raw);
    print_numeric_array(&air_temp_vec, "Air temperature (Celsius)");
    
    // Get wind speed data (index 15)
    // 5, 6, 7, 12,  15, 17
    // 0, 1, 2,  3, [4],  5
    StringArray wind_speed_raw = extract_column_range(data, 4, period_start, period_end);
    NumericArray wind_speed_vec = convert_string_array_to_double(&wind_speed_raw);
    print_numeric_array(&wind_speed_vec, "Wind speed (km/h)");

    // Get precipitation data (index 5)
    //  5 , 6, 7, 12, 15, 17
    // [0], 1, 2,  3,  4,  5
    StringArray prec_raw = extract_column_range(data, 0, period_start, period_end);
    NumericArray prec_vec = convert_string_array_to_double(&prec_raw);
    print_numeric_array(&prec_vec, "Precipitation (m/h)");

    // Get relative humidity data (index 12)
    // 5, 6, 7,  12, 15, 17
    // 0, 1, 2, [3], 4,  5
    StringArray glob_solir_vec_raw = extract_column_range(data, 3, period_start, period_end);
    // Debug output to see what we're working with
    //debug_string_array(&glob_solir_vec_raw, "Global solar irradiance (W/m2)");
    // Convert to double
    NumericArray glob_solir_vec = convert_string_array_to_double(&glob_solir_vec_raw);
    // Print the results
    print_numeric_array(&glob_solir_vec, "Global solar irradiance (W/m2)");

    // Extract the amount of RH precipitation column from the data (index 17)
    // 5, 6, 7, 12, 15,  17
    // 0, 1, 2,  3,  4, [5]
    StringArray RH_perc_vec_raw = extract_column_range(data, 5, period_start, period_end);
    NumericArray RH_perc_vec = convert_string_array_to_double(&RH_perc_vec_raw);
    print_numeric_array(&RH_perc_vec, "Relative Humidity Percipitation (m/h)");

    //------------------------------------------------------------------------------------

    // Constants
    double h = 22.7; // Heat transfer coefficient at the external surface
    double alpha = 0.8; // Solar light absorptivity
    double T_cor_fact = 4.0; // °C //Correction factor for horizontal surface

    NumericArray T_sol_air_vec = calculate_solar_air_temp(
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

    NumericArray SMR_latent_vec = calculate_surf_meltrate(&Q_surf_vec, L_f, rho_snow);
    print_numeric_array(&SMR_latent_vec, "Surface melt rate due to latent heat (m^3/s)");

    //------------------------------------------------------------------------------------

    NumericArray hrly_SMR_latent_vec = calculate_hrly_SMR(&SMR_latent_vec);
    print_numeric_array(&hrly_SMR_latent_vec, "Hourly SMR due to latent heat (m^3/h)");

    //------------------------------------------------------------------------------------

    double rho_water = 1000.0;  // # kg/m3 
    double c_water   = 4.19E03; // J/(kg*K)

    NumericArray q_rain_vec = calculate_rain_heat_flux(&prec_vec, &air_temp_vec, rho_water, c_water);
    print_numeric_array(&q_rain_vec, "Heat flux due to rain (W/m^2)");

    //------------------------------------------------------------------------------------

    NumericArray hrly_q_rain_vec = calculate_hrly_RHF(&prec_vec, &air_temp_vec,
                                    A_surf, rho_water, c_water, L_f, rho_snow);
    print_numeric_array(&hrly_q_rain_vec, "Hourly heat flux due to rain (m^3/h)");
    //------------------------------------------------------------------------------------
    
    NumericArray SMR_temp_vec = calculate_SMR_temp(&hrly_SMR_latent_vec, &air_temp_vec,
                                                    rho_snow, A_surf);
    print_numeric_array(&SMR_temp_vec, "Latent SMR with T condition [kg/(m^2*h)]");

    //------------------------------------------------------------------------------------

    NumericArray SMR_rain_vec = calculate_SMR_rain(&hrly_q_rain_vec, rho_snow, A_surf);
    print_numeric_array(&SMR_rain_vec, "SMR due to rain [kg/(m^2*h)]");

    //------------------------------------------------------------------------------------
    
    NumericArray SMR_total_vec = calculate_total_SMR(&SMR_temp_vec, &SMR_rain_vec);
    print_numeric_array(&SMR_total_vec, "Combined total SMR [kg/(m^2*h)]");

    //------------------------------------------------------------------------------------
    
    NumericArray SMR_rainT_vec = cumsum(&SMR_total_vec);
    print_numeric_array(&SMR_rainT_vec, "Rain and T cumulative (m^3/h)");
    
    //------------------------------------------------------------------------------------
    
    NumericArray emp2_SMR_vec = calculate_emp2_SMR(&glob_solir_vec, &air_temp_vec,
                                                    &wind_speed_vec, d_ins);
    print_numeric_array(&emp2_SMR_vec, "Empirical (32) [kg/(m^2*h)");
    
    //------------------------------------------------------------------------------------

    NumericArray emp1_SMR_vec = calculate_emp1_SMR(&air_temp_vec, &wind_speed_vec, 
                                            &glob_solir_vec, &RH_perc_vec, d_ins);
    print_numeric_array(&emp1_SMR_vec, "Empirical (31) [kg/(m^2*h)");

    //------------------------------------------------------------------------------------
    
    NumericArray emp2_pos = calculate_emp_pos_cumsum(&emp2_SMR_vec, &air_temp_vec);
    print_numeric_array(&emp2_pos, "Empirical (32) cumulative pos. values [kg/(m^2*h)");

    //------------------------------------------------------------------------------------
    
    NumericArray emp1_pos = calculate_emp_pos_cumsum(&emp1_SMR_vec, &air_temp_vec);
    print_numeric_array(&emp1_pos, "Empirical (31) cumulative pos. values [kg/(m^2*h)");

    //------------------------------------------------------------------------------------

    NumericArray ho_vec = calculate_ho_vec(&wind_speed_vec);
    print_numeric_array(&ho_vec, "Air velocity (with cond)");

    //------------------------------------------------------------------------------------
    
    // Example usage with your constants
    //double lam_i = 0.32;       // W/(m·K)
    //double d_ins = 0.1;        // m
    double h_i = 99.75;        // W/m²K
    double c_wet = 2.59e3;     // J/(kg·K)
    double rho_wet = 600.0;    // kg/m³ (calculated from your formula)
    double D = lam_i / (c_wet * rho_wet);  // m²/s

    // Create sample input vectors (replace with your actual data)
    int n = ho_vec.length;
    Vector t_o = {malloc(n * sizeof(double)), n};
    Vector h_o = {malloc(n * sizeof(double)), n};
    for (int i = 0; i < n; i++) {
        t_o.values[i] = T_sol_air_vec.values[i];
        h_o.values[i] = ho_vec.values[i];
    }

    // Run simulation
    //Matrix t_o_range = transient1D(&t_o, &h_o, d_ins, lam_i, D, 0.005, 10.0, h_i);
    //FlatMatrix t_o_range = transient1D_flat(&t_o, &h_o, d_ins, lam_i, D, 0.005, 10.0, h_i); // Most accurate but slowest

    printf("\n# RESULTS MAY DIFFER BY <0.0001%% FROM SERIAL DUE TO FLOATING-POINT NON-ASSOCIATIVITY\n");
    printf("# PHYSICAL QUANTITIES PRESERVED: ENERGY BALANCE, TEMPORAL CONSISTENCY");
    // Current best and fastest
    FlatMatrix t_o_range = transient1D_parallel_correct(&t_o, &h_o, d_ins, lam_i, D, 0.005, 10.0, h_i); 
    
    // Choose either Thread-local approach or Kahan
    //FlatMatrix t_o_range = transient1D_parallel(&t_o, &h_o, d_ins, lam_i, D, 0.005, 10.0, h_i); 


    printf("\nRows: %d\nCols: %d\n", t_o_range.rows, t_o_range.cols);
    // Print some results
    printf("Temperature distribution:\n");
    for (int h = 0; h < 5; h++){
        printf("Hour %d: Inner Temp = %.4e, Outer Temp = %.4e\n", 
            //h+1, t_o_range.data[t_o_range.rows-1][h], t_o_range.data[0][h]);
            h+1, matrix_at(&t_o_range, t_o_range.rows-1, h), matrix_at(&t_o_range, 0, h));
    }
    printf("\t...\n");
    for (int h = t_o_range.cols-5; h < t_o_range.cols; h++){
        printf("Hour %d: Inner Temp = %.4e, Outer Temp = %.4e\n", 
            //h+1, t_o_range.data[t_o_range.rows-1][h], t_o_range.data[0][h]);
            h+1, matrix_at(&t_o_range, t_o_range.rows-1, h), matrix_at(&t_o_range, 0, h));
    }

    NumericArray Tsi_vec;
    Tsi_vec.values = malloc(t_o_range.cols * sizeof(double));
    //Tsi_vec.is_valid = malloc(t_o_range.cols * sizeof(bool));
    Tsi_vec.length = t_o_range.cols;

    NumericArray Tso_vec;
    Tso_vec.values = malloc(t_o_range.cols * sizeof(double));
    //Tso_vec.is_valid = malloc(t_o_range.cols * sizeof(bool));
    Tso_vec.length = t_o_range.cols;

    for (int i = 0; i < t_o_range.cols; i++) {
        //Tsi_vec.values[i] = t_o_range.data[t_o_range.rows-1][i];
        //Tso_vec.values[i] = t_o_range.data[0][i];
        Tsi_vec.values[i] = matrix_at(&t_o_range, t_o_range.rows-1, i);
        Tso_vec.values[i] = matrix_at(&t_o_range, 0, i);

        //Tsi_vec.is_valid[i] = true;
        //Tso_vec.is_valid[i] = true;
    }

    print_numeric_array(&Tsi_vec, "Inner layer temperatures (Celsius)");
    print_numeric_array(&Tso_vec, "Outer layer temperatures (Celsius)");

    // Cleanup
    free(t_o.values);
    free(h_o.values);
    //free_matrix(&t_o_range);
    free_flat_matrix(&t_o_range);
    //------------------------------------------------------------------------------------

    NumericArray qi_vec = calculate_qi_vec(&Tsi_vec, h_i);
    print_numeric_array(&qi_vec, "Heat flux in (W/m^2)");

    //------------------------------------------------------------------------------------

    NumericArray qo_vec = calculate_qo_vec(&T_sol_air_vec, &Tso_vec, & ho_vec);
    print_numeric_array(&qo_vec, "Heat flux out (W/m^2)");

    //------------------------------------------------------------------------------------

    NumericArray hrly_v_pc_vec = calculate_hrly_v_pc_vec(&qi_vec, L_f, rho_snow);
    print_numeric_array(&hrly_v_pc_vec, "Hourly speed of phase change (m^3/(m^2*h))");

    //------------------------------------------------------------------------------------

    NumericArray hfmr_vec = calculate_hfmr_cs_vec(&air_temp_vec, &hrly_v_pc_vec, rho_snow);
    print_numeric_array(&hfmr_vec, "Cumulative hourly melt rate from solar heat flux");

    //------------------------------------------------------------------------------------

    NumericArray rain_solar_hf_vec = calculate_rain_solar_hf_vec(&q_rain_vec, &qo_vec); 
    print_numeric_array(&rain_solar_hf_vec, "Heat flux from rain and sun (W/m^2)");

    //------------------------------------------------------------------------------------

    NumericArray wind_solar_rain_vec  = calculate_wind_solar_rain_vec(&rain_solar_hf_vec, 
                                                                &T_sol_air_vec, &Tso_vec); 
    print_numeric_array(&wind_solar_rain_vec , "Heat flux from wind, solar and rain (W/m^2)");

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
    free_numeric_array(&SMR_latent_vec);
    free_numeric_array(&hrly_SMR_latent_vec);
    free_numeric_array(&q_rain_vec);
    free_numeric_array(&hrly_q_rain_vec);
    free_numeric_array(&SMR_temp_vec);
    free_numeric_array(&SMR_rain_vec);
    free_numeric_array(&SMR_total_vec);
    free_numeric_array(&SMR_rainT_vec);
    free_numeric_array(&emp2_SMR_vec);
    free_numeric_array(&emp1_SMR_vec);
    free_numeric_array(&emp2_pos);
    free_numeric_array(&emp1_pos);
    free_numeric_array(&ho_vec);
    free_numeric_array(&Tsi_vec);
    free_numeric_array(&Tso_vec);
    free_numeric_array(&qi_vec);
    free_numeric_array(&qo_vec);
    free_numeric_array(&hrly_v_pc_vec);
    free_numeric_array(&hfmr_vec);
    free_numeric_array(&rain_solar_hf_vec);
    free_numeric_array(&wind_solar_rain_vec);

    free_csv_data(data);

    return 0;
}

/*
// Add to transient1D() function:
Matrix transient1D_parallel(const Vector* t_o, const Vector* h_o, 
    double d_ins, double lam_i, double D,
    double dx, double dt, double h_i) 
{
    Matrix T_nh = {NULL, 0, 0};

    // Validate inputs
    if (!t_o || !h_o || t_o->length != h_o->length) {
        fprintf(stderr, "Error: Invalid temperature or h_o vectors\n");
        return T_nh;
    }

    const double t_i = 0.0;  // Inner temperature (°C)
    const int n_el = (int)(d_ins / dx);  // Number of elements
    const int nodes = n_el + 1;          // Number of nodes
    const int nr_hour = t_o->length;     // Number of hours
    const int nh = (int)(3600.0 / dt);   // Time steps per hour

    // Allocate result matrix
    T_nh.data = malloc(nodes * sizeof(double*));
    for (int i = 0; i < nodes; i++) {
        T_nh.data[i] = malloc(nr_hour * sizeof(double));
    }
    T_nh.rows = nodes;
    T_nh.cols = nr_hour;

    // Initialize temperature distribution (shared)
    Vector T_n_initial = {calloc(nodes, sizeof(double)), nodes};

    // Main loop over hours - PARALLELIZED
    #pragma omp parallel
    {
        // Each thread gets private workspace
        Vector T_n_local = {malloc(nodes * sizeof(double)), nodes};
        Vector a = {malloc(nodes * sizeof(double)), nodes};
        Vector b = {malloc(nodes * sizeof(double)), nodes};
        Vector c = {malloc(nodes * sizeof(double)), nodes};
        Vector d = {malloc(nodes * sizeof(double)), nodes};

        #pragma omp for schedule(dynamic)
        for (int h = 0; h < nr_hour; h++) 
        {
            // Start from initial condition each hour
            memcpy(T_n_local.values, T_n_initial.values, nodes * sizeof(double));

            const double dFo = D * dt / (dx * dx);
            const double dBio_i = h_i * dx / lam_i;
            const double dBio_o = h_o->values[h] * dx / lam_i;

            // Time steps within one hour
            for (int step = 0; step < nh; step++) 
            {
                // Set up tridiagonal system
                b.values[0] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o;
                c.values[0] = -2.0 * dFo;
                d.values[0] = T_n_local.values[0] + 2.0 * dFo * dBio_o * t_o->values[h];

                for (int j = 1; j < nodes-1; j++) 
                {
                    a.values[j] = -dFo;
                    b.values[j] = 1.0 + 2.0 * dFo;
                    c.values[j] = -dFo;
                    d.values[j] = T_n_local.values[j];
                }

                a.values[nodes-1] = -2.0 * dFo;
                b.values[nodes-1] = 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i;
                d.values[nodes-1] = T_n_local.values[nodes-1] + 2.0 * dFo * dBio_i * t_i;

                // Solve system
                Vector solution = solve_tdma(&a, &b, &c, &d);
                memcpy(T_n_local.values, solution.values, nodes * sizeof(double));
                free(solution.values);
            }

            // Store results for this hour
            #pragma omp critical
            {
                for (int j = 0; j < nodes; j++) {
                    T_nh.data[j][h] = T_n_local.values[j];
                }
            }
        }

        // Free thread-local memory
        free(T_n_local.values);
        free(a.values);
        free(b.values);
        free(c.values);
        free(d.values);
    }

    // Cleanup shared memory
    free(T_n_initial.values);

    return T_nh;
}
*/

/*
CSVData* read_csv_with_encoding_parallel(const char* file_path, int* columns_to_keep, int num_cols_to_keep) {
    const char* encoding = detect_encoding(file_path);
    printf("\nAssuming encoding: %s\n", encoding);

    FILE* file = fopen(file_path, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    // First pass to count rows and columns (this remains sequential)
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
    if (!csv_data->data) {
        free(cols_kept);
        free(csv_data);
        fclose(file);
        return NULL;
    }

    // Pre-allocate all row pointers to simplify parallel access
    for (int i = 0; i < rows; i++) {
        csv_data->data[i] = calloc(cols_kept_count, sizeof(char*));
        if (!csv_data->data[i]) {
            // Cleanup previously allocated rows
            for (int j = 0; j < i; j++) {
                free(csv_data->data[j]);
            }
            free(csv_data->data);
            free(cols_kept);
            free(csv_data);
            fclose(file);
            return NULL;
        }
    }

    // Read all lines into a buffer for parallel processing
    rewind(file);
    char** line_buffer = calloc(rows, sizeof(char*));
    if (!line_buffer) {
        for (int i = 0; i < rows; i++) {
            free(csv_data->data[i]);
        }
        free(csv_data->data);
        free(cols_kept);
        free(csv_data);
        fclose(file);
        return NULL;
    }

    // Read all lines sequentially first
    for (int i = 0; i < rows; i++) {
        line_buffer[i] = calloc(100000, sizeof(char));
        if (!line_buffer[i] || !fgets(line_buffer[i], 100000, file)) {
            // Clean up on error
            for (int j = 0; j <= i; j++) {
                free(line_buffer[j]);
            }
            free(line_buffer);
            for (int j = 0; j < rows; j++) {
                free(csv_data->data[j]);
            }
            free(csv_data->data);
            free(cols_kept);
            free(csv_data);
            fclose(file);
            return NULL;
        }
    }

    // Process lines in parallel
    bool error_occurred = false;
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < rows; i++) {
            if (error_occurred) continue; // Skip if an error has occurred

            char** all_columns = parse_csv_line(line_buffer[i], max_columns);
            if (!all_columns) {
                #pragma omp critical
                {
                    error_occurred = true;
                }
                continue;
            }

            // Process columns for this row
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
    }

    // Clean up the line buffer
    for (int i = 0; i < rows; i++) {
        free(line_buffer[i]);
    }
    free(line_buffer);

    fclose(file);

    // Check if an error occurred during parallel processing
    if (error_occurred) {
        for (int i = 0; i < rows; i++) {
            free_csv_row(csv_data->data[i], cols_kept_count);
        }
        free(csv_data->data);
        free(cols_kept);
        free(csv_data);
        return NULL;
    }

    return csv_data;
}
*/

/*
NumericArray convert_string_array_to_double_parallel(StringArray* string_arr) {
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
    
    #pragma omp parallel for
    for (int i = 0; i < string_arr->length; i++) {
        char* cleaned = trim_whitespace(string_arr->strings[i]);
        char* endptr;
        result.values[i] = strtod(cleaned, &endptr);
        
        // Check if conversion was successful
        result.is_valid[i] = (*cleaned != '\0' && *endptr == '\0');
        if (!result.is_valid[i]) {
            #pragma omp critical
            {
                fprintf(stderr, "Warning: Could not convert '%s' to double\n", string_arr->strings[i]);
            }
            result.values[i] = 0.0; // Default value for invalid
        }
    }
    
    return result;
}
*/

/*
Vector solve_tdma_gsl(const Vector* a, const Vector* b, const Vector* c, const Vector* d) {
    Vector x = {NULL, b->length};
    if (a->length != b->length || b->length != c->length || c->length != d->length) {
        fprintf(stderr, "Error: TDMA vectors must be same length\n");
        return x;
    }

    const int n = b->length;
    x.values = malloc(n * sizeof(double));

    gsl_vector *diag = gsl_vector_alloc(n);
    gsl_vector *superdiag = gsl_vector_alloc(n-1);
    gsl_vector *subdiag = gsl_vector_alloc(n-1);
    gsl_vector *rhs = gsl_vector_alloc(n);

    // Copy data into GSL vectors
    for (int i = 0; i < n-1; i++) {
        gsl_vector_set(subdiag, i, a->values[i+1]);  // a[1..n-1]
        gsl_vector_set(superdiag, i, c->values[i]);   // c[0..n-2]
    }
    for (int i = 0; i < n; i++) {
        gsl_vector_set(diag, i, b->values[i]);        // b[0..n-1]
        gsl_vector_set(rhs, i, d->values[i]);         // d[0..n-1]
    }

    // Solve the system
    gsl_linalg_solve_tridiag(diag, superdiag, subdiag, rhs, rhs);

    // Copy result back
    for (int i = 0; i < n; i++) {
        x.values[i] = gsl_vector_get(rhs, i);
    }

    // Free GSL vectors
    gsl_vector_free(diag);
    gsl_vector_free(superdiag);
    gsl_vector_free(subdiag);
    gsl_vector_free(rhs);

    return x;
}

Matrix transient1D_gsl(const Vector* t_o, const Vector* h_o, 
    double d_ins, double lam_i, double D,
    double dx, double dt, double h_i) 
{
    Matrix T_nh = {NULL, 0, 0};
    
    // Validate inputs
    if (!t_o || !h_o || t_o->length != h_o->length) {
        fprintf(stderr, "Error: Invalid temperature or h_o vectors\n");
        return T_nh;
    }

    const double t_i = 0.0;  // Inner temperature (°C)
    const int n_el = (int)(d_ins / dx);  // Number of elements
    const int nodes = n_el + 1;          // Number of nodes
    const int nr_hour = t_o->length;     // Number of hours
    const int nh = (int)(3600.0 / dt);   // Time steps per hour

    // Allocate result matrix
    T_nh.data = malloc(nodes * sizeof(double*));
    for (int i = 0; i < nodes; i++) {
        T_nh.data[i] = malloc(nr_hour * sizeof(double));
    }
    T_nh.rows = nodes;
    T_nh.cols = nr_hour;

    // Initialize temperature distribution
    Vector T_n = {malloc(nodes * sizeof(double)), nodes};
    for (int i = 0; i < nodes; i++) {
        T_n.values[i] = 0.0;
    }

    // Pre-allocate all GSL vectors once
    gsl_vector *diag = gsl_vector_alloc(nodes);
    gsl_vector *superdiag = gsl_vector_alloc(nodes-1);
    gsl_vector *subdiag = gsl_vector_alloc(nodes-1);
    gsl_vector *rhs = gsl_vector_alloc(nodes);
    gsl_vector *solution = gsl_vector_alloc(nodes);

    for (int h = 0; h < nr_hour; h++) {
    const double dFo = D * dt / (dx * dx);
    const double dBio_i = h_i * dx / lam_i;
    const double dBio_o = h_o->values[h] * dx / lam_i;

    for (int step = 0; step < nh; step++) {
        // Directly fill GSL vectors (no intermediate Vector structs)
        gsl_vector_set(diag, 0, 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o);
        gsl_vector_set(superdiag, 0, -2.0 * dFo);
        gsl_vector_set(rhs, 0, T_n.values[0] + 2.0 * dFo * dBio_o * t_o->values[h]);

        for (int j = 1; j < nodes-1; j++) {
        gsl_vector_set(subdiag, j-1, -dFo);
        gsl_vector_set(diag, j, 1.0 + 2.0 * dFo);
        gsl_vector_set(superdiag, j, -dFo);
        gsl_vector_set(rhs, j, T_n.values[j]);
        }

        gsl_vector_set(subdiag, nodes-2, -2.0 * dFo);
        gsl_vector_set(diag, nodes-1, 1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i);
        gsl_vector_set(rhs, nodes-1, T_n.values[nodes-1] + 2.0 * dFo * dBio_i * t_i);

        // Solve
        gsl_linalg_solve_tridiag(diag, superdiag, subdiag, rhs, solution);

        // Update T_n directly from GSL vector
        for (int j = 0; j < nodes; j++) {
            T_n.values[j] = gsl_vector_get(solution, j);
        }
    }

        // Store results for this hour
        for (int j = 0; j < nodes; j++) {
            T_nh.data[j][h] = T_n.values[j];
        }
    }

    // Free GSL vectors
    gsl_vector_free(diag);
    gsl_vector_free(superdiag);
    gsl_vector_free(subdiag);
    gsl_vector_free(rhs);
    gsl_vector_free(solution);

    // Cleanup
    free(T_n.values);

    return T_nh;
}
*/

/*
int find_index_parallel(CSVData* data, int column_index, const char* value) {
    if (!data || column_index < 0 || column_index >= data->cols_kept) {
        return -1;
    }
    
    int found_index = -1;
    
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < data->rows; i++) {
            // Skip if already found by another thread
            if (found_index != -1) continue;
            
            if (strcmp(data->data[i][column_index], value) == 0) {
                #pragma omp critical
                {
                    // Only update if this is the first occurrence found
                    if (found_index == -1 || i < found_index) {
                        found_index = i;
                    }
                }
                
                // Optional: Break out of the loop
                #pragma omp cancel for
            }
        }
    }
    
    return found_index;
}
*/