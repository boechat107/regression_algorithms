#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "aakr.h"

#define MALLOC(type,n) (type*)malloc(n*sizeof(type))

Signals_matrix *
read_data_file (FILE *filep)
{
    char var_value[64];
    int var_count = 0, 
        vec_count = 0,
        eof;
    /* Loop to count the number of variables. */
    while (1)
    {
        char line_end;
        eof = fscanf(filep, "%s%c", var_value, &line_end);
        if (eof == EOF)
            break;
        var_count++;
        if (line_end == '\n')
            break;
    }
    vec_count++;
    /* Loop to count the number of vectors and to verify if the input file is consistent.
     * */
    while (1)
    {
        char line_end;
        eof = fscanf(filep, "%s%c", var_value, &line_end);
        if (eof == EOF)
            break;
        for (int i = 1; i < var_count-1; i++)
        {
            eof = fscanf(filep, "%s%c", var_value, &line_end);
            if (eof == EOF || line_end == '\n')
            {
                printf("ERROR: Input file is not consistent! Line %d\n", vec_count);
                return NULL;
            }
        }
        eof = fscanf(filep, "%s%c", var_value, &line_end);
        if (eof == EOF || line_end != '\n')
        {
            printf("ERROR: Input file is not consistent! Line %d\n", vec_count);
            return NULL;
        }
        vec_count++;
    }
    /* Starts the file from the beginning again. */
    fseek(filep, 0, SEEK_SET);
    /* Output variable initialization. */
    Signals_matrix *output = MALLOC(Signals_matrix,1);
    output->nvars = var_count;
    output->nvectors = 0;
    output->vectors = MALLOC(Signals_node,vec_count);
    double *vector_data;
    for (int i = 0; i < vec_count; i++)
    {
        for (int var = 0; var < var_count; var++)
        {
            fscanf(filep, "%s", var_value);
            if (var == 0)
            {
                vector_data = MALLOC(double,var_count);
                output->vectors[i].index = i;
                output->vectors[i].data = vector_data;
                output->nvectors++;
            }
            vector_data[var] = atof(var_value);
        }
    }
    return output;
}


Signals_matrix *
scale_data (Signals_matrix *input, double *mean, double *std)
{
    Signals_matrix *output = signals_matrix_init(input->nvars, input->nvectors);
    output->nvectors = input->nvectors;
    for (int vec = 0; vec < input->nvectors; vec++)
    {
        output->vectors[vec].data = MALLOC(double,output->nvars);
        for (int var = 0; var < input->nvars; var++)
        {
            output->vectors[vec].data[var] = (input->vectors[vec].data[var] - mean[var])
                / std[var];
        }
    }
    return output;
}

double **
calculate_mean_std (Signals_matrix *input)
{
    double **output = MALLOC(double *,2);
    output[0] = MALLOC(double, input->nvars);
    output[1] = MALLOC(double, input->nvars);
    for (int i = 0; i < input->nvars; i++)
    {
        output[0][i] = 0.0;
        output[1][i] = 0.0;
    }
    for (int vec = 0; vec < input->nvectors; vec++)
    {
        for (int var = 0; var < input->nvars; var++)
        {
            output[0][var] += input->vectors[vec].data[var];
            output[1][var] += pow(input->vectors[vec].data[var], 2);
        }
    }
    for (int var = 0; var < input->nvars; var++)
    {
        /* Mean */
        output[0][var] = output[0][var] / (double)input->nvectors;
        /* Variance */
        output[1][var] = (output[1][var] / (double)input->nvectors) 
            - pow(output[0][var], 2);
        /* Standard deviation */
        output[1][var] = sqrt(output[1][var]);
    }
    return output;
}

void
print_signals_matrix (Signals_matrix *input)
{
    printf("%d\n", input->nvectors);
    for (int i = 0; i < input->nvectors; i++)
    {
        for (int var = 0; var < input->nvars; var++)
            printf("%f ",input->vectors[i].data[var]);
        printf("\n");
    }
}


int
main(int argc, const char *argv[])
{
    FILE *filep = fopen(argv[1], "r");
    if (filep == NULL) 
    {
        fprintf(stderr, "ERROR: It is not possible to open the file %s.\n\n", argv[1]);
        return 1;
    }
    Signals_matrix *matrix = read_data_file(filep);
    fclose(filep);
    assert(matrix != NULL);

    filep = fopen(argv[2], "r");
    if (filep == NULL) 
    {
        fprintf(stderr, "ERROR: It is not possible to open the file %s.\n\n", argv[2]);
        return 1;
    }
    Signals_matrix *testData = read_data_file(filep);
    fclose(filep);
    assert(testData != NULL);

    double **statistics = calculate_mean_std(matrix);
    /* Data scaling to zero mean and unit standard deviation. */
    Signals_matrix *trainNormData = scale_data(matrix, statistics[0], statistics[1]);
    Signals_matrix *testNormData = scale_data(testData, statistics[0], statistics[1]);

    Signals_matrix *memory = memory_vector_selection(trainNormData, 300);
//    print_signals_matrix(memory);

    filep = fopen("out.txt", "w+");
    double *prediction = MALLOC(double,matrix->nvars);
    for (int query = 0; query < testData->nvectors; query++)
    {
        aakr_prediction(memory, .6, testNormData->vectors[query].data, prediction);
        for (int var = 0; var < matrix->nvars; var++)
        {
            prediction[var] = prediction[var]*statistics[1][var] + statistics[0][var];
            fprintf(filep, "%.9f\t", prediction[var]);
        }
        fprintf(filep, "\n");
    }
    fclose(filep);

    free(prediction);
    free(statistics[0]);
    free(statistics[1]);
    free(statistics);
    free_signals_matrix(matrix);
    free_signals_matrix(memory);
    free_signals_matrix(testData);
    free_signals_matrix(trainNormData);
    free_signals_matrix(testNormData);
    return 0;
}
