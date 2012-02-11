#include <stdio.h>
#include <stdlib.h>
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
    if (matrix == NULL)
        return 0;
    fclose(filep);
    Signals_matrix *memory = memory_vector_selection(matrix, 72);
//    print_signals_matrix(memory);

    free_signals_matrix(matrix);
    free_signals_matrix(memory);
    return 0;
}
