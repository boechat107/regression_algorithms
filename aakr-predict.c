#include <stdio.h>
#include <stdlib.h>
#include "aakr.h"

#define MALLOC(type,n) (type*)malloc(n*sizeof(type))

Signals_matrix *
read_data_file (FILE *filep)
{
    char var_value[64];
    int var_count = 0;
    int eof;
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
    /* Starts the file from the beginning again. */
    fseek(filep, 0, SEEK_SET);
    /* Output variable initialization. */
    Signals_matrix *output = MALLOC(Signals_matrix,1);
    output->nvars = var_count;
    output->nvectors = 0;
    output->head = NULL;
    Signals_node *new_node,
                 *list_marker;
    double *vector_data;
    int var = 0;
    while (1)
    {
        eof = fscanf(filep, "%s", var_value);
        if (eof == EOF)
            break;
        var = var % var_count;
        if (var == 0)
        {
            vector_data = MALLOC(double,var_count);
            new_node = MALLOC(Signals_node,1);
            new_node->index = output->nvectors;
            new_node->data = vector_data;
            new_node->next = NULL;
            if (output->nvectors == 0)
                output->head = new_node;
            else
                list_marker->next = new_node;
            list_marker = new_node;
            output->nvectors++;
        }
        vector_data[var] = atof(var_value);
        var++;
    }
    return output;
}


void
print_signals_matrix (Signals_matrix *input)
{
    printf("%d\n", input->nvectors);
    Signals_node *marker = input->head;
    while (marker != NULL)
    {
        int var;
        for (var = 0; var < input->nvars; var++)
            printf("%f ", marker->data[var]);
        printf("\n");
        marker = marker->next;
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
    Signals_matrix *memory = memory_vector_selection(matrix, 72);
    print_signals_matrix(memory);

    free_signals_matrix(matrix);
    free_signals_matrix(memory);
    return 0;
}
