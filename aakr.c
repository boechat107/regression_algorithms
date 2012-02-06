#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "aakr.h"

#ifdef DEBUG
#define print_debug(msg); printf(msg);
#else
#define print_debug(msg); 
#endif

#define MALLOC(type,n) (type*)malloc(n*sizeof(type))
#define PI 3.14159265

/* ----------------------------------------------------------------------------------------- 
 * Auxiliar functions
 * ----------------------------------------------------------------------------------------- 
 */
static void
initMaxArrays (double* array, int size)
{
    int i;
    for (i = 0; i < size; i++) 
    {
        array[i] = 0;
    }
}


static void
initMinArrays (double* array, int size)
{
    int i;
    for (i = 0; i < size; i++) 
    {
        array[i] = 1e10;
    }
}


static Signals_node*
new_signals_node (double* array, int index)
{
    Signals_node* node = MALLOC(Signals_node,1);
    node->data = array;
    node->index = index;
    node->next = NULL;
    return node;
}


static void
copyDoubleValues (double* orig, double* dest, int size)
{
    int i;
    for (i = 0; i < size; i++) 
    {
        dest[i] = orig[i];
    }
}


/* ----------------------------------------------------------------------------------------- 
 * Main functions
 * ----------------------------------------------------------------------------------------- 
 */
void
free_signals_matrix (Signals_matrix *input)
{
    Signals_node *marker = input->head;
    Signals_node *liberator = input->head;
    while (marker != NULL)
    {
        marker = marker->next;
        free(liberator->data);
        free(liberator);
        liberator = marker;
    }
    free(input);
}


double*
aakr_prediction (Signals_matrix* memory, double bandwidth, double* query)
{
    double distances[memory->nvectors];
    double weights[memory->nvectors];
    double *prediction = MALLOC(double,memory->nvars);
    double weights_sum = 0;
    int var;
    Signals_node *memory_vec = memory->head;
    int vec = 0;
    /* Calculation of the vector of distances between the query and the memory vectors. */
    while (memory_vec != NULL)
    {
        distances[vec] = 0.0;
        for (var = 0; var < memory->nvars; var++)
            distances[vec] += pow(memory_vec->data[var] - query[var], 2.0);
        distances[vec] = sqrt(distances[vec]);
        /* Calculation of the weights using the Gaussian kernel. */
        weights[vec] = (1/sqrt(2 * PI * pow(bandwidth, 2)))
            * exp(- pow(distances[vec], 2) / (2*pow(bandwidth, 2)));
        weights_sum += weights[vec];
        vec++;
        memory_vec = memory_vec->next;
    }
    /* Calculation of the prediction vector, the function's output. */
    for (var = 0; var < memory->nvars; var++) 
    {
        prediction[var] = 0;
        vec = 0;
        memory_vec = memory->head;
        while (memory_vec != NULL)
        {
            prediction[var] += weights[vec] * memory_vec->data[var];
            vec++;
            memory_vec = memory_vec->next;
        }
        prediction[var] = prediction[var] / weights_sum;
    }
    return prediction;
}


struct _Norm
{
    int index;
    double norm;
};
typedef struct _Norm _Norm;

Signals_matrix*
memory_vector_selection (Signals_matrix* train_data, int num_vectors)
{
    if (num_vectors >= train_data->nvectors)
        return train_data;
    /* Output variable initialization. */
    Signals_matrix* memory = MALLOC(Signals_matrix,1);
    memory->nvars = train_data->nvars;
    memory->nvectors = 0;
    memory->head = NULL;
    /* Iterators on bands, vectors of the training data in the current band and the column
     * of the data. */
    int band;
    int vec;
    int var;
    /* Just a table to mark the training vectors that were selected to the memory list. */
    char selected_vectors[train_data->nvectors];
    for (vec = 0; vec < train_data->nvectors; vec++) 
        selected_vectors[vec] = 'n';
    /* The training data is divided in a num_bands parts. On each part, the vectors with
     * the minimum and maximum values of each variable are selected. */
    int num_bands = (int)floor(((float)num_vectors) / ((float)(2*train_data->nvars)) + .5);
    int band_size = (int)floor(((float)train_data->nvectors) / ((float)num_bands));
    /* Values of the variables where one is the maximum (minimum) of a variable. Although
     * only one variable is analysed, the entire data vector is stored for the ouput
     * structure, the memory. */
    double *minValues, *maxValues;
    /* The respective indexes of the array of maximum (minimum) values. */
    int minIndex = -1,
        maxIndex = -1;
    /* Pointer to iterate on the training data. */
    Signals_node *matrix_vec = train_data->head;
    /* A marker of the band's beginning. It's used to restar the matrix_vec for every
     * variable. */
    Signals_node *band_start = train_data->head;
    /* A marker to reduce the number of memory vector that must be analysed to remove
     * duplicated entries. Just the nodes starting at this pointer are analysed. */
    Signals_node *first_band_memory = NULL;
    /* A marker to insert new memory nodes to the memory list. */
    Signals_node *memory_marker = NULL;
    /* Iterator to analyse duplicated entries in the memory list. */
    Signals_node *mem_rep_marker;
    /* A new entry of the memory list. */
    Signals_node *new_node;
    /* Flags to indicate if the selected training vector is duplicated. */
    char flag_max = 'y';
    char flag_min = 'y';
    for (band = 0; band < num_bands; band++)
    {
        for (var = 0; var < train_data->nvars; var++)
        {
            minValues = MALLOC(double,train_data->nvars);
            maxValues = MALLOC(double,train_data->nvars);
            initMinArrays(minValues, train_data->nvars);
            initMaxArrays(maxValues, train_data->nvars);
            /* Searches the maximum (minimum) vector. */
            matrix_vec = band_start;
            for (vec = 0; vec < band_size; vec++)
            {
                if (matrix_vec->data[var] > maxValues[var])
                {
                    copyDoubleValues(matrix_vec->data, maxValues, train_data->nvars);
                    maxIndex = matrix_vec->index;
                }
                if (matrix_vec->data[var] < minValues[var])
                {
                    copyDoubleValues(matrix_vec->data, minValues, train_data->nvars);
                    minIndex = matrix_vec->index;
                }
                matrix_vec = matrix_vec->next;
            }
            /* Duplicated entries verification. */
            mem_rep_marker = first_band_memory;
            while(mem_rep_marker != NULL)
            {
                if (mem_rep_marker->index == maxIndex)
                    flag_max = 'n';
                if (mem_rep_marker->index == minIndex)
                    flag_min = 'n';
                mem_rep_marker = mem_rep_marker->next;
            }
            /* Vectors with maximum values are added to the memory list, or not. */
            if (flag_max == 'y')
            {
                selected_vectors[maxIndex] = 'y';
                new_node = new_signals_node(maxValues, maxIndex);
                if (memory_marker != NULL)
                {
                    memory_marker->next = new_node;
                    memory_marker = memory_marker->next;
                    memory->nvectors++;
                }
                else
                {
                    memory->head = new_node;
                    memory_marker = new_node;
                }
            }
            else
                free(maxValues);
            /* Vectors with minimum values are added to the memory list, or not. */
            if (flag_min == 'y' && minIndex != maxIndex)
            {
                selected_vectors[minIndex] = 'y';
                new_node = new_signals_node(minValues, minIndex);
                if (memory_marker != NULL)
                {
                    memory_marker->next = new_node;
                    memory_marker = memory_marker->next;
                    memory->nvectors++;
                }
                else
                {
                    memory->head = new_node;
                    memory_marker = new_node;
                }
            }
            else
                free(minValues);
            flag_min = 'y';
            flag_max = 'y';
        }
        band_start = matrix_vec;
        first_band_memory = memory_marker;
    }
    /* Normally, after the selection of the min-max vectors, the desired number of memory
     * vectors is not obtained. So, it is used the vector ordering method to obtain the
     * rest of memory vectors. */
    if (memory->nvectors < num_vectors)
    {
        double index_data[train_data->nvectors][train_data->nvars];
        _Norm norm_index[train_data->nvectors];
        matrix_vec = train_data->head;
        /* List of the vector norm values are created. After that, this list must be
         * sorted in a ascend order of norm. */
        while (matrix_vec != NULL)
        {
            double norm = 0.0;
            /* Calculates the norm of the vector. */
            for (var = 0; var < train_data->nvars; var++)
                norm += pow(matrix_vec->data[var], 2);
            norm = sqrt(norm);
            norm_index[matrix_vec->index].norm = norm;
            norm_index[matrix_vec->index].index = matrix_vec->index;
            matrix_vec = matrix_vec->next;
        }
    }
    return memory;
}



