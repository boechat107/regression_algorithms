#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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
    for (int i = 0; i < input->nvectors; i++)
    {
        free(input->vectors[i].data);
    }
    free(input->vectors);
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
    Signals_node *memory_vec;
    int vec = 0;
    /* Calculation of the vector of distances between the query and the memory vectors. */
    for (vec = 0; vec < memory->nvectors; vec++)
    {
        memory_vec = &memory->vectors[vec];
        distances[vec] = 0.0;
        for (var = 0; var < memory->nvars; var++)
            distances[vec] += pow(memory_vec->data[var] - query[var], 2.0);
        distances[vec] = sqrt(distances[vec]);
        /* Calculation of the weights using the Gaussian kernel. */
        weights[vec] = (1/sqrt(2 * PI * pow(bandwidth, 2)))
            * exp(- pow(distances[vec], 2) / (2*pow(bandwidth, 2)));
        weights_sum += weights[vec];
    }
    /* Calculation of the prediction vector, the function's output. */
    for (var = 0; var < memory->nvars; var++) 
    {
        prediction[var] = 0;
        vec = 0;
        for (vec = 0; vec < memory->nvectors; vec++)
        {
            memory_vec = &memory->vectors[vec];
            prediction[var] += weights[vec] * memory_vec->data[var];
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
    memory->vectors = MALLOC(Signals_node,num_vectors);
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
    int num_bands = (int)floor(((float)num_vectors) / ((float)(2*train_data->nvars)));
    int band_size = (int)floor(((float)train_data->nvectors) / ((float)num_bands));
    for (band = 0; band < num_bands; band++)
    {
        for (var = 0; var < train_data->nvars; var++)
        {
            double maxValue = -1e10,
                   minValue = 1e10;
            /* The respective indexes of the array of maximum (minimum) values. */
            int minIndex = -1,
                maxIndex = -1;
            /* Loop to search the maximum (minimum) vector. */
            for (vec = band*band_size; vec < (band+1)*band_size; vec++)
            {
                if (train_data->vectors[vec].data[var] > maxValue)
                {
                    maxValue = train_data->vectors[vec].data[var];
                    maxIndex = vec;
                }
                if (train_data->vectors[vec].data[var] < minValue)
                {
                    minValue = train_data->vectors[vec].data[var];
                    minIndex = vec;
                }
            }
            /* Vectors with maximum values are added to the memory list, or not. */
            if (selected_vectors[maxIndex] == 'n')
            {
                selected_vectors[maxIndex] = 'y';
                double *data = MALLOC(double,train_data->nvars);
                copyDoubleValues(
                        train_data->vectors[maxIndex].data,
                        data,
                        train_data->nvars);
                memory->vectors[memory->nvectors].data = data;
                memory->vectors[memory->nvectors].index = maxIndex;
                memory->nvectors++;
            }
            /* Vectors with minimum values are added to the memory list, or not. */
            if (selected_vectors[minIndex] == 'n')
            {
                selected_vectors[minIndex] = 'y';
                double *data = MALLOC(double,train_data->nvars);
                copyDoubleValues(
                        train_data->vectors[minIndex].data,
                        data,
                        train_data->nvars);
                memory->vectors[memory->nvectors].data = data;
                memory->vectors[memory->nvectors].index = minIndex;
                memory->nvectors++;
            }
        }
    }
    /* Normally, after the selection of the min-max vectors, the desired number of memory
     * vectors is not obtained. So, it is used the vector ordering method to obtain the
     * rest of memory vectors. */
    int nMissingVectors = memory->nvectors - num_vectors;
    if (nMissingVectors > 0)
    {
        int normArraySize = train_data->nvectors - memory->nvectors;
        _Norm norm_index[normArraySize];
        /* A array of the vector norm values are created. Just the vectors that were not
         * added to the memory are analysed. After that, this array must be
         * sorted in a ascend order of norm. */
        int i = 0;
        for (vec = 0; vec < train_data->nvectors; vec++)
        {
            if (selected_vectors[vec] == 'y')
                continue;
            double norm = 0.0;
            /* Calculates the norm of the vector. */
            for (var = 0; var < train_data->nvars; var++)
                norm += pow(train_data->vectors[vec].data[var], 2);
            norm = sqrt(norm);
            norm_index[i].norm = norm;
            norm_index[i].index = vec;
            i++;
        }
        assert(i == normArraySize);
        /* Sorts the array of norm values in ascend order. */
//        sort_by_norm(norm_index, normArraySize);
        int partitionSize = normArraySize / nMissingVectors;
        for (vec = 0; memory->nvectors < num_vectors; vec += partitionSize)
        {
            selected_vectors[norm_index[vec].index] = 'y';
            memory->vectors[memory->nvectors].index = norm_index[vec].index;
            copyDoubleValues(
                    train_data->vectors[norm_index[vec].index].data,
                    memory->vectors[memory->nvectors].data,
                    train_data->nvars);
            memory->nvectors++;
        }
    }
    return memory;
}



