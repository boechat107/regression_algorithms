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

/**
 * Partition function of the quicksort algorithm.
 */
static int
partition (_Norm *array, int p, int r)
{
    int index;
    double norm;
    _Norm *x = &array[r];
    int i = p - 1;
    for (int j = p; j < r; j++)
    {
        if (array[j].norm <= x->norm)
        {
            i++;
            /* Exchange A[i] <-> A[j]. */
            index = array[j].index;
            norm = array[j].norm;
            array[j].index = array[i].index;
            array[j].norm = array[i].norm;
            array[i].index = index;
            array[i].norm = norm;
        }
    }
    /* Exchange A[i+1] <-> A[r]. */
    i++;
    index = array[i].index;
    norm = array[i].norm;
    array[i].index = array[r].index;
    array[i].norm = array[r].norm;
    array[r].index = index;
    array[r].norm = norm;
    return i;
}


/**
 *    Quicksort.
 *    The array is sorted in the ascend order of norm.
 */
static void
sort_by_norm (_Norm *array, int p, int r)
{
    if (p < r)
    {
        int q = partition(array, p, r);
        sort_by_norm(array, p, q-1);
        sort_by_norm(array, q+1, r);
    }
}

void
printNorm (_Norm *array, int size)
{
    for (int i = 0; i < size; i++)
        printf("Index %d and Norm %f\n", array[i].index, array[i].norm);
}
/**
 *    \brief brief description
 *    \param Nome var description
 *    \return Nome var description    
 */
Signals_matrix*
memory_vector_selection (Signals_matrix* trainData, int num_vectors)
{
    if (num_vectors >= trainData->nvectors)
        return trainData;
    /* Output variable initialization. */
    Signals_matrix* memory = MALLOC(Signals_matrix,1);
    memory->nvars = trainData->nvars;
    memory->nvectors = 0;
    memory->vectors = MALLOC(Signals_node,num_vectors);
    /* Iterators on bands, vectors of the training data in the current band and the column
     * of the data. */
    int band;
    int vec;
    int var;
    /* Just a table to mark the training vectors that were selected to the memory list. */
    char selected_vectors[trainData->nvectors];
    for (vec = 0; vec < trainData->nvectors; vec++) 
        selected_vectors[vec] = 'n';
    /* The training data is divided in a num_bands parts. On each part, the vectors with
     * the minimum and maximum values of each variable are selected. */
    int num_bands = (int)floor(((float)num_vectors) / ((float)(2*trainData->nvars)));
    int band_size = (int)floor(((float)trainData->nvectors) / ((float)num_bands));
    for (band = 0; band < num_bands; band++)
    {
        for (var = 0; var < trainData->nvars; var++)
        {
            double maxValue = -1e10,
                   minValue = 1e10;
            /* The respective indexes of the array of maximum (minimum) values. */
            int minIndex = -1,
                maxIndex = -1;
            /* Loop to search the maximum (minimum) vector. The index of the maximal
             * (minimal) vector is saved. */
            for (vec = band*band_size; vec < (band+1)*band_size; vec++)
            {
                if (trainData->vectors[vec].data[var] > maxValue)
                {
                    maxValue = trainData->vectors[vec].data[var];
                    maxIndex = vec;
                }
                if (trainData->vectors[vec].data[var] < minValue)
                {
                    minValue = trainData->vectors[vec].data[var];
                    minIndex = vec;
                }
            }
            /* Here it's verified occurrences of duplicated entries in memory. Every time
             * that a vector is selected for memory, its index is marked as used. So,
             * before adding the vector to the memory, a list of marked vectors is
             * analysed. */
            if (selected_vectors[maxIndex] == 'n')
            {
                selected_vectors[maxIndex] = 'y';
                double *data = MALLOC(double,trainData->nvars);
                copyDoubleValues(
                        trainData->vectors[maxIndex].data,
                        data,
                        trainData->nvars);
                memory->vectors[memory->nvectors].data = data;
                memory->vectors[memory->nvectors].index = maxIndex;
                memory->nvectors++;
            }
            if (selected_vectors[minIndex] == 'n')
            {
                selected_vectors[minIndex] = 'y';
                double *data = MALLOC(double,trainData->nvars);
                copyDoubleValues(
                        trainData->vectors[minIndex].data,
                        data,
                        trainData->nvars);
                memory->vectors[memory->nvectors].data = data;
                memory->vectors[memory->nvectors].index = minIndex;
                memory->nvectors++;
            }
        }
    }
    /* Normally, after the selection of the min-max vectors, the desired number of memory
     * vectors is not obtained. So, it is used the vector ordering method to obtain the
     * rest of memory vectors. */
    int nMissingVectors = num_vectors - memory->nvectors;
    if (nMissingVectors > 0)
    {
        int normArraySize = trainData->nvectors - memory->nvectors;
        _Norm norm_index[normArraySize];
        /* A array of the vector norm values are created. Just the vectors that were not
         * added to the memory are analysed. After that, this array must be
         * sorted in a ascend order of norm. */
        int i = 0;
        for (vec = 0; vec < trainData->nvectors; vec++)
        {
            if (selected_vectors[vec] == 'y')
                continue;
            double norm = 0.0;
            /* Calculates the norm of the vector. */
            for (var = 0; var < trainData->nvars; var++)
                norm += pow(trainData->vectors[vec].data[var], 2);
            norm = sqrt(norm);
            norm_index[i].norm = norm;
            norm_index[i].index = vec;
            i++;
        }
        assert(i == normArraySize);
        /* Sorts the array of norm values in ascend order. */
//        printNorm(norm_index, normArraySize);
//        printf("\n\n");
        sort_by_norm(norm_index, 0, normArraySize-1);
//        printNorm(norm_index, normArraySize);
        int partitionSize = normArraySize / nMissingVectors;
        for (vec = 0; memory->nvectors < num_vectors; vec += partitionSize)
        {
            selected_vectors[norm_index[vec].index] = 'y';
            double *data = MALLOC(double,trainData->nvars);
            copyDoubleValues(
                    trainData->vectors[norm_index[vec].index].data,
                    data,
                    trainData->nvars);
            memory->vectors[memory->nvectors].data = data;
            memory->vectors[memory->nvectors].index = norm_index[vec].index;
            memory->nvectors++;
        }
    }
    return memory;
}



