/**
 *    @file       aakr.h
 *    @brief      Functions to make predictions using the AAKR model.
 *    @author     Andre A. Boechat <boechat@das.ufsc.br>
 *    @date       February 03, 2012
 *
 *    Reference: 
 *    "Technical Review of On-Line Monitoring Techniques for Performance Assessment Vol.
 *    2: Theoretical Issues"
 *
 *    TODO: verify the necessity of the index in the Signals_node structure.
 */

#ifndef AUTOASSOCIATIVE_KERNEL_REGRESSION_H
#define AUTOASSOCIATIVE_KERNEL_REGRESSION_H


struct Signals_node
{
    double* data;
    int index;
};

struct Signals_matrix
{
    int nvars;
    int nvectors;
    struct Signals_node* vectors;
};

typedef struct Signals_matrix Signals_matrix;
typedef struct Signals_node Signals_node;


/**
 *    \brief Creates and initialize a Signals_matrix structure.
 *    \param nvars Number of variables.
 *    \param nvectors Number of variables observations or samples.
 *    \return Returns a initialized Signals_matrix structure, where nvectors value is
 *    equal to zero (no one of Signals_node has data) and a array of Signals_node with
 *    size equal to nvectors (input parameter) is already alocated.    
 */
Signals_matrix *
signals_matrix_init (int nvars, int nvectors);

void
aakr_prediction (Signals_matrix* memory, double bandwidth, double* query, double *output);


Signals_matrix*
memory_vector_selection (Signals_matrix* train_data, int num_vectors);

void
free_signals_matrix (Signals_matrix *input);

#endif /* end of include guard: AUTOASSOCIATIVE_KERNEL_REGRESSION_H */
