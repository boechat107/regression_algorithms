/**
 *    @file       aakr.h
 *    @brief      Functions to make predictions using the AAKR model.
 *    @author     Andre A. Boechat <boechat@das.ufsc.br>
 *    @date       February 03, 2012
 *
 *    Reference: 
 *    "Technical Review of On-Line Monitoring Techniques for Performance Assessment Vol.
 *    2: Theoretical Issues"
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


double*
aakr_prediction (Signals_matrix* memory, double bandwidth, double* query);


Signals_matrix*
memory_vector_selection (Signals_matrix* train_data, int num_vectors);

void
free_signals_matrix (Signals_matrix *input);

#endif /* end of include guard: AUTOASSOCIATIVE_KERNEL_REGRESSION_H */
