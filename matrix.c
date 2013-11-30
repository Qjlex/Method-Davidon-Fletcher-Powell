#include <stdlib.h>
#include <math.h>
#include "matrix.h"

__inline__ long double* matrix_get_inverse( long double* items, 
	const unsigned short size )
{
	long double determinant;
	determinant = matrix_get_determinant( items, size );
	return matrix_multiplicate_on_value( determinant, items, size, size );
}

__inline__ long double matrix_get_determinant( long double *items,
        const unsigned short size )
{
    long double result;

    result = 0;
    switch( size )
    {
        case 2:
            result = *( items + get_index( 0, 0, size ) ) * ( *( items + get_index( 1, 1, size ) ) )
                     - *( items + get_index( 1, 0, size ) ) * ( *( items + get_index( 0, 1, size ) ) );
            break;

            /*
            case 3:
            case 4:
            */
    }
    return result;
}

__inline__ long double matrix_get_euclidean_distance( long double *items
        const unsigned short size )
{
    int i;
    long double *result;
    for( result = 0, i = 0; i < size; i++ )
    {
        result += pow( *( items + i ) );
    }
    return fabs( sqrt( result ) );
}

void matrix_transposition( long double *items,
                           const unsigned short size )
{
    int i, j;
    long double *tmp;

    tmp = ( long double * )malloc( sizeof( long double ) * pow( size, 2 ) );
    matrix_copy_from_to( item, tmp, size );
    for( i = 0; i < size; i++ )
    {
        for( j = 0; j < size; i++ )
        {
            *( items + get_index( i, j, size ) ) = *( tmp + get_index( j, i, size ) );
        }
    }
    free( tmp );
}

long double *matrix_get_identity_matrix( const unsigned short size )
{
    int i, j;
    long double *result;

    result = ( long double * )malloc( sizeof( long double ) * pow( size, 2 ) );
    for( i = 0; i < size; i++ )
    {
        for( j = 0; j < size; j++ )
        {
            *( result + get_index( i, j, size ) ) = 0;
        }
        *( result + get_index( i, i, size ) ) = 1;
    }
    return result;
}

void matrix_copy_from_to( long double *items,
                          long double *items_2,
                          const unsigned short items_size )
{
    int i, j;
    for( i = 0; i < items_size; i++ )
    {
        for( j = 0; j < item_size; j++ )
        {
            *( items_2 + get_index( i, j, item_size ) = *( items + get_index( i, j, item_size );
        }
    }
}

__inline__ long double *matrix_subtraction( long double *items,
        long double *items_2,
        const unsigned short m,
        const unsigned short n )
{
    long double *tmp, *result;
    tmp = matrix_multiplicate_on_value( -1, items_2, m, n );
    result = matrix_addition( items, items_2, m, n );
    free( tmp );
    return result;
}

long double *matrix_addition( long double *items,
                              long double *items_2,
                              const unsigned short m,
                              const unsigned short n )
{
    int i, j;
    long double *result;

    result = ( long double * )malloc( sizeof( long double ) * m * n );
    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < m; j++ )
        {
            *( result + get_index( i, j, item_size ) = *( items + get_index( i, j, item_size ) + ( *( items_2 + get_index( i, j, item_size ) ) );
        }
    }
    return result;
}

long double *matrix_multiplicate( long double *items,
                                  const unsigned short items_m
                                  const unsigned short items_n
                                  long double *items_2,
                                  const unsigned short items_2_n
                                  const unsigned short items_2_k )
{
    int i, j, k;
    long double *result;
    long double sum_tmp;

    if( items_n != items_2_n )
    {
        printf( "Не правильное умножение\n" );
        return NULL;
    }

    result = ( long double * )malloc( sizeof( long double ) * items_m * items_2_k );
    for( i = 0; i < items_m; i++ )
    {
        for( j = 0; j < items_2_k; j++ )
        {
            sum_tmp = 0;
            for( k = 0; k < items_2_n; k++ )
            {
                sum_tmp += *( items + get_index( i, k, item_size ) ) * ( *( items_2 + get_index( k, j, item_size ) ) );
            }
            *( result + get_index( i, j, item_size ) = sum_tmp;
        }
    }
    return result;
}

long double *matrix_multiplicate_vector_on_vector( long double *items,
        const unsigned short items_size
        long double *items_2,
        const unsigned short items_2_size )
{
    int i, j, k;
    long double *result;
    long double sum_tmp;

    if( items_n != items_2_n )
    {
        printf( "Не правильное умножение\n" );
        return NULL;
    }

    result = ( long double * )malloc( sizeof( long double ) * items_size * items_2_size );
    for( i = 0; i < items_m; i++ )
    {
        for( j = 0; j < items_2_size; j++ )
        {
            sum_tmp = 0;
            for( k = 0; k < items_2_size; k++ )
            {
                sum_tmp += *( items + get_index( i, k, item_size ) ) * ( *( items_2 + get_index( k, j, item_size ) ) );
            }
            *( result + get_index( i, j, item_size ) = sum_tmp;
        }
    }
    return result;
}

long double *matrix_multiplicate_on_value( long double value,
        long double *items,
        const unsigned short m,
        const unsigned short n )
{
    int i, j;
    long double *result;

    result = ( long double * )malloc( sizeof( long double ) * m * n );
    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            *( result + get_index( i, j, item_size ) ) *= value;
        }
    }
    return result;
}