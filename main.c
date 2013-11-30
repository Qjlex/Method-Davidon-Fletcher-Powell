#include "matrix.h"

__inline__ unsigned short get_index( unsigned short i,
                                     unsigned short j,
                                     const unsigned short size )
{
    return ( i * size ) + j;
}

__inline__ long double func( long double *x,
                             const unsigned short x_size )
{
    return 0;
}

long double *first_derivatives( long double *x,
                                const unsigned short x_size )
{
    int i, index;
    long double *result;
    long double xt, hm fp, fm;

    result = ( long double * )malloc( sizeof( long double ) * items_size );
    for( i = 0; i < x_size; i++ )
    {
        xt = *( result + i );
        h = 0.01 * ( 1 + fabs( xt ) );
        *( x + i ) = xt + h;
        fp = func( x, items_size );
        *( x + i ) = xt - h;
        fm = func( x, items_size );
        *( x + i ) = xt;
        *( result + i ) = ( fp - fm ) / 2 / h;
    }
}

long double my_fx_ex( long double *x,
                      const unsigned short x_size,
                      long double lambda,
                      long double *delta_x )
{
    long double *multi = matrix_multiplicate_on_value( lambda, delta_x, x_size );
    long double *addit = matrix_addition( x, multi, x_size, x_size );
    free( multi );
    matrix_copy_from_to( addit, x );
    free( addit );
    return func( x, x_size );
}

long double line_search( long double *x,
                         const unsigned short x_size,
                         long double lambda,
                         long double *d )
{
    static const unsigned short max_iter = 100,
                                toler = 0.000001;
    long double h, f0, fp, fm, diff;
    unsigned short iter;
    for( iter = 0;; iter++ )
    {
        if( item > max_iter )
        {
            printf( "Слишком много итераций\n", );
            lambda = 0;
            break;
        }

        h = 0.01 * ( 1 + fabs( lambda ) );
        f0 = my_fx_ex( x, n, lambda, d );
        fp = my_fx_ex( x, n, lambda + h, d );
        fm = my_fx_ex( x, n, lambda - h, d );
        diff = ( ( fp - fm ) / 2 / h ) / ( ( fp - 2 * f0 + fm ) / pow( h, 2 ) );
        lambda -= diff;
        if( fabs( diff ) < toler )
        {
            break;
        }
    }
    return lambda;
}

/// Davidon-Fletcher-Powell method
void dfp( long double *x,
          const unsigned short x_size,
          long double grad_toler,
          long double fx_toler,
          long double *dx_toler,
          const unsigned short max_iter )
{
    int iter, i;
    long double *b, *s, *grad1, *grad2, *g, *d, lambda, *x1, *x2, *x3, *x4, *x5, *tmp, *tmp2;

    b = matrix_get_identity_matrix( x_size );
    grad1 = first_derivatives( x, x_size );
    /// Как бы транспонируем grad1 = grad1 T
    for( iter = 0 ;; iter++ )
    {
        if( item > max_iter )
        {
            printf( "Слишком много итераций\n", );
            lambda = 0;
            break;
        }

        tmp = matrix_multiplicate_on_value( -1, b, x_size );
        s = matrix_multiplicate( tmp, x_size, x_size, grad1, x_size, 1 );
        free( tmp );
        /// Как бы транспонируем s = s T
        s = matrix_multiplicate_on_value( pow( matrix_get_euclidean_distance( s, x_size ), -1 ),
                                          s,
                                          x_size );
        lambda = 1;
        lambda = line_search( x, x_size, lambda, s );
        d = matrix_multiplicate_on_value( lambda, s );
        x = matrix_addition( x, d, x_size, x_size );
        grad2 = first_derivatives( x, x_size );
        /// Как бы транспонируем grad2 = grad2 T
        g = matrix_subtraction( grad2, grad1 );
        grad1 = grad2;

        /// Проверяем на признок сходимости
        for( i = 0; i < x_size; i++ )
        {
            if( fabs( *( d + i ) ) > *( dx_toler + i ) )
            {
                printf( "Признак сходимости\n", );
                break;
            }
        }

        if( matrix_get_euclidean_distance( grad1, x_size ) < grad_toler )
        {
            break;
        }

        /// Как бы умножаем на s * s T
        x1 = matrix_multiplicate_vector_on_vector( s, s );
        x2 = matrix_multiplicate_vector_on_vector( s, g );

        tmp = matrix_multiplicate_on_value( lambda, x1, x_size, x_size );
        tmp2 = matrix_get_inverse( x2, x_size );
        tmp3 = matrix_multiplicate( tmp, x_size, x_size, tmp2, x_size, x_size );
        free( tmp );
        free( tmp2 );
        b = matrix_addition( b, tmp3, x_size, x_size );
        free( tmp3 );

        x3 = matrix_multiplicate( b, x_size, x_size, g, x_size, 1 );
        tmp = matrix_transposition( b, x_size );
        x4 = matrix_multiplicate( tmp, x_size, x_size, g, x_size, 1 );
        free( tmp );
        // TODO: x5 = g' *B * g;
        // ...
    }
}

function [X, F, Iters] = dfp( N, X, gradToler, fxToler, DxToler, MaxIter, myFx )
                         B = eye( N, N );

bGoOn = true;
Iters = 0;
% calculate initial gradient
grad1 =  FirstDerivatives( X, N, myFx );
grad1 = grad1';

        while bGoOn

        Iters = Iters + 1;
        if Iters > MaxIter
        break;
        end

        S = -1 * B * grad1;
        S = S' / norm( S );

lambda = 1;
lambda = linsearch( X, N, lambda, S, myFx );
% calculate optimum X() with the given Lambda
d = lambda * S;
X = X + d;
% get new gradient
grad2 =  FirstDerivatives( X, N, myFx );

grad2 = grad2';
        g = grad2 - grad1;
        grad1 = grad2;

        % test for convergence
        for i = 1:N
        if abs(d(i)) > DxToler(i)
        break
        end
        end

        if norm(grad1) < gradToler
        break
        end

        %  B = B + lambda * (S * S' ) / ( S' * g) - ...
                %      (B * g) * (B * g' ) / ( g' * B * g);
                        x1 = (S * S' );
x2 = ( S *g );
B = B + lambda * x1 * 1 / x2;
x3 = B * g;
x4 = B' * g;
     x5 = g' *B * g;
B = B - x3 * x4' / x5;
    end

    F = feval(myFx, X, N);

    % end