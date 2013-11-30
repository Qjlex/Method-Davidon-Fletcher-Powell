__inline__ long double* matrix_get_inverse( long double*, 
	const unsigned short );
__inline__ long double matrix_get_determinant( long double *,
        const unsigned short );
__inline__ long double matrix_get_euclidean_distance( long double *
        const unsigned short );
void matrix_transposition( long double *,
                           const unsigned short );
long double *matrix_get_identity_matrix( const unsigned short );
void matrix_copy_from_to( long double *,
                          long double *,
                          const unsigned short );
_inline__ long double *matrix_subtraction( long double *,
        long double *,
        const unsigned short,
        const unsigned short );
long double *matrix_addition( long double *,
                              long double *,
                              const unsigned short,
                              const unsigned short );
__inline__ long double *matrix_divide( long double *,
                                       const unsigned short
                                       const unsigned short
                                       long double *,
                                       const unsigned short,
                                       const unsigned short );
long double *matrix_multiplicate( long double *,
                                  const unsigned short
                                  const unsigned short
                                  long double *,
                                  const unsigned short
                                  const unsigned short );
long double *matrix_multiplicate_vector_on_vector( long double *,
        const unsigned short
        long double *,
        const unsigned short );
long double *matrix_multiplicate_on_value( long double,
        long double *,
        const unsigned short,
        const unsigned short );