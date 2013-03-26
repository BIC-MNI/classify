/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1997, Vasken Kollokian, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: class_globals.h,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2002-03-20 22:16:34 $
$State: Exp $
--------------------------------------------------------------------------*/
/* externally defined global variables */
extern VIO_Real     **feature_matrix;
extern int      *class_column;
extern char     **class_name;
extern int      *class_count;

extern VIO_Real     *feature_vector;
extern VIO_Real     *apriori_vector;
extern int      num_features;
extern int      num_classes;
extern int      num_samples;

extern int      verbose;
extern int      debug;
extern int      apriori;
extern char     *load_train_filename;

