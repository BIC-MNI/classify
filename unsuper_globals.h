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
$RCSfile: unsuper_globals.h,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2002-03-20 22:16:34 $
$State: Exp $
--------------------------------------------------------------------------*/
/* externally defined global variables */

extern int      num_features;
extern int      num_classes;
extern int      num_samples;
extern int      *class_count;
extern char     **class_name;
extern int      *first_volume_sizes;
extern VIO_Real     *feature_vector;
extern int      verbose;
extern int      debug;
extern int      max_class_index;
extern int      num_train_vols;
extern char     **trainvol_filename;
extern char     *classname_buffer;
extern VIO_Volume   *train_volume;
extern VIO_Volume   *in_volume;
extern char     *load_train_filename;
extern int      *create_fuzzy_volume;
