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
$RCSfile: class_protos.h,v $
$Revision: 1.2 $
$Author: bert $
$Date: 2005-02-11 20:15:38 $
$State: Exp $
--------------------------------------------------------------------------*/

/* Function declarations and prototyping */
typedef void (*Init_Training)(char *filename);
typedef void (*Load_Training)(char *filename);
typedef void (*Save_Training)(char *filename);
typedef void (*Train)();
typedef void (*Classify)(int *class_num, VIO_Real *class_probs, int *class_labels);

/* FUNCTION PROTOTYPES */

void parse_arguments(int argc, char *argv[]);
void load_input_volumes(void);
void load_tag_file( char *);
void load_training_samples( char *);
void classify_volume(void);
void create_empty_classified_volume(void);
void write_classified_volume(void);
void create_feature_matrix_from_tagfile();
void create_feature_matrix_from_tagvolume();
void allocate_memory(void);
void cleanup_memory(void);
void initialize_fuzzy_volumes(void);
void write_fuzzy_volumes(void);
void set_classifier_functions( int classifier_index );
void load_mask_volume(char *mask_file);
int  voxel_is_in_volume( VIO_Real vox1, VIO_Real vox2, VIO_Real vox3);
void decide_fuzzy_volumes(void);
void convert_features_to_slice_caching(void);
void load_train_volumes(char **tagvolume_filename);
int  volume_size_is_ok( VIO_Volume vol );

