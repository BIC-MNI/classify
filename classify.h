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
$RCSfile: classify.h,v $
$Revision: 1.4 $
$Author: bert $
$Date: 2005-02-11 20:15:49 $
$State: Exp $
--------------------------------------------------------------------------*/
#define MD 0
#define KNN 1
#define ANN 2
#define C45 3
#define HCM 4
#define FCM 5
#define BAY 6

// JPL: values.h does not exist on OS X, so here I'm adding a hack and
// defining the necessary values myself
//#include <values.h>
#include <limits.h>
#include <float.h>
#define MAXSHORT   SHRT_MAX
#define MAXINT     INT_MAX 
#define MAXDOUBLE  DBL_MAX

  // JPL: end values.h hack

/* GLOBAL VARIABLES */  
Init_Training init_training = 0;
Load_Training load_training = 0;
Save_Training save_training = 0;
Train         train = 0;
Classify      classify = 0;

VIO_Status     status;            /* status of loading and saving */
int        num_features;      /* number of volumes to be processed */  
char       *pname;            /* the name of the invoked program */
char       *history;          /* command line added to volume's history */
int        verbose = FALSE;   /* inform user as things are processed */
int        clobber = FALSE;   /* overwrite existing volume file */
int        clobber_training  = FALSE;   /* overwrite existing training file */
int        clobber_fuzzy  = FALSE;   /* overwrite existing fuzzy volume files */
int        clobber_all = FALSE;   /* overwrite any existing volume files */
int        debug = 0;         /* print debugging info, relative to level */
int        dump_features = 0; /* print feature matrix, & class column, and exit */
int        classifier = 0;    /* default classifier is 0 - minimum distance */
int        train_only = FALSE;/* train only the classifier, don't classify */
int        test_classifier = FALSE; /* option to test classifier */
char       *fuzzy = NULL;           /* set fuzzy string for each class */
int        *create_fuzzy_volume;    /* int array to indicate fuzzy volume creation */
int        fuzzy_available = FALSE;     /* fuzzy classification available */
int        use_apriori;                 /* apriori use flag */
char       *fuzzy_prefix = NULL;        /* fuzzy filename prefix */ 
int        cache_set = TRUE;           /* set the caching option to true  */
int        apriori = FALSE;            /* set the apriori option to false */
int        max_cache_size = 0;          /* set the max cache size to 0 */

nc_type    type = NC_BYTE;    /* specify what type of data to read from volumes */
nc_type    fuzzy_type = NC_BYTE;    /* specify type of fuzzy volumes */
int        sign = FALSE;      /* specify whether it is signed or not */


VIO_Volume     *in_volume;                /* pointer to array of volumes */
VIO_Volume     *fuzzy_volume;             /* pointer to array of fuzzy volumes */
char       **fuzzy_filename;          /* pointer to array of fuzzy volume filenames */
char       *fuzzy_path = NULL;        /* pointer to array of fuzzy volume path */

int        *first_volume_sizes;       /* 1D array to hold sizes for 1st vol */
int        first_volume_num_dims;     /* to hold num of dimensions */
char       **first_volume_dim_names;  /* 1D array of char* to hold the dim names */

VIO_Volume     classified_volume;         /* classified volume pointer */

char       **input_filename;          /* 1D char array to hold volume filenames*/
char       *output_filename;          /* classified volume filename */
char       *tagfile_filename = NULL;  /* tagfile file filename */
char       **trainvol_filename = NULL; /* tagfile volumes filename pointer */
char       *tagvolume_buffer;         /* buffer to hold all tag volumes in */ 
char       *classname_buffer;         /* buffer to hold class names in tag volume*/ 
char       *load_train_filename;      /* custom input training point filename */
char       *save_train_filename;      /* custom output training point filename */
char       *param_filename;           /* classifier parameter filename */

int        num_samples;               /* number of samples in a tag file */

int        n_tag_volumes;             /* number of volumes in a tag file */
char       **labels;                  /* array of labels indicating classes */
char       **class_name = NULL;       /* array indicating class names */
int        *class_count;              /* array indicating each class count */
int        max_class_index = 0;       /* highest class number to set voxelmax */
VIO_Real       **tags;                    /* matrix to hold world coordinates */

int        v1_ptr, v2_ptr, v3_ptr;    /* pointer to voxels 1 2 and 3 */

VIO_Real       fuzzy_image_min = 0.0;    /* fuzzy volume image min & max */
VIO_Real       fuzzy_image_max = 1.0;
VIO_Real       fuzzy_voxel_min = 0.0;    /* fuzzy volume voxel min & max */
VIO_Real       fuzzy_voxel_max = 255.0;
  
VIO_Real       **feature_matrix;    /* matrix to hold feature vector of all samples */

VIO_Real       *feature_vector;     /* vector to hold feature vector of one voxel */
VIO_Real       *apriori_vector;     /* vector to hold apriori class probabilities */
int        *class_column;       /* array of integers indicating class labels */

int        num_classes = 0;     /* the number of classes to be processes */

VIO_Real       *class_probs;        /* array to indicate fuzzy class confidence */
int        *class_labels;       /* array to indicate fuzzy class labels */

char       *mask_filename;        /* filename of the mask volume */
VIO_Volume     mask_volume;           /* volume variable to hold the mask */
VIO_Volume     *train_volume;         /* 1d array variable to hold the tag volumes */

int        user_mask_class = 0;   /* default class value for masked voxels */
VIO_Real       user_mask_value = 1.0; /* default mask value for masked voxels */
int        block_sizes[3] = {1,1,1};        /* default block size, (1 voxel) */
int        user_block_sizes[3] = {1,-1,-1}; /* user_selected block size (1slice) */
VIO_Real       output_range[2] = { -MAXDOUBLE, -MAXDOUBLE }; /* range of output values */
int        supervised;                      /* denote type of classifier */
int        num_train_vols = 0;              /* number of training volumes specified */

/* ARGTABLE */

#ifdef __DATE__
#ifdef __TIME__
#define VERSIONSTR VERSION " built " __DATE__ " " __TIME__
#else
#define VERSIONSTR VERSION " built " __DATE__
#endif /* __TIME__ not defined */
#else
#define VERSIONSTR VERSION
#endif /* __DATE__ not defined */

ArgvInfo argTable[] = {
  { NULL, ARGV_VERINFO, VERSIONSTR, NULL, "" },

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions to specify verbosity and clobbering.\n"},

  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Show progress."},
 
  {"-debug", ARGV_INT, (char *) NULL, (char *) &debug,
     "Show debugging information depending on level specified."},

  {"-dump_features", ARGV_CONSTANT, (char *) TRUE, (char *) &dump_features,
     "Dump the feature matrix and the class column of training set and exit."},

  {"-clob_output", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
     "Overwrite classified output file."},

  {"-clob_training", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber_training,
     "Overwrite training output file."},

  {"-clob_fuzzy", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber_fuzzy,
     "Overwrite fuzzy output volume files."},

  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber_all,
     "Overwrite any existing files (equal to '-clob_out -clob_training -clob_fuzzy')."},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions to specify input data type.\n"},

  {"-byte", ARGV_CONSTANT, (char *) NC_BYTE, (char *) &type,
     "Read input volumes as Byte values."},
  
  {"-short", ARGV_CONSTANT, (char *) NC_SHORT, (char *) &type,
     "Read input volumes as Short integer values."},
  
  {"-long", ARGV_CONSTANT, (char *) NC_LONG, (char *) &type,
     "Read input volumes as Long integer values."},

  {"-float", ARGV_CONSTANT, (char *) NC_FLOAT, (char *) &type,
     "Read input volumes as Single-precision floating point values."},

  {"-double", ARGV_CONSTANT, (char *) NC_DOUBLE, (char *) &type,
     "Read input volumes as Double-precision floating point values."},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions for sign of input data.\n"  },
  
  {"-signed", ARGV_CONSTANT, (char *) TRUE, (char *) &sign,
     "Read input volumes as Signed values."},
  
  {"-unsigned", ARGV_CONSTANT, (char *) FALSE, (char *) &sign,
     "Read input volumes as Unsigned values."},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions to specify output data type.\n"},

  {"-output_range", ARGV_FLOAT, (char *) 2, (char *) &output_range,
     "Set the valid range of the output volume. The default is to use [0, num_classes]."},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions to specify training.\n"  },
  
  {"-tagfile", ARGV_STRING, (char *) NULL, (char *) &tagfile_filename,
     "Input training points as tag file (.tag)"},      

  {"-volume", ARGV_STRING, (char *) NULL, (char *) &tagvolume_buffer,
     "Input tag or apriori volumes for corresponding classifiers, comma-separated,"},      
  {"-class_names", ARGV_STRING, (char *) NULL, (char *) &classname_buffer,
     "Input class names for the specified tag volumes, comma-separated"},      

  {"-load_train", ARGV_STRING, (char *) NULL, (char *) &load_train_filename,
     "Input training points as a custom training file."},      

  {"-save_train", ARGV_STRING, (char *) NULL, (char *) &save_train_filename,
     "Write training points as a custom training file."},      

  {"-train_only", ARGV_CONSTANT, (char *) TRUE, (char *) &train_only,
     "Only train the classifier - do not classify the volume."},

  {"-parameter", ARGV_STRING, (char *) NULL, (char *) &param_filename,
     "Specify a parameter file to load a particular classifier parameters."},      

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions to create fuzzy classification volumes.\n"},

  {"-apriori", ARGV_CONSTANT, (char *) TRUE, (char *) &apriori, 
     "Specify that apriori probability volumes be used."},

  {"-fuzzy", ARGV_STRING, (char *) NULL, (char *) &fuzzy,
     "Specify a string ex. '011...' to indicate classes for fuzzy classification.\n\t\tNote: The sequence of switches reflects class order, not class labels.\n\t\tSpecifying 'all' creates fuzzy volumes for all classes."},

  {"-fprefix", ARGV_STRING, (char *) NULL, (char *) &fuzzy_prefix,
     "Specify a prefix for fuzzy volume filenames. ( default: fuzzy_volume"},      

  {"-fpath", ARGV_STRING, (char *) NULL, (char *) &fuzzy_path,
     "Specify a path to write fuzzy volumes in, (include ending /)."},      

  {"-fuzzy_voxel_min", ARGV_FLOAT, (char *) NULL, (char *) &fuzzy_voxel_min,
     "Set the minimum voxel value of the fuzzy classified volume."},

  {"-fuzzy_voxel_max", ARGV_FLOAT, (char *) NULL, (char *) &fuzzy_voxel_max,
     "Set the maximum voxel value of the fuzzy classified volume."},

  {"-fuzzy_image_min", ARGV_FLOAT, (char *) NULL, (char *) &fuzzy_image_min,
     "Set the minimum image value of the fuzzy classified volume"},

  {"-fuzzy_image_max", ARGV_FLOAT, (char *) NULL, (char *) &fuzzy_image_max,
     "Set the maximum image value of the fuzzy classified volume."},
  
  {"-fbyte", ARGV_CONSTANT, (char *) NC_BYTE, (char *) &fuzzy_type,
      "Create fuzzy volumes with Byte values."},
  
  {"-fshort", ARGV_CONSTANT, (char *) NC_SHORT, (char *) &fuzzy_type,
     "Create fuzzy volumes with Short Integer values."},
  
  {"-flong", ARGV_CONSTANT, (char *) NC_LONG, (char *) &fuzzy_type,
     "Create fuzzy volumes with Long Integer values."},

  {"-ffloat", ARGV_CONSTANT, (char *) NC_FLOAT, (char *) &fuzzy_type,
     "Create fuzzy volumes with Single-precision floating point values."},

  {"-fdouble", ARGV_CONSTANT, (char *) NC_DOUBLE, (char *) &fuzzy_type,
     "Create fuzzy volumes with Double-precision floating point values."},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "\nSpecify masking options.\n"  },
  
  {"-mask", ARGV_STRING, (char *) NULL, (char *) &mask_filename, 
     "Specify volume which will be used as a mask."},
 
  {"-user_mask_value", ARGV_FLOAT, (char *) NULL, (char *) &user_mask_value,
     "Specify the mask value. (If the mask is a classified volume)"},
 
  {"-user_mask_class", ARGV_INT, (char *) NULL, (char *) &user_mask_class,
     "Specify the class label of masked voxels (used in the classified volume)."},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nSpecify volume caching options.\n"  },
  
  {"-nocache", ARGV_CONSTANT, (char *) FALSE, (char *) &cache_set, 
     "Specify that NO volume caching is used."},
 
  {"-max_cache_size", ARGV_INT, (char *) NULL, (char *) &max_cache_size,
     "Set the default maximum bytes in the cache."},

  {"-block_sizes", ARGV_INT, (char *) 3, (char *) &user_block_sizes,
     "Set the size of the cache block for the feature volumes."},

    
  /************************  S T A R T  *************************************/

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions to specify a classifier.\n"  },

  {"-min", ARGV_CONSTANT, (char *) MD, (char *) &classifier,
     "Use the 'Minimum Distance' classifier."},

  {"-knn", ARGV_CONSTANT, (char *) KNN, (char *) &classifier,
     "Use the 'K Nearest Neighbour' classifier."},

  {"-ann", ARGV_CONSTANT, (char *) ANN, (char *) &classifier,
     "Use the 'Artificial Neural Network' classifier."},

  {"-hcm", ARGV_CONSTANT,  (char *) HCM, (char *) &classifier,
     "Use the unsupervised Hard C Means (K-means) classifier."},

  {"-fcm", ARGV_CONSTANT,  (char *) FCM, (char *) &classifier,
     "Use the unsupervised Fuzzy C Means classifier."},

  {"-bayes", ARGV_CONSTANT,  (char *) BAY, (char *) &classifier,
     "Use the 'Bayesian'  classifier."},

  /**************************  E N D  ***************************************/
  {NULL, ARGV_HELP, NULL, NULL,
     "\n"  },

  {NULL, ARGV_END, NULL, NULL, NULL}
};
  

