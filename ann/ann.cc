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
$RCSfile: ann.cc,v $
$Revision: 1.2 $
$Author: claude $
$Date: 2011-05-27 20:47:01 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include "EBTKS/Minc.h"
#include "EBTKS/FileIO.h"
#include "EBTKS/SimpleArray.h" 
//#include "TrainingSet.h"    
#include "EBTKS/ValueMap.h"
#include "ann.h"
#include "EBTKS/backProp.h"
extern "C" {
  #include <volume_io.h>
  #include "../class_globals.h"
}

//Global variables:
double           minInput  = -1.0;
double           maxInput  = 1.0;
double           minTarget = 0.1;
double           maxTarget = 0.9;
unsigned         nNodesInHiddenLayer = 10;
Array<LinearMap> inputMaps;
IntArray         nodeToClassMap;
BP_ANN          *BP = 0;

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_init_training
@INPUT         : parameter file name
@OUTPUT        :
@RETURNS       :
@DESCRIPTION   : initializes training for an artificial neural network
@METHOD        :
@GLOBALS       :
@CALLS         : 
@CREATED       : Jul 24, 1995 (Samson ANTEL)
@MODIFIED      : Aug 06, 1995 (Alex Zijdenbos)
------------------------------------------------------------------------------ */
void ann_init_training(char *param_filename)
{
  if (BP) {
    delete BP;
    BP = 0;
  }

  InputFile annFile;
  if (!param_filename || !annFile.attach(param_filename)) {
    if (param_filename) {
      cerr << "Couldn't open parameter file " << param_filename << endl;
      exit(EXIT_FAILURE);
    }

    if (verbose)
      cout << "No parameters specified - using defaults" << endl;

    UnsignedArray topology(3);
    topology[(unsigned int)0] = num_features;
    topology[(unsigned int)1] = nNodesInHiddenLayer;
    topology[(unsigned int)2] = num_classes;
    BP = new BP_ANN(topology, verbose);

    if (!BP) {
      cerr << "Couldn't create BP_ANN" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {
    if (verbose)
      cout << "Reading network parameters..." << endl;

    BP = new BP_ANN(annFile.stream(), verbose);

    annFile.attach(param_filename); // Re-attach to rewind input stream
    if (!BP) {
      cerr << "Couldn't create BP_ANN" << endl;
      exit(EXIT_FAILURE);
    }

    ann_load_train_param(annFile.stream());
  }

  if (num_features != BP->nInputNodes()) {
    cerr << "Warning: num_features (" << num_features 
	 << ") doesn't match # ann input nodes ("
	 << BP->nInputNodes() << "). Forcing " << num_features << endl;
    BP->nInputNodes(num_features);
  }
  if (num_classes != BP->nOutputNodes()) {
    cerr << "Warning: num_classes (" << num_classes 
	 << ") doesn't match # ann output nodes ("
	 << BP->nOutputNodes() << "). Forcing " << num_classes << endl;
    BP->nOutputNodes(num_classes);
  }

  /*
  if (num_features != BP->nInputNodes()) {
    cerr << "num_features (" << num_features << ") doesn't match # ann input nodes ("
	 << BP->nInputNodes() << ")" << endl;
    exit(EXIT_FAILURE);
  }
  if (num_classes != BP->nOutputNodes()) {
    cerr << "num_classes (" << num_classes << ") doesn't match # ann output nodes ("
	 << BP->nOutputNodes() << ")" << endl;
    exit(EXIT_FAILURE);
  }
  */
}

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_train_samples
@INPUT         :
@OUTPUT        :
@RETURNS       :
@DESCRIPTION   : takes a feature matrix and trains an ANN classifier on it.
@METHOD        :
@GLOBALS       :
@CALLS         : 
@CREATED       : Jul 24, 1995 (Samson ANTEL)
@MODIFIED      : Aug 06, 1995 (Alex Zijdenbos)
------------------------------------------------------------------------------ */
void ann_train_samples(void)
{
  if (!BP) {
    fprintf(stderr, "BP_ANN not initialized!\n");
    exit(EXIT_FAILURE);
  }

  unsigned i, j;

  inputMaps.newSize(num_features);

  for (j=0; j<num_features; j++){
    double mmax = -MAXDOUBLE;
    double mmin = MAXDOUBLE;
    for (unsigned i=0; i<num_samples; i++){
      double value = feature_matrix[i][j];
      if (value > mmax)
	mmax=value;
      if (value < mmin)
	mmin=value;
    }

    inputMaps[j](mmin, mmax, minInput, maxInput);
  }
       
  IntArray temp(class_column, num_samples);
  int minClass = min(temp);
  int maxClass = max(temp);

  UnsignedArray classToNodeMap(maxClass - minClass + 1);
  classToNodeMap.clear(0);

  nodeToClassMap.newSize(num_classes);
  nodeToClassMap.clear(minClass);

  unsigned nodeCtr = 1;
  for (i = 0; i < num_samples; i++) {
    unsigned classID = class_column[i] - minClass;

    if ((classID != 0) && !classToNodeMap[classID]) {
      classToNodeMap[classID] = nodeCtr;
      nodeToClassMap[nodeCtr] = classID + minClass;
      nodeCtr++;
    }
  }

  if (verbose)
    ann_save_train_param(cout);

  // Remap feature_matrix
  if (verbose)
    cout << "Remapping feature matrix..." << flush;
  for (i = 0; i < num_samples; i++)
    for (j = 0; j < num_features; j++)
      feature_matrix[i][j] = inputMaps[j](feature_matrix[i][j]);
  if (verbose)
    cout << "Done" << endl;

  BP->initTraining(num_samples);
  BP->shuffleInterval(25);
  
  DblArray  targetValues(minTarget, num_classes);
  double   *targetValuePtr = targetValues.contents(); // For efficiency only

  int sample = 0;
  while (sample >= 0 && sample < num_samples) {
    // Create target values array
    double *node = targetValuePtr + classToNodeMap[class_column[sample] - minClass];

    *node = maxTarget;
    sample = BP->train(feature_matrix[sample], targetValuePtr, ann_error_monitor);
    *node = minTarget;
  }

  /*
  if (verbose)
    cout << "Creating training set..." << flush;
  
  TrainingSet TS(num_samples, num_features, num_classes, minTarget, maxTarget);

  for (i = 0; i < num_samples; i++){
    DblArray inputValues(num_features);
    for (unsigned j = 0; j < num_features; j++)
      inputValues[j] = inputMaps[j](feature_matrix[i][j]);
    TS.add(classToNodeMap[class_column[i] - minClass], inputValues);
  }

  if (verbose)
    cout << "Done" << endl << "Training network..." << endl;

  BP->train(TS, ann_error_monitor);
  */
}

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_save_training;
@INPUT         :
@OUTPUT        :
@RETURNS       :
@DESCRIPTION   : saves the results of an ANN training session.
@METHOD        :
@GLOBALS       :
@CALLS         : ann_save_train_param;
@CREATED       : Jul 24, 1995 (Samson ANTEL)
@MODIFIED      : Aug 05, 1995 (Alex ZIJDENBOS)
------------------------------------------------------------------------------ */
void ann_save_training(char *save_train_filename)
{
  if (!BP) {
    cerr << "BP_ANN not trained!" << endl;
    exit(EXIT_FAILURE);
  }

  if (verbose)
    cout << "Saving trained network..." << flush;

  OutputFile save_train_file(save_train_filename, ios::out, OutputFile::NO_COMPRESS);
 
  if (!save_train_file) {
    cerr << "Cannot open " << save_train_filename << endl;
    exit(EXIT_FAILURE);
  }
 
  BP->save(save_train_file);
  ann_save_train_param(save_train_file);
 
  if (verbose)
    cout << "Done" << endl;
}

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_load_training
@INPUT         :
@OUTPUT        :
@RETURNS       : 
@DESCRIPTION   : loads the results of an ANN training session
@METHOD        :
@GLOBALS       :
@CALLS         : ann_load_train_param(FILE *)
@CREATED       : Jul 24, 1995 (Samson ANTEL)
@MODIFIED      : Aug 06, 1995 (Alex ZIJDENBOS)
------------------------------------------------------------------------------ */
void ann_load_training(char *load_train_filename)
{
  InputFile load_train_file(load_train_filename);
  if (!load_train_file) {
    cerr << "Cannot open " << load_train_filename << endl;
    exit(EXIT_FAILURE);
  }

  if (BP)
    delete BP;
  
  if (verbose)
    cout << "Reading network..." << endl;

  BP = new BP_ANN(load_train_file.stream(), verbose);
  if (!BP) {
    cerr << "Cannot create BP_ANN" << endl;
    exit(EXIT_FAILURE);
  }

  num_features = BP->nInputNodes();
  num_classes = BP->nOutputNodes();

  load_train_file.attach(load_train_filename); // Re-attach to rewind input stream
  ann_load_train_param(load_train_file.stream());

  if (verbose)
    ann_save_train_param(cout);
}

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_classify_sample
@INPUT         :
@OUTPUT        : class_number, class_prob, class_labels
@RETURNS       : 
@DESCRIPTION   : given a feature vector, return a class, class probability values,
                 and the corresponding class numbers.
@METHOD        :
@GLOBALS       :
@CALLS         :
@CREATED       : Jul 24, 1995 (Samson ANTEL)
@MODIFIED      : Aug 07, 1995 (Alex ZIJDENBOS)
------------------------------------------------------------------------------ */

void ann_classify_sample(int *class_num, double *class_prob, int *class_labels)
{  
  if (!BP) {
    fprintf(stderr, "BP_ANN doesn't exist!\n");
    exit(EXIT_FAILURE);
  }

  DblArray inputValues(feature_vector, num_features);
  DblArray outputValues(num_classes);

  BP->evaluate(map(inputValues, inputMaps), outputValues);

  unsigned node;
  outputValues.max(&node);
  *class_num = nodeToClassMap[node];

  if (class_prob) {
    outputValues.asCarray(class_prob);
    if (class_labels)
      nodeToClassMap.asCarray(class_labels);
  }

//    cout << "Classify: " << inputValues << " => " << *class_num << endl;
} 

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_save_train_param
@INPUT         :
@OUTPUT        :
@RETURNS       : 
@DESCRIPTION   : appends nCycles, mmax, mmin to saved network file
@METHOD        :
@GLOBALS       :
@CALLS         :
@CREATED       : Aug 01, 1995 (Samson ANTEL)
@MODIFIED      : Aug 06, 1995 (Alex ZIJDENBOS)
------------------------------------------------------------------------------ */

ostream& ann_save_train_param(ostream& OS)
{
  OS << "min_input:         " << minInput << endl
     << "max_input:         " << maxInput << endl
     << "min_target:        " << minTarget << endl
     << "max_target:        " << maxTarget << endl
     << "input_maps:        " << num_features << endl;

  for (unsigned i = 0; i < num_features; i++)
    OS << inputMaps[i].factor() << " " << inputMaps[i].offset() << endl;

  //  OS << "node_to_class_map: " << nodeToClassMap << endl;
  
  return OS;
}

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_load_train_param
@INPUT         :
@OUTPUT        :
@RETURNS       : 
@DESCRIPTION   : loads nCycles, mmax, mmin from saved network file
@METHOD        :
@GLOBALS       :
@CALLS         :
@CREATED       : Aug 01, 1995 (Samson ANTEL)
@MODIFIED      : Aug 06, 1995 (Alex ZIJDENBOS)
------------------------------------------------------------------------------ */

istream& ann_load_train_param(istream& IS)
{
// Set defaults
  minInput  = -1.0;
  maxInput  = 1.0;
  minTarget = 0.1;
  maxTarget = 0.9;
  nodeToClassMap = IntArray(0, 1, num_classes - 1);

// Load values from stream
  MString key;
  while (IS >> key) {
    if (key.contains("min_input"))
      IS >> minInput;
    else if (key.contains("max_input"))
      IS >> maxInput;
    else if (key.contains("min_target"))
      IS >> minTarget;
    else if (key.contains("max_target"))
      IS >> maxTarget;
    else if (key.contains("input_maps")) {
      unsigned nMaps;
      IS >> nMaps;
      inputMaps.newSize(nMaps);
      for (unsigned i = 0; i < nMaps; i++)
	IS >> inputMaps[i].factor() >> inputMaps[i].offset();
    }
    else if (key.contains("node_to_class_map")) {
      for (unsigned i = 0; i < num_classes; i++)
	IS >> nodeToClassMap[i];

      if (class_name == NULL) {
	typedef char *charPtr;
	class_name = new charPtr[num_classes];
	assert(class_name);
	for (unsigned i = 0; i < num_classes; i++) {
	  MString intString;
	  intString += nodeToClassMap[i];
	  class_name[i] = new char[intString.length() + 1];
	  strcpy(class_name[i], intString);
	}
      }
    }
  }

  return IS;
}

/* ----------------------------- MNI Header ----------------------------------
@NAME          : ann_error_monitor
@INPUT         : cycle #, cycle MSE
@OUTPUT        :
@RETURNS       : 
@DESCRIPTION   : ANN output error monitor
@METHOD        :
@GLOBALS       :
@CALLS         :
@CREATED       : Aug 06, 1995 (Alex ZIJDENBOS)
@MODIFIED      : 
------------------------------------------------------------------------------ */

void
ann_error_monitor(unsigned cycle, double error)
{
  if (verbose)
    cout << "Cycle: " << cycle << " error: " << error << endl;
}
