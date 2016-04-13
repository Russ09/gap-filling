/*=========================================================================
 Author: Laurent Risser, Francois-Xavier Vialard
 
 Disclaimer: This software has been developed for research purposes only, and hence should 
 not be used as a diagnostic tool. In no event shall the authors or distributors
 be liable to any direct, indirect, special, incidental, or consequential 
 damages arising of the use of this software, its documentation, or any 
 derivatives thereof, even if the authors have been advised of the possibility 
 of such damage. 
 =========================================================================*/


#ifndef _LD_SCICALCPACK_H
#define _LD_SCICALCPACK_H


//Parallelization with openmp
#ifndef __APPLE__
  #define COMPILE_WITH_OPENMP
  #include <omp.h>
#endif

//related to the i/o nifti and added instead of 'irtkImage.h' in irtk
#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352


// C++ header files
#include <iostream>
using namespace std;
#include <iomanip>
#include <fstream>
#include <complex>
#include <algorithm>
#include <string>
#include <limits>
#include <ctime>


// C header files
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>



#include "nifti1.h"
#include "nifti1_io.h"


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           1:   FUNCTIONS FOR THE CLASS "ScalarField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///1.1) conventional scalar field

class ScalarField{
  /// ******************************************************************************
private:
  //size of the fields
  int NXtY;      //NX*NY
  int NXtYtZ;    //NX*NY*NZ
  
  //scalar field
  float * ScalField;
  
  /// ******************************************************************************
public:
  //size of the fields
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  int NT;        //image size on time
  
  float Image2World[4][4];  //quaternion to go from the image coordinates to the world coordinates
  float World2Image[4][4];  //quaternion to go from the world coordinates to the image coordinates
  
  /// Constructor and destructor
  ScalarField();
  ~ScalarField();
  
  /// functions associated to the field
  //put a value in the scalar field
  virtual inline void P(float value,int x,int y,int z=0,int t=0){
    this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x]=value;
  }
  
  //Add a value to the scalar field
  virtual inline void Add(float value,int x,int y, int z=0,int t=0){
    this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x]+=value;
  }
  
  // put a the same value at every points of the scalar field
  virtual void PutToAllVoxels(float cste,int t=0);
  
  //get a value from the scalar field
  virtual inline float G(int x,int y,int z=0,int t=0){
    return this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x];
  }
  
  //get a value from the scalar field by linear interpolation
  virtual float G(double x,double y,double z=0.,int t=0);
  virtual float G(float x, float y, float z=0., int t=0);
  
  //get a value using linear interpolation or nearest neigbhour
  //-> Here the coordinate not defined in the image space but in 'coordSpace' 
  //-> The 4*4 matrix 'coordSpace2imageSpace' is a quaternion that allows to transform the coordinates from 'coordSpace' to the image space
  virtual float G(float coordSpace2imageSpace[4][4],double x,double y,double z=0.,int t=0,int NN=0);
  virtual float G(float coordSpace2imageSpace[4][4],float x,float y,float z=0.,int t=0,int NN=0);
  virtual float G(float coordSpace2imageSpace[4][4],int x,int y,int z=0.,int t=0,int NN=0);
  
  virtual float G_NoExtrapo(float coordSpace2imageSpace[4][4],float x,float y,float z=0.,int t=0);
  virtual float G_NoExtrapo(float x,float y,float z=0.,int t=0);
  
  // get the maximum absolute values out of the scalar field
  virtual float GetMaxAbsVal(int t=0);
  
  //read a scalar field (in a nifti image)
  virtual void Read(char *);
  
  //read a scalar field in a ROI (in a nifti image)
  // -> advanced memory managment for large images: only allocate memory for the ROI
  // -> Inputs are different than in Read_ROI_Given_ImageToWorld:
  //     We give here the min and max {X,Y,Z,T} in the input image defining the outputed ROI
  virtual void Read_only_ROI(char *,int,int,int,int,int,int,int,int);

  //read a scalar field in a ROI (in a nifti image)
  // -> advanced memory managment for large images: only allocate memory for the ROI
  // -> Inputs are different than in Read_only_ROI:
  //     We give here the ImageToWorld coordinates and the size (in voxels) of the outputed ROI
  // -> remark that nearest ngbh intensity is considered to sample the output image
  virtual void Read_ROI_Given_ImageToWorld_and_Size(float ROI_Image2World[4][4],int NBX,int NBY,int NBZ,char * RefImageName);

  //read a scalar field in a ROI (in a nifti image)
  // -> advanced memory managment for large images: undersample the image direcly by averaging blocks of voxels of size BlockS*BlockS*BlockS
  //Use 'Read_and_Undersample' or 'Read_and_Interpolate' to have finer resamplings requiring more memory
  virtual void Read_directly_Undersampled(char * ImageName,int BlockS);
  
  //read a scalar field and perform linear interpolation to give it a specific size
  virtual void Read_and_Interpolate(char *,int,int,int);
  
  ///read a scalar field and perform undersampling by a factor 'factor'
  virtual void Read_and_Undersample(char * ImageName,float factor,float NN=0);

  ///read a scalar field and expend its domain
  virtual void ReadAndExpend(char * ,int ,int ,int ,int ,int ,int);
    
  //create a void scalar field. All values are initialized to 'cste' which is null by default. No message is printed if Verbose==1.
  virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1,float cste=0.0,int Verbose=1);
  
  //Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
  virtual void SlashFieldSize(int verbative=1);

  //write a scalar field in a nifti image
  //The optional 2nd file is an input file containing the headers 
  virtual void Write(char *);
  
  //write a scalar field in a nifti image
  //The optional 2nd file is an input file containing the headers 
  virtual void Write(char *,char *);

  //return the number of voxels in a ScalarField
  virtual inline int GetNbVoxels(void){
    return this->NX *this->NY*this->NZ*this->NT;
  }
  
  
  ///General function used to align the cumulative histogram of *this to the input cumulative histogram.
  ///-> called by GreyLevAlignment, GreyLevAlignment_UsingRefROIs, HistoEqualization, HistoEqualizationInROI
  virtual void HistoMatch(int NbBinsCumHisto,float * TrgCumHisto_x_axis,float * TrgCumHisto_y_axis,float * LocCumHisto_x_axis,float * LocCumHisto_y_axis);
  
  ///grey levels alignment of the grey levels using optimal transport 
  ///-> optimal transportation minmizes the wassertein 1 distance between the linearly aligned histograms
  virtual void GreyLevAlignment(ScalarField * RefImage);

  ///Same as "GreyLevAlignment" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI.
  ///Note that "LocImageROI" must have the same size as "this" and that "RefImageROI" must have the same size as "RefImage"
  ///ROIs are defined by non-null values
  virtual void GreyLevAlignment_UsingRefROIs(ScalarField * RefImage,ScalarField * RefImageROI,ScalarField * LocImageROI);

  ///Grey levels histogram equalization using optimal transport 
  virtual void HistoEqualization();

  ///Same as "HistoEqualization" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI. 
  ///Note that "LocImageROI" must have the same size as "this"
  ///ROIs are defined by non-null values
  virtual void HistoEqualizationInROI(ScalarField * LocImageROI);
  
  ///Compute the histogram. The histogram is normalized (sum of values equal to 1). 
  /// -> Input_BinsNb is the number of bins in the histogram
  /// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
  void CptHistogram(int Input_BinsNb,float * Output_Histo_x_axis,float * Output_Histo_y_axis,int useLogHisto=0);

  ///Same as "CptHistogramCptHistogram" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI.
  ///Note that "ImageROI" must have the same size as "this" and the ROI is defined by non-null values
  void CptHistogram_InROI(ScalarField * ImageROI,int Input_BinsNb,float * Output_Histo_x_axis,float * Output_Histo_y_axis,int useLogHisto=0);

  ///Compute the cumulative histogram. The cumulative histogram is normalized (last value equal to 1). 
  /// -> Input_BinsNb is the number of bins in the histogram
  /// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
  void CptCumulativeHistogram(int Input_BinsNb,float * Output_CumHisto_x_axis,float * Output_CumHisto_y_axis,int useLogHisto=0);

  ///Same as "CptCumulativeHistogram" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI.
  ///Note that "ImageROI" must have the same size as "this" and the ROI is defined by non-null values
  void CptCumulativeHistogram_InROI(ScalarField * ImageROI,int Input_BinsNb,float * Output_CumHisto_x_axis,float * Output_CumHisto_y_axis,int useLogHisto=0);

};

///get the 'Image 2 World matrix' of an image without loading the image  (only its header)
void Get_Image2World(char *,float LocI2W[4][4]);



///1.2) scalar field of integers


class IntScalarField{
  /// ******************************************************************************
private:
  //size of the fields
  int NXtY;      //NX*NY
  
  //scalar field
  int * ScalField;
  
  /// ******************************************************************************
public:
  //size of the fields
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  
  
  /// Constructor and destructor
  IntScalarField();
  ~IntScalarField();


  /// functions associated to the field
  
  //allocate memory for the field and set all values to zero
  virtual void CreateVoidField(int NBX,int NBY,int NBZ=1);
  
  //write a scalar field in a nifti image
  virtual void Write(char *);
  
  // put a the same value at every points of the scalar field
  virtual void PutToAllVoxels(int cste);
  
  
  //put a value in the scalar field
  virtual inline void P(int value,int x,int y,int z=0){
    this->ScalField[z*this->NXtY+y*this->NX+x]=value;
  }
  
  //get a value from the scalar field
  virtual inline int G(int x,int y,int z=0){
    return this->ScalField[z*this->NXtY+y*this->NX+x];
  }
};


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           2:   FUNCTIONS FOR THE CLASS "VectorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class VectorField{
  /// ******************************************************************************
private:
  //size of the fields
  int NXtY;      //NX*NY
  int NXtYtZ;    //NX*NY*NZ
  int NXtYtZtT;  //NX*NY*NZ*NT
  
  //scalar field
  float * VecField;
  
  /// ******************************************************************************
public:
  //size of the fields
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  int NT;        //image size on time
  
  float Image2World[4][4];  //quaternion to go from the image coordinates to the world coordinates
  float World2Image[4][4];  //quaternion to go from the world coordinates to the image coordinates

  
  /// Constructor
  VectorField();
  
  /// Destructor
  ~VectorField();
  
  //put a value in the scalar field
  virtual inline void P(float value,int IdDirec,int x,int y,int z=0,int t=0){
    this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x]=value;
  }
  
  //Add a value to the scalar field
  virtual inline void Add(float value,int IdDirec,int x,int y, int z=0,int t=0){
    this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x]+=value;
  }
  
  //put the same value at all entries of the vector field
  virtual void PutToAllVoxels(float cste,int t=0);
  
  //get a value from the scalar field
  virtual inline float G(int IdDirec,int x,int y,int z=0,int t=0){
    return this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x];
  }
  
  //get a value from the scalar field by linear interpolation
  virtual float G(int IdDirec,double x,double y,double z=0.,int t=0);
  virtual float G(int IdDirec,float x,float y,float z=0.,int t=0);
  
  
  // get the maximum of the absolute values of the vector field
  virtual float GetMaxAbsVal(int t=0);
  
  //read a vector field (from 3 nifti images -> X, Y, Z)
  virtual void Read(char *,char *,char *);
  
  //read a vector field and perform linear interpolation to give it a specific size
  //  If rescaleVF!=0, the values of the vector field are rescaled proportionally to
  //  the re-sizing (usefull for velocity fields)
  virtual void Read_and_Interpolate(char *,char *,char *,int,int,int,int rescaleVF=0);
  
  //create a void vector field. No message is printed if  Verbose!=1
  virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1,int Verbose=1);

  //Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
  virtual void SlashFieldSize(int verbative=1);
  
  
  //for each point (x,y,z) and direction direc: this(direc,x,y,z,0) <- \sum_t TimeDepVF(direc,x,y,z,t)
  //TimeDepVF must obviously have the same size as this. If yes, 1 is returned / 0 is returned otherwise 
  virtual int ProjectTimeDepVF_to_StationaryVF(VectorField * TimeDepVF);
  
  //write a vector field (in 3 nifti images -> X, Y, Z)
  //The optional 4th file is an input file containing the headers 
  virtual void Write(char *,char *,char *);
  virtual void Write(char *,char *,char *,char *);
  
  //Copy by reference the vectorized velocity field -> Access to a value: IdDirec*NX*NY*NZ*NT + t*NX*NY*NZ + z*NX*NY + y*NX + x
  //The whole vector size is also given in VectorSize
  virtual float * CopyByRefVectorizedVF(int * VectorSize);

};


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                            3:   FUNCTIONS FOR THE CLASS "TensorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class TensorField{
  /// ******************************************************************************
private:
  //size of the fields
  //int NXtY;      //NX*NY
  //int NXtYtZ;    //NX*NY*NZ
  //int NXtYtZtT;  //NX*NY*NZ*NT
  //int NXtYtZtTt3;  //NX*NY*NZ*NT*3
  
  int NXt9;         //NX*9
  int NXtYt9;       //NX*NY*9
  int NXtYtZt9;     //NX*NY*NZ*9  
  
  
  //scalar field
  float * TField;
  
  /// ******************************************************************************
public:
  //size of the fields
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  int NT;        //image size on time
  
  /// Constructor
  TensorField();
  
  /// Destructor
  ~TensorField();
  
  //create a void tensor field
  virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1);
  
  //put a value in the tensor field
  virtual void P(float value,int IdDirec1,int IdDirec2,int x,int y,int z=0,int t=0);
  
  //add a value in the tensor field
  virtual void Add(float value,int IdDirec1,int IdDirec2,int x,int y,int z=0,int t=0);
  
  //get a value from the tensor field
  virtual float G(int IdDirec1,int IdDirec2,int x,int y,int z=0,int t=0);
  
  //put all values to 0
  virtual void PutAllValuesToZero();
  
  //Add a tensorised vector to the existing tensor
  virtual void AddTensorisedVector(float vec[3],int x,int y,int z=0,int t=0,float weight=1);
  
  //Perform a principal component analysis of the 3*3 tensor
  //The outputs are:
  // lambda1,lambda2,lambda3: the eigenvalues in decrasing order
  // vec1: 1st eigenvector  (must be initialised as vec1[3])
  // vec2: 2nd eigenvector  (must be initialised as vec2[3])
  // vec3: 3rd eigenvector  (must be initialised as vec3[3])
  virtual void PCA(float vec1[3],float vec2[3],float vec3[3],float * lambda1, float * lambda2, float * lambda3, int x,int y,int z=0,int t=0);
  
};





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           4: CLASSES TO PERFORM CONVOLUTION AND DECONVOLUTION USING FFT
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///4.1: Main class (sufficiently general to have kernels which are not separable but memory consuming)
class FFTconvolver3D{
  /// ******************************************************************************
private:
  //size of the inputs (that will be copied in images having sizes = 2^{...})
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  
  //size of the fields transformed by fft
  int NXfft;
  int NYfft;
  int NZfft;
  
  //fields transformed by fft
  ScalarField RealSignalForFFT;
  ScalarField ImagSignalForFFT;
  ScalarField RealFilterForFFT;
  ScalarField ImagFilterForFFT;
  
  //temporary scalar or vector field
  ScalarField ImageTemp;
  VectorField VecFieldTemp;
  
  //design a kernel that is the sum of up to 7 Gaussians
  //... if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
  void MakeSumOf7AnisotropicGaussianFilters(float weight1=100,float sigmaX1=1,float sigmaY1=1,float sigmaZ1=1,float weight2=0,float sigmaX2=1,float sigmaY2=1,float sigmaZ2=1,float weight3=0,float sigmaX3=1,float sigmaY3=1,float sigmaZ3=1,float weight4=0,float sigmaX4=1,float sigmaY4=1,float sigmaZ4=1,float weight5=0,float sigmaX5=1,float sigmaY5=1,float sigmaZ5=1,float weight6=0,float sigmaX6=1,float sigmaY6=1,float sigmaZ6=1,float weight7=0,float sigmaX7=1,float sigmaY7=1,float sigmaZ7=1,int NormalizeWeights=1);

  //Design a kernel representing a 3D second order B-Spline
  //StepSize controls the extent of the kernel  (step between two nodes in voxels)
  void MakeBSplineFilter(int StepSize);

  
  //Fast Fourier Transform
  void DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal);
  
  //Inverse Fast Fourier Transform
  void InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal);
  
  //Fast Fourier Transform of numerical recipies (slighly modified)
  void four1NR(float data[], unsigned long nn, int isign);
  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  FFTconvolver3D();
  
  //Destructor
  ~FFTconvolver3D();
  
  //Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
  //7 Gaussians (set some weights to 0 if less Gaussians are required)
  //* NX, NY, NZ: is the size of the input image
  //* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //* w2,sX2,sY2,sZ2,: weight of the 2nd Gaussian kernel and std. dev. in direction X, Y, Z
  //* w3,sX3,sY3,sZ3,: weight of the 3rd Gaussian kernel and std. dev. in direction X, Y, Z
  //* w4,sX4,sY4,sZ4,: weight of the 4th Gaussian kernel and std. dev. in direction X, Y, Z
  //* w5,sX5,sY5,sZ5,: weight of the 5th Gaussian kernel and std. dev. in direction X, Y, Z
  //* w6,sX6,sY6,sZ6,: weight of the 6th Gaussian kernel and std. dev. in direction X, Y, Z
  //* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolver(int NBX,int NBY, int NBZ, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
  
  //Allocate memory to perform convolutions with matrices  (for now: with the function Convolve_divFree)
  virtual void InitiateMatrixValuedConvolver(int NBX,int NBY, int NBZ);
  
  
  //Initiate the complex fields for the FFT and the smoothing kernels being parametrized by B-splines
  virtual void InitiateSplineConvolver(int NBX,int NBY, int NBZ,int GridStep);
  
  //change the kernel of the convolver (same notations as the constructor)
  //... here the new kernel is normalized (w1+...+w7 is normalized to 1)
  virtual void ChangeKernel(float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., 
                            float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., 
                            float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., 
                            float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1.,
                            float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1.,
                            float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1.,
                            float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
  
  //change the kernel of the convolver (same notations as the constructor)
  //... here the new kernel is not normalized
  virtual void ChangeKernel_SingleScale(float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1.);
  
  
  //convolution of a 3D scalar field using the predifined kernel
  virtual void Convolution(ScalarField *);
  
  //convolution of a 3D vector field using a divergence free (matrix valued) kernel, based on a Gaussian kernel of std dev sigma
  //   ->    \phi_{\alpha}(x) = \{ [ 2 (s-1) \alpha -  4 \alpha^2 ||x||^2  ] Id  + 4 \alpha^2 xx' \}  e^{\alpha ||x||^2}  ... where s is the image dimension  
  //   ->    F. J. Narcowich and J. D. Ward,   Generalized Hermite interpolation via matrix-valued conditionally positive definite functions, Math. Comp. 63 (1994), 661-687
  virtual void Convolution_DivergenceFreeFilter(VectorField *, float sigma);

  //convolution of the real scalar field defined inside of the class using the predifined kernel
  virtual void Convolution();
  
  //deconvolution of a 3D scalar field using the predifined kernel
  // !!! NOT VALIDATED !!!
  virtual void Deconvolution(ScalarField *);
  
  //put a value in the real part of the field that is transformed by the class
  virtual void P(float value,int x,int y, int z=0);
  
  //put a value in the real part of the field that is transformed by the class
  virtual float G(int x,int y, int z=0);
};


void SmoothVFUsingFFT(VectorField * SmoothedField,FFTconvolver3D * FFTconvolver_loc);



///4.2: class where we consider two regions

class MultiRegionFFTConvolver{
  /// ******************************************************************************
private:
  FFTconvolver3D Region0_convolver;
  FFTconvolver3D Region1_convolver;
  ScalarField PartitionOfUnity;
  int TypeOfConvolver1,TypeOfConvolver0; //(==-1 if non-initialised convolver / ==0 if standard convolver / ==1 if divergence free convolver)
  VectorField TempVF0;
  VectorField TempVF1;
  float SigmaReg0_IfDivFree,SigmaReg1_IfDivFree;  //div. free kernels only use a Gaussian kernel (and not the sum of kernel) with stddiv=sX1 when calling 'InitiateConvolver_RegN'
  float TypicalScaleVF0;  //Typical scale of the update VF0
  float TypicalScaleVF1;  //Typical scale of the update VF1
  
  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  MultiRegionFFTConvolver();
  
  //Destructor
  ~MultiRegionFFTConvolver();
  
  //Initiate the convolver of region 0
  //7 Gaussians (set some weights to 0 if less Gaussians are required)
  //* NX, NY, NZ: is the size of the input image
  //* If DivFree==1 then the kernel is divergence free
  //* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //* ...
  //* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  //IMPORTANT REMARK: If divergence free kernel, we only define a Gaussian kernels of stddev sX1
  virtual void InitiateConvolver_Reg0(int NBX,int NBY, int NBZ, int DivFree=0, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
  
  //Initiate the convolver of region 1
  //7 Gaussians (set some weights to 0 if less Gaussians are required)
  //* NX, NY, NZ: is the size of the input image
  //* If DivFree==1 then the kernel is divergence free
  //* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //* ...
  //* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  //IMPORTANT REMARK: If divergence free kernel, we only define a Gaussian kernels of stddev sX1
  virtual void InitiateConvolver_Reg1(int NBX,int NBY, int NBZ, int DivFree=0, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);

  //Set the partition of unity (it contains values between 0 and 1 -> if==0 then region 0 / if==1 then region 1 / otherwise something inbetween)
  virtual void SetPartitionOfUnity(ScalarField * RefPartitionOfUnity);
  
  
  //convolution of a 3D vector field using the predifined kernel
  virtual void Convolution(VectorField * VF);
};


///4.3: light weight convolver (requires little memory but can only have separable kernels. Contains the most common functions of FFTconvolver3D, plus a direct smoothing of velocity fields)
class LightFFTconvolver3D{
  /// ******************************************************************************
private:
  //size of the inputs (that will be copied in images having sizes = 2^{...})
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  
  //size of the fields transformed by fft
  int NXfft;
  int NYfft;
  int NZfft;
  
  //fields transformed by fft
  ScalarField RealSignalForFFT_X;
  ScalarField ImagSignalForFFT_X;
  ScalarField RealSignalForFFT_Y;
  ScalarField ImagSignalForFFT_Y;
  ScalarField RealSignalForFFT_Z;
  ScalarField ImagSignalForFFT_Z;
  
  ScalarField RealFilterForFFT_X;
  ScalarField ImagFilterForFFT_X;
  ScalarField RealFilterForFFT_Y;
  ScalarField ImagFilterForFFT_Y;
  ScalarField RealFilterForFFT_Z;
  ScalarField ImagFilterForFFT_Z;
  
  //design a kernel that is the sum of up to 7 Gaussians
  //... if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
  void MakeSumOf7AnisotropicGaussianFilters(float weight1=100,float sigmaX1=1,float sigmaY1=1,float sigmaZ1=1,float weight2=0,float sigmaX2=1,float sigmaY2=1,float sigmaZ2=1,float weight3=0,float sigmaX3=1,float sigmaY3=1,float sigmaZ3=1,float weight4=0,float sigmaX4=1,float sigmaY4=1,float sigmaZ4=1,float weight5=0,float sigmaX5=1,float sigmaY5=1,float sigmaZ5=1,float weight6=0,float sigmaX6=1,float sigmaY6=1,float sigmaZ6=1,float weight7=0,float sigmaX7=1,float sigmaY7=1,float sigmaZ7=1,int NormalizeWeights=1);

  //Fast Fourier Transform
  // -> if axis == 0 -> FFT on X axis
  // -> if axis == 1 -> FFT on Y axis
  // -> if axis == 2 -> FFT on Z axis
  void DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis);
  
  //Inverse Fast Fourier Transform
  // -> if axis == 0 -> IFFT on X axis
  // -> if axis == 1 -> IFFT on Y axis
  // -> if axis == 2 -> IFFT on Z axis
  void InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis);
  
  //Fast Fourier Transform of numerical recipies (slighly modified)
  void four1NR(float data[], unsigned long nn, int isign);
  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  LightFFTconvolver3D();
  
  //Destructor
  ~LightFFTconvolver3D();
  
  //Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
  //7 Gaussians (set some weights to 0 if less Gaussians are required)
  //* NX, NY, NZ: is the size of the input image
  //* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //...
  //* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolver(int NBX,int NBY, int NBZ, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.,int NormalizeWeights=1);
  
  //change the kernel of the convolver (same notations as the constructor)
  //... here the new kernel is normalized (w1+...+w7 is normalized to 1)
  virtual void ChangeKernel(float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1.,float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1.,float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1.,float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1.,float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1.,float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1.,float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.,int NormalizeWeights=1);

  //convolution of a 3D scalar field using the predifined kernel
  //If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
  virtual void Convolution(ScalarField *,int TimeFrame=-1);
  
  //convolution of a 3D vector field using the predifined kernel
  //If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
  virtual void Convolution(VectorField *,int TimeFrame=-1);

  //convolution of a 3D vector field using the predifined kernel. Convolution is performed in the ROI defined by (xmin, xmax, ymin, ymax, zmin, zmax) only.
  //If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
  virtual void ConvolutionInROI(VectorField *,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax,int TimeFrame=-1);

  //Hack to perform convolution of the 3D vector field 'VF' in a masked region with mirror conditions.
  //  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
  virtual void Convolution_Mask_Mirror(VectorField * VF,ScalarField * Mask, int MaskId=1);
  
  //Hack to perform convolution of the 3D scalar field 'SF' in a masked region with mirror conditions.
  //  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
  virtual void Convolution_Mask_Mirror(ScalarField * SF,ScalarField * Mask, int MaskId=1);
};



///4.4: class where we consider N regions to smooth

class MultiRegionFFTConvolver2{
  /// ******************************************************************************
private:
  LightFFTconvolver3D * Region_convolver;
  ScalarField PartitionOfUnity; //3D + channels
  int * xmin;
  int * xmax;
  int * ymin;
  int * ymax;
  int * zmin;
  int * zmax;
  int NumberOfRegions;
  VectorField TempVF1;
  VectorField TempVF2;

  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  MultiRegionFFTConvolver2();
  
  //Destructor
  ~MultiRegionFFTConvolver2();
  //Initiate the convolver in all regions using same kernel
  //-> 'Part_Of_Unity' is a 3D mask which define different ROIs. It only contains integer values each of them associated to a ROI. Partition of unity is first 
  //   defined by splitting this mask into several channels, each of them having an intensity equal to 1 in the corresponding ROI and 0 otherwise. Each channel
  //   is then smoothed with a Gaussian kernel of stddev 'sigmaPOI'.
  //-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
  //     * NX, NY, NZ: is the size of the input image
  //     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //     * ...
  //     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolver(ScalarField * Part_Of_Unity, float sigmaPOI=5., float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);

  //Initiate the convolver in all regions using same kernel
  //-> Part_Of_Unity is a '3D + channels' scalar field which encodes the partition of unity in the different channels.
  //   * Its size and number of channels (NBT actually) defines the size and the number of regions of the convolver 
  //   * The maximum point-wise sum of the probabilities may be different to 1: normalisation will be automatically performed
  //   * Point-wise sum of the probabilities may vary in space. If so, a background region will be automatically defined
  //-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
  //     * NX, NY, NZ: is the size of the input image
  //     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //     * ...
  //     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolverWithActualPOI(ScalarField * Part_Of_Unity, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);

  //save the actual partition of unity (after potential treatments in InitiateConvolver or undersampling)
  //-> 1st char* is the name in which the image is saved
  //-> 2nd char* is the name of the image that will be used in the header of the saved image
  virtual void SaveActualParitionOfUnity(char *,char *);

  //change the smoothing kernel in one region
  virtual void ChangeKernelInOneRegion(int IdRegion, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
  
  //Update the parition of unity 
  //-> it must have the same size and number of layers/times as in the current POU (tested)
  //-> to make sense, the new POU must have a sum of intensities across layers/times equals to 1 at each voxel (not tested)
  virtual void UpdatePartitionOfUnity(ScalarField * Part_Of_Unity);

  //convolution of a 3D vector field using the predifined kernel
  virtual void Convolution(VectorField * VF);
  virtual void ConvolutionOnROI(VectorField *,int idRegion);  // Hack FX
  
  //return the number of regions considered
  virtual int GetRegionNb();
};



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                  5: CLASS TO MANAGE THE MUTUAL INFORMATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MImanager{
  /// ******************************************************************************
private:
  //size of the treated images
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  
  //nbBins
  int NumberOfBinsS;
  int NumberOfBinsT;
  float MinGreyLevelsS;
  float MinGreyLevelsT;
  float SizeStepsGreyLevelsS; 
  float SizeStepsGreyLevelsT; 
  
  //margin on which the mutual information is not computed, plus mask indicator
  int Margin;
  int indicatorMaskDefined;
  
  //registered images and mask (ROI is where the mask intensities are 1)
  ScalarField * SrcImage;
  ScalarField * TrgImage;
  ScalarField * Mask;
  
  //joint Histo
  float ** JointHistogram;
  float * MarginalHistogramS;
  float * MarginalHistogramT;
  
  //joint and marginal entropies
  float JointEntropy;
  float MarginalEntropyS;
  float MarginalEntropyT;
  
  //normalised mutual informations
  float MI;
  
  //indicator saying that the saved computations are made on up-to-date data
  int indicatorUpdatedSrcHisto;
  int indicatorUpdatedTrgHisto;
  int indicatorUpdatedJointHisto;
  int indicatorUpdatedMI;
  
  
  //normalized the intensities adapted to the bins of the histograms
  float GiveFloatBinT(float intensity);
  float GiveFloatBinS(float intensity);
  
  
  //when computing the histograms give the contribution of 'intensity' in the bin 'IdBin' for the source and target image
  float GiveValueParzenWindowS(float intensity,int IdBin);
  float GiveValueParzenWindowT(float intensity,int IdBin);
  
  //compute the histograms
  void ComputeJointHistogramAndEntropy();
  void ComputeMarginalHistogramAndEntropyS();
  void ComputeMarginalHistogramAndEntropyT();
  
  
  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  MImanager();
  
  //Destructor
  ~MImanager();
  
  //Initiate the MI manager without any mask
  //Important: the current MImanager will always point to the source and target images defined here (with intensities that can change!!!)
  virtual void Initiate(ScalarField * SourceImage,ScalarField * TargetImage,int NbBinsSrc=30,int NbBinsTrg=30, int LocMargin=3);

  //Initiate the MI manager with a mask. The MI and MI gradients will only be computed where the mask equals 1
  virtual void Initiate(ScalarField * SourceImage,ScalarField * TargetImage,ScalarField * ROI_Mask,int NbBinsSrc=30,int NbBinsTrg=30, int LocMargin=3);

  //returns the normalized mutual information (plus update all the histograms)
  virtual float EvaluateMI();
  
  //returns the estimated gradient of normalized mutual information
  virtual void EvaluateGradMI(VectorField * Gradient);
  
  //Indicate to the MI manager that the Source image has changed
  virtual void IndicateSrcHasChanged();
  
  //Indicate to the MI manager that the Target image has changed
  virtual void IndicateTrgHasChanged();
};







///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           6: LOW LEVEL FUNCTIONS MAKING USE OF THE CLASSES ScalarField AND VectorField 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///+++++++++++++++++++++++++++++++++      6.1: Diffusion stuffs       +++++++++++++++++++++++++++++++++

///Isotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
void Diffusion_3D(ScalarField * SField, float alpha, float dTau,int ITERATIONS_NB, float dx=1,float dy=1, float dz=1);
void Diffusion_3D(VectorField * VField, float alpha, float dTau,int ITERATIONS_NB, float dx=1,float dy=1, float dz=1);

///Isotropic diffusion of 'SField' or 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
///If optNoMaskNoDef==1 then the values for which the mask is 0 are set to 0
void Diffusion_3D(ScalarField * SField,ScalarField * Mask, int MaskId, float alpha, float dTau,int ITERATIONS_NB, int optNoMaskNoDef=0, float dx=1,float dy=1, float dz=1);
void Diffusion_3D(VectorField * VField,ScalarField * Mask, int MaskId, float alpha, float dTau,int ITERATIONS_NB,  int optNoMaskNoDef=0,int direction=-1, float dx=1,float dy=1, float dz=1);


///Isotropic diffusion of 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///!!! -> Dirichlet conditions are considered at the domain boundaries !!!
void Diffusion_3Dbis(VectorField * VField,ScalarField * Mask,int MaskId, float alpha, float dTau,int ITERATIONS_NB, int optNoMaskNoDef=0, float dx=1,float dy=1, float dz=1);


///Anisotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Directions and intensities of diffusion are defined using ax, ay, az.
///The size of each voxel is defined by dx, dy, dz.
///Remark: very similar to the classic paper of Perona-Malik
void anisoDiff_3D(ScalarField * SField,float ax, float ay, float az, float dTau,int ITERATIONS_NB, float dx=1,float dy=1, float dz=1);


///Smooth the field 'SmoothedField' for fluid like regularisation with sliding conditions at th boundaries of 'LocMask'.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> 'NeaBoun' must contain a vector field of the same size as 'SmoothedField' (treated as a temporary field)
///-> 'NormalCmp'  must contain a scalar field of the same size as 'SmoothedField' (treated as a temporary field)
///-> 'TempSF'  must contain a scalar field of the same size as 'SmoothedField' (treated as a temporary field)
void SmoothNormalAndTangentContributions(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB,int NbRegions, float * IdRegions,VectorField * NeaBoun,ScalarField * NormalCmp,ScalarField * TempSF, float dx=1,float dy=1, float dz=1);




///Smooth the field 'SmoothedField' for diffusion like regularisation with sliding conditions at the boundaries of 'LocMask'.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> If 'NoDistMapEstim' == 1 the distance map (NeaBoun) is not estimated here. In that case 'TempSF' is also not reestimated.
void SmoothWithinRegions(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB, int NbRegions,float * IdRegions, float dx=1,float dy=1, float dz=1,int DirchletNeumann=1);


///Smooth the field 'SmoothedField' for diffusion like regularisation with sliding conditions at the boundaries of 'LocMask'. 
///Contrary to what is done in 'SmoothWithinRegions', we put to 0 the contributions too close to the bounaries.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary in 'LocMask'
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which the SmoothedField is set to zero before smoothing
void SmoothWithinRegions2(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB, int NbRegions,float * IdRegions,VectorField * NeaBoun,ScalarField * TempSF, float dx=1,float dy=1, float dz=1,int DistMapToEstim=0);



///Remove the normal contributions of a velocity field that are too close to a boundary
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which SmoothedField has reduced normal contirbutions
void RemoveNormalContributions(VectorField * SmoothedField,VectorField * NeaBoun,ScalarField * TempSF,float BoundaMargin, float dx,float dy, float dz,int SetBoundaryToZero=0);



///Remove the normal contributions of a velocity field that are too close to a boundary
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which SmoothedField has reduced normal contirbutions
void RemoveNormalContributionsOLD(VectorField * SmoothedField,VectorField * NeaBoun,ScalarField * TempSF,float BoundaMargin, float dx=1,float dy=1, float dz=1);



///Smooth the field 'SmoothedField' for diffusion like regularisation with sliding conditions at the boundaries of 'LocMask'. 
///Contrary to what is done in SmoothWithinRegions, we reduce the contribution of the normals close to the bounary.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> 'NeaBoun' must contain a vector field of the same size as 'SmoothedField' (treated as a temporary field)
///-> 'TempSF'  must contain a scalar field of the same size as 'SmoothedField' (treated as a temporary field)
///-> BoundaMargin is the margin around the bounadrie (in voxels) in which the SmoothedField is set to zero before smoothing
///DEPRECATED FUNCTION - SHOULD NOT BE USED
void SmoothWithinRegionsWithSliding(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB, int NbRegions,float * IdRegions,VectorField * NeaBoun,ScalarField * TempSF, float dx=1,float dy=1, float dz=1,float BoundaMargin=1,int NoDistMapEstim=0,int BoundaryCompensation=1);


///+++++++++++++++++++++++++++++++++      6.2: derivatives and diamonds      +++++++++++++++++++++++++++++++++

///6.2.1: derivatives

//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_Grad_ScalarField(ScalarField * SField,VectorField * Gradient,int SpecificTimeFrame=-1,float DeltaX=1);
void Cpt_Grad_MaskedScalarField(ScalarField * SField,VectorField * Gradient,ScalarField * Mask,int SpecificTimeFrame=-1,float DeltaX=1);


//Compute (d VField(X) / d x) + (d VField(Y) / d y) + (d VField(Z) / d z) and put the result in 'GradScalVF'
//where 'VField' is a vector field and 'GradScalVF' a scalar field
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_Grad_Scal_VectorField(VectorField * VField,ScalarField * GradScalVF,int SpecificTimeFrame=-1,float DeltaX=1);

//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame=-1,float DeltaX=1);




///6.2.2: diamonds

//VF = SF1 \diamond_{Diff} SF2
//which means: VF = - SF2 \nabla SF1                    
void Diamond_Diff(ScalarField * SF1,ScalarField * SF2,VectorField * VF);


//V  = SF1 \diamond_{Trans} SF2
//which means: V =  - \int_{R^3} ( SF2 \nabla SF1  )
void Diamond_Trans(ScalarField * SF1,ScalarField * SF2,float V[3]);


//T = V1 \diamond_{R^3} V2 
//which means:  T =   -0.5  (V1 \otimes V2 - V2 \otimes V1 )   
void Diamond_R3(float V1[3],float V2[3],float T[3][3]);


//T = SF1 \diamond_{SO(3)} SF2 
//which means:  T =   -0.5  \int_{R^3} \left(    SF2(x) (  (\nabla SF1)(x) \otimes x - x \otimes (\nabla SF1)(x) )  \right) dx
//Remark: 
//  x represents where we are in the image. 
//  Here it can have two interpretations: (1) image coordinates and (2) where we are depending on the center of rotation
//  If you don't want the center of rotation to be the image coordinate (0,0,0), fill the image coordinate you want in {Cx,Cy,Cz}
void Diamond_SO3(ScalarField * SF1,ScalarField * SF2,float T[3][3],float Cx=0,float Cy=0,float Cz=0);




///+++++++++++++++++++++++++++++++++      6.3: image deformation and mappings      +++++++++++++++++++++++++++++++++

///6.3.0: affine transformation

///Image transformation using the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfo(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,ScalarField * TransformedImage);


///6.3.1: displacement fields

///Compose RefField with UpdateField. The result is saved in RefField  (we consider that a disp. f. points from the source to the target)
void DisplacementFieldCompose(VectorField * RefField,VectorField * UpdateField);

///Compose InvUpdateField with InvRefField. The result is saved in InvRefField  (we consider that an inv. disp. f. points from the target source to the source)
void InvDisplacementFieldCompose(VectorField * InvUpdateField,VectorField * InvRefField);

///We consider here that the vectors of DispField point from 'StaticImage' to 'DeformedImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
///remark: 'DefoField' has to be sufficiently smooth to allow an accurate inversion
void ProjectImageUsingDispField(VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0);


///We consider here that the vectors of InvDispField point from 'DeformedImage' to 'StaticImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
void ProjectImageUsingInvDispField(VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0);



///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndInvDispField(float ProjectCS_2_OriginCS[4][4],VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN=0);
void ProjectImageUsingDispFieldAndInvDispField(VectorField * DispField,VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage, int NN=0);


void ProjectImageUsingAffineTransfoAndDispField(float ProjectCS_2_OriginCS[4][4],VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage);


///6.3.2: stationary (steady) velocity fields


///Integrate a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptDefFromSteadyVeloField(VectorField * VeloField,VectorField * DeformationField,int log2TimeStepNb,float MultFactor=1);

///Integrate the inverse of a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptInvDefFromSteadyVeloField(VectorField * VeloField,VectorField * InvDeformationField,int log2TimeStepNb,float MultFactor=1);


///compute the Lie Braket of VF1 and VF2. Put the result in VF3.   (subfunction of ComposeTwoLogFieldsUsingBCH)
void LieBracket(VectorField * VF1,VectorField * VF2,VectorField * VF3);

///RefVeloField and UpdateVeloField are two steady velocity fields. Their exponentials are respectively the current deformation
///and the update of the deformation. This function approximates the velocity field which is the log of the composition 
///between the current deformation and the update deformation (Baker-Campbell-Hausdorff formula) (cf Vercauteren MICCAI 2008)
///In output, RefVeloField is the updated velocity field. UpdateVeloField is also modified for algorithmic reasons but it
///represents nothing pertinent as an output.
void ComposeTwoLogFieldsUsingBCH(VectorField * RefVeloField,VectorField * UpdateVeloField);
void ComposeTwoLogFieldsUsingSum(VectorField * RefVeloField,VectorField * UpdateVeloField);

///... explicit name ....
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
void ProjectImageUsingInverseSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0,float factor=1);

void ProjectImageUsingSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0,float factor=1);
void ProjectImageUsingAffineTransfoAndSteadyVeloField(float ProjectCS_2_OriginCS[4][4],VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN=0,float factor=1);
void ProjectImageUsingDispFieldAndSteadyVeloField(VectorField * DispField,VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN=0,float factor=1);

///6.3.3: time-dependent velocity fields (no momenta)


//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'.
//'MappingAtRefTimeStep' is the mapping at refTimeStep (possibly not the identity).
//* An iterative leap-frog like technique is performed to compute the backward mapping. The more
//  iterations (='ConvergenceSteps') the more accurate the mapping. For most of the deformations
//  1 iteration is far enough but more iterations are suitable if the deformations have large 
//  Jacobians.
//* 'DeltaX' is the spatial step between two voxels.
void CptMappingFromVeloField(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps=1,float DeltaX=1);
void CptMappingFromVeloField_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps=1,float DeltaX=1);

//Do the same thing as the initialisation of CptMappingFromVeloField -> load the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
//The function 'CptMappingFromVeloField2_Increment' then allows to increment the field backward or forward according to 'VeloField'
void CptMappingFromVeloField2_Init(VectorField * MappingAtRefTimeStep,VectorField * Map);
void CptMappingFromVeloField2_Init_IniIdMap(VectorField * Map);

//Do the same thing as an incrementation of CptMappingFromVeloField from 'CurrentTimeStep' and backward (BackwardOrForward==-1) or forward (BackwardOrForward==1) 
//The function 'CptMappingFromVeloField2_Init' is then supposed to have loaded the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
void CptMappingFromVeloField2_Increment(VectorField * VeloField,VectorField * Map,int CurrentTimeStep,int BackwardOrForward,int ConvergenceSteps=1,float DeltaX=1);

//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialMapping' which is the partial mapping of 'MappingAtRefTimeStep' from the time 
//subdivision 'refTimeStep' due to the contribution of 'PartialVeloField'. Note, that an Identity mapping 'MappingId' is //also defined in the inputs (to avoid defining it each time the function is used)
void CptPartialMappingFromVeloFields(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * MappingId,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps=1,float DeltaX=1);
void CptPartialMappingFromVeloFields_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps=1,float DeltaX=1);

//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialLocMap' which is the partial mapping ONLY AT 'TargetSubdiv' FROM 'SourceSubdiv' due to the contribution of PartialVeloField.
//-> PartialLocMap therefore represents where are the coordinates of the points of time subdivision 'SourceSubdiv' when transported on time subdivision 'TargetSubdiv'
void ComputeLagrangianPartialMapping(int SourceSubdiv,int TargetSubdiv,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialLocBmap,float DeltaX=1);

//Compute the projection of a 3D image 'ImagToPropag' using Mapping 'Map'.
//The image is projected at the time step 'TimeStepProj' of 'Map' and stored in 'ImageTimeT'.
//
//Importantly, we consider that 'Map' has an identity transformation at the time step (say 't') where 'ImagToPropag' is.
//'Project3Dimage' can therefore perform a forward projection (t<TimeStepProj) or a backward projection (TimeStepProj<t).
void Project3Dimage(ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj);

void Project3DImageUsingAffineTransfoAndTimeDepVF(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj);

void Project3DImageUsingDispFieldAndTimeDepVF(VectorField * DispField,ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj);

void Project3Dimage(VectorField * VFToPropag,VectorField * Map,VectorField * VFTimeT,int TimeStepProj);

void EulerScheme(int refTimeStep,VectorField * VeloField,VectorField * Map,float DeltaX=1);


//Not defined yet
void SimpleEulerStep(VectorField *VectorField1,VectorField *VectorField2, int time, float dt);


//By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
//'VeloField4Measure' in the length of the flow from each point of the field. The length of flow
//is projected AT T=0 and returned in the 3D scalar field 'LengthOfFlow'
// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
//   is computed.
void CptLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps=3,float DeltaX=1);

//By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
//'VeloField4Measure' AT THE CURRENT TIME in the length of the flow from each point of the field. The length of flow
//is returned in the 3D+t scalar field 'LengthOfFlow'
// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
//   is computed.
void CptEvoLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps=3,float DeltaX=1);





///6.3.4: time-dependent velocity fields (with momenta)

// Computes the transport of the initial momentum by the diffeo and stores it in Momentum
//void TransportMomentum(ScalarField *InitialMomentum, VectorField *InvDiffeo, ScalarField *Momentum,float DeltaX,int t=0);
void TransportMomentum(ScalarField *InitialMomentum, VectorField *InvDiffeo, ScalarField *Momentum,int t=0,float DeltaX=1);

// Computes cste * transport momentum from the initial momentum by the diffeo and add it in Image.
void AddTransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,float cste=1.0, int t=0);

// Computes the transport image from the initial image by the diffeo and stores it in Image.
void TransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image,int t=0);

// Computes cste * transport image from the initial image by the diffeo and add it in Image.
void AddTransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image,float cste=1.0,int t=0);



///+++++++++++++++++++++++++++++++++      6.4: linear algebra       +++++++++++++++++++++++++++++++++


// Copies the values of a VectorField1(t=0) in VectorField2(t)
void DeepCopy(VectorField *VectorField1,VectorField *VectorField2,int t=0);
void DeepCopy(VectorField *VectorField,ScalarField* ScalarField,int direc,int t=0);

// Copies the values of a ScalarField1(t=0) in ScalarField2(t)
void DeepCopy(ScalarField *ScalarField1,ScalarField *ScalarField2,int t=0);

//copy and multiplies by 'cste" the values of VectorField1 in VectorField2
void DeepCopyMultiplyVectorField(VectorField *VectorField1,VectorField *VectorField2, float cste);

// Compute the L^2 scalar product and store it in ScalarField0
void ScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t=0,float cste = 1.0);
void ScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t=0, float cste = 1.0);

// Compute the L^2 scalar product between two vectorfields at time t and add it to ScalarField0
void AddScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t=0);
void AddScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t=0);

// Add  cste * ScalarField1 at time t1 to ScalarField2 at time t2
void AddScalarField(ScalarField *ScalarField1, ScalarField *ScalarField2,float cste,int t1 = 0,int t2=0);
// Add  cste * VectorField1 at time t1 to VectorField2 at time t2
void AddVectorField(VectorField *VectorField1, VectorField *VectorField2,float cste,int t1 = 0,int t2=0);
// Multiply a vector field by the cste
void MultiplyVectorField(VectorField *VectorField1, float cste,int t=0);
// 
void SumVectorField(VectorField *VectorField1, VectorField *VectorField2, VectorField *Output, int t1=0,int t2=0,int t3=0, float cste1 = 1.0,float cste2 =1.0);

// Compute the product element by element of ScalarField and VectorField and store it in VectorField2
void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2);
void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2,int time);
void AddProduct(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2);

// Compute the dot product
float DotProduct(ScalarField *ScalarField1, ScalarField *ScalarField2,int t1=0,int t2=0);
float DotProduct(VectorField *ScalarField1, VectorField *ScalarField2, int t1=0,int t2=0);



///+++++++++++++++++++++++++++++++++      6.5: misc       +++++++++++++++++++++++++++++++++


///compute the L_2 norm of the difference between two scalar fields
float CalcSqrtSumOfSquaredDif(ScalarField * I1,ScalarField * I2);

///Shrink regions contained in an image/mask to single points which can be considered as the regions center. These points are estimated using iterative erosions.
void ShrinkRegionsToSinglePoints(ScalarField * Mask, float MaskId,ScalarField * CenterPointsImage);

///Compute the map of the nearsest boundaries in the 3D mask 'Mask' (propagate the nearest boundaries identification until at a maximum distance of MaxDist)
///Outputs:
///  * NeaBoun: vector field pointing the nearest boundary
///  * TempSF: In the end the points for which the nearest boundary is identified have a value equal to 1
///Option:
///If NbRegions>0 then masks with partial volume effects are treated. In that case NbRegions is the nb of regions 
///considered in the mask and IdRegions are the corresponding regions 
///If 'LowCostPropagation==1' a wise algorithm propagates the distance map at a low algorithmic cost with a small approximation (otherwise no approximation but much slower algorithm)
void Cpt_NearestBoundaryOLD(ScalarField * Mask,VectorField * NeaBoun,ScalarField * TempSF,float MaxDist,int NbRegions=-1, float * IdRegions=NULL,float dx=1,float dy=1, float dz=1,int LowCostPropagation=1);

void Cpt_NearestBoundary2D(ScalarField * Mask,VectorField * NeaBoun,ScalarField * TempSF,int NbRegions=-1, float * IdRegions=NULL,float dx=1,float dy=1);

///compute the distance map of the edges in a 3D Mask. Greedy algorithm updating the distances in a band.
///-> 'NeaBoun' contains a vector field pointing to the nearest boundaries
///-> 'Distance' contains the distance to the nearest boundary
void Cpt_NearestBoundary(ScalarField * Mask,VectorField * NeaBoun,ScalarField * Distance);


///Stochastically propagate the non-null values of ImageIn in a ROI.
/// -> The ROI is defined where we have the non-null values of Mask.
/// -> Intensity propagation takes into account the distance between neighbors in world coordinates
void PropagateNonNullValues_Sto(ScalarField * ImageIn,ScalarField * Mask);

///Propagate the non-null values of ImageIn in a ROI. The value given to each null value of ImageIn in the ROI is the nearest non-null value of ImageIn in world coordinates
/// -> The ROI is defined where we have the non-null values of Mask.
/// -> Intensity propagation takes into account the distance between neighbors in world coordinates
void PropagateNonNullValues_NN(ScalarField * ImageIn,ScalarField * Mask);


void Cpt_DistMap(ScalarField * SegImage,float IDregion,ScalarField * DistMap);


void ImageDilation(ScalarField * TreatedImage,int iterationNb);


void ImageErosion(ScalarField * TreatedImage,int iterationNb);



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                   7: OTHER FUNCTIONS OF SCIENTIFIC COMPUTATION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// 7.1) ++++++++++++++++++ 3*3 and 4*4 matrices ++++++++++++++++++ 

//Solve the problem: MX=D where D is a known vector, M a tridiagonal matrix and X the unknown vector.
// Inputs are a,b,c,d,n where  M(i,i-1)=a(i), M(i,i)=b(i), M(i,i+1)=c(i), D(i)=d(i), D in R^n and M in R^n*R^n.
// Output is X where X in R^n.  Warning: will modify c and d!
void TridiagonalSolveFloat(const float *, const float *, float *, float *, float *, int);


//Perform the eigenvalue decomposition of a 3*3 matrix
//Adapted from the algorithm having the same name in 'numerical recipes'.
//Input:  The 3*3 matrix 'MatIni' that has to be symetric.
//Ouputs: 'ValP' is a vector of size 3 which contains the eigen values (in decreasing order). 
//        'VecP' is a 3*3 matrix containg the eigen vectors in columns.
void jacobi3(float **MatIni,float *ValP, float **VecP);


///compute two orthogonal vectors tvec1 and tvec2 in R^3 which are orthogonal to nvec
///the norm of tvec1 and tvec2 is defined as equal to the one of nvec
void CptVecsTangentPlane(float nvec[3],float tvec1[3],float tvec2[3]);

///normalize a vector
void VecNormalize(float vec[3],float norm);



///compute the determinant of a 3*3 matrix
float determinant_3t3matrix(float m[3][3]);

///compute the comatrix of a 3*3 matrix
void comatrix_3t3matrix(float m1[3][3],float m2[3][3]);

///Estimate the exponential of a 3*3 matrix
void Exponential_3t3matrix(float m1[3][3],float m2[3][3]);

///transpose a 3*3 matrix
void transpose_3t3matrix(float m1[3][3],float m2[3][3]);

///inverse of a 3*3 matrix
void invert_3t3matrix(float m1[3][3],float m2[3][3]);

///multiply two 3*3 matrices
void mult_3t3mat_3t3mat(float m1[3][3], float m2[3][3], float MatRes[3][3]);

///multiply a vector of size 3 and a 3*3 matrix
void mult_3t3mat_3vec(float mat[3][3], float vectIni[3], float vectRes[3]);

///In a 3D domain, project a point "Pt2proj" to a line defined by a point "LinePt" and a vector "LineVec". Result is saved in "ProjPt"
float Project_point_to_a_line(float Pt2proj[3],float LinePt[3],float LineVec[3],float ProjPt[3]);

///scalar product between two vectors of size 3
float scalarProd_3vec(float v1[3],float v2[3]);

///cross product between two vectors of size 3
void crossProd_3vec(float v1[3],float v2[3],float OutV[3]);

///project a vector 'vec' to a plan defined by the vectors (Pvec1,Pvec2)
void Project_3vec_plan(float Pvec1[3],float Pvec2[3],float vec[3]);

///inverse of a quaternion
void invert_4t4quaternion(float q1[4][4],float q2[4][4]);

///multiply a vector of size 4 and a 4*4 matrix
void mult_4t4mat_4vec(float mat[4][4], float vectIni[4], float vectRes[4]);

///multiply two 4*4 matrix representing quaternions: mat_i1 * mat_i2 -> mat_o
void mult_quat4t4mat_quat4t4mat(float mat_i1[4][4], float mat_i2[4][4], float mat_o[4][4]);

///read a 4*4 matrix in a text file
void Read_quat4t4mat(char *,float locmat[4][4]);

///write a 4*4 matrix in a text file
void Write_quat4t4mat(char *,float locmat[4][4]);


/// 7.2) ++++++++++++++++++ Square matrices ++++++++++++++++++ 

class SquareMatrix{
private:

  int CopiedByRef;
  
public:
  //contains the matrix
  float * LocMatrix;
  
    //matrix size
  int NI;
  int NJ;
  
  /// Constructors
  SquareMatrix();
  
  /// Destructor
  ~SquareMatrix();
  
  //read and write the matrix in CSV files with space delimiters
  virtual void Read(char *);
  virtual void Write(char *);

  //allocate memory for a n*n matrix
  virtual void CreateVoidMatrix(int n); 
  
  //put all values to zero or make an eye matrix
  virtual void SetToZero();
  virtual void Eye();
  
  //show matrix
  virtual void Show();
  
  //copy by reference another square matrix
  virtual void CopyByRef(SquareMatrix * copiedMatrix);
  
  //copy another square matrix.
  virtual void Copy(SquareMatrix * copiedMatrix);
  
  //get a matrix value
  virtual inline float G(int i,int j){
    return this->LocMatrix[j+this->NI*i];
  }
  
  //put/add a matrix value
  virtual inline void P(float value,int i,int j){
    this->LocMatrix[j+this->NI*i]=value;
  }
  
  virtual inline void Add(float value,int i,int j){
    this->LocMatrix[j+this->NI*i]+=value;
  }
  
  //transpose matrix
  virtual void Transpose();

  //symmetrize a matrix -> Mat <- (Mat + Mat')/2
  virtual void symmetrize();
  
  //Put/Add the outer-product of two vectors in the square matrix
  //The size of the vector must be coherent with the matrix
  virtual void PutOuterProduct(float * v1,float * v2);
  virtual void AddOuterProduct(float * v1,float * v2);
  
  //Matrix/Matrix multiplication and sum: (this M2)   and  (this + M2)
  //The result is stored in this
  virtual void Mult(SquareMatrix * M2);
  virtual void Add(SquareMatrix * M2);
  
  //weighted sum of two matrices: this <- (weight_this * this) + (weight_M2 * M2)
  virtual void WeightedSum(float weight_this,SquareMatrix * M2,float weight_M2);
  
  //weighted sum of this with identity matrix: this <- (weight_this * this) + (weight_M2 * Id)
  virtual void WeightedSumWithId(float weight_this,float weight_M2);
  
  //multiplication by a scalar / addition with a scalar
  virtual void Mult(float coefMult);
  virtual void Add(float b);
  
  //Matrix/Vector multiplication. Return the product in V
  virtual void Mult(float * V);
  
  //sum and max of the values or absolute values
  virtual float SumValues();
  virtual float MaxValue();
  virtual float SumAbsValues();
  virtual float MaxAbsValue();
  
  //normalize the matrix so that the average sum of the columns is 1
  virtual void MatrixNormalization();

  //+=, -=, /=, and *= operators
  virtual inline void PlusEq(int i,int j,float value){ this->LocMatrix[j+this->NI*i]+=value;}  //same as Add except the order of the inputs
  virtual inline void MinusEq(int i,int j,float value){ this->LocMatrix[j+this->NI*i]-=value;}
  virtual inline void TimesEq(int i,int j,float value){ this->LocMatrix[j+this->NI*i]*=value;}
  virtual inline void DivEq(int i,int j,float value){ this->LocMatrix[j+this->NI*i]/=value;}
 
};


///Perform a Housholder reduction of a real symmetric matrix z[][]
///Outputs:
///  -> z[][] is replaced by the orthogonal matrix effecting the transformation. 
///  -> d[] returns the diagonal elements of the tri-diagonal matrix
///  -> e[] the off-diagonal elements with e[0] = 0.
///
///Adapted from www.physics.sdsu.edu/~johnson/phys580/tred2.c which was very slightly adapted from tred2() in numerical recipies in C++ 3rd edition (differences where checked)
///-> Validated on small matrices - results differ from those obtained using the hess function in matlab but make sense when reconstructing the original matrix
void tred2(SquareMatrix * z, float *d, float *e);


///The function tqli() determine eigenvalues and eigenvectors of a real 
///symmetric tri-diagonal matrix, or a real, symmetric matrix previously
///reduced by function tred2[] to tri-diagonal form. 
///Input: 
///  d[] contains the diagonal element 
///  e[] the sub-diagonal of the tri-diagonal matrix. 
///  z[][] is the matrix output from tred2() or identity otherwise
///Outputs:
///  d[] contains the eigenvalues 
///  e[] is destroyed
///  z[][] its k'th column is the normalized eigenvector corresponding to d[k]. 
///
///Adapted from www.physics.sdsu.edu/~johnson/phys580/tqli.c which was very slightly adapted from tred2() in numerical recipies in C++ 3rd edition (differences where checked)
///-> Validated on small matrices - results differ from those obtained using pca in matlab but make sense when reconstructing the original matrix
void tqli(float *d, float *e, SquareMatrix * z);


///Construct the QR decomposition of a. The upper triangluar matrix R and the transpose of the orthogonal matrix Q are stored.
///1 is returned if a singularity is encountered during the decomposition, but the decomposition is still completed in this case; 0 is returned otherwise
///Adapted from the algorithm of numerical recipes in c++ - third edition
int QRdecomp(SquareMatrix * a,SquareMatrix * qt, SquareMatrix * r);

///rank eigenvalues and eigenvetors with decreasing eigenvalues 
void rankEigenvalues(float *eigenvalues, SquareMatrix * eigenvectors);

///show recomposed matrix after a PCA
void ShowRecomposedMatrixAfterPCA(float *eigenvalues, SquareMatrix * eigenvectors);

///perform the PCA of a reasonably large SYMMETRIC matrix
void LargeMatrixPCA(SquareMatrix * eigenvectors, float *eigenvalues);


/// 7.3) ++++++++++++++++++ Other scientific computation stuffs ++++++++++++++++++

///minmod function
float minmodfunc(float a,float b);

///computes sqrt(a^2+b^2) without destructive underflow or overflow
///From numerical recipies (2nd edition)
float pythag(float a, float b);


///inline functions from numerical recipies (3rd edition)
inline float nr_SIGN(const float &a, const double &b){
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
  }

inline float nr_MAX(const float &a, const float &b){
  return b > a ? (b) : (a);
  }

///Generate a random number following a normal distribution  (mean=0 std=1)
float sampleNormal();

///Generate a random number following an uniform distribution  (in [0,1])
float sampleUniform();


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                        8: LANDMARKS (points and curves)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// 8.1) ++++++++++++++++++ point landmarks ++++++++++++++++++ 
class LDMK_Points{
  /// ******************************************************************************
private:
  //number of LDMK_Points
  int LDMK_Points_Nb;
  
  //landmark coordinates
  float * Lx;
  float * Ly;
  float * Lz;
  float * Val;
  
  /// ******************************************************************************
public:
  
  /// Constructor and destructor
  LDMK_Points();
  ~LDMK_Points();
  
  ///read LDMK_Points in a CSV file
  ///If (withValues!=1) the file should have this format:
  ///[Point 1: x]\t[Point 1: y]\t[Point 1:  z]
  ///[Point 2: x]\t[Point 2: y]\t[Point 2:  z]
  ///...
  ///or, if (withValues==1):
  ///[Point 1: x]\t[Point 1: y]\t[Point 1:  z]\t[Point 1:  Val]
  ///[Point 2: x]\t[Point 2: y]\t[Point 2:  z]\t[Point 2:  Val]
  ///...
  virtual void Read(char *,int withValues=0);
  
  ///read the LDMK_Points as the non-null values of a scalar field
  virtual void ReadInScalarField(ScalarField * img3d);
  
  ///Write LDMK_Points in a CSV file
  virtual void Write(char *);
  
  ///Show LDMK_Points
  virtual void Show();
  
  ///translate the points
  virtual void Translate(float DecX,float DecY,float DecZ);

  ///transform the LDMK_Points coordinates and diameters from voxels to mm according to the image 2 world properties of RefSF 
  virtual void VoxelsToMillimeters(ScalarField * RefSF);

  ///transform the LDMK_Points coordinates and diameters from mm to voxels according to the image 2 world properties of RefSF 
  virtual void MillimetersToVoxels(ScalarField * RefSF);

  ///Get the X, Y, Z coordinates of the LDMK_Points
  virtual float GetX(int Id=0);
  virtual float GetY(int Id=0);
  virtual float GetZ(int Id=0);
  virtual float GetVal(int Id=0);
  
  ///Get the number of LDMK_Points
  virtual int Get_LDMK_PointsNumber(void);
};
  
  
  
  
  
/// 8.2) ++++++++++++++++++ curve landmarks ++++++++++++++++++ 
class LDMK_Curves{
/// ******************************************************************************
private:
  int NbSeg;       //segments number (amount of curves)
  int *NbEl;       //elements number in each segment (amount of nodes in each segment). NbEl[i]=0 means no memory is allocated at segment i
  float **x;      //x[i][j] coordinate of segment j / element i on the x axis
  float **y;      //y[i][j] coordinate of segment j / element i on the y axis
  float **z;      //z[i][j] coordinate of segment j / element i on the z axis
  float **d;      //d[i][j] diametre of segment j / element i
  
  //functions used for skeletonisation
  virtual void NgbhTransform(ScalarField * CubeIn,ScalarField * CubeTemp,ScalarField * CubeOut,int opt1,int opt2);
  virtual int IsNotCenterlinePt(ScalarField * CubeIn,int FilterID);
  
  
/// ******************************************************************************
public:
  
  /// Constructors and destructor
  LDMK_Curves();
  LDMK_Curves(int SegNb, int ElNb);
  ~LDMK_Curves();
  
  ///read LDMK_Curves in a mv3d file
  virtual void Read(char *);
  
  ///read a LDMK_Curves in a raw ascii file -> The file only contains the coordinates
  ///If SubSamplingFactor is higher than 1, one point of RawFileName is loaded every 'SubSamplingFactor' points
  virtual void ReadInRawFile(char * RawFileName,int SubSamplingFactor=1);
  
  ///generate the LDMK_Curves by skeletonizing the img3d (shape=1 / backgroud=0)
  ///The 18-neighborhood is considered for the shape and the 6-neighborhood is considered for the background
  ///Remark 1: Radii represent the nearest boundary in mm
  ///Remark 2: Coordinates are in mm
  ///Remark 3: Post-treatments are usually necessary to merge the segments (often subdivided into several segments) and useless segments
  virtual void Skeletonize(ScalarField * img3d);
  
  ///Reconstruct the 3D volume form the network (shape=1 / backgroud=0)
  ///warning: all intensities of RefImageDomain will be modified so that they are the reconstructed volume in the end of the computations
  virtual void Generate3DVolume(char * OutputImageName,ScalarField * RefImageDomain);
  
  ///create a LDMK_Curves containing one curve with NbElem elements
  virtual void Create_LDMK_1_Curve(int NbElem);
  
  ///create a LDMK_Curves containing NbCurves void curves (no elements)
  virtual void Create_LDMK_N_VoidCurves(int N);
  
  ///write LDMK_Curves in a mv3d file
  ///Set Preserve_IDs to 1 to preserve the original segment identifiers (default). They are optimaly resampled otherwise.
  virtual void Write(char *, int Preserve_IDs=1);

  ///Export LDMK_Curves in a vtk file
  virtual void ExportAsVtk(char *,int ShowVolume=1);

  ///Get the coordinate value of x, y, z or get the diameter of the element 'IdEl' of segment 'IdSeg'
  virtual inline float GetX(int IdSeg, int IdEl){return x[IdSeg][IdEl];}
  virtual inline float GetY(int IdSeg, int IdEl){return y[IdSeg][IdEl];}
  virtual inline float GetZ(int IdSeg, int IdEl){return z[IdSeg][IdEl];}
  virtual inline float GetD(int IdSeg, int IdEl){return d[IdSeg][IdEl];}
  
  ///Put a coordinate value for x, y, z or put a diameter to the element 'IdEl' of segment 'IdSeg'
  virtual inline void PutX(float value, int IdSeg, int IdEl){x[IdSeg][IdEl]=value;}
  virtual inline void PutY(float value, int IdSeg, int IdEl){y[IdSeg][IdEl]=value;}
  virtual inline void PutZ(float value, int IdSeg, int IdEl){z[IdSeg][IdEl]=value;}
  virtual inline void PutD(float value, int IdSeg, int IdEl){d[IdSeg][IdEl]=value;}
  
  ///get the number of segments
  virtual inline int GetSegNumber(){return NbSeg;}
  
  ///get the number of elements in segment 'IdSeg'
  virtual inline int GetElNumber(int IdSeg){return NbEl[IdSeg];}
  
  ///Reduce the number of elements in a segment - consider the elements between FirstEl and LastEl (included)
  virtual void ReduceElNumber(int IdSeg, int FirstEl, int LastEl);
  
  ///Merge two segments at their nearest extremity. The new segment is saved in Seg1.
  virtual void MergeSegments(int Seg1, int Seg2);
  
  ///Resample the number of elements in a segment. The resampled segment will be (roughly) homogeneously resamped in space.
  virtual void ResampleSegment(int IdSeg, int NewElNb);

  ///Delete a segment
  virtual void DeleteSegment(int Seg);

  ///Generate a segment of size 'NbElements' with only null values
  ///The ID of the segment in 'this' is returned
  virtual int GenerateVoidSegment(int NbElements);

  ///count the number of segments related to the segment end SegEnd (= 0 pr 1) of the segment Seg
  virtual int CountLinkedSeg(int Seg,int SegEnd,float epsilon=0.0001);

  ///Clean-up a moderately large LDMK_Curves structure
  /// -> two segments ends are supposed linked if their distance is less than epsilon
  /// -> merge the segments linked by a node with only two segments
  /// -> remove the segments linked to only one node and for which the node has more than two segments related to other nodes
  virtual void CleanUp(float epsilon=0.0001);
  
  /// -> Smooth the segments
  virtual void Smooth(int ItNb);
  virtual void SmoothOneSegment(int SegmentID,int ItNb);

  ///transform the LDMK_Curves coordinates and diameters from voxels to mm according to the image 2 world properties of RefSF 
  virtual void VoxelsToMillimeters(ScalarField * RefSF);

  ///transform the LDMK_Curves coordinates and diameters from mm to voxels according to the image 2 world properties of RefSF 
  virtual void MillimetersToVoxels(ScalarField * RefSF);

  ///Affine transformation (TransfoMat is from the source/template to the target)
  virtual void AffineTransfo(float TransfoMat[4][4]);

  ///translate a network
  virtual void Translate(float DecX,float DecY,float DecZ);

  ///estimate network ROI
  virtual void EstimateROI(float * LowerX,float * LowerY,float * LowerZ,float * WidthX,float * WidthY,float * WidthZ);
  
};



///return a weight for BSpline interpolation. x is in [-1.5,1.5] and the weights are centered in 0   
///-> it will then be necessary to use the values at pts: [pt-2] [pt-1] [pt] [pt+1] [pt+2] for interpolation
float BSplineWeight(float x);

///return the value of a curve 'densified' using B-spline interpolation at a given location
///location should have a value in [0 , 1]
float BSplineCurveSampler(float * Values, int NbValues,float location);
  
  

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// 9: VectorField convolver taking into accound the whole image domain (projected on a B-spline basis)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//BSplineBasis_of_VecField
class BS_GlobalConvolver{
  /// ******************************************************************************
  private:
    VectorField * RefVecField;               //reference vector field  (defined on all points)
    int StepSize;                            //step between the control points in voxels  (separated steps could be used on X, Y and Z if the voxel size is anisotropic)
    int NbCtPtsX,NbCtPtsY,NbCtPtsZ,NbCtPtsT; // number of control points on the X, Y, Z and T axis.
    FFTconvolver3D SplineConvolver;          //used to perform quick projection of the VF on the basis
    
    SquareMatrix NodesWeigther_MatrixStyle;  //matrix M used to weight the terms of the basis with a matrix-vector multiplication
    int NodesWeigther_NbNodes;               //nb of nodes in the basis (M is of size (3*NodesWeigther_NbNodes)*(3*NodesWeigther_NbNodes))
    int NodesWeigther_Type;                  //If NodesWeigther_Type==1 then matrix style / If NodesWeigther_Type==2 then vector style / no weighting action otherwise
    
    VectorField ResidualField;
    
    VectorField TmpVecAtControlPts;
    VectorField VecAtControlPts;

    ///add the values of other control points to the local ones
    virtual void  Add_VecAtControlPts(BS_GlobalConvolver * OtherBasis);
    
    ///multiply the control points by the coef MultCoef
    virtual void  MultCoef(float MultCoef);
    
    ///get the value at a node of the grid
    virtual inline float G(int IdDirec,int x,int y,int z=0,int t=0){
      return this->VecAtControlPts.G(IdDirec,x,y,z,t);
    }

    ///link a bspline basis with vector fields and define the grid step size
    virtual void  LinkWithVecField(VectorField * LinkedVF,int GridStepSize);
    
    ///load a matrix which will weight the terms of the basis with a matrix-vector multiplication
    ///Weigths of the basis are considered as a vector with (first x) / (then y) / (then z) / (then t) / (then direction)
    ///Basis terms are weighted when using 'WeightBasisTerms'
    virtual void Load_MatrixNodesWeigther(char * Mfile); 
    
    ///same as Load_MatrixNodesWeigther but the matrix does not come from a file but is copied by reference
    virtual void Load_MatrixNodesWeigther_usingCopyByRef(SquareMatrix * RefMatrix, int Verbose=0); 
    
    ///compress the information of RefVectorField at the control points VecAtControlPts
    virtual void  RefVecField_2_VecAtControlPts();

    ///uncompress the information at the control points VecAtControlPts in RefVectorField 
    virtual void  VecAtControlPts_2_RefVecField();
    
    ///weight the terms of the basis with the strategy predefined in 'Load_WeightingMatrix' or 'Load_WeightingVector'
    ///Nothing is made if no strategy is defined
    virtual void NodesWeight();
    
  /// ******************************************************************************
  public:
    
    /// Constructor
    BS_GlobalConvolver();
  
    /// Destructor
    ~BS_GlobalConvolver();
    
    
    ///initiate the convolver
    ///-> link a bspline basis with the VECTOR FIELD 'LinkedVF' which will be smoothed by the convolver
    ///-> define the grid step size of a bspline basis on which the VF will be projected
    ///-> Either copy by reference or load in an ascii file the square matrix (linear operator) which smoothes 
    ///     the weights of the basis nodes (must be coherent with the VF and grid step size). If no matrix
    ///     is defined, the convolver is equivalent to identity.
    virtual void  InitiateConvolver(VectorField * LinkedVF,int GridStepSize,char * Mfile);
    virtual void  InitiateConvolver(VectorField * LinkedVF,int GridStepSize,SquareMatrix * RefMatrix);
    virtual void  InitiateConvolver(VectorField * LinkedVF,int GridStepSize);
    
    ///perform the global smoothing
    virtual void  SmoothLinkedVF();
    
    ///Write an identity template matrix which could be used with 'Load_WeightingMatrix'
    virtual void WriteTemplateNodesWeighterIdMat(VectorField * LinkedVF,int GridStepSize,char * Mfile); 

    ///Write the weights on the control points
    virtual void WriteControlPointWeights(char * FileNameX,char * FileNameY,char * FileNameZ); 

    
};



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                     10: SURFACE MESHES (with triangles)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// 10.1) ++++++++++++++++ Class nodes list  (simplified compared to UtilzGraph) ++++++++++++++++ 

class NodesList{
private:
  /// parameters
  int * Node;       //Identifier of a node
  float * Weight;       //Corresponding weight
  int NbNodes;      //number of nodes in the list
  int WithWeight;
  
protected:
  
public:
  
  /// Constructor
  NodesList();
  
  /// Destructor
  ~NodesList();
  
  ///allocate memory for the list
  ///size; size of the list
  ///with_weight: allocate memory for weights if ==1 (in addition to the node identifiers)
  virtual void allocate(int size,int with_weight=0);
  
  ///Get and put functions
  virtual inline int GetNbNodes();
  virtual inline void ReduceNbNodes(int size);
  virtual inline int GetNode(int entry);
  virtual inline void PutNode(int entry, int ID_node);
  virtual inline float GetWeight(int entry);
  virtual inline void PutWeight(int entry, float locWeight);

};

/// 10.2) ++++++++++++++++ Class Graph (simplified compared to UtilzGraph) +++++++++++++++++


class Graph{
private:
  
protected:
  
public:
  
  /// Constructor
  Graph();
  
  /// Destructor
  ~Graph();
  
  /// public functions
  
  
  /// public Parameters
  int NbNodes;       //number of nodes in the graph  (node have ID from 0 to NbNodes-1)
  int NbEdges;       //number of directed edges in the graph
  int isSymmetric;   //==1 if the directed graph is symmetric / ==0 otherwise
  NodesList * NgbhsOfNode;   //for each node i in the graph, NgbhsOfNode[i] contains the list of neighbor nodes and the weights in the corresponding edges
};


/// 10.3) ++++++++++++++++ Class ValuedGraph +++++++++++++++++
/// -> simplified compared to UtilzGraph but with AffineTransfo and NonRigidTransfo

class ValuedGraph{
private:
  int NbNodes;
  Graph Connections;
  float * NodeValues;
  
  int KnownNodeCoordinates;
  float *  NodeCoordX;
  float *  NodeCoordY;
  float *  NodeCoordZ;
  
  int KnownTriangularFaces;
  int ** TriangularFaces;
  int NbTriangularFaces;
  
protected:
  
public:
  
  /// Constructor
  ValuedGraph();
  
  /// Destructor
  ~ValuedGraph();
  
  /// public functions
  
  ///Read the graph in an ascci POLYDATA vtk file containing a mesh with values at each node. 
  ///The imported mesh is a graph with weigths which are the distances between the nodes
  ///Returns 1 if the vtk file seems to be read normally and 0 otherwise.
  virtual int ReadMeshInVtkFile(char *);
  
  ///Save the graph as a mesh in an ascci POLYDATA vtk file. 
  ///The graph will only be saved if 'KnownNodeCoordinates' and 'KnownTriangularFaces' are equal to 1.
  ///Returns 1 if the graph is saved and 0 otherwise.
  virtual int WriteMeshInVtkFile(char *);
  
  ///Affine transformation of the mesh coordinates  (if KnownNodeCoordinates==1)
  ///DefQuaternion is a 4*4 matrix: [R_xx R_xy R_xz T_x ; R_yx R_yy R_yz T_y ; R_zx R_zy R_zz T_z ; 0 0 0 1]
  virtual void AffineTransfo(float DefQuaternion[4][4]);
  
  ///Non rigid deformation with a DisplacementField
  virtual void NonRigidTransfo(VectorField * DisplacementField);

};


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                     11: EMPIRICAL MODE DECOMPOSITION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///Compute the local minima and maxima from the image
void ExtractLocalExtrema(ScalarField * InputImag,ScalarField * TmpImag,ScalarField * LocalMax, ScalarField * LocalMin);

///Compute the local minima and maxima from the image
void ExtractLocalExtrema2(ScalarField * InputImag,ScalarField * LocalMax, ScalarField * LocalMin);

///Add simulated Local Extrema to a map of Local extrema to numerically stabilize the normalized convolution
void AddSimulatedLocalExtrema(ScalarField * LocalExtrema,ScalarField * TmpImag,ScalarField * PointsToAdd);

///Estimate a good standard deviation of the Jacobian to perform normalized convolution of the input local extra map
float EstimateGoodSigma(ScalarField * LocalMax, ScalarField * LocalMin);

///Compute the normalized convolution of local minima/maxima map
void NormalizedConvolution(ScalarField * LocalExtrema,ScalarField * InputImag,LightFFTconvolver3D * TmpConvolver,ScalarField * TmpImag,float sigma,ScalarField * OutputConvol);

///Estimate the IMF from the input image
void EstimateIMF(ScalarField * InputImag,int NbIt,ScalarField * TmpImag1,ScalarField * LocalMax,ScalarField * LocalMin,ScalarField * TmpImag4,LightFFTconvolver3D * TmpConvolver,ScalarField * OutputIMF);

///Perform an EMD
void PerformEMD(char InputImagFile[256],char OutputPrefix[256]);




///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                     12: metamorphoses
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///update CurrentMetamo using a pseudo-gradient method to manage the non-differentiability at 0 -- here the regularisation is everywhere the same (lambdaMetamo)
void PseudoGradMetamoUpdate(ScalarField * CurrentMetamo,ScalarField * UpdateForcesMetamo,float lambdaMetamo,float refMultFactor=0.1,float ThreshZero=0.01);

///update CurrentMetamo using a pseudo-gradient method to manage the non-differentiability at 0 -- here the regularisation depends on space
void PseudoGradMetamoUpdate(ScalarField * CurrentMetamo,ScalarField * UpdateForcesMetamo,ScalarField *  lambdaMetamoField,float refMultFactor=0.1,float ThreshZero=0.01);

void CptMetamoReg1stOrderDerivatives(ScalarField * CurrentMetamo,ScalarField *  lambdaMetamoField,ScalarField * tempFl,ScalarField * tempFl2,float locLambda);






#endif
