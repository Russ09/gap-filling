/*=========================================================================
 
 Author: Laurent Risser, Francois-Xavier Vialard
 
 Disclaimer: This software has been developed for research purposes only, and hence should 
 not be used as a diagnostic tool. In no event shall the authors or distributors
 be liable to any direct, indirect, special, incidental, or consequential 
 damages arising of the use of this software, its documentation, or any 
 derivatives thereof, even if the authors have been advised of the possibility 
 of such damage. 
 
 
 =========================================================================*/

#include <SciCalcPack.h>



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           1:   FUNCTIONS FOR THE CLASS "ScalarField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///1.1) conventional scalar field

///constructor
ScalarField::ScalarField(void){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  this->NT=0;
}

///destructor
ScalarField::~ScalarField(void){
  if ((this->ScalField!=NULL)&&(this->NX>0)) delete this->ScalField;
}


///put a value
//-> inline function in the .h file

/// add a value
//-> inline function in the .h file

///put a the same value at every points of the scalar field
void ScalarField::PutToAllVoxels(float cste,int t)
{
  int x,y,z;
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) { this->P(cste,x,y,z,t); }
}

///get a value
//-> inline function in the .h file

///get a value using linear interpolation
float ScalarField::G(float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  //values out of the image
  if (x<0.) x=0.0001;
  if (x>=this->NX-1.) x=this->NX-1.0001;
  if (y<0.) y=0.0001;
  if (y>=this->NY-1.) y=this->NY-1.0001;
  if (z<0.) z=0.0001;
  if (z>=this->NZ-1.) z=this->NZ-1.0001;
  if (t<0) t=0;
  if (t>this->NT-1) t=this->NT-1;
  
  //closest entire value
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  //interpolation
  if (this->NZ==1){ //2D IMAGE
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
    
    InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
  }
  else{//3D IMAGE
    wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
    wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
    
    InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
  }
  
  return InterpoGreyLevel;
}


///same as above
float ScalarField::G(double x,double y,double z,int t){
  return this->G((float)x,(float)y,(float)z,t);
}




///get a value using linear interpolation or nearest neigbhour
float ScalarField::G(float coordSpace2imageSpace[4][4],float x,float y,float z,int t,int NN){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xt,yt,zt;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  
  //transform x,y,z to the image space
  xt=x;
  yt=y;
  zt=z;
  
  x=xt*coordSpace2imageSpace[0][0]+yt*coordSpace2imageSpace[0][1]+zt*coordSpace2imageSpace[0][2]+coordSpace2imageSpace[0][3];
  y=xt*coordSpace2imageSpace[1][0]+yt*coordSpace2imageSpace[1][1]+zt*coordSpace2imageSpace[1][2]+coordSpace2imageSpace[1][3];
  z=xt*coordSpace2imageSpace[2][0]+yt*coordSpace2imageSpace[2][1]+zt*coordSpace2imageSpace[2][2]+coordSpace2imageSpace[2][3];
  
  if (NN==0){ //A) trilinear interpolation
    //values out of the image
    if (x<0.) x=0.0001;
    if (x>=this->NX-1.) x=this->NX-1.0001;
    if (y<0.) y=0.0001;
    if (y>=this->NY-1.) y=this->NY-1.0001;
    if (z<0.) z=0.0001;
    if (z>=this->NZ-1.) z=this->NZ-1.0001;
    if (t<0) t=0;
    if (t>this->NT-1) t=this->NT-1;
    
    //closest entire value
    xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
    yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
    zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
    
    //interpolation
    if (this->NZ==1){ //2D IMAGE
      wmm=xwm*ywm;
      wmp=xwm*ywp;
      wpm=xwp*ywm;
      wpp=xwp*ywp;
      
      InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
    }
    else{//3D IMAGE
      wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
      wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
      
      InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    }
  }
  else{//B) nearest neighbor interpolation
      //values out of the image
      xi=static_cast<int>(x+0.5);
      yi=static_cast<int>(y+0.5);
      zi=static_cast<int>(z+0.5);
      
      if (xi<0) xi=0;  if (xi>=this->NX) xi=this->NX-1;  
      if (yi<0) yi=0;  if (yi>=this->NY) yi=this->NY-1;  
      if (zi<0) zi=0;  if (zi>=this->NZ) zi=this->NZ-1;  
      
      InterpoGreyLevel=this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
  }
  
  return InterpoGreyLevel;
}


///same as above
float ScalarField::G(float coordSpace2imageSpace[4][4],double x,double y,double z,int t,int NN){
  return this->G(coordSpace2imageSpace,(float)x,(float)y,(float)z,t,NN);
}

///same as above
float ScalarField::G(float coordSpace2imageSpace[4][4],int x,int y,int z,int t,int NN){
  return this->G(coordSpace2imageSpace,(float)x,(float)y,(float)z,t,NN);
}





///get a value using linear interpolation -- here 0 is returned if we're out of the image
float ScalarField::G_NoExtrapo(float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  //values out of the image
  if (x<0.) return 0;
  if (x>=this->NX-1.) return 0;
  if (y<0.) return 0;
  if (y>=this->NY-1.) return 0;
  if (z<0.) return 0;
  if (z>=this->NZ-1.) return 0;
  if (t<0) return 0;
  if (t>this->NT-1) return 0;

  
  //closest entire value
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  //interpolation
  if (this->NZ==1){ //2D IMAGE
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
    
    InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
  }
  else{//3D IMAGE
    wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
    wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
    
    InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
  }
  
  return InterpoGreyLevel;
}




///get a value using linear interpolation or nearest neigbhour -- here 0 is returned if we're out of the image
float ScalarField::G_NoExtrapo(float coordSpace2imageSpace[4][4],float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xt,yt,zt;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  
  //transform x,y,z to the image space
  xt=x;
  yt=y;
  zt=z;
  
  x=xt*coordSpace2imageSpace[0][0]+yt*coordSpace2imageSpace[0][1]+zt*coordSpace2imageSpace[0][2]+coordSpace2imageSpace[0][3];
  y=xt*coordSpace2imageSpace[1][0]+yt*coordSpace2imageSpace[1][1]+zt*coordSpace2imageSpace[1][2]+coordSpace2imageSpace[1][3];
  z=xt*coordSpace2imageSpace[2][0]+yt*coordSpace2imageSpace[2][1]+zt*coordSpace2imageSpace[2][2]+coordSpace2imageSpace[2][3];
  
    //values out of the image
    if (x<0.) return 0;
    if (x>=this->NX-1.) return 0;
    if (y<0.) return 0;
    if (y>=this->NY-1.) return 0;
    if (z<0.) return 0;
    if (z>=this->NZ-1.) return 0;
    if (t<0) return 0;
    if (t>this->NT-1) return 0;
    
    //closest entire value
    xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
    yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
    zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
    
    //interpolation
    if (this->NZ==1){ //2D IMAGE
      wmm=xwm*ywm;
      wmp=xwm*ywp;
      wpm=xwp*ywm;
      wpp=xwp*ywp;
      
      InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
    }
    else{//3D IMAGE
      wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
      wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
      
      InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    }
  
  return InterpoGreyLevel;
}




///get the maximum absolute values out of the scalar field
float ScalarField::GetMaxAbsVal(int t){
  float max=0.0;
  int x,y,z;
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
  {
    if(max<abs(this->G(x,y,z,t))){max = abs(this->G(x,y,z,t));}
  }
  return max;
}


///read a scalar field (in a nifti image)
void ScalarField::Read(char * ImageName){
  nifti_1_header hdr;
  FILE *fp;
  int ret,i;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z,t;
  float a,b,c,d,qfac;

/*BEGIN TEST*/
//nifti_image * toto;
//toto = new nifti_image [1];
//toto=nifti_image_read( ImageName,1);
/*END TEST*/


  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  //to open some multichannel 3D images which have the channels in the 5th dimension instead of the 4th dimension
  if ((hdr.dim[4]==1)&&(hdr.dim[5]>1)){
    hdr.dim[4]=hdr.dim[5];
    hdr.dim[5]=1;
  }
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=hdr.dim[1])||(this->NY!=hdr.dim[2])||(this->NZ!=hdr.dim[3])||(this->NT!=hdr.dim[4]))
      cout << "The size of a scalar field is changed\n";
  
  
  //fill the parameters of the class
  this->NX=static_cast<int>(hdr.dim[1]);
  this->NY=static_cast<int>(hdr.dim[2]);
  this->NZ=static_cast<int>(hdr.dim[3]);
  this->NT=static_cast<int>(hdr.dim[4]);
  if (this->NT<1) {this->NT=1; hdr.dim[4]=1;}
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;

  //cout << "Size of " << ImageName << ": " << this->NX << " " <<  this->NY << " " <<  this->NZ << endl;
  
  //Image to world matrix
  if (hdr.sform_code>0){ 
    //METHOD 3 of nifti1.h
    //cout << "Orientation of image " << ImageName << " opened using method 3 (alternative - normal)" << endl;
    this->Image2World[0][0]=hdr.srow_x[0];  this->Image2World[0][1]=hdr.srow_x[1];  this->Image2World[0][2]=hdr.srow_x[2];  this->Image2World[0][3]=hdr.srow_x[3];
    this->Image2World[1][0]=hdr.srow_y[0];  this->Image2World[1][1]=hdr.srow_y[1];  this->Image2World[1][2]=hdr.srow_y[2];  this->Image2World[1][3]=hdr.srow_y[3];
    this->Image2World[2][0]=hdr.srow_z[0];  this->Image2World[2][1]=hdr.srow_z[1];  this->Image2World[2][2]=hdr.srow_z[2];  this->Image2World[2][3]=hdr.srow_z[3];
    this->Image2World[3][0]=0;              this->Image2World[3][1]=0;              this->Image2World[3][2]=0;              this->Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      //cout << "Orientation of image " << ImageName << " opened using method 2 (normal)" << endl;
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      this->Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); this->Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       this->Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     this->Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   this->Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     this->Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       this->Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;                               this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                    this->Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      //set the dixel dimensions to 1 if they are absolutely unrealistic 
      //if ((hdr.pixdim[1]<0.000000001)&&(hdr.pixdim[1]>1000000000)) {hdr.pixdim[1]=1;  hdr.qoffset_x=0;}
      //if ((hdr.pixdim[2]<0.000000001)&&(hdr.pixdim[2]>1000000000)) {hdr.pixdim[2]=1;  hdr.qoffset_y=0;}
      //if ((hdr.pixdim[3]<0.000000001)&&(hdr.pixdim[3]>1000000000)) {hdr.pixdim[3]=1;  hdr.qoffset_z=0;}
      
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      this->Image2World[0][0]=hdr.pixdim[1];  this->Image2World[0][1]=0;             this->Image2World[0][2]=0;             this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=0;              this->Image2World[1][1]=hdr.pixdim[2]; this->Image2World[1][2]=0;             this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=0;              this->Image2World[2][1]=0;             this->Image2World[2][2]=hdr.pixdim[3]; this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;              this->Image2World[3][1]=0;             this->Image2World[3][2]=0;             this->Image2World[3][3]=1;
    }
  }

//  cout << "Image to world matrix of " <<  ImageName << ":" << endl;
//  int j;
//  for (i=0;i<4;i++){
//    for (j=0;j<4;j++){
//      cout << this->Image2World[i][j] << " ";
//    }
//    cout << endl;
//  }
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
//  cout << "World to image matrix of " <<  ImageName << ":" << endl;
//  for (i=0;i<4;i++){
//    for (j=0;j<4;j++){
//      cout << this->World2Image[i][j] << " ";
//    }
//    cout << endl;
//  }
  
  //print a little header information
  //fprintf(stderr, "\n%s header information:",ImageName);
  //fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
  //fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
  //fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
  //fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
  //fprintf(stderr, "\n");
  
  //message about the type of grey levels
//   cout << endl;
//   cout << ImageName << " contains";
//   if (hdr.datatype==2) cout << " unsigned char ";
//   if (hdr.datatype==4) cout << " signed short ";
//   if (hdr.datatype==8) cout << " signed int ";
//   if (hdr.datatype==16) cout << " float ";
//   if (hdr.datatype==64) cout << " double ";
//   if (hdr.datatype==256) cout << " signed char ";
//   if (hdr.datatype==512) cout << " unsigned short ";
//   if (hdr.datatype==768) cout << " unsigned int ";
//   cout << "pixels and has a resolution of " << hdr.dim[1] << "*"  << hdr.dim[2] << "*"  << hdr.dim[3] << "*"  << hdr.dim[4] << " voxels" << endl;
//   cout << endl;

  //2) read the image
  
  //allocate the memory for the image
  this->ScalField = new float [this->NXtYtZ*this->NT];

  // jump to data offset
  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
  
  //test the multiplicatory value of hdr
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    //cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << ImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }

  //load the image
  if (hdr.datatype==2){ //2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data2, sizeof(unsigned char), 1, fp);
      this->P(static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==4){ //4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data4, sizeof(signed short), 1, fp);
      this->P(static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==8){  //8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data8, sizeof(signed int), 1, fp);
      this->P(static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==16){  //16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data16, sizeof(float), 1, fp);
      this->P(static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==64){  //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data64, sizeof(double), 1, fp);
      this->P(static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==256){  //256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data256, sizeof(signed char), 1, fp);
      this->P(static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==512){  //512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data512, sizeof(unsigned short), 1, fp);
      this->P(static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==768){  //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data768, sizeof(unsigned int), 1, fp);
      this->P(static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else{
    cout << "I can't open an image with grey levels of this type" << endl;
  }
  
  fclose(fp);

}


#define HEADER_N_VOX_OFFSET(x,y,z,t,s) (long)(hdr.vox_offset+((x)+(y*OrigNBX)+(z*OrigNBY*OrigNBX)+(t*OrigNBZ*OrigNBY*OrigNBX))*sizeof(s))

#define HEADER_N_VOX_OFFSET2(x,y,z,t,sizeof_s) (long)(hdr.vox_offset+((x)+(y*OrigNBX)+(z*OrigNBY*OrigNBX)+(t*OrigNBZ*OrigNBY*OrigNBX))*sizeof_s)

///read a scalar field in a ROI (in a nifti image)
/// -> advanced memory managment for large images: only allocate memory for the ROI
/// -> Inputs are different than in Read_ROI_Given_ImageToWorld_and_Size:
///     We give here the min and max {X,Y,Z,T} in the input image defining the outputed ROI
void ScalarField::Read_only_ROI(char * ImageName,int xMin,int xMax,int yMin,int yMax,int zMin,int zMax,int tMin,int tMax){
  nifti_1_header hdr;
  int tempInt;
  FILE *fp;
  int ret,i;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z,t;
  float a,b,c,d,qfac;
  int OrigNBX,OrigNBY,OrigNBZ,OrigNBT;
  float RealOriginX,RealOriginY,RealOriginZ;
  float VoxOriginX,VoxOriginY,VoxOriginZ;
  int sizeof_currentType;

  
  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) open the file and read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  
      //to open some multichannel 3D images which have the channels in the 5th dimension instead of the 4th dimension
  if ((hdr.dim[4]==1)&&(hdr.dim[5]>1)){
    hdr.dim[4]=hdr.dim[5];
    hdr.dim[5]=1;
  }

  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=hdr.dim[1])||(this->NY!=hdr.dim[2])||(this->NZ!=hdr.dim[3])||(this->NT!=hdr.dim[4]))
      cout << "The size of a scalar field is changed\n";
  
  
  //fill the parameters of the class
  if (hdr.dim[4]<1) {hdr.dim[4]=1;}
  
  if (xMax<xMin) {tempInt=xMin; xMin=xMax; xMax=tempInt;}
  if (yMax<yMin) {tempInt=yMin; yMin=yMax; yMax=tempInt;}
  if (zMax<zMin) {tempInt=zMin; zMin=zMax; zMax=tempInt;}
  if (tMax<tMin) {tempInt=tMin; tMin=tMax; tMax=tempInt;}
  
  if (xMin<0) xMin=0;
  if (yMin<0) yMin=0;
  if (zMin<0) zMin=0;
  if (tMin<0) tMin=0;
   
  if (xMax>(static_cast<int>(hdr.dim[1]))) xMax=(static_cast<int>(hdr.dim[1]));
  if (yMax>(static_cast<int>(hdr.dim[2]))) yMax=(static_cast<int>(hdr.dim[2]));
  if (zMax>(static_cast<int>(hdr.dim[3]))) zMax=(static_cast<int>(hdr.dim[3]));
  if (tMax>(static_cast<int>(hdr.dim[4]))) tMax=(static_cast<int>(hdr.dim[4]));

  this->NX=xMax-xMin;
  this->NY=yMax-yMin;
  this->NZ=zMax-zMin;
  this->NT=tMax-tMin;
  
  OrigNBX=static_cast<int>(hdr.dim[1]);
  OrigNBY=static_cast<int>(hdr.dim[2]);
  OrigNBZ=static_cast<int>(hdr.dim[3]);
  OrigNBT=static_cast<int>(hdr.dim[4]);

  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;

  //Image to world matrix
  if (hdr.sform_code>0){ 
    this->Image2World[0][0]=hdr.srow_x[0];  this->Image2World[0][1]=hdr.srow_x[1];  this->Image2World[0][2]=hdr.srow_x[2];  this->Image2World[0][3]=hdr.srow_x[3];
    this->Image2World[1][0]=hdr.srow_y[0];  this->Image2World[1][1]=hdr.srow_y[1];  this->Image2World[1][2]=hdr.srow_y[2];  this->Image2World[1][3]=hdr.srow_y[3];
    this->Image2World[2][0]=hdr.srow_z[0];  this->Image2World[2][1]=hdr.srow_z[1];  this->Image2World[2][2]=hdr.srow_z[2];  this->Image2World[2][3]=hdr.srow_z[3];
    this->Image2World[3][0]=0;              this->Image2World[3][1]=0;              this->Image2World[3][2]=0;              this->Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      this->Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); this->Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       this->Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     this->Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   this->Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     this->Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       this->Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;                               this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                    this->Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      this->Image2World[0][0]=hdr.pixdim[1];  this->Image2World[0][1]=0;             this->Image2World[0][2]=0;             this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=0;              this->Image2World[1][1]=hdr.pixdim[2]; this->Image2World[1][2]=0;             this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=0;              this->Image2World[2][1]=0;             this->Image2World[2][2]=hdr.pixdim[3]; this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;              this->Image2World[3][1]=0;             this->Image2World[3][2]=0;             this->Image2World[3][3]=1;
    }
  }

  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //message about the type of grey levels
  cout << endl;
  cout << ImageName << " contains";
  if (hdr.datatype==2) {cout << " unsigned char ";         sizeof_currentType=static_cast<int>(sizeof(unsigned char));}
  else if (hdr.datatype==4) {cout << " signed short ";     sizeof_currentType=static_cast<int>(sizeof(signed short));}
  else if (hdr.datatype==8) {cout << " signed int ";       sizeof_currentType=static_cast<int>(sizeof(signed int));}
  else if (hdr.datatype==16) {cout << " float ";           sizeof_currentType=static_cast<int>(sizeof(float));}
  else if (hdr.datatype==64) {cout << " double ";          sizeof_currentType=static_cast<int>(sizeof(double));}
  else if (hdr.datatype==256) {cout << " signed char ";    sizeof_currentType=static_cast<int>(sizeof(signed char));}
  else if (hdr.datatype==512) {cout << " unsigned short "; sizeof_currentType=static_cast<int>(sizeof(unsigned short));}
  else if (hdr.datatype==768) {cout << " unsigned int ";   sizeof_currentType=static_cast<int>(sizeof(unsigned int));}
  else{
    cout << endl;
    cout << "I can't open an image with grey levels of this type" << endl;
    exit(0);
  }
  cout << "pixels and has a resolution of " << OrigNBX << "*"  << OrigNBY << "*"  << OrigNBZ << "*"  << OrigNBT << " voxels" << endl;
  cout << "The ROI has a resolution " << this->NX << "*" << this->NY << "*" << this->NZ << "*" << this->NT << endl;
  cout << endl;
  
  //test the multiplicatory value of hdr
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    //cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << ImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }

  
  //2) read the image -- OLD MEMORY EFFICIENT VERSION (WHICH SEEMS TO CONTAIN A BUG)
  /*
  //allocate the memory for the image
  this->ScalField = new float [this->NXtYtZ*this->NT];

  // jump to data offset
  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);

  
  //load the image  (MODIFIED FILE BUT STILL SEEMS TO CONTAIN A BUG)
  t=tMin;
  while (t<tMax){
    z=zMin;
    while (z<zMax){
      y=yMin;
      while (y<yMax){
        x=xMin;
        ret = fseek(fp, HEADER_N_VOX_OFFSET2(x,y,z,t,sizeof_currentType), SEEK_SET);
        while (x<xMax){
          if (hdr.datatype==16){//16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data16, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==2){//2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data2, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==8){//8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data8, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==4){//4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data4, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==64){ //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data64, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==256){//256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data256, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==512){//512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data512, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==768){ //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data768, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          x++;
        }
        y++;
      }
      z++;
    }
    t++;
  }
   
  fclose(fp);
  */
  
  //2) read the image -- HIGHLY MEMORY CONSUMING (BUT BUG FREE) VERSION
  
  cerr << "The original algorithm which only allocates memory for the ROI seems to contain a bug. " << endl;
  cerr << "An alternative version consuming more memory but bug free is used instead." << endl;

  //allocate the memory for the image
  this->ScalField = new float [this->NXtYtZ*this->NT];

  //read the orginal image
  ScalarField tmpImag;
  tmpImag.Read(ImageName);
  
  //copy the ROI
  for (t=0;t<this->NT;t++) for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++) 
    this->P(tmpImag.G(x+xMin,y+yMin,z+zMin,t+tMin),x,y,z,t);
  
  
  //3) update the image to world matrix
  VoxOriginX=static_cast<float>(xMin);
  VoxOriginY=static_cast<float>(yMin);
  VoxOriginZ=static_cast<float>(zMin);
  
  
  RealOriginX=this->Image2World[0][0]*VoxOriginX+this->Image2World[0][1]*VoxOriginY+this->Image2World[0][2]*VoxOriginZ+this->Image2World[0][3];
  RealOriginY=this->Image2World[1][0]*VoxOriginX+this->Image2World[1][1]*VoxOriginY+this->Image2World[1][2]*VoxOriginZ+this->Image2World[1][3];
  RealOriginZ=this->Image2World[2][0]*VoxOriginX+this->Image2World[2][1]*VoxOriginY+this->Image2World[2][2]*VoxOriginZ+this->Image2World[2][3];
  
  this->Image2World[0][3]=RealOriginX;
  this->Image2World[1][3]=RealOriginY;
  this->Image2World[2][3]=RealOriginZ;

  invert_4t4quaternion(this->Image2World,this->World2Image);
  
}




///read a scalar field in a ROI (in a nifti image)
/// -> advanced memory managment for large images: only allocate memory for the ROI
/// -> Inputs are different than in Read_only_ROI:
///     We give here the min and max {X,Y,Z,T} in the input image defining the outputed ROI
void ScalarField::Read_ROI_Given_ImageToWorld_and_Size(float ROI_Image2World[4][4],int NBX,int NBY,int NBZ,char * RefImageName){
  nifti_1_header hdr;
  FILE *fp;
  int ret;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z;
  float a,b,c,d,qfac;
  int OrigNBX,OrigNBY,OrigNBZ;
  int sizeof_currentType;
  float RefIma_Image2World[4][4];
  float RefIma_World2Image[4][4];
  float ImaCoord_ROI2RefIma[4][4];
  float Orig_coord[4];
  float ROI_coord[4];
  int x_ref,y_ref,z_ref;
  
  
  //1) generate a void image in this and set its image 2 world properties
  
  this->CreateVoidField(NBX,NBY,NBZ);
  
  this->Image2World[0][0]=ROI_Image2World[0][0]; this->Image2World[0][1]=ROI_Image2World[0][1]; this->Image2World[0][2]=ROI_Image2World[0][2]; this->Image2World[0][3]=ROI_Image2World[0][3];
  this->Image2World[1][0]=ROI_Image2World[1][0]; this->Image2World[1][1]=ROI_Image2World[1][1]; this->Image2World[1][2]=ROI_Image2World[1][2]; this->Image2World[1][3]=ROI_Image2World[1][3];
  this->Image2World[2][0]=ROI_Image2World[2][0]; this->Image2World[2][1]=ROI_Image2World[2][1]; this->Image2World[2][2]=ROI_Image2World[2][2]; this->Image2World[2][3]=ROI_Image2World[2][3];
  this->Image2World[3][0]=ROI_Image2World[3][0]; this->Image2World[3][1]=ROI_Image2World[3][1]; this->Image2World[3][2]=ROI_Image2World[3][2]; this->Image2World[3][3]=ROI_Image2World[3][3];
  
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  cout << this->Image2World[0][0] << " " << this->Image2World[0][1] <<  " " << this->Image2World[0][2] <<  " " << this->Image2World[0][3] << endl;
  cout << this->Image2World[1][0] << " " << this->Image2World[1][1] <<  " " << this->Image2World[1][2] <<  " " << this->Image2World[1][3] << endl;
  cout << this->Image2World[2][0] << " " << this->Image2World[2][1] <<  " " << this->Image2World[2][2] <<  " " << this->Image2World[2][3] << endl;
  cout << this->Image2World[3][0] << " " << this->Image2World[3][1] <<  " " << this->Image2World[3][2] <<  " " << this->Image2World[3][3] << endl;


  //2) open the reference file and read its properties

  //2.1) open the file, read the header and make some security treatments on the header
  fp = fopen(RefImageName,"rb");
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  if (hdr.dim[4]<1) {hdr.dim[4]=1;}
  
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    //cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << RefImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }

  
  //2.2) Image to world matrix of the reference image
  if (hdr.sform_code>0){ 
    RefIma_Image2World[0][0]=hdr.srow_x[0];  RefIma_Image2World[0][1]=hdr.srow_x[1];  RefIma_Image2World[0][2]=hdr.srow_x[2];  RefIma_Image2World[0][3]=hdr.srow_x[3];
    RefIma_Image2World[1][0]=hdr.srow_y[0];  RefIma_Image2World[1][1]=hdr.srow_y[1];  RefIma_Image2World[1][2]=hdr.srow_y[2];  RefIma_Image2World[1][3]=hdr.srow_y[3];
    RefIma_Image2World[2][0]=hdr.srow_z[0];  RefIma_Image2World[2][1]=hdr.srow_z[1];  RefIma_Image2World[2][2]=hdr.srow_z[2];  RefIma_Image2World[2][3]=hdr.srow_z[3];
    RefIma_Image2World[3][0]=0;              RefIma_Image2World[3][1]=0;              RefIma_Image2World[3][2]=0;              RefIma_Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      RefIma_Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); RefIma_Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       RefIma_Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     RefIma_Image2World[0][3]=hdr.qoffset_x;
      RefIma_Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     RefIma_Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   RefIma_Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     RefIma_Image2World[1][3]=hdr.qoffset_y;
      RefIma_Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     RefIma_Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       RefIma_Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); RefIma_Image2World[2][3]=hdr.qoffset_z;
      RefIma_Image2World[3][0]=0;                               RefIma_Image2World[3][1]=0;                                 RefIma_Image2World[3][2]=0;                                    RefIma_Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      
      //put the voxel dimensions in image to world
      cout << "Orientations of " << RefImageName << " were basically estimated..." << endl;
      RefIma_Image2World[0][0]=hdr.pixdim[1];  RefIma_Image2World[0][1]=0;             RefIma_Image2World[0][2]=0;             RefIma_Image2World[0][3]=hdr.qoffset_x;
      RefIma_Image2World[1][0]=0;              RefIma_Image2World[1][1]=hdr.pixdim[2]; RefIma_Image2World[1][2]=0;             RefIma_Image2World[1][3]=hdr.qoffset_y;
      RefIma_Image2World[2][0]=0;              RefIma_Image2World[2][1]=0;             RefIma_Image2World[2][2]=hdr.pixdim[3]; RefIma_Image2World[2][3]=hdr.qoffset_z;
      RefIma_Image2World[3][0]=0;              RefIma_Image2World[3][1]=0;             RefIma_Image2World[3][2]=0;             RefIma_Image2World[3][3]=1;
    }
  }

  invert_4t4quaternion(RefIma_Image2World,RefIma_World2Image);
  
  //2.3) Size of the reference image
  OrigNBX=static_cast<int>(hdr.dim[1]);
  OrigNBY=static_cast<int>(hdr.dim[2]);
  OrigNBZ=static_cast<int>(hdr.dim[3]);

  
  //2.4) compute the matrix to transform the voxel coordinates from the output ROI to the input reference image
  mult_quat4t4mat_quat4t4mat(RefIma_World2Image,this->Image2World,ImaCoord_ROI2RefIma);
  
  cout << "ROI to Ref. image in image coordinates:" << endl;
  cout << ImaCoord_ROI2RefIma[0][0] << " " << ImaCoord_ROI2RefIma[0][1] <<  " " << ImaCoord_ROI2RefIma[0][2] <<  " " << ImaCoord_ROI2RefIma[0][3] << endl;
  cout << ImaCoord_ROI2RefIma[1][0] << " " << ImaCoord_ROI2RefIma[1][1] <<  " " << ImaCoord_ROI2RefIma[1][2] <<  " " << ImaCoord_ROI2RefIma[1][3] << endl;
  cout << ImaCoord_ROI2RefIma[2][0] << " " << ImaCoord_ROI2RefIma[2][1] <<  " " << ImaCoord_ROI2RefIma[2][2] <<  " " << ImaCoord_ROI2RefIma[2][3] << endl;
  cout << ImaCoord_ROI2RefIma[3][0] << " " << ImaCoord_ROI2RefIma[3][1] <<  " " << ImaCoord_ROI2RefIma[3][2] <<  " " << ImaCoord_ROI2RefIma[3][3] << endl;

  
  //2.5) define sizeof(type in the ref image)
  if (hdr.datatype==2) {        sizeof_currentType=static_cast<int>(sizeof(unsigned char));}
  else if (hdr.datatype==4) {   sizeof_currentType=static_cast<int>(sizeof(signed short));}
  else if (hdr.datatype==8) {   sizeof_currentType=static_cast<int>(sizeof(signed int));}
  else if (hdr.datatype==16) {  sizeof_currentType=static_cast<int>(sizeof(float));}
  else if (hdr.datatype==64) {  sizeof_currentType=static_cast<int>(sizeof(double));}
  else if (hdr.datatype==256) { sizeof_currentType=static_cast<int>(sizeof(signed char));}
  else if (hdr.datatype==512) { sizeof_currentType=static_cast<int>(sizeof(unsigned short));}
  else if (hdr.datatype==768) { sizeof_currentType=static_cast<int>(sizeof(unsigned int));}
  else{
    cout << endl;
    cout << "I can't open an image with grey levels of this type" << endl;
    exit(0);
  }
  
  
  //3) read the image
  ROI_coord[3]=1;
  for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++) {
    //3.1) go to where the point is in the reference image
    ROI_coord[0]=x;
    ROI_coord[1]=y;
    ROI_coord[2]=z;
    
    mult_4t4mat_4vec(ImaCoord_ROI2RefIma,ROI_coord,Orig_coord);
  
    x_ref=static_cast<int>(Orig_coord[0]+0.5);
    y_ref=static_cast<int>(Orig_coord[1]+0.5);
    z_ref=static_cast<int>(Orig_coord[2]+0.5);
    
    if (x_ref<0) x_ref=0;
    if (y_ref<0) y_ref=0;
    if (z_ref<0) z_ref=0;
    if (x_ref>=OrigNBX) x_ref=OrigNBX-1;
    if (y_ref>=OrigNBY) y_ref=OrigNBY-1;
    if (z_ref>=OrigNBZ) z_ref=OrigNBZ-1;
  
    ret = fseek(fp, HEADER_N_VOX_OFFSET2(x_ref,y_ref,z_ref,0,sizeof_currentType), SEEK_SET);
    
    //3.2) pick up the grey level
    if (hdr.datatype==16){//16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data16, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==2){//2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data2, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==8){//8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data8, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==4){//4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data4, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==64){ //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data64, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==256){//256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data256, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==512){//512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data512, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==768){ //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data768, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
  }
  


  fclose(fp);

}




///read a scalar field in a ROI (in a nifti image)
/// -> advanced memory managment for large images: undersample the image direcly by averaging blocks of voxels of size BlockS*BlockS*BlockS
///Use 'Read_and_Undersample' or 'Read_and_Interpolate' to have finer resamplings requiring more memory
void ScalarField::Read_directly_Undersampled(char * ImageName,int BlockS){
  nifti_1_header hdr;
  int tempInt;
  FILE *fp;
  int ret,i;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z,t;
  float a,b,c,d,qfac;
  float BlockVolume;
  float TempValue;
  int LocTmpX,LocTmpY,LocTmpZ;
  
  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) open the file and read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=hdr.dim[1])||(this->NY!=hdr.dim[2])||(this->NZ!=hdr.dim[3])||(this->NT!=hdr.dim[4]))
      cout << "The size of a scalar field is changed\n";
  
  
  //fill the parameters of the class
  if (hdr.dim[4]<1) {hdr.dim[4]=1;}
  
  this->NX=static_cast<int>(hdr.dim[1]/BlockS);
  this->NY=static_cast<int>(hdr.dim[2]/BlockS);
  this->NZ=static_cast<int>(hdr.dim[3]/BlockS);
  this->NT=hdr.dim[4];
  
  BlockVolume=static_cast<float>(BlockS*BlockS*BlockS);
  
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;

  //Image to world matrix
  if (hdr.sform_code>0){ 
    this->Image2World[0][0]=hdr.srow_x[0];  this->Image2World[0][1]=hdr.srow_x[1];  this->Image2World[0][2]=hdr.srow_x[2];  this->Image2World[0][3]=hdr.srow_x[3];
    this->Image2World[1][0]=hdr.srow_y[0];  this->Image2World[1][1]=hdr.srow_y[1];  this->Image2World[1][2]=hdr.srow_y[2];  this->Image2World[1][3]=hdr.srow_y[3];
    this->Image2World[2][0]=hdr.srow_z[0];  this->Image2World[2][1]=hdr.srow_z[1];  this->Image2World[2][2]=hdr.srow_z[2];  this->Image2World[2][3]=hdr.srow_z[3];
    this->Image2World[3][0]=0;              this->Image2World[3][1]=0;              this->Image2World[3][2]=0;              this->Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      this->Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); this->Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       this->Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     this->Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   this->Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     this->Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       this->Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;                               this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                    this->Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      this->Image2World[0][0]=hdr.pixdim[1];  this->Image2World[0][1]=0;             this->Image2World[0][2]=0;             this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=0;              this->Image2World[1][1]=hdr.pixdim[2]; this->Image2World[1][2]=0;             this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=0;              this->Image2World[2][1]=0;             this->Image2World[2][2]=hdr.pixdim[3]; this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;              this->Image2World[3][1]=0;             this->Image2World[3][2]=0;             this->Image2World[3][3]=1;
    }
  }

  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //message about the type of grey levels
  cout << endl;
  cout << ImageName << " contains";
  if (hdr.datatype==2) cout << " unsigned char ";
  if (hdr.datatype==4) cout << " signed short ";
  if (hdr.datatype==8) cout << " signed int ";
  if (hdr.datatype==16) cout << " float ";
  if (hdr.datatype==64) cout << " double ";
  if (hdr.datatype==256) cout << " signed char ";
  if (hdr.datatype==512) cout << " unsigned short ";
  if (hdr.datatype==768) cout << " unsigned int ";
  cout << "pixels and has a resolution of " << hdr.dim[1] << "*"  << hdr.dim[2] << "*"  << hdr.dim[3] << "*"  << hdr.dim[4] << " voxels" << endl;
  cout << "The undersampled image has a resolution " << this->NX << "*" << this->NY << "*" << this->NZ << "*" << this->NT << endl;
  cout << endl;
  
  //2) read the image
  
  //allocate the memory for the image and fill it with zeros
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(0,x,y,z,t);
  

  // jump to data offset
  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);

  //test the multiplicatory value of hdr
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    //cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << ImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }
  
  //load the image
  if (hdr.datatype==2){ //2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data2, sizeof(unsigned char), 1, fp);
      
      TempValue=static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==4){ //4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++    
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data4, sizeof(signed short), 1, fp);
      
      TempValue=static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==8){  //8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data8, sizeof(signed int), 1, fp);
      
      TempValue=static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==16){  //16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data16, sizeof(float), 1, fp);
      
      TempValue=static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==64){  //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data64, sizeof(double), 1, fp);
      
      TempValue=static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==256){  //256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data256, sizeof(signed char), 1, fp);
      
      TempValue=static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==512){  //512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data512, sizeof(unsigned short), 1, fp);
      
      TempValue=static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==768){  //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data768, sizeof(unsigned int), 1, fp);
      
      TempValue=static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else{
    cout << "I can't open an image with grey levels of this type" << endl;
  }
  
  fclose(fp);
  
  //3) update the image to world matrix
  this->Image2World[0][0]*=BlockS; this->Image2World[0][1]*=BlockS; this->Image2World[0][2]*=BlockS; 
  this->Image2World[1][0]*=BlockS; this->Image2World[1][1]*=BlockS; this->Image2World[1][2]*=BlockS; 
  this->Image2World[2][0]*=BlockS; this->Image2World[2][1]*=BlockS; this->Image2World[2][2]*=BlockS; 

  invert_4t4quaternion(this->Image2World,this->World2Image);
  
}


///read a scalar field and expend its domain
void ScalarField::ReadAndExpend(char * ImageName,int addX1,int addX2,int addY1,int addY2,int addZ1,int addZ2){
  ScalarField OrigSF;
  int x,y,z,t;
  float x2,y2,z2;
  float minImag;
  
  
  //read the scalar field at the original format
  OrigSF.Read(ImageName);
  
  //fill the parameters of the class and allocate the memory for the scalar field
  this->NX=OrigSF.NX+addX1+addX2;
  this->NY=OrigSF.NY+addY1+addY2;
  this->NZ=OrigSF.NZ+addZ1+addZ2;
  this->NT=OrigSF.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  
  this->Image2World[0][0]=OrigSF.Image2World[0][0];  this->Image2World[0][1]=OrigSF.Image2World[0][1];  this->Image2World[0][2]=OrigSF.Image2World[0][2];  this->Image2World[0][3]=OrigSF.Image2World[0][3]-addX1*OrigSF.Image2World[0][0];
  this->Image2World[1][0]=OrigSF.Image2World[1][0];  this->Image2World[1][1]=OrigSF.Image2World[1][1];  this->Image2World[1][2]=OrigSF.Image2World[1][2];  this->Image2World[1][3]=OrigSF.Image2World[1][3]-addY1*OrigSF.Image2World[1][1];
  this->Image2World[2][0]=OrigSF.Image2World[3][0];  this->Image2World[2][1]=OrigSF.Image2World[2][1];  this->Image2World[2][2]=OrigSF.Image2World[2][2];  this->Image2World[2][3]=OrigSF.Image2World[2][3]-addZ1*OrigSF.Image2World[2][2];
  this->Image2World[3][0]=0;                         this->Image2World[3][1]=0;                         this->Image2World[3][2]=0;                         this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //fill the image image...
  
  //...init
  minImag=OrigSF.G(0,0,0,0);
  
  for(t=0;t<OrigSF.NT;t++) for(z=0;z<OrigSF.NZ;z++) for(y=0;y<OrigSF.NY;y++) for(x=0;x<OrigSF.NX;x++) if (minImag>OrigSF.G(x,y,z,t)) minImag=OrigSF.G(x,y,z,t);

  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) this->P(minImag-1,x,y,z,t);
  
  //... original image
  for(t=0;t<OrigSF.NT;t++) for(z=0;z<OrigSF.NZ;z++) for(y=0;y<OrigSF.NY;y++) for(x=0;x<OrigSF.NX;x++)
    if ((x+addX1>=0)&&(x+addX1<this->NX)&& (y+addY1>=0)&&(y+addY1<this->NY)&& (z+addZ1>=0)&&(z+addZ1<this->NZ))
      this->P(OrigSF.G(x,y,z,t),x+addX1,y+addY1,z+addZ1,t);

  //... image extension
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) if (this->G(x,y,z,t)<minImag){
    x2=x;
    if (x<addX1) x2=addX1;
    if (x>=OrigSF.NX+addX1) x2=OrigSF.NX+addX1-1;
    y2=y;
    if (y<addY1) y2=addY1;
    if (y>=OrigSF.NY+addY1) y2=OrigSF.NY+addY1-1;
    z2=z;
    if (z<addZ1) z2=addZ1;
    if (z>=OrigSF.NZ+addZ1) z2=OrigSF.NZ+addZ1-1;
    
    this->P(this->G(x2,y2,z2,t),x,y,z,t);
  
  }
}



///read a scalar field and perform linear interpolation to give it a specific size
void ScalarField::Read_and_Interpolate(char * ImageName,int NBX,int NBY,int NBZ){
  ScalarField OrigSF;
  int x,y,z,t;
  float x2,y2,z2;
  float factorX,factorY,factorZ;
  
  
  //read the scalar field at the original format
  OrigSF.Read(ImageName);
  
  //fill the parameters of the class and allocate the memory for the scalar field
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=OrigSF.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  factorX=(static_cast<float>(OrigSF.NX)/static_cast<float>(this->NX));
  factorY=(static_cast<float>(OrigSF.NY)/static_cast<float>(this->NY));
  factorZ=(static_cast<float>(OrigSF.NZ)/static_cast<float>(this->NZ));
  
  this->Image2World[0][0]=OrigSF.Image2World[0][0]*factorX;  this->Image2World[0][1]=OrigSF.Image2World[0][1]*factorX;  this->Image2World[0][2]=OrigSF.Image2World[0][2]*factorX;  this->Image2World[0][3]=OrigSF.Image2World[0][3];
  this->Image2World[1][0]=OrigSF.Image2World[1][0]*factorY;  this->Image2World[1][1]=OrigSF.Image2World[1][1]*factorY;  this->Image2World[1][2]=OrigSF.Image2World[1][2]*factorY;  this->Image2World[1][3]=OrigSF.Image2World[1][3];
  this->Image2World[2][0]=OrigSF.Image2World[3][0]*factorZ;  this->Image2World[2][1]=OrigSF.Image2World[2][1]*factorZ;  this->Image2World[2][2]=OrigSF.Image2World[2][2]*factorZ;  this->Image2World[2][3]=OrigSF.Image2World[2][3];
  this->Image2World[3][0]=0;                                 this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                 this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //interpolate the original image
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t),x,y,z,t);
        }
      }
    }
  }
}

///read a scalar field and perform undersampling by a factor 'factor'
void ScalarField::Read_and_Undersample(char * ImageName,float factor, float NN){
  ScalarField OrigSF;
  int x,y,z,t;
  float x2,y2,z2;
  float factorX,factorY,factorZ;
  
  
  //read the scalar field at the original format
  OrigSF.Read(ImageName);
  
  //fill the parameters of the class and allocate the memory for the scalar field
  this->NX=static_cast<int>(OrigSF.NX/factor);
  this->NY=static_cast<int>(OrigSF.NY/factor);
  if (OrigSF.NZ>1)
    this->NZ=static_cast<int>(OrigSF.NZ/factor);
  else
    this->NZ=static_cast<int>(OrigSF.NZ/1);
  this->NT=OrigSF.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  factorX=(static_cast<float>(OrigSF.NX)/static_cast<float>(this->NX));
  factorY=(static_cast<float>(OrigSF.NY)/static_cast<float>(this->NY));
  factorZ=(static_cast<float>(OrigSF.NZ)/static_cast<float>(this->NZ));
  
  this->Image2World[0][0]=OrigSF.Image2World[0][0]*factorX;  this->Image2World[0][1]=OrigSF.Image2World[0][1]*factorX;  this->Image2World[0][2]=OrigSF.Image2World[0][2]*factorX;  this->Image2World[0][3]=OrigSF.Image2World[0][3];
  this->Image2World[1][0]=OrigSF.Image2World[1][0]*factorY;  this->Image2World[1][1]=OrigSF.Image2World[1][1]*factorY;  this->Image2World[1][2]=OrigSF.Image2World[1][2]*factorY;  this->Image2World[1][3]=OrigSF.Image2World[1][3];
  this->Image2World[2][0]=OrigSF.Image2World[3][0]*factorZ;  this->Image2World[2][1]=OrigSF.Image2World[2][1]*factorZ;  this->Image2World[2][2]=OrigSF.Image2World[2][2]*factorZ;  this->Image2World[2][3]=OrigSF.Image2World[2][3];
  this->Image2World[3][0]=0;                                 this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                 this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //interpolate the original image
  if (NN==0){  //triliear interpolation
    for(t=0;t<this->NT;t++){
      for(z=0;z<this->NZ;z++){ 
        z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
        for(y=0;y<this->NY;y++){ 
          y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
          for(x=0;x<this->NX;x++){
            x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
            this->P(OrigSF.G(x2,y2,z2,t),x,y,z,t);
          }
        }
      }
    }
  }
  else{   //nearest neighbor
    for(t=0;t<this->NT;t++){
      for(z=0;z<this->NZ;z++){ 
        z2=floor(0.5+static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1));
        for(y=0;y<this->NY;y++){ 
          y2=floor(0.5+static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1));
          for(x=0;x<this->NX;x++){
            x2=floor(0.5+static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1));
            this->P(OrigSF.G(x2,y2,z2,t),x,y,z,t);
          }
        }
      }
    }
  }
}

///create a void scalar field. All values are initialized to 'cste' which is null by default. No message is printed if Verbose!=1.
void ScalarField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT,float cste,int Verbose){
  int x,y,z,t;
  int HasSameSize;

  //check if the size is changed
  if ((this->NX==NBX)&&(this->NY==NBY)&&(this->NZ==NBZ)&&(this->NT==NBT))
    HasSameSize=1;
  else 
    HasSameSize=0;
  
  if ((this->NX!=0)&&(HasSameSize==0)){
    if (Verbose==1) cout << "The size of a non-null scalar field is changed\n";
    delete this->ScalField;
  }
  
  //image size
  if (HasSameSize==0){
    this->NX=NBX;
    this->NY=NBY;
    this->NZ=NBZ;
    this->NT=NBT;
    this->NXtY=this->NX*this->NY;
    this->NXtYtZ=this->NXtY*this->NZ;
  }  
  
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //allocate memory to cast (and eventually transform) the original template and target images
  //    -->  ScalarField[ptSF(x,y,z)]= gray level at (x,y,z)
  if (HasSameSize==0){
    this->ScalField = new float [this->NXtYtZ*this->NT];
  }
  
  //set all entries of the field at 0.
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(cste,x,y,z,t);
}



///Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
void ScalarField::SlashFieldSize(int verbative){
  
  //image size
  this->NX=1;
  this->NY=1;
  this->NZ=1;
  this->NT=1;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  if (verbative==1) cout << "Slash some memory" << endl;
  
  //re-allocate memory
  if (this->NX!=0) delete this->ScalField;
  this->ScalField = new float [1];
}



///write a scalar field in a nifti image
void ScalarField::Write(char * OutputImageName){
  int i;
  int x,y,z,t;
  nifti_1_header hdr;
  nifti1_extender pad={0,0,0,0};
  FILE *fp;
  int ret;
  float *data=NULL;
  
  
  //1) create the header
  memset((void *)&hdr,0, sizeof(hdr));
  hdr.sizeof_hdr = MIN_HEADER_SIZE;
  hdr.dim[0] = 4;
  hdr.dim[1] = this->NX;
  hdr.dim[2] = this->NY;
  hdr.dim[3] = this->NZ;
  hdr.dim[4] = this->NT;
  hdr.datatype = NIFTI_TYPE_FLOAT32;
  hdr.bitpix = 32; 
  hdr.qform_code=0; // should ideally be set to 1 but I don't set the values of 'quatern_b', 'quatern_c' and 'quatern_d'
  hdr.pixdim[1] = sqrt(this->Image2World[0][0]*this->Image2World[0][0]+this->Image2World[0][1]*this->Image2World[0][1]+this->Image2World[0][2]*this->Image2World[0][2]);
  hdr.pixdim[2] = sqrt(this->Image2World[1][0]*this->Image2World[1][0]+this->Image2World[1][1]*this->Image2World[1][1]+this->Image2World[1][2]*this->Image2World[1][2]);
  hdr.pixdim[3] = sqrt(this->Image2World[2][0]*this->Image2World[2][0]+this->Image2World[2][1]*this->Image2World[2][1]+this->Image2World[2][2]*this->Image2World[2][2]);
  hdr.qoffset_x=this->Image2World[0][3];
  hdr.qoffset_y=this->Image2World[1][3];
  hdr.qoffset_z=this->Image2World[2][3];
  hdr.pixdim[4] = 1.0;
  hdr.sform_code=1;
  hdr.srow_x[0]=this->Image2World[0][0];  hdr.srow_x[1]=this->Image2World[0][1];  hdr.srow_x[2]=this->Image2World[0][2];  hdr.srow_x[3]=this->Image2World[0][3];
  hdr.srow_y[0]=this->Image2World[1][0];  hdr.srow_y[1]=this->Image2World[1][1];  hdr.srow_y[2]=this->Image2World[1][2];  hdr.srow_y[3]=this->Image2World[1][3];
  hdr.srow_z[0]=this->Image2World[2][0];  hdr.srow_z[1]=this->Image2World[2][1];  hdr.srow_z[2]=this->Image2World[2][2];  hdr.srow_z[3]=this->Image2World[2][3];
  hdr.vox_offset = (float) NII_HEADER_SIZE;
  hdr.scl_inter = 0.0;
  hdr.scl_slope = 1.0;
  hdr.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
  strncpy(hdr.magic, "n+1\0", 4);

  //2) save the image OutputImageName
  //allocate and fill the buffer 
  data = new float [hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]];
  
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(x,y,z,t);
    i++;
  }
  
  // write first 348 bytes of header  
  fp = fopen(OutputImageName,"wb");
  ret = fwrite(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  // write extender pad and image data  
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4], fp);
  
  fclose(fp);
  
}



///write a scalar field in a nifti image
//The 2nd file is an input file containing the headers 
void ScalarField::Write(char * OutputImageName, char * ImageForHeaderName){
  nifti_1_header hdr_ref;
  FILE *fp_header;
  int i;
  int x,y,z,t;
  nifti1_extender pad={0,0,0,0};
  FILE *fp;
  int ret;
  float *data=NULL;
  
  
  //1) read the header of ImageForHeaderName
  fp_header = fopen(ImageForHeaderName,"rb");
  fread(&hdr_ref, MIN_HEADER_SIZE, 1, fp_header);
  fclose(fp_header);
  

  
  //print a little header information
  //fprintf(stderr, "\n%s header information:",ImageForHeaderName);
  //fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr_ref.dim[1],hdr_ref.dim[2],hdr_ref.dim[3],hdr_ref.dim[4]);
  //fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr_ref.datatype,hdr_ref.bitpix);
  //fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr_ref.scl_slope,hdr_ref.scl_inter);
  //fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr_ref.vox_offset));
  //fprintf(stderr, "\n");
  

  //2) save the image OutputImageName
  hdr_ref.dim[0] = 4;
  hdr_ref.dim[1] = this->NX;
  hdr_ref.dim[2] = this->NY;
  hdr_ref.dim[3] = this->NZ;
  hdr_ref.dim[4] = this->NT;
  hdr_ref.datatype = NIFTI_TYPE_FLOAT32;
  hdr_ref.bitpix = 32; 
  hdr_ref.scl_inter = 0.0;
  hdr_ref.scl_slope = 1.0;
  hdr_ref.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;

  
  // allocate and fill the buffer 
  //data = (float *) malloc(sizeof(float) * hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]);
  data = new  float [hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]];
  
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  // write first 348 bytes of header  
  fp = fopen(OutputImageName,"wb");

  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  
  // write extender pad and image data  
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  
  fclose(fp);
}


///General function used to align the cumulative histogram of *this to the input cumulative histogram.
///-> called by GreyLevAlignment, GreyLevAlignment_UsingRefROIs, HistoEqualization, HistoEqualizationInROI
void ScalarField::HistoMatch(int NbBinsCumHisto,float * TrgCumHisto_x_axis,float * TrgCumHisto_y_axis,float * LocCumHisto_x_axis,float * LocCumHisto_y_axis){
  int i,j,k,t;
  float * LUT_x_axis;
  float * LUT_y_axis;
  float tmpFl,tmpFl2;
  int CorrespBin;
  float weight_m,weight_p;
  float CorrespGL;
  

  //1) allocate and compute the cumulative histograms
  LUT_x_axis = new float [NbBinsCumHisto];
  LUT_y_axis = new float [NbBinsCumHisto];
  
  
  //2) compute the look-up table (LUT) from the intensities of the Local image to those of the reference image
  for (i=0;i<NbBinsCumHisto;i++) LUT_x_axis[i]=LocCumHisto_x_axis[i];
  
  for (i=0;i<NbBinsCumHisto;i++){
    //... tmpFl is between 0 and 1
    tmpFl=LocCumHisto_y_axis[i]; 
    
    //... find the corresponding bin of tmpFl in the inverse of CumHisto_TrgImage
    CorrespBin=-2;
    for (j=0;j<NbBinsCumHisto;j++) if (CorrespBin==-2) if (TrgCumHisto_y_axis[j]>tmpFl){
      CorrespBin=j-1;
      weight_p=(tmpFl-TrgCumHisto_y_axis[CorrespBin])/(TrgCumHisto_y_axis[CorrespBin+1]-TrgCumHisto_y_axis[CorrespBin]);
      weight_m=1-weight_p;
      }
    
    //... fill the LUT
    if (CorrespBin<0){
      LUT_y_axis[i]=TrgCumHisto_x_axis[0];
      }
    else  if ((CorrespBin==-2)||(CorrespBin>=NbBinsCumHisto-1)){
      LUT_y_axis[i]=TrgCumHisto_x_axis[NbBinsCumHisto-1];
      }
    else LUT_y_axis[i]=weight_m*TrgCumHisto_x_axis[CorrespBin]+weight_p*TrgCumHisto_x_axis[CorrespBin+1];
  }
  
  //... final refinement to exactely match the min and max grey levels  (x-axis of the histgrams is the center of the areas where the grey levels are considered)
  LUT_y_axis[0]=TrgCumHisto_x_axis[0]-((TrgCumHisto_x_axis[1]-TrgCumHisto_x_axis[0])/2);
  LUT_y_axis[NbBinsCumHisto-1]=TrgCumHisto_x_axis[NbBinsCumHisto-1]+((TrgCumHisto_x_axis[1]-TrgCumHisto_x_axis[0])/2);
  
  //for (i=0;i<NbBinsCumHisto;i++) 
  //  cout << LUT_x_axis[i] << " -> " << LUT_y_axis[i] << endl;
  //cout << LUT_x_axis[NbBinsCumHisto-1] << " -> " << LUT_y_axis[NbBinsCumHisto-1] << endl;
  
  //3) resample the image grey levels
  for (t=0;t<this->NT;t++) for (i=0;i<this->NZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++){
    //current grey level
    tmpFl=this->G(k,j,i,t);
    
    //corresponding bin in the LUT
    if (tmpFl<LUT_x_axis[0]) CorrespGL=LUT_y_axis[0];
    else if (tmpFl>=LUT_x_axis[NbBinsCumHisto-1]) CorrespGL=LUT_y_axis[NbBinsCumHisto-1];
    else{
      tmpFl2=static_cast<float>(NbBinsCumHisto)*(tmpFl-LUT_x_axis[0])/(LUT_x_axis[NbBinsCumHisto-1]-LUT_x_axis[0]);
      CorrespBin=static_cast<int>(tmpFl2);
      weight_p=tmpFl2-floor(tmpFl2);
      weight_m=1-weight_p;
      
      if (CorrespBin>=NbBinsCumHisto-1){
        CorrespBin=NbBinsCumHisto-2;
        weight_p=1;
        weight_m=0;
      }
      
      CorrespGL=weight_m*LUT_y_axis[CorrespBin]+weight_p*LUT_y_axis[CorrespBin+1];
    }
  //put the new grey level value
  this->P(CorrespGL,k,j,i,t);
  }
}
 
///grey levels alignment of the grey levels using optimal transport 
///-> optimal transportation minmizes the wassertein 1 distance between the linearly aligned histograms
void ScalarField::GreyLevAlignment(ScalarField * RefImage){
  int NbBinsHisto;
  float * CumHisto_TrgImage_x_axis;
  float * CumHisto_TrgImage_y_axis;
  float * CumHisto_LocImage_x_axis;
  float * CumHisto_LocImage_y_axis;
  

  //1) allocate and compute the cumulative histograms
  NbBinsHisto=256;
  CumHisto_TrgImage_x_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_y_axis = new float [NbBinsHisto];
  CumHisto_LocImage_x_axis = new float [NbBinsHisto];
  CumHisto_LocImage_y_axis = new float [NbBinsHisto];
  
  /*+++++ begin hack +++++*/
  //for (i=0;i<NbBinsHisto;i++) CumHisto_TrgImage_x_axis[i]=CumHisto_LocImage_x_axis[NbBinsHisto-1]*static_cast<float>(i)/static_cast<float>(NbBinsHisto-1);
  //for (i=0;i<NbBinsHisto;i++) CumHisto_TrgImage_y_axis[i]=static_cast<float>(i)/static_cast<float>(NbBinsHisto-1);
  /*+++++ end hack +++++*/

  RefImage->CptCumulativeHistogram(NbBinsHisto,CumHisto_TrgImage_x_axis,CumHisto_TrgImage_y_axis,1);   //remark: the log of the histogram is actually used
  this->CptCumulativeHistogram(NbBinsHisto,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis,1);   //remark: the log of the histogram is actually used
  
  //2) do the job
  this->HistoMatch(NbBinsHisto,CumHisto_TrgImage_x_axis,CumHisto_TrgImage_y_axis,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis);

}

 
///Grey levels histogram equalization using optimal transport 
void ScalarField::HistoEqualization(){
  int NbBinsHisto;
  float * CumHisto_TrgImage_x_axis;
  float * CumHisto_TrgImage_y_axis;
  float * CumHisto_LocImage_x_axis;
  float * CumHisto_LocImage_y_axis;
  int i;

  //1) allocate and compute the cumulative histograms
  NbBinsHisto=256;
  CumHisto_TrgImage_x_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_y_axis = new float [NbBinsHisto];
  CumHisto_LocImage_x_axis = new float [NbBinsHisto];
  CumHisto_LocImage_y_axis = new float [NbBinsHisto];
  
  this->CptCumulativeHistogram(NbBinsHisto,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis,1);   //remark: the log of the histogram is actually used
  
  for (i=0;i<NbBinsHisto;i++) CumHisto_TrgImage_x_axis[i]=CumHisto_LocImage_x_axis[NbBinsHisto-1]*static_cast<float>(i)/static_cast<float>(NbBinsHisto-1);
  for (i=0;i<NbBinsHisto;i++) CumHisto_TrgImage_y_axis[i]=static_cast<float>(i)/static_cast<float>(NbBinsHisto-1);

  //2) do the job
  this->HistoMatch(NbBinsHisto,CumHisto_TrgImage_x_axis,CumHisto_TrgImage_y_axis,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis);
}




///Same as "GreyLevAlignment" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI.
///Note that "LocImageROI" must have the same size as "this" and that "RefImageROI" must have the same size as "RefImage"
///ROIs are defined by non-null values
void ScalarField::GreyLevAlignment_UsingRefROIs(ScalarField * RefImage,ScalarField * RefImageROI,ScalarField * LocImageROI){
  int NbBinsHisto;
  float * CumHisto_LocImage_x_axis;
  float * CumHisto_LocImage_y_axis;
  float * CumHisto_TrgImage_x_axis;
  float * CumHisto_TrgImage_y_axis;
   

  //1) allocate and compute the cumulative histograms
  NbBinsHisto=256;
  CumHisto_LocImage_x_axis = new float [NbBinsHisto];
  CumHisto_LocImage_y_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_x_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_y_axis = new float [NbBinsHisto];
  
  this->CptCumulativeHistogram_InROI(LocImageROI,NbBinsHisto,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis,1);   //remark: the log of the histogram is actually used
  RefImage->CptCumulativeHistogram_InROI(RefImageROI,NbBinsHisto,CumHisto_TrgImage_x_axis,CumHisto_TrgImage_y_axis,1);
  
  //2) do the job
  this->HistoMatch(NbBinsHisto,CumHisto_TrgImage_x_axis,CumHisto_TrgImage_y_axis,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis);

}

///Same as "HistoEqualization" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI. 
///Note that "LocImageROI" must have the same size as "this"
///ROIs are defined by non-null values
void ScalarField::HistoEqualizationInROI(ScalarField * LocImageROI){
  int NbBinsHisto;
  float * CumHisto_LocImage_x_axis;
  float * CumHisto_LocImage_y_axis;
  float * CumHisto_TrgImage_x_axis;
  float * CumHisto_TrgImage_y_axis;
  int i;
   
  //1) allocate and compute the cumulative histograms
  NbBinsHisto=256;
  CumHisto_LocImage_x_axis = new float [NbBinsHisto];
  CumHisto_LocImage_y_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_x_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_y_axis = new float [NbBinsHisto];
  
  this->CptCumulativeHistogram_InROI(LocImageROI,NbBinsHisto,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis,1);   //remark: the log of the histogram is actually used
  
  for (i=0;i<NbBinsHisto;i++) CumHisto_TrgImage_x_axis[i]=CumHisto_LocImage_x_axis[NbBinsHisto-1]*static_cast<float>(i)/static_cast<float>(NbBinsHisto-1);
  for (i=0;i<NbBinsHisto;i++) CumHisto_TrgImage_y_axis[i]=static_cast<float>(i)/static_cast<float>(NbBinsHisto-1);
  
  //2) do the job
  this->HistoMatch(NbBinsHisto,CumHisto_TrgImage_x_axis,CumHisto_TrgImage_y_axis,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis);

}




///Compute the histogram. The histogram is normalized (sum of values equal to 1). 
/// -> Input_BinsNb is the number of bins in the histogram
/// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
void ScalarField::CptHistogram(int Input_BinsNb,float * Output_Histo_x_axis,float * Output_Histo_y_axis,int useLogHisto){
  int i,j,k,t;
  float GLmin,GLmax;
  int locNZ,tmpBin;
  float Input_BinsNbFl,halfDeltaBin,VoxelsNb,tmpFl;
  
  
  //1) initiate
  if (this->NZ<=1) locNZ=1;
  else locNZ=this->NZ;
  
  for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]=0;
  
  Input_BinsNbFl=static_cast<float>(Input_BinsNb);
  
  
  //2) get the min and max grey level
  GLmin=this->G(0,0,0);
  GLmax=this->G(0,0,0);
  
  for (t=0;t<this->NT;t++)  for (i=0;i<locNZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++){
    if (GLmin>this->G(k,j,i,t)) GLmin=this->G(k,j,i,t);
    if (GLmax<this->G(k,j,i,t)) GLmax=this->G(k,j,i,t);
    }
    
  //3) fill the histogram x-axis
  halfDeltaBin=((GLmax-GLmin)/Input_BinsNbFl)/2;
  for (i=0;i<Input_BinsNb;i++) Output_Histo_x_axis[i]=halfDeltaBin+GLmin+(GLmax-GLmin)*(static_cast<float>(i)/Input_BinsNbFl);  
  
  //4) fill the histogram y-axis
  for (t=0;t<this->NT;t++)   for (i=0;i<locNZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++){
    tmpBin=static_cast<int>(Input_BinsNbFl*(this->G(k,j,i,t)-GLmin)/(GLmax-GLmin));
    if (tmpBin<0) tmpBin=0; //in case of numerical problems
    if (tmpBin>=Input_BinsNb) tmpBin=Input_BinsNb-1;  //in case of numerical problems
    Output_Histo_y_axis[tmpBin]++;
    }
  
  //5) normalize the histogram y-axis and eventually cpt its log before normalizing
  //5.1) standard case
  
  if (useLogHisto==0) {
    VoxelsNb=static_cast<float>(locNZ*this->NY*this->NX);
    for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]/=VoxelsNb;
  }
  
  //5.2) case where the log of the histogram is considered
  if (useLogHisto!=0) {
    tmpFl=0;
    for (i=0;i<Input_BinsNb;i++){
      Output_Histo_y_axis[i]=log(Output_Histo_y_axis[i]+1);
      tmpFl+=Output_Histo_y_axis[i];
    }
    
    for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]/=tmpFl;
  }
}



///Same as "CptHistogramCptHistogram" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI.
///Note that "ImageROI" must have the same size as "this" and the ROI is defined by non-null values
void ScalarField::CptHistogram_InROI(ScalarField * ImageROI,int Input_BinsNb,float * Output_Histo_x_axis,float * Output_Histo_y_axis,int useLogHisto){
  int i,j,k,t;
  float GLmin,GLmax;
  int locNZ,tmpBin;
  float Input_BinsNbFl,halfDeltaBin,VoxelsNb,tmpFl;
  int ROIexists;
  float epsilon;
  
  //1) initiate
  if (this->NZ<=1) locNZ=1;
  else locNZ=this->NZ;
  
  for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]=0;
  
  Input_BinsNbFl=static_cast<float>(Input_BinsNb);
  
  epsilon=0.0001;
  
  //2) Is there at least one voxel in the ROI
  ROIexists=0;
  
  for (t=0;t<this->NT;t++)  for (i=0;i<locNZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++)
    if (fabs(ImageROI->G(k,j,i,t))>epsilon){
      ROIexists=1;
      GLmin=this->G(k,j,i,t);
      GLmax=this->G(k,j,i,t);
      }
      
  if (ROIexists==0){
    cout << "No ROI defined in the ROI file. A ROI is defined here by non-null values" << endl;
    CptHistogram(Input_BinsNb,Output_Histo_x_axis,Output_Histo_y_axis,useLogHisto);
    return;
    }
    
  
  //2) get the min and max grey level
  for (t=0;t<this->NT;t++)  for (i=0;i<locNZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++) if (fabs(ImageROI->G(k,j,i,t))>epsilon){
    if (GLmin>this->G(k,j,i,t)) GLmin=this->G(k,j,i,t);
    if (GLmax<this->G(k,j,i,t)) GLmax=this->G(k,j,i,t);
    }
    
  //3) fill the histogram x-axis
  halfDeltaBin=((GLmax-GLmin)/Input_BinsNbFl)/2;
  for (i=0;i<Input_BinsNb;i++) Output_Histo_x_axis[i]=halfDeltaBin+GLmin+(GLmax-GLmin)*(static_cast<float>(i)/Input_BinsNbFl);  
  
  //4) fill the histogram y-axis
  for (t=0;t<this->NT;t++)   for (i=0;i<locNZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++) if (fabs(ImageROI->G(k,j,i,t))>epsilon){
    tmpBin=static_cast<int>(Input_BinsNbFl*(this->G(k,j,i,t)-GLmin)/(GLmax-GLmin));
    if (tmpBin<0) tmpBin=0; //in case of numerical problems
    if (tmpBin>=Input_BinsNb) tmpBin=Input_BinsNb-1;  //in case of numerical problems
    Output_Histo_y_axis[tmpBin]++;
    }
  
  //5) normalize the histogram y-axis and eventually cpt its log before normalizing
  //5.1) standard case
  
  if (useLogHisto==0) {
    VoxelsNb=static_cast<float>(locNZ*this->NY*this->NX);
    for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]/=VoxelsNb;
  }
  
  //5.2) case where the log of the histogram is considered
  if (useLogHisto!=0) {
    tmpFl=0;
    for (i=0;i<Input_BinsNb;i++){
      Output_Histo_y_axis[i]=log(Output_Histo_y_axis[i]+1);
      tmpFl+=Output_Histo_y_axis[i];
    }
    
    for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]/=tmpFl;
  }
  
  return;
}





///Compute the cumulative histogram. The cumulative histogram is normalized (last value equal to 1). 
/// -> Input_BinsNb is the number of bins in the histogram
/// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
void ScalarField::CptCumulativeHistogram(int Input_BinsNb,float * Output_CumHisto_x_axis,float * Output_CumHisto_y_axis,int useLogHisto){
  int i;
  float cumulatedvalues;
  
  //compute the histogram
  this->CptHistogram(Input_BinsNb,Output_CumHisto_x_axis,Output_CumHisto_y_axis,useLogHisto);
  
  //compute the cumulative histogram
  cumulatedvalues=0;
  for (i=0;i<Input_BinsNb;i++){
    cumulatedvalues+=Output_CumHisto_y_axis[i];
    Output_CumHisto_y_axis[i]=cumulatedvalues;
  }
}


///Same as "CptCumulativeHistogram" but the voxels taken into account to compute the Look-Up Table are only in the ROI defined by ImageROI.
///Note that "ImageROI" must have the same size as "this" and the ROI is defined by non-null values
void ScalarField::CptCumulativeHistogram_InROI(ScalarField * ImageROI,int Input_BinsNb,float * Output_CumHisto_x_axis,float * Output_CumHisto_y_axis,int useLogHisto){
  int i;
  float cumulatedvalues;
  
  //compute the histogram
  this->CptHistogram_InROI(ImageROI,Input_BinsNb,Output_CumHisto_x_axis,Output_CumHisto_y_axis,useLogHisto);
  
  //compute the cumulative histogram
  cumulatedvalues=0;
  for (i=0;i<Input_BinsNb;i++){
    cumulatedvalues+=Output_CumHisto_y_axis[i];
    Output_CumHisto_y_axis[i]=cumulatedvalues;
  }
}










///get the 'Image 2 World matrix' of an image without loading the image  (load its header only)
void Get_Image2World(char * ImageName,float LocI2W[4][4]){
  nifti_1_header hdr;
  FILE *fp;
  int ret,i;
  int x,y,z,t;
  float a,b,c,d,qfac;

  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  //2) Define the Image to world matrix
  if (hdr.sform_code>0){ //METHOD 3 of nifti1.h
    //cout << "Orientation of image " << ImageName << " opened using method 3 (alternative - normal)" << endl;
    LocI2W[0][0]=hdr.srow_x[0];  LocI2W[0][1]=hdr.srow_x[1];  LocI2W[0][2]=hdr.srow_x[2];  LocI2W[0][3]=hdr.srow_x[3];
    LocI2W[1][0]=hdr.srow_y[0];  LocI2W[1][1]=hdr.srow_y[1];  LocI2W[1][2]=hdr.srow_y[2];  LocI2W[1][3]=hdr.srow_y[3];
    LocI2W[2][0]=hdr.srow_z[0];  LocI2W[2][1]=hdr.srow_z[1];  LocI2W[2][2]=hdr.srow_z[2];  LocI2W[2][3]=hdr.srow_z[3];
    LocI2W[3][0]=0;              LocI2W[3][1]=0;              LocI2W[3][2]=0;              LocI2W[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      //cout << "Orientation of image " << ImageName << " opened using method 2 (normal)" << endl;
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      LocI2W[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); LocI2W[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       LocI2W[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     LocI2W[0][3]=hdr.qoffset_x;
      LocI2W[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     LocI2W[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   LocI2W[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     LocI2W[1][3]=hdr.qoffset_y;
      LocI2W[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     LocI2W[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       LocI2W[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); LocI2W[2][3]=hdr.qoffset_z;
      LocI2W[3][0]=0;                               LocI2W[3][1]=0;                                 LocI2W[3][2]=0;                                    LocI2W[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      LocI2W[0][0]=hdr.pixdim[1];  LocI2W[0][1]=0;             LocI2W[0][2]=0;             LocI2W[0][3]=hdr.qoffset_x;
      LocI2W[1][0]=0;              LocI2W[1][1]=hdr.pixdim[2]; LocI2W[1][2]=0;             LocI2W[1][3]=hdr.qoffset_y;
      LocI2W[2][0]=0;              LocI2W[2][1]=0;             LocI2W[2][2]=hdr.pixdim[3]; LocI2W[2][3]=hdr.qoffset_z;
      LocI2W[3][0]=0;              LocI2W[3][1]=0;             LocI2W[3][2]=0;             LocI2W[3][3]=1;
    }
  }
}




///1.2) scalar field of integers

///constructor
IntScalarField::IntScalarField(void){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
}

///destructor
IntScalarField::~IntScalarField(void){
  if ((this->ScalField!=NULL)&&(this->NX>0)) delete this->ScalField;
}


///put a value
//-> inline function in the .h file

/// add a value
//-> inline function in the .h file


//allocate memory for the field and set all values to zero
void IntScalarField::CreateVoidField(int NBX,int NBY,int NBZ){
  int x,y,z;
  int HasSameSize;

  //check if the size is changed
  if ((this->NX==NBX)&&(this->NY==NBY)&&(this->NZ==NBZ))
    HasSameSize=1;
  else 
    HasSameSize=0;
  
  if ((this->NX!=0)&&(HasSameSize==0)){
    cout << "The size of a non-null scalar field is changed\n";
    delete this->ScalField;
  }
  
  //image size
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NXtY=this->NX*this->NY;
  
  //allocate memory to cast (and eventually transform) the original template and target images
  //    -->  ScalarField[ptSF(x,y,z)]= gray level at (x,y,z)
  if (HasSameSize==0) this->ScalField = new int [this->NXtY*this->NZ];
  
  //set all entries of the field at 0.
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(0,x,y,z);
}




///write a scalar field in a nifti image
void IntScalarField::Write(char * OutputImageName){
  int i;
  int x,y,z;
  nifti_1_header hdr;
  nifti1_extender pad={0,0,0,0};
  FILE *fp;
  int ret;
  int *data=NULL;
  
  
  //1) create the header
  memset((void *)&hdr,0, sizeof(hdr));
  hdr.sizeof_hdr = MIN_HEADER_SIZE;
  hdr.dim[0] = 4;
  hdr.dim[1] = this->NX;
  hdr.dim[2] = this->NY;
  hdr.dim[3] = this->NZ;
  hdr.dim[4] = 1;
  hdr.datatype = NIFTI_TYPE_INT32;  //used to be NIFTI_TYPE_FLOAT32
  hdr.bitpix = 8*sizeof(int);                  //used to be 32
  hdr.qform_code=0; // should ideally be set to 1 but I don't set the values of 'quatern_b', 'quatern_c' and 'quatern_d'
  hdr.pixdim[1] = 1;
  hdr.pixdim[2] = 1;
  hdr.pixdim[3] = 1;
  hdr.pixdim[4] = 1.0;
  hdr.qoffset_x=0;
  hdr.qoffset_y=0;
  hdr.qoffset_z=0;
  hdr.sform_code=1;
  hdr.srow_x[0]=1;  hdr.srow_x[1]=0;  hdr.srow_x[2]=0;  hdr.srow_x[3]=0;
  hdr.srow_y[0]=0;  hdr.srow_y[1]=1;  hdr.srow_y[2]=0;  hdr.srow_y[3]=0;
  hdr.srow_z[0]=0;  hdr.srow_z[1]=0;  hdr.srow_z[2]=1;  hdr.srow_z[3]=0;
  hdr.vox_offset = (float) NII_HEADER_SIZE;
  hdr.scl_inter = 0.0;
  hdr.scl_slope = 1.0;
  hdr.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
  strncpy(hdr.magic, "n+1\0", 4);

  //2) save the image OutputImageName
  //allocate and fill the buffer 
  data = new int [hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]];
  
  i=0;
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(x,y,z);
    i++;
  }
  
  // write first 348 bytes of header  
  fp = fopen(OutputImageName,"wb");
  ret = fwrite(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  // write extender pad and image data  
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4], fp);
  
  fclose(fp);
  
}

///put a the same value at every points of the scalar field
void IntScalarField::PutToAllVoxels(int cste){
  int x,y,z;
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) this->P(cste,x,y,z);
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           2:   FUNCTIONS FOR THE CLASS "VectorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///constructor
VectorField::VectorField(void){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  this->NT=0;
}

///destructor
VectorField::~VectorField(void){
  if ((this->VecField!=NULL)&&(this->NX>0)) delete this->VecField;
}

///put a value
//-> inline function in the .h file

/// add a value
//-> inline function in the .h file

///put the same value at all entries of the vector field
void VectorField::PutToAllVoxels(float cste,int t){
  int i,x,y,z;
  for (i=0;i<3;i++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) { this->P(cste,i,x,y,z,t); }
}


///get a value
//-> inline function in the .h file

///get a value using linear interpolation
float VectorField::G(int IdDirec,float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  
  //values out of the image
  if (x<0.) x=0.0001;
  if (x>=this->NX-1.) x=this->NX-1.0001;
  if (y<0.) y=0.0001;
  if (y>=this->NY-1.) y=this->NY-1.0001;
  if (z<0.) z=0.0001;
  if (z>=this->NZ-1.) z=this->NZ-1.0001;
  if (t<0) t=0;
  if (t>this->NT-1) t=this->NT-1;
  
  //closest entire value
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  //interpolation
  if (this->NZ==1){ //2D IMAGE
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
    
    InterpoGreyLevel= wmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
  }
  else{//3D IMAGE
    wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
    wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
    
    InterpoGreyLevel= wmmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
  }
  
  return InterpoGreyLevel;
}

///same as above
float VectorField::G(int IdDirec,double x,double y,double z,int t){
  return this->G(IdDirec,(float) x,(float) y,(float) z,t);
}




///get the maximum of the absolute values of the vector field
float VectorField::GetMaxAbsVal(int t)
{
  float max=0.0;
  int direc,x,y,z;
  for(direc=0;direc<3;direc++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
  {
    if(max<abs(this->G(direc,x,y,z,t))){max = abs(this->G(direc,x,y,z,t));}
  }
  return max;
}


///read a vector field (in 3 nifti images -> X, Y, Z)
void VectorField::Read(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z){
  ScalarField VecField_X;
  ScalarField VecField_Y;
  ScalarField VecField_Z;
  int x,y,z,t;
  int DifferentSize;
  
  //read the vector fields
  VecField_X.Read(NameVecField_X);
  VecField_Y.Read(NameVecField_Y);
  VecField_Z.Read(NameVecField_Z);
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX==0){
    DifferentSize=1;
    }
  else{
    if ((this->NX!=VecField_X.NX)||(this->NY!=VecField_X.NY)||(this->NZ!=VecField_X.NZ)||(this->NT!=VecField_X.NT)){
      DifferentSize=1;
      delete this->VecField;
      }
    else{
      DifferentSize=0;
      }
    }
  
  //fill the parameters of the class and allocate the memory for the image
  this->NX=VecField_X.NX;
  this->NY=VecField_X.NY;
  this->NZ=VecField_X.NZ;
  this->NT=VecField_X.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXtYtZtT=this->NXtYtZ*this->NT;
  if (DifferentSize==1) this->VecField = new float [this->NXtYtZtT*3];
  this->Image2World[0][0]=VecField_X.Image2World[0][0];  this->Image2World[0][1]=VecField_X.Image2World[0][1];  this->Image2World[0][2]=VecField_X.Image2World[0][2];  this->Image2World[0][3]=VecField_X.Image2World[0][3];
  this->Image2World[1][0]=VecField_X.Image2World[1][0];  this->Image2World[1][1]=VecField_X.Image2World[1][1];  this->Image2World[1][2]=VecField_X.Image2World[1][2];  this->Image2World[1][3]=VecField_X.Image2World[1][3];
  this->Image2World[2][0]=VecField_X.Image2World[2][0];  this->Image2World[2][1]=VecField_X.Image2World[2][1];  this->Image2World[2][2]=VecField_X.Image2World[2][2];  this->Image2World[2][3]=VecField_X.Image2World[2][3];
  this->Image2World[3][0]=VecField_X.Image2World[3][0];  this->Image2World[3][1]=VecField_X.Image2World[3][1];  this->Image2World[3][2]=VecField_X.Image2World[3][2];  this->Image2World[3][3]=VecField_X.Image2World[3][3];

  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //cast the image to the format used in the class
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(VecField_X.G(x,y,z,t),0,x,y,z,t);
  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(VecField_Y.G(x,y,z,t),1,x,y,z,t);
  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(VecField_Z.G(x,y,z,t),2,x,y,z,t);
}


///read a scalar vector and perform linear interpolation to give it a specific size
void VectorField::Read_and_Interpolate(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z,int NBX,int NBY,int NBZ,int rescaleVF){
  ScalarField OrigSF;
  int x,y,z,t;
  float x2,y2,z2;
  float scaleFactor;
  
  //0) init
  scaleFactor=1.;
  
  //1) X DIRECTION
  //1.1) read the scalar field in direction X at the original format
  OrigSF.Read(NameVecField_X);
  
  //1.2) fill the parameters of the class and allocate the memory for the vector field
  //(the directions Y and Z are supposed in the same format)
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=OrigSF.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXtYtZtT=this->NXtYtZ*this->NT;
  this->VecField = new float [this->NXtYtZtT*3];
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;

  
  //1.3) manage the scale factor
  if (rescaleVF!=0) scaleFactor=((float)this->NX)/((float)OrigSF.NX);
  this->Image2World[0][0]=OrigSF.Image2World[0][0]/scaleFactor;  this->Image2World[0][1]=OrigSF.Image2World[0][1]/scaleFactor;  this->Image2World[0][2]=OrigSF.Image2World[0][2]/scaleFactor;  this->Image2World[0][3]=OrigSF.Image2World[0][3];
  
  //1.4) interpolate the original image to compute the vector field in direction X
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,0,x,y,z,t);
        }
      }
    }
  }

  //2) Y DIRECTION
  //2.1) read the scalar field in direction Y
  OrigSF.Read(NameVecField_Y);
  
  //2.2) manage the scale factor
  if (rescaleVF!=0) scaleFactor=((float)this->NY)/((float)OrigSF.NY);
  this->Image2World[1][0]=OrigSF.Image2World[1][0]/scaleFactor;  this->Image2World[1][1]=OrigSF.Image2World[1][1]/scaleFactor;  this->Image2World[1][2]=OrigSF.Image2World[1][2]/scaleFactor;  this->Image2World[1][3]=OrigSF.Image2World[1][3];
  
  //2.3) interpolate the original image to compute the vector field in direction Y
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,1,x,y,z,t);
        }
      }
    }
  }
  
  
  //3) Z DIRECTION
  //3.1) read the scalar field in direction Z
  OrigSF.Read(NameVecField_Z);
  
  //3.2) manage the scale factor
  if (rescaleVF!=0) scaleFactor=((float)this->NZ)/((float)OrigSF.NZ);
  this->Image2World[2][0]=OrigSF.Image2World[2][0]/scaleFactor;  this->Image2World[2][1]=OrigSF.Image2World[2][1]/scaleFactor;  this->Image2World[2][2]=OrigSF.Image2World[2][2]/scaleFactor;  this->Image2World[2][3]=OrigSF.Image2World[2][3];
  
  //3.3) interpolate the original image to compute the vector field in direction Z
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,2,x,y,z,t);
        }
      }
    }
  }
  
  //compute the invert World2Image matrix
  invert_4t4quaternion(this->Image2World,this->World2Image);

  
}



///create a void vector field. No message is printed if  Verbose!=1
void VectorField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT,int Verbose){
  int x,y,z,t,direc;
  int HasSameSize;
  
  //check if the size is changed
  if ((this->NX==NBX)&&(this->NY==NBY)&&(this->NZ==NBZ)&&(this->NT==NBT))
    HasSameSize=1;
  else 
    HasSameSize=0;
  
  if ((this->NX!=0)&&(HasSameSize==0)){
    if (Verbose==1) cout << "The size of a non-null vector field is changed\n";
    delete this->VecField;
  }
  
  //image size
  if (HasSameSize==0){
    this->NX=NBX;
    this->NY=NBY;
    this->NZ=NBZ;
    this->NT=NBT;
    this->NXtY=this->NX*this->NY;
    this->NXtYtZ=this->NXtY*this->NZ;
    this->NXtYtZtT=this->NXtYtZ*this->NT;
  }
  
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //allocate memory to cast (and eventually transform) the original template and target images
  if (HasSameSize==0){
    this->VecField = new float [this->NXtYtZtT*3];
    }
  
  //fill the vector field with zeros
  for(direc=0;direc<3;direc++) for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(0.,direc,x,y,z,t);
}



///for each point (x,y,z) and direction direc: this(direc,x,y,z,0) <- \sum_t TimeDepVF(direc,x,y,z,t)
///TimeDepVF must be allocated and have the same size as this. If yes, 1 is returned / 0 is returned otherwise 
int VectorField::ProjectTimeDepVF_to_StationaryVF(VectorField * TimeDepVF){
  int x,y,z,t,direc;
  int doTheJob;
  
  if ((this->NX==TimeDepVF->NX)&&(this->NY==TimeDepVF->NY)&&(this->NZ==TimeDepVF->NZ)) doTheJob=1;
  else doTheJob=0;
  
  if (doTheJob==1){
    for(direc=0;direc<3;direc++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) 
      this->P(0.,direc,x,y,z);
    
    for(direc=0;direc<3;direc++) for(t=0;t<TimeDepVF->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) 
      this->Add(TimeDepVF->G(direc,x,y,z,t),direc,x,y,z);
  }
  return doTheJob;
}


///Copy by reference the vectorized velocity field -> Access to a value: IdDirec*NX*NY*NZ*NT + t*NX*NY*NZ + z*NX*NY + y*NX + x
///The whole vector size is also given in VectorSize
float * VectorField::CopyByRefVectorizedVF(int * VectorSize){
  *VectorSize=this->NXtYtZtT*3;
  return this->VecField;
}



///Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
void VectorField::SlashFieldSize(int verbative){

  //image size
  this->NX=1;
  this->NY=1;
  this->NZ=1;
  this->NT=1;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXtYtZtT=this->NXtYtZ*this->NT;
  
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  if (verbative==1) cout << "Slash some memory" << endl;

  //re-allocate memory
  if (this->NX!=0) delete this->VecField;
  
  this->VecField = new float [3];
}



///write a vector field (from 3 nifti images -> X, Y, Z)
void VectorField::Write(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z){
  ScalarField OutputImage;
  int x,y,z,t;
  
  //image to cast
  OutputImage.CreateVoidField(this->NX, this->NY, this->NZ,this->NT);
  
  //write the vector field in X direction
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    OutputImage.P(this->G(0,x,y,z,t),x,y,z,t);
  
  OutputImage.Write(NameVecField_X);
  
  //write the vector field in Y direction
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    OutputImage.P(this->G(1,x,y,z,t),x,y,z,t);
  
  OutputImage.Write(NameVecField_Y);
  
  //write the vector field in Z direction
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    OutputImage.P(this->G(2,x,y,z,t),x,y,z,t);
  
  OutputImage.Write(NameVecField_Z);
}



///write a vector field (in 3 nifti images -> X, Y, Z)
//The 4th file is an input file containing the headers 
void VectorField::Write(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z, char * ImageForHeaderName){
  nifti_1_header hdr_ref;
  FILE *fp_header;
  int i;
  int x,y,z,t;
  nifti1_extender pad={0,0,0,0};
  FILE *fp;
  int ret;
  float *data=NULL;
  
  
  //1) read the header of ImageForHeaderName
  fp_header = fopen(ImageForHeaderName,"rb");
  fread(&hdr_ref, MIN_HEADER_SIZE, 1, fp_header);
  fclose(fp_header);
  
  //2)give image parameters to the header
  hdr_ref.dim[0] = 4;
  hdr_ref.dim[1] = this->NX;
  hdr_ref.dim[2] = this->NY;
  hdr_ref.dim[3] = this->NZ;
  hdr_ref.dim[4] = this->NT;
  hdr_ref.datatype = NIFTI_TYPE_FLOAT32;
  hdr_ref.bitpix = 32; 
  hdr_ref.scl_inter = 0.0;
  hdr_ref.scl_slope = 1.0;
  hdr_ref.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
  
  
  //3) allocate memory for the buffer
  //data = (float *) malloc(sizeof(float) * hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]);
  data = new float  [hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]];
  
  //4.1) write the vector field in X direction
  //... fill the buffer
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(0,x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  //... save the field
  fp = fopen(NameVecField_X,"wb");
  
  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  fclose(fp);
  

  //4.2) write the vector field in Y direction
  //... fill the buffer
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(1,x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  //... save the field
  fp = fopen(NameVecField_Y,"wb");
  
  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  fclose(fp);
  
  //4.3) write the vector field in Z direction
  //... fill the buffer
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(2,x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  //... save the field
  fp = fopen(NameVecField_Z,"wb");
  
  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  fclose(fp);
  
  //5) free the memory
  //free(data);
  delete data;

}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           3:   FUNCTIONS FOR THE CLASS "TensorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



///constructor
TensorField::TensorField(void){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  this->NT=0;
}

///destructor
TensorField::~TensorField(void){
  if ((this->TField!=NULL)&&(this->NX>0)) delete this->TField;
}


///constructor
void TensorField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT){
  int x,y,z,t,direc1,direc2;
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=NBX)||(this->NY!=NBY)||(this->NZ!=NBZ)||(this->NT!=NBT))
      cout << "The size of a tensor field is changed\n";
  
  //image size
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=NBT;
  this->NXt9=this->NX*9;
  this->NXtYt9=this->NX*NY*9;
  this->NXtYtZt9=this->NX*NY*NZ*9;
  
  
  //allocate memory to cast (and eventually transform) the original template and target images
  this->TField = new float [this->NXtYtZt9*NT];
  
  //fill the image with 0.
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) for(direc2=0;direc2<3;direc2++) for(direc1=0;direc1<3;direc1++) 
    this->P(0.,direc1,direc2,x,y,z,t);
}

///put a value
void TensorField::P(float value,int IdDirec1,int IdDirec2,int x,int y,int z,int t){
  this->TField[IdDirec1+3*IdDirec2+9*x+this->NXt9*y+this->NXtYt9*z+this->NXtYtZt9*t]=value;
}

///add a value
void TensorField::Add(float value,int IdDirec1,int IdDirec2,int x,int y,int z,int t){
  this->TField[IdDirec1+3*IdDirec2+9*x+this->NXt9*y+this->NXtYt9*z+this->NXtYtZt9*t]+=value;
}


///get a value
float TensorField::G(int IdDirec1,int IdDirec2,int x,int y,int z,int t){
  return this->TField[IdDirec1+3*IdDirec2+9*x+this->NXt9*y+this->NXtYt9*z+this->NXtYtZt9*t];
}

///Add a tensorised vector to the existing tensor
void TensorField::AddTensorisedVector(float vec[3],int x,int y,int z,int t,float weight){
  this->Add(weight*vec[0]*vec[0],0,0,x,y,z,t);
  this->Add(weight*vec[1]*vec[0],1,0,x,y,z,t);
  this->Add(weight*vec[2]*vec[0],2,0,x,y,z,t);
  this->Add(weight*vec[0]*vec[1],0,1,x,y,z,t);
  this->Add(weight*vec[1]*vec[1],1,1,x,y,z,t);
  this->Add(weight*vec[2]*vec[1],2,1,x,y,z,t);
  this->Add(weight*vec[0]*vec[2],0,2,x,y,z,t);
  this->Add(weight*vec[1]*vec[2],1,2,x,y,z,t);
  this->Add(weight*vec[2]*vec[2],2,2,x,y,z,t);
}

///put all values to zero
void TensorField::PutAllValuesToZero(){
  int x,y,z,t,direc1,direc2;
  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) for(direc2=0;direc2<3;direc2++) for(direc1=0;direc1<3;direc1++) 
    this->P(0.,direc1,direc2,x,y,z,t);
}

///Perform a principal component analysis of the 3*3 tensor
///The outputs are:
/// lambda1,lambda2,lambda3: the eigenvalues in decrasing order
/// vec1: 1st eigenvector  (must be initialised as vec1[3])
/// vec2: 2nd eigenvector  (must be initialised as vec2[3])
/// vec3: 3rd eigenvector  (must be initialised as vec3[3])
void TensorField::PCA(float * vec1,float * vec2,float * vec3,float * lambda1, float * lambda2, float * lambda3, int x,int y,int z,int t){
  float ** a;
  float ** q;
  float *  d;
  int i;
  
  //alloc temporary variables
  a = new float * [3];
  for(i=0;i<3;i++) a[i]=new float [3];
  q = new float * [3];
  for(i=0;i<3;i++) q[i]=new float [3];
  d = new float [3];
  
  //copy input values
  a[0][0]=this->G(0,0,x,y,z,t);  a[0][1]=this->G(0,1,x,y,z,t);  a[0][2]=this->G(0,2,x,y,z,t);
  a[1][0]=this->G(1,0,x,y,z,t);  a[1][1]=this->G(1,1,x,y,z,t);  a[1][2]=this->G(1,2,x,y,z,t);
  a[2][0]=this->G(2,0,x,y,z,t);  a[2][1]=this->G(2,1,x,y,z,t);  a[2][2]=this->G(2,2,x,y,z,t);
  
  
  //eigenvalue decomposition
  jacobi3(a,d,q);
  
  //copy output values
  *lambda1=d[0];  *lambda2=d[1];  *lambda3=d[2];
  
  vec1[0]=q[0][0];  vec1[1]=q[1][0];  vec1[2]=q[2][0];
  vec2[0]=q[0][1];  vec2[1]=q[1][1];  vec2[2]=q[2][1];
  vec3[0]=q[0][2];  vec3[1]=q[1][2];  vec3[2]=q[2][2];
  
  //delete allocated variables
  for(i=0;i<3;i++) delete a[i];
  delete a;
  for(i=0;i<3;i++) delete q[i];
  delete q;
  delete d;
  
}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           4: CLASS TO PERFORM CONVOLUTION AND DECONVOLUTION USING FFT
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///4.1: Main class

///Constructor
FFTconvolver3D::FFTconvolver3D(){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  
  this->NXfft=0;
  this->NYfft=0;
  this->NZfft=0;
}

///destructor
FFTconvolver3D::~FFTconvolver3D(void){}


///Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
///4 Gaussians (set some weights to 0 if less Gaussians are required)
/// * NX, NY, NZ: is the size of the input image
/// * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
/// * w2,sX2,sY2,sZ2,: weight of the 2nd Gaussian kernel and std. dev. in direction X, Y, Z
/// * w3,sX3,sY3,sZ3,: weight of the 3rd Gaussian kernel and std. dev. in direction X, Y, Z
/// * w4,sX4,sY4,sZ4,: weight of the 4th Gaussian kernel and std. dev. in direction X, Y, Z
/// * w5,sX5,sY5,sZ5,: weight of the 5th Gaussian kernel and std. dev. in direction X, Y, Z
/// * w6,sX6,sY6,sZ6,: weight of the 6th Gaussian kernel and std. dev. in direction X, Y, Z
/// * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void FFTconvolver3D::InitiateConvolver(int NBX,int NBY, int NBZ,float w1,float sX1,float sY1,float sZ1,float w2,float sX2,float sY2,float sZ2,float w3,float sX3,float sY3,float sZ3,float w4,float sX4,float sY4,float sZ4,float w5,float sX5,float sY5,float sZ5,float w6,float sX6,float sY6,float sZ6,float w7,float sX7,float sY7,float sZ7){
  
  //set the size of the original image
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  
  //set the size of images for the FFT
  this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
  this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
  this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
  

  //cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << "\n";
  
  //allocate memory for the images for the FFT
  this->RealSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - real part
  this->ImagSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - imaginary part
  this->RealFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - real part
  this->ImagFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - imaginary part
  
  //allocate memory for the temporary image
  this->ImageTemp.CreateVoidField(this->NXfft,this->NYfft,this->NZfft);
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf7AnisotropicGaussianFilters(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7);
}



void FFTconvolver3D::InitiateMatrixValuedConvolver(int NBX,int NBY, int NBZ){
  
  //set the size of the original image
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  
  //set the size of images for the FFT
  this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
  this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
  this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
  

  //cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << "\n";
  
  //allocate memory for the images for the FFT
  this->RealSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - real part
  this->ImagSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - imaginary part
  this->RealFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - real part
  this->ImagFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - imaginary part
  
  //allocate memory for the temporary vector field
  this->VecFieldTemp.CreateVoidField(this->NXfft,this->NYfft,this->NZfft);
}


///Initiate the complex fields for the FFT and the smoothing kernels being parametrized by B-splines
void FFTconvolver3D::InitiateSplineConvolver(int NBX,int NBY, int NBZ,int GridStep){
  
  //set the size of the original image
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  
  //set the size of images for the FFT
  this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
  this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
  this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
  

  //cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << "\n";
  
  //allocate memory for the images for the FFT
  this->RealSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - real part
  this->ImagSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - imaginary part
  this->RealFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - real part
  this->ImagFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - imaginary part
  
  //allocate memory for the temporary image
  this->ImageTemp.CreateVoidField(this->NXfft,this->NYfft,this->NZfft);
  
  //define the kernel and transform it in Fourier spaces
  this->MakeBSplineFilter(GridStep);
}



///change the kernel of the convolver (same notations as the constructor)
///... here the new kernel is normalized (w1+...+w7 is normalized to 1)
void FFTconvolver3D::ChangeKernel(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,
                                  float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,
                                  float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,
                                  float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,
                                  float weight5,float sigmaX5,float sigmaY5,float sigmaZ5,
                                  float weight6,float sigmaX6,float sigmaY6,float sigmaZ6,
                                  float weight7,float sigmaX7,float sigmaY7,float sigmaZ7){
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf7AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,weight2,sigmaX2,sigmaY2,sigmaZ2,weight3,sigmaX3,sigmaY3,sigmaZ3,weight4,sigmaX4,sigmaY4,sigmaZ4,weight5,sigmaX5,sigmaY5,sigmaZ5,weight6,sigmaX6,sigmaY6,sigmaZ6,weight7,sigmaX7,sigmaY7,sigmaZ7);
}

///change the kernel of the convolver (same notations as the constructor)
///... here the new kernel is not normalized
void FFTconvolver3D::ChangeKernel_SingleScale(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1){
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf7AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0);
}



///design a kernel that is the sum of up to 7 Gaussians and transform it in Fourier spaces
//if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
void FFTconvolver3D::MakeSumOf7AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,float weight5,float sigmaX5,float sigmaY5,float sigmaZ5,float weight6,float sigmaX6,float sigmaY6,float sigmaZ6,float weight7,float sigmaX7,float sigmaY7,float sigmaZ7,int NormalizeWeights){
  int k,x,y,z;
  float SumLoc;
  float weight,sigmaX,sigmaY,sigmaZ,sumWeight;
  
  //set RealFilterForFFT and ImagFilterForFFT to 0 in case it contains something
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(0.,x,y,z);
  
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT.P(0.,x,y,z);

  sumWeight=0;
  
  //compute and save the 7 kernels
  for (k=0;k<7;k++){
    //parameters of the current kernel
    if (k==6){weight=weight7; sigmaX=sigmaX7; sigmaY=sigmaY7; sigmaZ=sigmaZ7;}
    if (k==5){weight=weight6; sigmaX=sigmaX6; sigmaY=sigmaY6; sigmaZ=sigmaZ6;}
    if (k==4){weight=weight5; sigmaX=sigmaX5; sigmaY=sigmaY5; sigmaZ=sigmaZ5;}
    if (k==3){weight=weight4; sigmaX=sigmaX4; sigmaY=sigmaY4; sigmaZ=sigmaZ4;}
    if (k==2){weight=weight3; sigmaX=sigmaX3; sigmaY=sigmaY3; sigmaZ=sigmaZ3;}
    if (k==1){weight=weight2; sigmaX=sigmaX2; sigmaY=sigmaY2; sigmaZ=sigmaZ2;}
    if (k==0){weight=weight1; sigmaX=sigmaX1; sigmaY=sigmaY1; sigmaZ=sigmaZ1;}
    
    sumWeight+=weight;
    
    for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImageTemp.P(0,x,y,z);
    
    //design the current kernel with no influence of the weight
    for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp( -(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp( -(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp(-(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp(-(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    
    //normalization of the current filter
    SumLoc=0.;
    for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) SumLoc+=this->ImageTemp.G(x,y,z);
    
    //copy in RealFilterForFFT
    if (fabs(SumLoc)>0.01){
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.Add(weight*this->ImageTemp.G(x,y,z)/SumLoc,x,y,z);
    }
    else{
      cout << "One kernel appears to be null or almost null" << endl;
    }
  }
  
  //normalize the weights
  if (NormalizeWeights!=0) for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(this->RealFilterForFFT.G(x,y,z)/sumWeight,x,y,z);

  //Transform RealFilterForFFT and ImagFilterForFFT in Fourier spaces
  this->DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
}


///Design a kernel representing a 3D second order B-Spline
///StepSize controls the extent of the kernel  (step between two nodes in voxels)
void FFTconvolver3D::MakeBSplineFilter(int StepSize){
  int k,x,y,z;
  float SumLoc;
  float distX,distY,distZ,FlStepSize;
  
  //1) set RealFilterForFFT, ImagFilterForFFTn and ImageTemp to 0 in case it contains something
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImageTemp.P(0,x,y,z);
  
  //float BSplineWeight(float x);
  
  //2) design the kernel
  FlStepSize=static_cast<float>(StepSize);
  for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++){
    distX=(x)/FlStepSize; distY=(y)/FlStepSize; distZ=(z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
    
  for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
    distX=(this->NXfft-x)/FlStepSize; distY=(y)/FlStepSize; distZ=(z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
  
  for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++){
    distX=(x)/FlStepSize; distY=(this->NYfft-y)/FlStepSize; distZ=(z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
  
  for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
    distX=(this->NXfft-x)/FlStepSize; distY=(this->NYfft-y)/FlStepSize; distZ=(z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
  
  for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++){
    distX=(x)/FlStepSize; distY=(y)/FlStepSize; distZ=(this->NZfft-z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
  
  for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
    distX=(this->NXfft-x)/FlStepSize; distY=(y)/FlStepSize; distZ=(this->NZfft-z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
  
  for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++){
    distX=(x)/FlStepSize; distY=(this->NYfft-y)/FlStepSize; distZ=(this->NZfft-z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
  
  for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
    distX=(this->NXfft-x)/FlStepSize; distY=(this->NYfft-y)/FlStepSize; distZ=(this->NZfft-z)/FlStepSize;
    this->ImageTemp.P(BSplineWeight(distX)*BSplineWeight(distY)*BSplineWeight(distZ),x,y,z);
  }
  
  
  //3) normalization of the current filter
  SumLoc=0.;
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) SumLoc+=this->ImageTemp.G(x,y,z);
  
  //4) copy in RealFilterForFFT
  if (fabs(SumLoc)>0.001){
    for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(this->ImageTemp.G(x,y,z)/SumLoc,x,y,z);
  }
  else{
    cout << "The kernel appears to be null or almost null" << endl;
  }
  

  //5) Transform RealFilterForFFT and ImagFilterForFFT in Fourier spaces
  this->DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
}




///Fast Fourier Transform of numerical recipies (slighly modified)
void FFTconvolver3D::four1NR(float data[], unsigned long nn, int isign){
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2){
    if (j>i){
      tempr=data[j]; data[j]=data[i]; data[i]=tempr;
      tempr=data[j+1]; data[j+1]=data[i+1]; data[i+1]=tempr;
    }
    m=n >> 1;
    while ((m>=2) && (j>m)){
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}



#ifdef COMPILE_WITH_OPENMP

///Fast Fourier Transform
void FFTconvolver3D::DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  int MaxSizeXSizeYSizeZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  MaxSizeXSizeYSizeZ=SizeX;
  if (SizeY>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeY;
  if (SizeZ>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,dataX) 
  {
    dataX = new float [MaxSizeXSizeYSizeZ*2+1];
  
    //2) perform the fft along x axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){
      for (z = 0; z < SizeZ; z++){
        for (x = 0; x < SizeX; x++){
          dataX[2*x+1]=RealSignal->G(x, y, z);
          dataX[2*x+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeX, 1);
        for (x = 0; x < SizeX; x++){
          RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
          ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX, x, y, z);
        }  
      }
    }

    //3) perform the fft along y axis
    #pragma omp for
    for (x = 0; x < SizeX; x++){ 
      for (z = 0; z < SizeZ; z++) {
        for (y = 0; y < SizeY; y++){
          dataX[2*y+1]=RealSignal->G(x, y, z);
          dataX[2*y+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeY, 1);
        for (y = 0; y < SizeY; y++){
          RealSignal->P(dataX[2*y+1]/SqrtSizeY,x, y, z);
          ImaginarySignal->P(dataX[2*y+2]/SqrtSizeY, x, y, z);
        }
      }
    }
  
    //4) perform the fft along z axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){
      for (x = 0; x < SizeX; x++){
        for (z = 0; z < SizeZ; z++){
          dataX[2*z+1]=RealSignal->G(x, y, z);
          dataX[2*z+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeZ, 1);
        for (z = 0; z < SizeZ; z++){
          RealSignal->P(dataX[2*z+1]/SqrtSizeZ,x, y, z);
          ImaginarySignal->P(dataX[2*z+2]/SqrtSizeZ, x, y, z);
        }
      }
    }

    delete dataX;
  //END FORK FOR THREADS
  }

}

#else

///Fast Fourier Transform
void FFTconvolver3D::DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  
  //2) perform the fft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->G(x, y, z);
      dataX[2*x+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataX, (unsigned long)SizeX, 1);
    for (x = 0; x < SizeX; x++){
      RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
      ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX, x, y, z);
    }
  }
  delete dataX;
  
  //3) perform the fft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->G(x, y, z);
      dataY[2*y+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataY, (unsigned long)SizeY, 1);
    for (y = 0; y < SizeY; y++){
      RealSignal->P(dataY[2*y+1]/SqrtSizeY,x, y, z);
      ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
    }
  }
  delete dataY;
  
  
  //4) perform the fft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->G(x, y, z);
      dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataZ, (unsigned long)SizeZ, 1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->P(dataZ[2*z+1]/SqrtSizeZ,x, y, z);
      ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ, x, y, z);
    }
  }
  delete dataZ;
}

#endif



#ifdef COMPILE_WITH_OPENMP


///Inverse Fast Fourier Transform
void FFTconvolver3D::InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  int MaxSizeXSizeYSizeZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  MaxSizeXSizeYSizeZ=SizeX;
  if (SizeY>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeY;
  if (SizeZ>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,dataX) 
  {
    dataX = new float [MaxSizeXSizeYSizeZ*2+1];
  
    //2) perform the ifft along z axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){ 
      for (x = 0; x < SizeX; x++){
        for (z = 0; z < SizeZ; z++){
          dataX[2*z+1]=RealSignal->G(x, y, z);
          dataX[2*z+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeZ, -1);
        for (z = 0; z < SizeZ; z++){
          RealSignal->P(dataX[2*z+1]/SqrtSizeZ, x, y, z);
          ImaginarySignal->P(dataX[2*z+2]/SqrtSizeZ,x, y, z);
        }
      }
    }
    
    //3) perform the ifft along y axis
    #pragma omp for
    for (x = 0; x < SizeX; x++){
      for (z = 0; z < SizeZ; z++){
        for (y = 0; y < SizeY; y++){
          dataX[2*y+1]=RealSignal->G(x, y, z);
          dataX[2*y+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeY, -1);
        for (y = 0; y < SizeY; y++){
          RealSignal->P(dataX[2*y+1]/SqrtSizeY, x, y, z);
          ImaginarySignal->P(dataX[2*y+2]/SqrtSizeY, x, y, z);
        }
      }
    }
    
    //4) perform the ifft along x axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){
      for (z = 0; z < SizeZ; z++){
        for (x = 0; x < SizeX; x++){
          dataX[2*x+1]=RealSignal->G(x, y, z);
          dataX[2*x+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeX, -1);
        for (x = 0; x < SizeX; x++){
          RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
          ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX,x, y, z);
        }
      }
    }
    
    delete dataX;
  //END FORK FOR THREADS
  }
}


#else

///Inverse Fast Fourier Transform
void FFTconvolver3D::InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  
  //2) perform the ifft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->G(x, y, z);
      dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataZ, (unsigned long)SizeZ, -1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->P(dataZ[2*z+1]/SqrtSizeZ, x, y, z);
      ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ,x, y, z);
    }
  }
  delete dataZ;
  
  //3) perform the ifft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->G(x, y, z);
      dataY[2*y+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataY, (unsigned long)SizeY, -1);
    for (y = 0; y < SizeY; y++){
      RealSignal->P(dataY[2*y+1]/SqrtSizeY, x, y, z);
      ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
    }
  }
  delete dataY;
  
  //4) perform the ifft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->G(x, y, z);
      dataX[2*x+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataX, (unsigned long)SizeX, -1);
    for (x = 0; x < SizeX; x++){
      RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
      ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX,x, y, z);
    }
  }
  delete dataX;
}

#endif


///convolution of a 3D scalar field using the predifined kernel
void FFTconvolver3D::Convolution(ScalarField * SF){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  //1) Copy the orginal image in the image that will be transformed
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
  
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++) this->RealSignalForFFT.P(SF->G(x,y,z),x,y,z);


  //2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
  this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //3) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
  
  for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
    a=this->RealSignalForFFT.G(x, y, z);
    b=this->ImagSignalForFFT.G(x, y, z);
    c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
    d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
    
    this->RealSignalForFFT.P(a*c-b*d, x, y, z);
    this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
  }
  
  //4) IFFT
  this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //5) Copy the image that has been convolved in the orginal image
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
    SF->P(this->RealSignalForFFT.G(x,y,z),x,y,z);
  }
}




///convolution of a 3D vector field using a divergence free (matrix valued) kernel, based on a Gaussian kernel of std dev sigma
///   ->    \phi_{\alpha}(x) = \{ [ 2 (s-1) \alpha -  4 \alpha^2 ||x||^2  ] Id  + 4 \alpha^2 xx' \}  e^{\alpha ||x||^2}  ... where s is the image dimension  
///   ->    F. J. Narcowich and J. D. Ward,   Generalized Hermite interpolation via matrix-valued conditionally positive definite functions, Math. Comp. 63 (1994), 661-687
void FFTconvolver3D::Convolution_DivergenceFreeFilter(VectorField * VF, float sigma){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  float normSq;
  float Xsq,Ysq,Zsq,XtY,XtZ,YtZ;
  float t1,t2,alpha,mt1;
  int ImDim;
  int xtp,ytp,ztp;
  
  //0) init
  if (this->ImageTemp.NZ>1) ImDim=3;
  else  ImDim=2;
  
  alpha=1/(2*sigma*sigma);
  t1=2*(ImDim-1)*alpha;
  t2=4*alpha*alpha;
  
  for (i=0;i<ImDim;i++) for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->VecFieldTemp.P(0,i,x,y,z);
  
  //2) matrix/vector convolution
  for (i=0;i<ImDim;i++){ // loop on the dimension of the vector field VF  (and the corresponding dimension in the matrix valued kernel)
    for (j=0;j<ImDim;j++){//loop on the second dimension of the matrix valued kernel

      //2.1) Copy the orginal image in the image that will be transformed
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
      
      for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++) this->RealSignalForFFT.P(VF->G(i,x,y,z),x,y,z);

      //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
      
      //2.3) Copy the filter in the images that will be transformed
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(0.,x,y,z);
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT.P(0.,x,y,z);
      
      //\phi_{\alpha}(x) = \{ [ 2 (s-1) \alpha -  4 \alpha^2 ||x||^2  ] Id  + 4 \alpha^2 xx' \}  e^{- \alpha ||x||^2}
      for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++){
  xtp=x; ytp=y; ztp=z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      
      for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
  xtp=this->NXfft-x; ytp=y; ztp=z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++){
  xtp=x; ytp=this->NYfft-y; ztp=z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
  xtp=this->NXfft-x; ytp=this->NYfft-y; ztp=z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++){
  xtp=x; ytp=y; ztp=this->NZfft-z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
  xtp=this->NXfft-x; ytp=y; ztp=this->NZfft-z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++){
  xtp=x; ytp=this->NYfft-y; ztp=this->NZfft-z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++){
  xtp=this->NXfft-x; ytp=this->NYfft-y; ztp=this->NZfft-z;
        Xsq=(float)(xtp*xtp);  Ysq=(float)(ytp*ytp);  Zsq=(float)(ztp*ztp);  XtY=(float)(xtp*ytp); XtZ=(float)(xtp*ztp); YtZ=(float)(ytp*ztp); 
  normSq=Xsq+Ysq+Ysq;
        if (i==0) { if (j==0){mt1=t2*(Xsq-normSq)+t1;}  else if (j==1){mt1=t2*XtY;}                else if (j==2){mt1=t2*XtZ;} }
        if (i==1) { if (j==0){mt1=t2*XtY;}                else if (j==1){mt1=t2*(Ysq-normSq)+t1;}  else if (j==2){mt1=t2*YtZ;} }
        if (i==2) { if (j==0){mt1=t2*XtZ;}                else if (j==1){mt1=t2*YtZ;}                else if (j==2){mt1=t2*(Zsq-normSq)+t1;} }
        this->RealFilterForFFT.P(mt1*exp(-alpha*normSq),x,y,z);
        //if (i==j) this->RealFilterForFFT.P(exp(-alpha*normSq),x,y,z); else this->RealFilterForFFT.P(0,x,y,z);
      }
      
      
      //2.4) Transform RealFilterForFFT and ImagFilterForFFT in Fourier spaces
      this->DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
      
      
      //2.5) filtering in Fourier spaces
      CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
      
      for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
        a=this->RealSignalForFFT.G(x, y, z);
        b=this->ImagSignalForFFT.G(x, y, z);
        c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
        d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
        
        this->RealSignalForFFT.P(a*c-b*d, x, y, z);
        this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
      }
      
      //2.6) IFFT
      this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
      
      //2.7) Update the velocity in the current direction
      for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++)
        this->VecFieldTemp.Add(this->RealSignalForFFT.G(x,y,z),i,x,y,z);
    }
  }
  
  //3) Copy the velocity (in direction i) that has been treated in the orginal velocity field
  for (i=0;i<ImDim;i++) for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++){
    VF->P(this->VecFieldTemp.G(i,x,y,z),i,x,y,z);
  }
}





///convolution of the real scalar field defined inside of the class
void FFTconvolver3D::Convolution(){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  //1) Set to 0. all values that cannot be accessed by outside of the class
  for (z=this->NZ;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++)        for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZ;z++)           for (y=this->NY;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZ;z++)           for (y=0;y<this->NY;y++)           for (x=this->NX;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) 
    this->ImagSignalForFFT.P(0.,x,y,z);
  
  //2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
  this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //3) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
  
  for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
    a=this->RealSignalForFFT.G(x, y, z);
    b=this->ImagSignalForFFT.G(x, y, z);
    c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
    d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
    
    this->RealSignalForFFT.P(a*c-b*d, x, y, z);
    this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
  }
  
  //4) IFFT
  this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
}

///put a value in the real part of the field that is transformed by the class
void FFTconvolver3D::P(float value,int x,int y, int z){
  this->RealSignalForFFT.P(value,x,y,z);
}

///put a value in the real part of the field that is transformed by the class
float FFTconvolver3D::G(int x,int y, int z){
  return this->RealSignalForFFT.G(x,y,z);
}

///deconvolution of a 3D scalar field using the predifined kernel
/// !!! NOT VALIDATED !!!
void FFTconvolver3D::Deconvolution(ScalarField * SF){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  cout << "DECONVOLUTION SHOULD BE USED CARREFULLY HERE - NOT VALIDATED!!!\n";
  
  //1) Copy the orginal image in the image that will be transformed
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
  
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
    this->RealSignalForFFT.P(SF->G(x,y,z),x,y,z);
  }
  //2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
  this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //3) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
  
  for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
    a=this->RealSignalForFFT.G(x, y, z);
    b=this->ImagSignalForFFT.G(x, y, z);
    c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
    d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
    
    this->RealSignalForFFT.P((a*c+b*d)/(c*c+d*d), x, y, z);
    this->ImagSignalForFFT.P((c*b-a*d)/(c*c+d*d),x, y, z);
  }
  //4) IFFT
  this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //5) Copy the image that has been deconvolved in the orginal image
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
    SF->P(this->RealSignalForFFT.G(x,y,z),x,y,z);
  }
}




///compute the update vector field using FFT
void SmoothVFUsingFFT(VectorField * SmoothedField,FFTconvolver3D * FFTconvolver_loc){
  int x,y,z;
  int direction;
  
  //smooth the information in the 3 directions X, Y, Z  (3D IMAGE)
  if (SmoothedField->NZ>1) for (direction=0;direction<3;direction++){
    //copy the region in the convolver
    for (z = 2; z < SmoothedField->NZ-2; z++) for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++){
      FFTconvolver_loc->P(SmoothedField->G(direction,x,y,z),x,y,z);
    }
    
    //convolve the region
    FFTconvolver_loc->Convolution();
    
    //copy the smoothed region in SmoothedField
    for (z = 0; z < SmoothedField->NZ-0; z++) for (y = 0; y < SmoothedField->NY-0; y++)  for (x = 0; x < SmoothedField->NX-0; x++){
      SmoothedField->P(FFTconvolver_loc->G(x,y,z),direction,x,y,z);
    }
  }

  //smooth the information in the 2 directions X, Y  (2D IMAGE)
  if (SmoothedField->NZ==1)  for (direction=0;direction<2;direction++){
    //copy the region in the convolver
    for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++){
      FFTconvolver_loc->P(SmoothedField->G(direction,x,y),x,y);
    }
    
    //convolve the region
    FFTconvolver_loc->Convolution();
    
    //copy the smoothed region in SmoothedField
    for (y = 0; y < SmoothedField->NY-0; y++)  for (x = 0; x < SmoothedField->NX-0; x++){
      SmoothedField->P(FFTconvolver_loc->G(x,y),direction,x,y);
    }
  }
  
  
}

///4.2: class where we consider two regions


///Constructor
MultiRegionFFTConvolver::MultiRegionFFTConvolver(){
  this->TypeOfConvolver0=-1;
  this->TypeOfConvolver1=-1;
  this->PartitionOfUnity.CreateVoidField(0,0,0);
  this->TypicalScaleVF0=-1;
  this->TypicalScaleVF1=-1;
}

///destructor
MultiRegionFFTConvolver::~MultiRegionFFTConvolver(){
}


///Initiate the convolver of region 0
//7 Gaussians (set some weights to 0 if less Gaussians are required)
// * NX, NY, NZ: is the size of the input image
// * If DivFree==1 then the kernel is divergence free
// * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
// * ...
// * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void MultiRegionFFTConvolver::InitiateConvolver_Reg0(int NBX,int NBY, int NBZ, int DivFree, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){
  if (DivFree!=1){
    this->TypeOfConvolver0=0;
    this->Region0_convolver.InitiateConvolver(NBX,NBY,NBZ,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7);
  }
  else{
    this->TypeOfConvolver0=1;
    this->Region0_convolver.InitiateMatrixValuedConvolver(NBX,NBY,NBZ);
    this->SigmaReg0_IfDivFree=sX1;
  }
}

///Initiate the convolver of region 1
//7 Gaussians (set some weights to 0 if less Gaussians are required)
// * NX, NY, NZ: is the size of the input image
// * If DivFree==1 then the kernel is divergence free
// * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
// * ...
// * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void MultiRegionFFTConvolver::InitiateConvolver_Reg1(int NBX,int NBY, int NBZ, int DivFree, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){
  if (DivFree!=1){
    this->TypeOfConvolver1=0;
    this->Region1_convolver.InitiateConvolver(NBX,NBY,NBZ,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7);
  }
  else{
    this->TypeOfConvolver1=1;
    this->Region1_convolver.InitiateMatrixValuedConvolver(NBX,NBY,NBZ);
    this->SigmaReg1_IfDivFree=sX1;
  }
}



///Set the partition of unity (it contains values between 0 and 1 -> if==0 then region 0 / if==1 then region 1 / otherwise something inbetween)
void MultiRegionFFTConvolver::SetPartitionOfUnity(ScalarField * RefPartitionOfUnity){
  int i,j,k;
  float maxV,minV,epsilon;
  char fileName[256];
            

  //create the mask
  this->PartitionOfUnity.CreateVoidField(RefPartitionOfUnity->NX,RefPartitionOfUnity->NY,RefPartitionOfUnity->NZ);
  
  //fill the mask and check the maximum and minumum values
  maxV=RefPartitionOfUnity->G(0,0,0);
  minV=RefPartitionOfUnity->G(0,0,0);
  
  for (i=0;i<RefPartitionOfUnity->NZ;i++) for (j=0;j<RefPartitionOfUnity->NY;j++) for (k=0;k<RefPartitionOfUnity->NX;k++){
    this->PartitionOfUnity.P(RefPartitionOfUnity->G(k,j,i),k,j,i);
    if (RefPartitionOfUnity->G(k,j,i)>maxV) maxV=RefPartitionOfUnity->G(k,j,i);
    if (RefPartitionOfUnity->G(k,j,i)<minV) minV=RefPartitionOfUnity->G(k,j,i);
  }
  
  //post-treatments to make sure all values are between 0 and 1  (or equal to 0 if uniform image)
   ///COMMENT TO REMOVE
  //epsilon=0.0001;
  //if ((fabs(minV)>epsilon)||(fabs(maxV-1)>epsilon)){
//    cout << "Partition of unity is not contained between 0 and 1 -> modified to respect this condition" << endl;
    //
    //if (fabs(minV-maxV)<epsilon){
//      for (i=0;i<RefPartitionOfUnity->NZ;i++) for (j=0;j<RefPartitionOfUnity->NY;j++) for (k=0;k<RefPartitionOfUnity->NX;k++)
//  this->PartitionOfUnity.P(0,k,j,i);
    //}
    //else{
//      for (i=0;i<RefPartitionOfUnity->NZ;i++) for (j=0;j<RefPartitionOfUnity->NY;j++) for (k=0;k<RefPartitionOfUnity->NX;k++)
//  this->PartitionOfUnity.P((this->PartitionOfUnity.G(k,j,i)-minV)/(maxV-minV),k,j,i);
    //}
//    
    //
    //strcpy(fileName,"ActualParitionOfUnity.nii");
    //cout << "The actual partition of unity is saved in " << fileName << " (with image 2 world matrix equals to identity)" << endl;
    //this->PartitionOfUnity.Write(fileName);
  //} 
  ///END COMMENT TO REMOVE


  //In addition: initiate the temporary VectorFields
  TempVF0.CreateVoidField(RefPartitionOfUnity->NX,RefPartitionOfUnity->NY,RefPartitionOfUnity->NZ);
  TempVF1.CreateVoidField(RefPartitionOfUnity->NX,RefPartitionOfUnity->NY,RefPartitionOfUnity->NZ);
}


///convolution of a 3D vector field using the predifined kernel
void MultiRegionFFTConvolver::Convolution(VectorField * VF){ 
  int i,j,k,l;
  float tmpFl;
  
  //1) check that everything is initialized
  if (this->TypeOfConvolver0==-1){
    cout << "Convolver 0 is not initialized -> no convolution" << endl;
    return;
  }
  if (this->TypeOfConvolver1==-1){
    cout << "Convolver 1 is not initialized -> no convolution" << endl;
    return;
  }
  if ((this->PartitionOfUnity.NZ==0)||(this->PartitionOfUnity.NY==0)||(this->PartitionOfUnity.NX==0)){
    cout << "Partition of unity is not defined -> no convolution" << endl;
    return;
  }
  
  //2) Convolution
  //2.1) fill the temporary vector fields
  DeepCopy(VF,&this->TempVF0,0);
  for (l=0;l<3;l++) for (i=0;i<this->PartitionOfUnity.NZ;i++) for (j=0;j<this->PartitionOfUnity.NY;j++) for (k=0;k<this->PartitionOfUnity.NX;k++)
    this->TempVF0.P(this->TempVF0.G(l,k,j,i)*(1-this->PartitionOfUnity.G(k,j,i)),l,k,j,i);
    
  DeepCopy(VF,&this->TempVF1,0);
  for (l=0;l<3;l++) for (i=0;i<this->PartitionOfUnity.NZ;i++) for (j=0;j<this->PartitionOfUnity.NY;j++) for (k=0;k<this->PartitionOfUnity.NX;k++)
    this->TempVF1.P(this->TempVF1.G(l,k,j,i)*this->PartitionOfUnity.G(k,j,i),l,k,j,i);
  
  //2.2) smooth the temporary vector fields
  if (this->TypeOfConvolver0==0)
    SmoothVFUsingFFT(&this->TempVF0,&this->Region0_convolver);
  else
    this->Region0_convolver.Convolution_DivergenceFreeFilter(&this->TempVF0,this->SigmaReg0_IfDivFree);
  
  if (this->TypeOfConvolver1==0)
    SmoothVFUsingFFT(&this->TempVF1,&this->Region1_convolver);
  else
    this->Region1_convolver.Convolution_DivergenceFreeFilter(&this->TempVF1,this->SigmaReg1_IfDivFree);
  
  //2.3) scaling stuffs at the first time the class is used after allocation
  if (this->TypicalScaleVF0<0){
    for (i=0;i<this->PartitionOfUnity.NZ;i++) for (j=0;j<this->PartitionOfUnity.NY;j++) for (k=0;k<this->PartitionOfUnity.NX;k++){
       tmpFl=(this->TempVF0.G(0,k,j,i)*this->TempVF0.G(0,k,j,i))+(this->TempVF0.G(1,k,j,i)*this->TempVF0.G(1,k,j,i))+(this->TempVF0.G(2,k,j,i)*this->TempVF0.G(2,k,j,i));
       if (tmpFl>this->TypicalScaleVF0) this->TypicalScaleVF0=tmpFl;
    }
    this->TypicalScaleVF0=sqrt(this->TypicalScaleVF0);
    
    for (i=0;i<this->PartitionOfUnity.NZ;i++) for (j=0;j<this->PartitionOfUnity.NY;j++) for (k=0;k<this->PartitionOfUnity.NX;k++){
       tmpFl=(this->TempVF1.G(0,k,j,i)*this->TempVF1.G(0,k,j,i))+(this->TempVF1.G(1,k,j,i)*this->TempVF1.G(1,k,j,i))+(this->TempVF1.G(2,k,j,i)*this->TempVF1.G(2,k,j,i));
       if (tmpFl>this->TypicalScaleVF1) this->TypicalScaleVF1=tmpFl;
    }
    this->TypicalScaleVF1=sqrt(this->TypicalScaleVF1);
  }
  
  //2.4) compute the smoothed vector field
  for (l=0;l<3;l++) for (i=0;i<this->PartitionOfUnity.NZ;i++) for (j=0;j<this->PartitionOfUnity.NY;j++) for (k=0;k<this->PartitionOfUnity.NX;k++)
    VF->P((this->TempVF0.G(l,k,j,i)*(1-this->PartitionOfUnity.G(k,j,i))/this->TypicalScaleVF0)+(this->TempVF1.G(l,k,j,i)*this->PartitionOfUnity.G(k,j,i)/this->TypicalScaleVF1),l,k,j,i);
}


///4.3: light weight convolver

///Constructor
LightFFTconvolver3D::LightFFTconvolver3D(){
  this->NX=0;    this->NY=0;    this->NZ=0;
  this->NXfft=0; this->NYfft=0; this->NZfft=0;
}

///destructor
LightFFTconvolver3D::~LightFFTconvolver3D(void){
  this->NX=0;    this->NY=0;    this->NZ=0;
  this->NXfft=0; this->NYfft=0; this->NZfft=0;
}

///Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 4 Gaussians (set some weights to 0 if less Gaussians are required)
/// * NX, NY, NZ: is the size of the input image
/// * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
///   ...
/// * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void LightFFTconvolver3D::InitiateConvolver(int NBX,int NBY, int NBZ,float w1,float sX1,float sY1,float sZ1,float w2,float sX2,float sY2,float sZ2,float w3,float sX3,float sY3,float sZ3,float w4,float sX4,float sY4,float sZ4,float w5,float sX5,float sY5,float sZ5,float w6,float sX6,float sY6,float sZ6,float w7,float sX7,float sY7,float sZ7,int NormalizeWeights){
  //set the size of the original image
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  
  //set the size of images for the FFT
  this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
  this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
  this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
  

  //cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << "\n";
  
  //allocate memory for the images for the FFT
  this->RealSignalForFFT_X.CreateVoidField(this->NXfft, 1, 1); //image  - real part
  this->ImagSignalForFFT_X.CreateVoidField(this->NXfft, 1, 1); //image  - imaginary part
  this->RealSignalForFFT_Y.CreateVoidField(1, this->NYfft, 1); //image  - real part
  this->ImagSignalForFFT_Y.CreateVoidField(1, this->NYfft, 1); //image  - imaginary part
  this->RealSignalForFFT_Z.CreateVoidField(1, 1, this->NZfft); //image  - real part
  this->ImagSignalForFFT_Z.CreateVoidField(1, 1, this->NZfft); //image  - imaginary part

  this->RealFilterForFFT_X.CreateVoidField(this->NXfft, 1, 1); //filter - real part
  this->ImagFilterForFFT_X.CreateVoidField(this->NXfft, 1, 1); //filter - imaginary part
  this->RealFilterForFFT_Y.CreateVoidField(1, this->NYfft, 1); //filter - real part
  this->ImagFilterForFFT_Y.CreateVoidField(1, this->NYfft, 1); //filter - imaginary part
  this->RealFilterForFFT_Z.CreateVoidField(1, 1, this->NZfft); //filter - real part
  this->ImagFilterForFFT_Z.CreateVoidField(1, 1, this->NZfft); //filter - imaginary part
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf7AnisotropicGaussianFilters(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,NormalizeWeights);
}

///change the kernel of the convolver (same notations as the constructor)
void LightFFTconvolver3D::ChangeKernel(float w1,float sX1,float sY1,float sZ1,float w2,float sX2,float sY2,float sZ2,float w3,float sX3,float sY3,float sZ3,float w4,float sX4,float sY4,float sZ4,float w5,float sX5,float sY5,float sZ5,float w6,float sX6,float sY6,float sZ6,float w7,float sX7,float sY7,float sZ7,int NormalizeWeights){
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf7AnisotropicGaussianFilters(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,NormalizeWeights);
}


#ifdef COMPILE_WITH_OPENMP

///convolution of a 3D scalar field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(ScalarField * SF,int TimeFrame){
  int x,y,z,i;
  float a,b,c,d;
  float CoefMult;
  int MinTimeFrame;
  int MaxTimeFrame;
  ScalarField loc_RealSignalForFFT;     //for openmp
  ScalarField loc_ImagSignalForFFT;     //for openmp
  
  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=SF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,a,b,c,d,i,CoefMult,loc_RealSignalForFFT,loc_ImagSignalForFFT) 
  {
    for (i=MinTimeFrame;i<=MaxTimeFrame;i++){ //loop on the time frames
      //1) convolution on X axis
      loc_RealSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);     //added
      loc_ImagSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);     //added
  
      #pragma omp for
      for (y = 0; y < SF->NY; y++){
        for (z = 0; z < SF->NZ; z++){
          //1.1) Copy the orginal image in the image that will be transformed
          for (x=0;x<this->NXfft;x++) loc_RealSignalForFFT.P(0.,x,0,0);
          for (x=0;x<this->NXfft;x++) loc_ImagSignalForFFT.P(0.,x,0,0);
          
          for (x = 0; x < SF->NX; x++) loc_RealSignalForFFT.P(SF->G(x,y,z,i),x,0,0);
    
    
          //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NX);
          
          for (x = 0; x < loc_RealSignalForFFT.NX; x++){
            a=loc_RealSignalForFFT.G(x, 0, 0);
            b=loc_ImagSignalForFFT.G(x, 0, 0);
            c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
            d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,x,0,0);
            loc_ImagSignalForFFT.P(c*b+a*d,x,0,0);
          }
          
          //1.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.5) Copy the image that has been convolved in the orginal image
          for (x = 0; x < SF->NX; x++){
            SF->P(loc_RealSignalForFFT.G(x,0,0),x,y,z,i);
          }
        }
      }
    
      //2) convolution on Y axis
      loc_RealSignalForFFT.CreateVoidField(1, this->NYfft, 1,1,0,0);     //added
      loc_ImagSignalForFFT.CreateVoidField(1, this->NYfft, 1,1,0,0);     //added
  
      #pragma omp for
      for (x = 0; x < SF->NX; x++){
        for (z = 0; z < SF->NZ; z++) {
          //2.1) Copy the orginal image in the image that will be transformed
          for (y=0;y<this->NYfft;y++) loc_RealSignalForFFT.P(0.,0,y,0);
          for (y=0;y<this->NYfft;y++) loc_ImagSignalForFFT.P(0.,0,y,0);
          
          for (y=0;y<SF->NY;y++) loc_RealSignalForFFT.P(SF->G(x,y,z,i),0,y,0);
    
    
          //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NY);
          
          for (y = 0; y < loc_RealSignalForFFT.NY; y++){
            a=loc_RealSignalForFFT.G(0, y, 0);
            b=loc_ImagSignalForFFT.G(0, y, 0);
            c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
            d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,y,0);
            loc_ImagSignalForFFT.P(c*b+a*d,0,y,0);
          }
          
          //2.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.5) Copy the image that has been convolved in the orginal image
          for (y = 0; y < SF->NY; y++){
            SF->P(loc_RealSignalForFFT.G(0,y,0),x,y,z,i);
          }
        }
      }


      //3) convolution on Z axis
      loc_RealSignalForFFT.CreateVoidField(1, 1, this->NZfft,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1, 1, this->NZfft,1,0,0);    //added
  
      #pragma omp for
      for (y = 0; y < SF->NY; y++){
        for (x = 0; x < SF->NX; x++){
          //3.1) Copy the orginal image in the image that will be transformed
          for (z=0;z<this->NZfft;z++) loc_RealSignalForFFT.P(0.,0,0,z);
          for (z=0;z<this->NZfft;z++) loc_ImagSignalForFFT.P(0.,0,0,z);
          
          for (z = 0; z < SF->NZ; z++) loc_RealSignalForFFT.P(SF->G(x,y,z,i),0,0,z);
    
    
          //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NZ);
          
          for (z = 0; z < loc_RealSignalForFFT.NZ; z++){
            a=loc_RealSignalForFFT.G(0, 0, z);
            b=loc_ImagSignalForFFT.G(0, 0, z);
            c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
            d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,0,z);
            loc_ImagSignalForFFT.P(c*b+a*d,0,0,z);
          }
          
          //3.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.5) Copy the image that has been convolved in the orginal image
          for (z = 0; z < SF->NZ; z++){
            SF->P(loc_RealSignalForFFT.G(0,0,z),x,y,z,i);
          }
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///convolution of a 3D scalar field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(ScalarField * SF,int TimeFrame){
  int x,y,z,i;
  float a,b,c,d;
  float CoefMult;
  int MinTimeFrame;
  int MaxTimeFrame;

  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=SF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  for (i=MinTimeFrame;i<=MaxTimeFrame;i++){ //loop on the time frames
    //1) convolution on X axis
    for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      for (x = 0; x < SF->NX; x++) this->RealSignalForFFT_X.P(SF->G(x,y,z,i),x,0,0);


      //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
      
      for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
        a=this->RealSignalForFFT_X.G(x, 0, 0);
        b=this->ImagSignalForFFT_X.G(x, 0, 0);
        c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
        d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
        
        this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
        this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
      }
      
      //1.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.5) Copy the image that has been convolved in the orginal image
      for (x = 0; x < SF->NX; x++){
        SF->P(this->RealSignalForFFT_X.G(x,0,0),x,y,z,i);
      }
    }
    
    //2) convolution on Y axis
    for (z = 0; z < SF->NZ; z++) for (x = 0; x < SF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      for (y=0;y<SF->NY;y++) this->RealSignalForFFT_Y.P(SF->G(x,y,z,i),0,y,0);


      //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
      
      for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
        a=this->RealSignalForFFT_Y.G(0, y, 0);
        b=this->ImagSignalForFFT_Y.G(0, y, 0);
        c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
        d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
        
        this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
        this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
      }
      
      //2.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.5) Copy the image that has been convolved in the orginal image
      for (y = 0; y < SF->NY; y++){
        SF->P(this->RealSignalForFFT_Y.G(0,y,0),x,y,z,i);
      }
    }

    //3) convolution on Z axis
    for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      for (z = 0; z < SF->NZ; z++) this->RealSignalForFFT_Z.P(SF->G(x,y,z,i),0,0,z);


      //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
      
      for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
        a=this->RealSignalForFFT_Z.G(0, 0, z);
        b=this->ImagSignalForFFT_Z.G(0, 0, z);
        c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
        d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
        
        this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
        this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
      }
      
      //3.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.5) Copy the image that has been convolved in the orginal image
      for (z = 0; z < SF->NZ; z++){
        SF->P(this->RealSignalForFFT_Z.G(0,0,z),x,y,z,i);
      }
    }
  }
}

#endif

#ifdef COMPILE_WITH_OPENMP

///convolution of a 3D vector field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(VectorField * VF,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;
  ScalarField loc_RealSignalForFFT;     //for openmp
  ScalarField loc_ImagSignalForFFT;     //for openmp

  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,a,b,c,d,i,j,CoefMult,loc_RealSignalForFFT,loc_ImagSignalForFFT) 
  {
    for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
      //1) convolution on X axis
      loc_RealSignalForFFT.CreateVoidField(this->NXfft, 1,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(this->NXfft, 1,1,1,0,0);    //added

      #pragma omp for
      for (y = 0; y < VF->NY; y++){
        for (z = 0; z < VF->NZ; z++){
          //1.1) Copy the orginal image in the image that will be transformed
          for (x=0;x<this->NXfft;x++) loc_RealSignalForFFT.P(0.,x,0,0);
          for (x=0;x<this->NXfft;x++) loc_ImagSignalForFFT.P(0.,x,0,0);
          
          for (x = 0; x < VF->NX; x++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),x,0,0);
    
    
          //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NX);
          
          for (x = 0; x < loc_RealSignalForFFT.NX; x++){
            a=loc_RealSignalForFFT.G(x, 0, 0);
            b=loc_ImagSignalForFFT.G(x, 0, 0);
            c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
            d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,x,0,0);
            loc_ImagSignalForFFT.P(c*b+a*d,x,0,0);
          }
          
          //1.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.5) Copy the image that has been convolved in the orginal image
          for (x = 0; x < VF->NX; x++){
            VF->P(loc_RealSignalForFFT.G(x,0,0),i,x,y,z,j);
          }
        }
      }
      
      //2) convolution on Y axis
      loc_RealSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      
      #pragma omp for
      for (x = 0; x < VF->NX; x++){ 
        for (z = 0; z < VF->NZ; z++){
          //2.1) Copy the orginal image in the image that will be transformed
          for (y=0;y<this->NYfft;y++) loc_RealSignalForFFT.P(0.,0,y,0);
          for (y=0;y<this->NYfft;y++) loc_ImagSignalForFFT.P(0.,0,y,0);
          
          for (y = 0; y < VF->NY; y++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,y,0);
    
    
          //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NY);
          
          for (y = 0; y < loc_RealSignalForFFT.NY; y++){
            a=loc_RealSignalForFFT.G(0, y, 0);
            b=loc_ImagSignalForFFT.G(0, y, 0);
            c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
            d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,y,0);
            loc_ImagSignalForFFT.P(c*b+a*d,0,y,0);
          }
          
          //2.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.5) Copy the image that has been convolved in the orginal image
          for (y = 0; y < VF->NY; y++){
            VF->P(loc_RealSignalForFFT.G(0,y,0),i,x,y,z,j);
          }
        }
      }
      
      //3) convolution on Z axis
      loc_RealSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      
      #pragma omp for
      for (y = 0; y < VF->NY; y++){ 
        for (x = 0; x < VF->NX; x++){
          //3.1) Copy the orginal image in the image that will be transformed
          for (z=0;z<this->NZfft;z++) loc_RealSignalForFFT.P(0.,0,0,z);
          for (z=0;z<this->NZfft;z++) loc_ImagSignalForFFT.P(0.,0,0,z);
          
          for (z = 0; z < VF->NZ; z++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,0,z);
    
    
          //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NZ);
          
          for (z = 0; z < loc_RealSignalForFFT.NZ; z++){
            a=loc_RealSignalForFFT.G(0, 0, z);
            b=loc_ImagSignalForFFT.G(0, 0, z);
            c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
            d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,0,z);
            loc_ImagSignalForFFT.P(c*b+a*d,0,0,z);
          }
          
          //3.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.5) Copy the image that has been convolved in the orginal image
          for (z = 0; z < VF->NZ; z++){
            VF->P(loc_RealSignalForFFT.G(0,0,z),i,x,y,z,j);
          }
        }
      }
    } 
  //END FORK FOR THREADS
  }
}

#else

///convolution of a 3D vector field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(VectorField * VF,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;

  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }

  
  for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
    //1) convolution on X axis
    for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      for (x = 0; x < VF->NX; x++) this->RealSignalForFFT_X.P(VF->G(i,x,y,z,j),x,0,0);


      //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
      
      for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
        a=this->RealSignalForFFT_X.G(x, 0, 0);
        b=this->ImagSignalForFFT_X.G(x, 0, 0);
        c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
        d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
        
        this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
        this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
      }
      
      //1.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.5) Copy the image that has been convolved in the orginal image
      for (x = 0; x < VF->NX; x++){
        VF->P(this->RealSignalForFFT_X.G(x,0,0),i,x,y,z,j);
      }
    }
    
    //2) convolution on Y axis
    for (z = 0; z < VF->NZ; z++) for (x = 0; x < VF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      for (y = 0; y < VF->NY; y++) this->RealSignalForFFT_Y.P(VF->G(i,x,y,z,j),0,y,0);


      //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
      
      for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
        a=this->RealSignalForFFT_Y.G(0, y, 0);
        b=this->ImagSignalForFFT_Y.G(0, y, 0);
        c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
        d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
        
        this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
        this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
      }
      
      //2.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.5) Copy the image that has been convolved in the orginal image
      for (y = 0; y < VF->NY; y++){
        VF->P(this->RealSignalForFFT_Y.G(0,y,0),i,x,y,z,j);
      }
    }
    
    //3) convolution on Z axis
    for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      for (z = 0; z < VF->NZ; z++) this->RealSignalForFFT_Z.P(VF->G(i,x,y,z,j),0,0,z);


      //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
      
      for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
        a=this->RealSignalForFFT_Z.G(0, 0, z);
        b=this->ImagSignalForFFT_Z.G(0, 0, z);
        c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
        d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
        
        this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
        this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
      }
      
      //3.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.5) Copy the image that has been convolved in the orginal image
      for (z = 0; z < VF->NZ; z++){
        VF->P(this->RealSignalForFFT_Z.G(0,0,z),i,x,y,z,j);
      }
    }
  } 
}

#endif

#ifdef COMPILE_WITH_OPENMP

///convolution of a 3D vector field using the predifined kernel. Convolution is performed in the ROI defined by (xmin, xmax, ymin, ymax, zmin, zmax) only.
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::ConvolutionInROI(VectorField * VF,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;
  ScalarField loc_RealSignalForFFT;     //for openmp
  ScalarField loc_ImagSignalForFFT;     //for openmp

  //0) Define the time frames to smooth and check the boundaries
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  if (xmin<0)      {xmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (ymin<0)      {ymin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (zmin<0)      {zmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (xmax>VF->NX) {xmax=VF->NX; cout << "Region boundaries are not well defined" << endl;}
  if (ymax>VF->NY) {ymax=VF->NY; cout << "Region boundaries are not well defined" << endl;}
  if (zmax>VF->NZ) {zmax=VF->NZ; cout << "Region boundaries are not well defined" << endl;}

  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,a,b,c,d,i,j,CoefMult,loc_RealSignalForFFT,loc_ImagSignalForFFT) 
  {
    for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
      //1) convolution on X axis
      loc_RealSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);    //added
      
      #pragma omp for
      for (y = ymin; y < ymax; y++){ 
        for (z = zmin; z < zmax; z++){
          //1.1) Copy the orginal image in the image that will be transformed
          for (x=0;x<this->NXfft;x++) loc_RealSignalForFFT.P(0.,x,0,0);
          for (x=0;x<this->NXfft;x++) loc_ImagSignalForFFT.P(0.,x,0,0);
          
          for (x = 0; x < VF->NX; x++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),x,0,0);
  
          //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NX);
          
          for (x = 0; x < loc_RealSignalForFFT.NX; x++){
            a=loc_RealSignalForFFT.G(x, 0, 0);
            b=loc_ImagSignalForFFT.G(x, 0, 0);
            c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
            d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,x,0,0);
            loc_ImagSignalForFFT.P(c*b+a*d,x,0,0);
          }
          
          //1.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.5) Copy the image that has been convolved in the orginal image
          for (x = xmin; x < xmax; x++){
            VF->P(loc_RealSignalForFFT.G(x,0,0),i,x,y,z,j);
          }
        }
      }
      
      
      //2) convolution on Y axis
      loc_RealSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      
      #pragma omp for
      for (x = xmin; x < xmax; x++){ 
        for (z = zmin; z < zmax; z++){
          //2.1) Copy the orginal image in the image that will be transformed
          for (y=0;y<this->NYfft;y++) loc_RealSignalForFFT.P(0.,0,y,0);
          for (y=0;y<this->NYfft;y++) loc_ImagSignalForFFT.P(0.,0,y,0);
          
          for (y = 0; y < VF->NY; y++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,y,0);
    
    
          //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NY);
          
          for (y = 0; y < loc_RealSignalForFFT.NY; y++){
            a=loc_RealSignalForFFT.G(0, y, 0);
            b=loc_ImagSignalForFFT.G(0, y, 0);
            c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
            d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,y,0);
            loc_ImagSignalForFFT.P(c*b+a*d,0,y,0);
          }
          
          //2.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.5) Copy the image that has been convolved in the orginal image
          for (y = ymin; y < ymax; y++){
            VF->P(loc_RealSignalForFFT.G(0,y,0),i,x,y,z,j);
          }
        }
      }
      
      //3) convolution on Z axis
      loc_RealSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      
      #pragma omp for
      for (y = ymin; y < ymax; y++){
        for (x = xmin; x < xmax; x++){
          //3.1) Copy the orginal image in the image that will be transformed
          for (z=0;z<this->NZfft;z++) loc_RealSignalForFFT.P(0.,0,0,z);
          for (z=0;z<this->NZfft;z++) loc_ImagSignalForFFT.P(0.,0,0,z);
          
          for (z = 0; z < VF->NZ; z++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,0,z);
    
    
          //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NZ);
          
          for (z = 0; z < loc_RealSignalForFFT.NZ; z++){
            a=loc_RealSignalForFFT.G(0, 0, z);
            b=loc_ImagSignalForFFT.G(0, 0, z);
            c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
            d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,0,z);
            loc_ImagSignalForFFT.P(c*b+a*d,0,0,z);
          }
          
          //3.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.5) Copy the image that has been convolved in the orginal image
          for (z = zmin; z < zmax; z++){
            VF->P(loc_RealSignalForFFT.G(0,0,z),i,x,y,z,j);
          }
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///convolution of a 3D vector field using the predifined kernel. Convolution is performed in the ROI defined by (xmin, xmax, ymin, ymax, zmin, zmax) only.
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::ConvolutionInROI(VectorField * VF,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;

  //0) Define the time frames to smooth and check the boundaries
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  if (xmin<0)      {xmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (ymin<0)      {ymin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (zmin<0)      {zmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (xmax>VF->NX) {xmax=VF->NX; cout << "Region boundaries are not well defined" << endl;}
  if (ymax>VF->NY) {ymax=VF->NY; cout << "Region boundaries are not well defined" << endl;}
  if (zmax>VF->NZ) {zmax=VF->NZ; cout << "Region boundaries are not well defined" << endl;}

  
  for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
    //1) convolution on X axis
    for (z = zmin; z < zmax; z++) for (y = ymin; y < ymax; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      for (x = 0; x < VF->NX; x++) this->RealSignalForFFT_X.P(VF->G(i,x,y,z,j),x,0,0);


      //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
      
      for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
        a=this->RealSignalForFFT_X.G(x, 0, 0);
        b=this->ImagSignalForFFT_X.G(x, 0, 0);
        c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
        d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
        
        this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
        this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
      }
      
      //1.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.5) Copy the image that has been convolved in the orginal image
      for (x = xmin; x < xmax; x++){
        VF->P(this->RealSignalForFFT_X.G(x,0,0),i,x,y,z,j);
      }
    }
    
    //2) convolution on Y axis
    for (z = zmin; z < zmax; z++) for (x = xmin; x < xmax; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      for (y = 0; y < VF->NY; y++) this->RealSignalForFFT_Y.P(VF->G(i,x,y,z,j),0,y,0);


      //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
      
      for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
        a=this->RealSignalForFFT_Y.G(0, y, 0);
        b=this->ImagSignalForFFT_Y.G(0, y, 0);
        c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
        d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
        
        this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
        this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
      }
      
      //2.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.5) Copy the image that has been convolved in the orginal image
      for (y = ymin; y < ymax; y++){
        VF->P(this->RealSignalForFFT_Y.G(0,y,0),i,x,y,z,j);
      }
    }
    
    //3) convolution on Z axis
    for (y = ymin; y < ymax; y++) for (x = xmin; x < xmax; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      for (z = 0; z < VF->NZ; z++) this->RealSignalForFFT_Z.P(VF->G(i,x,y,z,j),0,0,z);


      //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
      
      for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
        a=this->RealSignalForFFT_Z.G(0, 0, z);
        b=this->ImagSignalForFFT_Z.G(0, 0, z);
        c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
        d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
        
        this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
        this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
      }
      
      //3.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.5) Copy the image that has been convolved in the orginal image
      for (z = zmin; z < zmax; z++){
        VF->P(this->RealSignalForFFT_Z.G(0,0,z),i,x,y,z,j);
      }
    }
  } 
}
#endif

///Hack to perform convolution of the 3D vector field 'VF' in a masked region with mirror conditions.
///  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
///  -> sX1,sY1,sZ1: size of the smoothing kernel
void LightFFTconvolver3D::Convolution_Mask_Mirror(VectorField * VF,ScalarField * Mask, int MaskId){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i;
  int PreviousInMask,MirrorInMask;
  int Something;
  float MaskId_fl;
  
  MaskId_fl=static_cast<float>(MaskId);
  
  for (i=0;i<3;i++){ // loop on the time frames and directions
    //1) convolution on X axis
    for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      //1.1.1) Init
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      Something=0;
      for (x = 0; x < VF->NX; x++)  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_X.P(VF->G(i,x,y,z),x,0,0);
	Something=1;
      }
      
      if (Something==1){
	//1.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (x = 0; x < VF->NX; x++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}  

	//1.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (x = VF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_X.G(x,0,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_X.G(x,0,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}
	
	//1.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (x = 0; x < VF->NX; x++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)<0)  this->RealSignalForFFT_X.P(VF->G(i,MirrorInMask,y,z),x,0,0);
	  }
	}

	//1.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (x = VF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>VF->NX-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)>0)  this->RealSignalForFFT_X.P(VF->G(i,MirrorInMask,y,z),x,0,0);
	  }
	}	

	//1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
	
	for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
	  a=this->RealSignalForFFT_X.G(x, 0, 0);
	  b=this->ImagSignalForFFT_X.G(x, 0, 0);
	  c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  
	  this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
	  this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
	}
	
	//1.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.5) Copy the image that has been convolved in the orginal image
	for (x = 0; x < VF->NX; x++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  VF->P(this->RealSignalForFFT_X.G(x,0,0),i,x,y,z);
	}
      }
    }
  
    
    //2) convolution on Y axis
    for (z = 0; z < VF->NZ; z++) for (x = 0; x < VF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      //2.1.1) Init
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      Something=0;
      for (y = 0; y < VF->NY; y++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_Y.P(VF->G(i,x,y,z),0,y,0);
	Something=1;
      }
      
      
      if (Something==1){
	//2.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (y = 0; y < VF->NY; y++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (y = VF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Y.G(0,y,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Y.G(0,y,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (y = 0; y < VF->NY; y++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)<0) this->RealSignalForFFT_Y.P(VF->G(i,x,MirrorInMask,z),0,y,0);
	  }
	}

	//2.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (y = VF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>VF->NY-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)>0)  this->RealSignalForFFT_Y.P(VF->G(i,x,MirrorInMask,z),0,y,0);
	  }
	}


	//2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
	
	for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
	  a=this->RealSignalForFFT_Y.G(0, y, 0);
	  b=this->ImagSignalForFFT_Y.G(0, y, 0);
	  c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  
	  this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
	  this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
	}
	
	//2.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.5) Copy the image that has been convolved in the orginal image
	for (y = 0; y < VF->NY; y++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  VF->P(this->RealSignalForFFT_Y.G(0,y,0),i,x,y,z);
	}
      }
    }
    
    //3) convolution on Z axis
    for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      
      //3.1.1) Init
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      Something=0;
      for (z = 0; z < VF->NZ; z++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){ 
	this->RealSignalForFFT_Z.P(VF->G(i,x,y,z),0,0,z);
	Something=1;
      }
      
      
            
      if (Something==1){
	//3.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (z = VF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Z.G(0,0,z)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Z.G(0,0,z))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)<0) this->RealSignalForFFT_Z.P(VF->G(i,x,y,MirrorInMask),0,0,z);
	  }
	}

	//3.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (z = VF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>VF->NZ-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)>0) this->RealSignalForFFT_Z.P(VF->G(i,x,y,MirrorInMask),0,0,z);
	  }
	}


	
	//3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
	
	for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
	  a=this->RealSignalForFFT_Z.G(0, 0, z);
	  b=this->ImagSignalForFFT_Z.G(0, 0, z);
	  c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  
	  this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
	  this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
	}
	
	//3.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.5) Copy the image that has been convolved in the orginal image
	for (z = 0; z < VF->NZ; z++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  VF->P(this->RealSignalForFFT_Z.G(0,0,z),i,x,y,z);
	}
      }
    }
  } 
}




///Hack to perform convolution of the 3D scalar field 'SF' in a masked region with mirror conditions.
///  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
void LightFFTconvolver3D::Convolution_Mask_Mirror(ScalarField * SF,ScalarField * Mask, int MaskId){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i;
  int PreviousInMask,MirrorInMask;
  int Something;
  float MaskId_fl;
  
  MaskId_fl=static_cast<float>(MaskId);
  
    //1) convolution on X axis
    for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      //1.1.1) Init
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      Something=0;
      for (x = 0; x < SF->NX; x++)  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_X.P(SF->G(x,y,z),x,0,0);
	Something=1;
      }
      
      if (Something==1){
	//1.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (x = 0; x < SF->NX; x++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}  

	//1.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (x = SF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_X.G(x,0,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_X.G(x,0,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}
	
	//1.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (x = 0; x < SF->NX; x++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)<0)  this->RealSignalForFFT_X.P(SF->G(MirrorInMask,y,z),x,0,0);
	  }
	}

	//1.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (x = SF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>SF->NX-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)>0)  this->RealSignalForFFT_X.P(SF->G(MirrorInMask,y,z),x,0,0);
	  }
	}	

	//1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
	
	for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
	  a=this->RealSignalForFFT_X.G(x, 0, 0);
	  b=this->ImagSignalForFFT_X.G(x, 0, 0);
	  c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  
	  this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
	  this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
	}
	
	//1.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.5) Copy the image that has been convolved in the orginal image
	for (x = 0; x < SF->NX; x++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  SF->P(this->RealSignalForFFT_X.G(x,0,0),x,y,z);
	}
      }
    }
  
    
    //2) convolution on Y axis
    for (z = 0; z < SF->NZ; z++) for (x = 0; x < SF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      //2.1.1) Init
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      Something=0;
      for (y = 0; y < SF->NY; y++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_Y.P(SF->G(x,y,z),0,y,0);
	Something=1;
      }
      
      
      if (Something==1){
	//2.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (y = 0; y < SF->NY; y++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (y = SF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Y.G(0,y,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Y.G(0,y,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (y = 0; y < SF->NY; y++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)<0) this->RealSignalForFFT_Y.P(SF->G(x,MirrorInMask,z),0,y,0);
	  }
	}

	//2.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (y = SF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>SF->NY-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)>0)  this->RealSignalForFFT_Y.P(SF->G(x,MirrorInMask,z),0,y,0);
	  }
	}


	//2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
	
	for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
	  a=this->RealSignalForFFT_Y.G(0, y, 0);
	  b=this->ImagSignalForFFT_Y.G(0, y, 0);
	  c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  
	  this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
	  this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
	}
	
	//2.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.5) Copy the image that has been convolved in the orginal image
	for (y = 0; y < SF->NY; y++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  SF->P(this->RealSignalForFFT_Y.G(0,y,0),x,y,z);
	}
      }
    }
    
    //3) convolution on Z axis
    for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      
      //3.1.1) Init
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      Something=0;
      for (z = 0; z < SF->NZ; z++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){ 
	this->RealSignalForFFT_Z.P(SF->G(x,y,z),0,0,z);
	Something=1;
      }
      
      
            
      if (Something==1){
	//3.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (z = SF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Z.G(0,0,z)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Z.G(0,0,z))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)<0) this->RealSignalForFFT_Z.P(SF->G(x,y,MirrorInMask),0,0,z);
	  }
	}

	//3.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (z = SF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>SF->NZ-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)>0) this->RealSignalForFFT_Z.P(SF->G(x,y,MirrorInMask),0,0,z);
	  }
	}


	
	//3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
	
	for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
	  a=this->RealSignalForFFT_Z.G(0, 0, z);
	  b=this->ImagSignalForFFT_Z.G(0, 0, z);
	  c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  
	  this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
	  this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
	}
	
	//3.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.5) Copy the image that has been convolved in the orginal image
	for (z = 0; z < SF->NZ; z++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  SF->P(this->RealSignalForFFT_Z.G(0,0,z),x,y,z);
	}
      }
    }
}









///design a kernel that is the sum of up to 7 Gaussians and transform it in Fourier spaces
//if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
void LightFFTconvolver3D::MakeSumOf7AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,float weight5,float sigmaX5,float sigmaY5,float sigmaZ5,float weight6,float sigmaX6,float sigmaY6,float sigmaZ6,float weight7,float sigmaX7,float sigmaY7,float sigmaZ7,int NormalizeWeights){
  int k,x,y,z;
  float SumLoc;
  float weight,sigmaX,sigmaY,sigmaZ,sumWeight,CubicRootOfWeight;
  
  //1) set RealFilterForFFT and ImagFilterForFFT to 0 in case it contains something
  for (x=0;x<this->NXfft;x++) this->RealFilterForFFT_X.P(0.,x,0,0);
  for (y=0;y<this->NYfft;y++) this->RealFilterForFFT_Y.P(0.,0,y,0);
  for (z=0;z<this->NZfft;z++) this->RealFilterForFFT_Z.P(0.,0,0,z);
  
  for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT_X.P(0.,x,0,0);
  for (y=0;y<this->NYfft;y++) this->ImagFilterForFFT_Y.P(0.,0,y,0);
  for (z=0;z<this->NZfft;z++) this->ImagFilterForFFT_Z.P(0.,0,0,z);
  
  sumWeight=0;
  
  //2) compute and save the 7 kernels
  for (k=0;k<7;k++){
    //parameters of the current kernel
    if (k==6){weight=weight7; sigmaX=sigmaX7; sigmaY=sigmaY7; sigmaZ=sigmaZ7;}
    if (k==5){weight=weight6; sigmaX=sigmaX6; sigmaY=sigmaY6; sigmaZ=sigmaZ6;}
    if (k==4){weight=weight5; sigmaX=sigmaX5; sigmaY=sigmaY5; sigmaZ=sigmaZ5;}
    if (k==3){weight=weight4; sigmaX=sigmaX4; sigmaY=sigmaY4; sigmaZ=sigmaZ4;}
    if (k==2){weight=weight3; sigmaX=sigmaX3; sigmaY=sigmaY3; sigmaZ=sigmaZ3;}
    if (k==1){weight=weight2; sigmaX=sigmaX2; sigmaY=sigmaY2; sigmaZ=sigmaZ2;}
    if (k==0){weight=weight1; sigmaX=sigmaX1; sigmaY=sigmaY1; sigmaZ=sigmaZ1;}
    
    CubicRootOfWeight=static_cast<float>(pow(static_cast<double>(weight),0.3333));
    sumWeight+=weight;
    
    //2.1) design the current kernel in Z direction...
    //...normalisation factor
    SumLoc=0.;
    for (z=0;z<this->NZfft/2;z++)
      SumLoc+=exp( -((float)(z*z))/(2.*sigmaZ*sigmaZ));
    
    for (z=this->NZfft/2;z<this->NZfft;z++)
      SumLoc+=exp( -((float)((this->NZfft-z)*(this->NZfft-z)))/(2.*sigmaZ*sigmaZ));
    
    //...fill the normalised values
    if (SumLoc>=0.0001){
      for (z=0;z<this->NZfft/2;z++)
        this->RealFilterForFFT_Z.Add(CubicRootOfWeight*exp( -((float)(z*z)/(2.*sigmaZ*sigmaZ)))/SumLoc,0,0,z);
    
      for (z=this->NZfft/2;z<this->NZfft;z++)
        this->RealFilterForFFT_Z.Add(CubicRootOfWeight*exp( -((float)((this->NZfft-z)*(this->NZfft-z)))/(2.*sigmaZ*sigmaZ))/SumLoc,0,0,z);
    }
    else{
        cout << "Kernel on z axis has a problem" << endl;
      }
    
    
    //2.2) design the current kernel in Y direction
    //...normalisation factor
    SumLoc=0.;
    for (y=0;y<this->NYfft/2;y++)
      SumLoc+=exp( -((float)(y*y))/(2.*sigmaY*sigmaY));
    
    for (y=this->NYfft/2;y<this->NYfft;y++)
      SumLoc+=exp( -((float)((this->NYfft-y)*(this->NYfft-y)))/(2.*sigmaY*sigmaY));
    
    //...fill the normalised values
    if (SumLoc>=0.0001){
      for (y=0;y<this->NYfft/2;y++)
	      this->RealFilterForFFT_Y.Add(CubicRootOfWeight*exp( -((float)(y*y))/(2.*sigmaY*sigmaY))/SumLoc,0,y,0);
      
      for (y=this->NYfft/2;y<this->NYfft;y++)
	      this->RealFilterForFFT_Y.Add(CubicRootOfWeight*exp( -((float)((this->NYfft-y)*(this->NYfft-y)))/(2.*sigmaY*sigmaY))/SumLoc,0,y,0);
    }
    else{
        cout << "Kernel on y axis has a problem" << endl;
      }

    
    //2.3) design the current kernel in X direction
    //...normalisation factor
    SumLoc=0.;
    for (x=0;x<this->NXfft/2;x++)
      SumLoc+=exp( -((float)(x*x))/(2.*sigmaX*sigmaX));
    
    for (x=this->NXfft/2;x<this->NXfft;x++)
      SumLoc+=exp( -((float)((this->NXfft-x)*(this->NXfft-x)))/(2.*sigmaX*sigmaX));

    //...fill the normalised values
    if (SumLoc>=0.0001){
      for (x=0;x<this->NXfft/2;x++)
	this->RealFilterForFFT_X.Add(CubicRootOfWeight*exp( -((float)(x*x))/(2.*sigmaX*sigmaX))/SumLoc,x,0,0);
      
      for (x=this->NXfft/2;x<this->NXfft;x++)
	this->RealFilterForFFT_X.Add(CubicRootOfWeight*exp( -((float)((this->NXfft-x)*(this->NXfft-x)))/(2.*sigmaX*sigmaX))/SumLoc,x,0,0);
        }
    else{
        cout << "Kernel on x axis has a problem" << endl;
      }

  }
  
  //3) normalize the weights
  if (NormalizeWeights!=0){
    for (x=0;x<this->NXfft;x++) this->RealFilterForFFT_X.P(this->RealFilterForFFT_X.G(x,0,0)/sumWeight,x,0,0);
    for (y=0;y<this->NYfft;y++) this->RealFilterForFFT_Y.P(this->RealFilterForFFT_Y.G(0,y,0)/sumWeight,0,y,0);
    for (z=0;z<this->NZfft;z++) this->RealFilterForFFT_Z.P(this->RealFilterForFFT_Z.G(0,0,z)/sumWeight,0,0,z);
    
    //cout << "Note that the weights of the kernel are normalized -> play with the -VFpenalizer option to put more or less weight on the regularisation energy" << endl;
  }

  //4) Transform the RealFilterForFFT_. and ImagFilterForFFT_. in Fourier spaces
  this->DirectFFT(&this->RealFilterForFFT_X,&this->ImagFilterForFFT_X,0);
  this->DirectFFT(&this->RealFilterForFFT_Y,&this->ImagFilterForFFT_Y,1);
  this->DirectFFT(&this->RealFilterForFFT_Z,&this->ImagFilterForFFT_Z,2);

}


///Fast Fourier Transform
/// -> if axis == 0 -> FFT on X axis
/// -> if axis == 1 -> FFT on Y axis
/// -> if axis == 2 -> FFT on Z axis
void LightFFTconvolver3D::DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;

  SqrtSizeX=sqrt(SizeX);
  SqrtSizeY=sqrt(SizeY);
  SqrtSizeZ=sqrt(SizeZ);
  
  
  //2) perform the fft along x axis
  if (axis==0){
    dataX = new float [SizeX*2+1];
    for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
      for (x = 0; x < SizeX; x++){
        dataX[2*x+1]=RealSignal->G(x, y, z);
        dataX[2*x+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataX, (unsigned long)SizeX, 1);
      for (x = 0; x < SizeX; x++){
        RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
        ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX, x, y, z);
      }
    }
    delete dataX;
  }
  
  //3) perform the fft along y axis
  if (axis==1){
    dataY = new float [SizeY*2+1];
    for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
      for (y = 0; y < SizeY; y++){
        dataY[2*y+1]=RealSignal->G(x, y, z);
        dataY[2*y+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataY, (unsigned long)SizeY, 1);
      for (y = 0; y < SizeY; y++){
        RealSignal->P(dataY[2*y+1]/SqrtSizeY,x, y, z);
        ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
      }
    }
    delete dataY;
  }
  
  //4) perform the fft along z axis
  if (axis==2){
    dataZ = new float [SizeZ*2+1];
    for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
      for (z = 0; z < SizeZ; z++){
        dataZ[2*z+1]=RealSignal->G(x, y, z);
        dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataZ, (unsigned long)SizeZ, 1);
      for (z = 0; z < SizeZ; z++){
        RealSignal->P(dataZ[2*z+1]/SqrtSizeZ,x, y, z);
        ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ, x, y, z);
      }
    }
    delete dataZ;
  }
}


///Inverse Fast Fourier Transform
/// -> if axis == 0 -> IFFT on X axis
/// -> if axis == 1 -> IFFT on Y axis
/// -> if axis == 2 -> IFFT on Z axis
void LightFFTconvolver3D::InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SqrtSizeX=sqrt(SizeX);
  
  SizeY=RealSignal->NY;
  SqrtSizeY=sqrt(SizeY);
  
  SizeZ=RealSignal->NZ;
  SqrtSizeZ=sqrt(SizeZ);
  
  
  //2) perform the ifft along z axis
  if (axis==2){
    dataZ = new float [SizeZ*2+1];
    for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
      for (z = 0; z < SizeZ; z++){
        dataZ[2*z+1]=RealSignal->G(x, y, z);
        dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataZ, (unsigned long)SizeZ, -1);
      for (z = 0; z < SizeZ; z++){
        RealSignal->P(dataZ[2*z+1]/SqrtSizeZ, x, y, z);
        ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ,x, y, z);
      }
    }
    delete dataZ;
  }
  
  //3) perform the ifft along y axis
  if (axis==1){
    dataY = new float [SizeY*2+1];
    for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
      for (y = 0; y < SizeY; y++){
        dataY[2*y+1]=RealSignal->G(x, y, z);
        dataY[2*y+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataY, (unsigned long)SizeY, -1);
      for (y = 0; y < SizeY; y++){
        RealSignal->P(dataY[2*y+1]/SqrtSizeY, x, y, z);
        ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
      }
    }
    delete dataY;
  }
  
  //4) perform the ifft along x axis
  if (axis==0){
    dataX = new float [SizeX*2+1];
    for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
      for (x = 0; x < SizeX; x++){
        dataX[2*x+1]=RealSignal->G(x, y, z);
        dataX[2*x+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataX, (unsigned long)SizeX, -1);
      for (x = 0; x < SizeX; x++){
        RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
        ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX,x, y, z);
      }
    }
    delete dataX;
  }
}  


///Fast Fourier Transform of numerical recipies (slighly modified)
void LightFFTconvolver3D::four1NR(float data[], unsigned long nn, int isign){
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2){
    if (j>i){
      tempr=data[j]; data[j]=data[i]; data[i]=tempr;
      tempr=data[j+1]; data[j+1]=data[i+1]; data[i+1]=tempr;
    }
    m=n >> 1;
    while ((m>=2) && (j>m)){
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}




















///4.4: class where we consider N regions to smooth


///Constructor
MultiRegionFFTConvolver2::MultiRegionFFTConvolver2(){
  NumberOfRegions=0;
}

///destructor
MultiRegionFFTConvolver2::~MultiRegionFFTConvolver2(){
}



//Initiate the convolver in all regions using same kernel
//-> 'Part_Of_Unity' is a 3D mask which define different ROIs. It only contains integer values each of them associated to a ROI. Partition of unity is first 
//   defined by splitting this mask into several channels, each of them having an intensity equal to 1 in the corresponding ROI and 0 otherwise. Each channel
//   is then smoothed with a Gaussian kernel of stddev 'sigmaPOI'.
//-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
//     * NX, NY, NZ: is the size of the input image
//     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
//     * ...
//     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void MultiRegionFFTConvolver2::InitiateConvolver(ScalarField * Part_Of_Unity, float sigmaPOI, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){
  int i,j;
  float tmpFl,maxSum,minSum,maxValue;
  int NBX,NBY,NBZ;
  int z,y,x;
  float * RegionsFound;
  int RegionsNb;
  float epsilon;
  int tmpInt;
  LightFFTconvolver3D TmpLightFFTconvolver;
  float x_mm,y_mm,z_mm;
  
  NBX=Part_Of_Unity->NX;
  NBY=Part_Of_Unity->NY;
  NBZ=Part_Of_Unity->NZ;
  
  //1) find the number of regions and check that all values are integers
  RegionsFound = new float [100];
  RegionsNb=0;
  tmpInt=1;
  epsilon=0.000001;
  
  //1.1) detect if the mask is an actual mask
  for(z=0;z<NBZ;z++) for(y=0;y<NBY;y++) for(x=0;x<NBX;x++) if (fabs(Part_Of_Unity->G(x,y,z))>epsilon)
      if (fabs(1-(floor(Part_Of_Unity->G(x,y,z))/Part_Of_Unity->G(x,y,z)))>epsilon)
        tmpInt=0;
  
  if ((Part_Of_Unity->NT>1)||(tmpInt==0)){
    cout << endl;
    cout << "WARNING: The mask used to generate the partition of unity has not only integer values and/or has several channels -> considered as the actual partition of unity" << endl;
    cout << endl;
    this->InitiateConvolverWithActualPOI(Part_Of_Unity,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7);
    return;
    }
  
  //1.2) find the number of regions and their identifier
  for(z=0;z<NBZ;z++) for(y=0;y<NBY;y++) for(x=0;x<NBX;x++) {
    tmpInt=0;
    for (i=0;i<RegionsNb;i++) if (fabs(Part_Of_Unity->G(x,y,z)-RegionsFound[i])<epsilon) tmpInt=1;
    
    if ((tmpInt==0)&&(RegionsNb<99)){
      RegionsFound[RegionsNb]=Part_Of_Unity->G(x,y,z);
      RegionsNb++;
      }
    }
    
    if (RegionsNb==99)
      cout << "Too much regions found (more than 99!). Further results may be false." << endl;
  
  this->NumberOfRegions=RegionsNb;
  
  
  for (i=0;i<RegionsNb-1;i++) for (j=i+1;j<RegionsNb;j++) if (RegionsFound[i]>RegionsFound[j]){
      tmpFl=RegionsFound[i];
      RegionsFound[i]=RegionsFound[j];
      RegionsFound[j]=tmpFl;
    }
  
  cout << endl;
  cout << RegionsNb << " regions found: "; for (i=0;i<RegionsNb;i++) cout << RegionsFound[i] << " "; cout << endl; 
  cout << endl;
  
  
  //2) initiate and treat the partition of unity
  
  this->PartitionOfUnity.CreateVoidField(NBX,NBY,NBZ,this->NumberOfRegions);
  
  for(z=0;z<NBZ;z++) for(y=0;y<NBY;y++) for(x=0;x<NBX;x++){
    //2.1) find the local ROI
    tmpInt=0;
    for (i=0;i<RegionsNb;i++) if (fabs(Part_Of_Unity->G(x,y,z)-RegionsFound[i])<epsilon) tmpInt=i;
    
    //2.2) set it in the partition of unity
    this->PartitionOfUnity.P(1,x,y,z,tmpInt);
  }
  
  
  //2.3) smooth the partition of unity
  x_mm=sqrt(this->PartitionOfUnity.Image2World[0][0]*this->PartitionOfUnity.Image2World[0][0]+this->PartitionOfUnity.Image2World[0][1]*this->PartitionOfUnity.Image2World[0][1]+this->PartitionOfUnity.Image2World[0][2]*this->PartitionOfUnity.Image2World[0][2]);
  y_mm=sqrt(this->PartitionOfUnity.Image2World[1][0]*this->PartitionOfUnity.Image2World[1][0]+this->PartitionOfUnity.Image2World[1][1]*this->PartitionOfUnity.Image2World[1][1]+this->PartitionOfUnity.Image2World[1][2]*this->PartitionOfUnity.Image2World[1][2]);
  z_mm=sqrt(this->PartitionOfUnity.Image2World[2][0]*this->PartitionOfUnity.Image2World[2][0]+this->PartitionOfUnity.Image2World[2][1]*this->PartitionOfUnity.Image2World[2][1]+this->PartitionOfUnity.Image2World[2][2]*this->PartitionOfUnity.Image2World[2][2]);
 
  TmpLightFFTconvolver.InitiateConvolver(NBX,NBY,NBZ,1,sigmaPOI/x_mm,sigmaPOI/y_mm,sigmaPOI/z_mm);
  TmpLightFFTconvolver.Convolution(&this->PartitionOfUnity);
  
  
  //3) initiate the region convolvers
  this->Region_convolver = new LightFFTconvolver3D [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++)
    this->Region_convolver[i].InitiateConvolver(NBX,NBY,NBZ,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,0);
  
  //4) initiate the temporary vector fields
  this->TempVF1.CreateVoidField(NBX,NBY,NBZ);
  this->TempVF2.CreateVoidField(NBX,NBY,NBZ);

  //5) initiate the ROI containing information in each region
  xmin= new int [this->NumberOfRegions];
  xmax= new int [this->NumberOfRegions];
  ymin= new int [this->NumberOfRegions];
  ymax= new int [this->NumberOfRegions];
  zmin= new int [this->NumberOfRegions];
  zmax= new int [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++) xmin[i]=NBX-1;
  for (i=0;i<this->NumberOfRegions;i++) ymin[i]=NBY-1;
  for (i=0;i<this->NumberOfRegions;i++) zmin[i]=NBZ-1;
  for (i=0;i<this->NumberOfRegions;i++) xmax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) ymax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) zmax[i]=1;
  
  for (i=0;i<this->NumberOfRegions;i++){
    maxValue=0;
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++)
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue) maxValue=this->PartitionOfUnity.G(x,y,z,i);
    
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue/20){
	      if (xmin[i]>x) xmin[i]=x;
	      if (ymin[i]>y) ymin[i]=y;
	      if (zmin[i]>z) zmin[i]=z;
      	
	      if (xmax[i]<x) xmax[i]=x+1;
	      if (ymax[i]<y) ymax[i]=y+1;
	      if (zmax[i]<z) zmax[i]=z+1;
      }
    }
    
    cout << "ROI containing region " << i << ": ";
    cout << "x=["<< xmin[i] << "," << xmax[i] << "], ";
    cout << "y=["<< ymin[i] << "," << ymax[i] << "], ";
    cout << "z=["<< zmin[i] << "," << zmax[i] << "]" << endl;
  }
}


///Initiate the convolver in all regions using same kernel
///-> Part_Of_Unity is a '3D + channels' scalar field which encodes the partition of unity in the different channels.
///   * Its size and number of channels (NBT actually) defines the size and the number of regions of the convolver 
///   * The maximum point-wise sum of the probabilities may be different to 1: normalisation will be automatically performed
///   * Point-wise sum of the probabilities may vary in space. If so, a background region will be automatically defined
///-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
///     * NX, NY, NZ: is the size of the input image
///     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
///     * ...
///     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void MultiRegionFFTConvolver2::InitiateConvolverWithActualPOI(ScalarField * Part_Of_Unity, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){
  int i;
  float tmpFl,maxSum,minSum,maxValue;
  int NBX,NBY,NBZ;
  int z,y,x;
  
  //1) initiate the size of the convolver and the number of regions
  NBX=Part_Of_Unity->NX;
  NBY=Part_Of_Unity->NY;
  NBZ=Part_Of_Unity->NZ;
  this->NumberOfRegions=Part_Of_Unity->NT;
  
  //2) initiate and treat the partition of unity
  
  //2.1) check whether an additional 'background channel' is necessary
  for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
    
    tmpFl=0;
    
    for (i=0;i<this->NumberOfRegions;i++) tmpFl+=Part_Of_Unity->G(x,y,z,i);
    
    if ((x==0)&&(y==0)&&(z==0)){ maxSum=tmpFl; minSum=tmpFl; }
    if (maxSum<tmpFl) maxSum=tmpFl;
    if (minSum>tmpFl) minSum=tmpFl;
  }
  
  if ((maxSum-minSum)/maxSum>0.9){ //more than 10% difference between the point-wise partions of unity
    this->NumberOfRegions++;  // => a 'background channel' will be generated
    cout << "+++ A background channel is created +++" << endl;
  }
  
  //2.2) generate and fill the partition of unity
  this->PartitionOfUnity.CreateVoidField(NBX,NBY,NBZ,this->NumberOfRegions);
  
  for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
    tmpFl=0;
    for (i=0;i<Part_Of_Unity->NT;i++){
      this->PartitionOfUnity.P(Part_Of_Unity->G(x,y,z,i)/maxSum,x,y,z,i);
      tmpFl+=Part_Of_Unity->G(x,y,z,i)/maxSum;
    }
    
    if (Part_Of_Unity->NT<this->NumberOfRegions){
       this->PartitionOfUnity.P(1-tmpFl,x,y,z,this->NumberOfRegions-1);
    }
  }
  
  //2) initiate the region convolvers
  this->Region_convolver = new LightFFTconvolver3D [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++)
    this->Region_convolver[i].InitiateConvolver(NBX,NBY,NBZ,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,0);
  
  //3) initiate the temporary vector fields
  this->TempVF1.CreateVoidField(NBX,NBY,NBZ);
  this->TempVF2.CreateVoidField(NBX,NBY,NBZ);

  //4) initiate the ROI containing information in each region
  xmin= new int [this->NumberOfRegions];
  xmax= new int [this->NumberOfRegions];
  ymin= new int [this->NumberOfRegions];
  ymax= new int [this->NumberOfRegions];
  zmin= new int [this->NumberOfRegions];
  zmax= new int [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++) xmin[i]=NBX-1;
  for (i=0;i<this->NumberOfRegions;i++) ymin[i]=NBY-1;
  for (i=0;i<this->NumberOfRegions;i++) zmin[i]=NBZ-1;
  for (i=0;i<this->NumberOfRegions;i++) xmax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) ymax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) zmax[i]=1;
  
  for (i=0;i<this->NumberOfRegions;i++){
    maxValue=0;
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++)
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue) maxValue=this->PartitionOfUnity.G(x,y,z,i);
    
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue/20){
	      if (xmin[i]>x) xmin[i]=x;
	      if (ymin[i]>y) ymin[i]=y;
	      if (zmin[i]>z) zmin[i]=z;
      	
	      if (xmax[i]<x) xmax[i]=x+1;
	      if (ymax[i]<y) ymax[i]=y+1;
	      if (zmax[i]<z) zmax[i]=z+1;
      }
    }
    
    cout << "ROI containing region " << i << ": ";
    cout << "x=["<< xmin[i] << "," << xmax[i] << "], ";
    cout << "y=["<< ymin[i] << "," << ymax[i] << "], ";
    cout << "z=["<< zmin[i] << "," << zmax[i] << "]" << endl;
  }
}

///save the actual partition of unity (after potential treatments in InitiateConvolver or undersampling)
///-> 1st char* is the name in which the image is saved
///-> 2nd char* is the name of the image that will be used in the header of the saved image
void MultiRegionFFTConvolver2::SaveActualParitionOfUnity(char * OutputImageName, char * ImageForHeaderName){
  this->PartitionOfUnity.Write(OutputImageName,ImageForHeaderName);
  }


///Update the parition of unity 
///-> it must have the same size and number of layers/times as in the current POU (tested)
///-> to make sense, the new POU must have a sum of intensities across layers/times equals to 1 at each voxel (not tested)
void MultiRegionFFTConvolver2::UpdatePartitionOfUnity(ScalarField * Part_Of_Unity){
  int z,y,x,t,i;
  
  //1) Check the size of the POU and nb of regions
  if (this->PartitionOfUnity.NX!=Part_Of_Unity->NX){cout << "Partition of unity not changed" << endl; return;}
  if (this->PartitionOfUnity.NY!=Part_Of_Unity->NY){cout << "Partition of unity not changed" << endl; return;}
  if (this->PartitionOfUnity.NZ!=Part_Of_Unity->NZ){cout << "Partition of unity not changed" << endl; return;}
  if (this->PartitionOfUnity.NT!=Part_Of_Unity->NT){cout << "Partition of unity not changed" << endl; return;}
  
  //2) fill the partition of unity
  for(t=0;t<Part_Of_Unity->NT;t++) for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++)
      this->PartitionOfUnity.P(Part_Of_Unity->G(x,y,z,t),x,y,z,t);
  
  //3) change region-based ROIs  (maybe not optimal but it is not re-estimated, so some time is gained here for sure)
  for (i=0;i<this->NumberOfRegions;i++){
    this->xmin[i]=0;   this->ymin[i]=0;  this->zmin[i]=0;
    this->xmax[i]=Part_Of_Unity->NX;  this->ymax[i]=Part_Of_Unity->NY;  this->zmax[i]=Part_Of_Unity->NZ;
  }
}


///change the smoothing kernel in one region
void MultiRegionFFTConvolver2::ChangeKernelInOneRegion(int IdRegion, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){

  this->Region_convolver[IdRegion].ChangeKernel(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,0);

}



///convolution of a 3D vector field using the predifined kernel
void MultiRegionFFTConvolver2::Convolution(VectorField * VF){
  int z,y,x,l;
  int IdRegion;
  
  //1) init
  this->TempVF2.PutToAllVoxels(0);
  
  //2) multi-region convolution...
  for (IdRegion=0;IdRegion<this->NumberOfRegions;IdRegion++){
    //... init temp VF
    DeepCopy(VF,&this->TempVF1);
    
    //... convolve with local convolver
    this->Region_convolver[IdRegion].ConvolutionInROI(&this->TempVF1,this->xmin[IdRegion],this->xmax[IdRegion],this->ymin[IdRegion],this->ymax[IdRegion],this->zmin[IdRegion],this->zmax[IdRegion]);
    
    //... update the final velocity field 
    for(z=this->zmin[IdRegion];z<this->zmax[IdRegion];z++) for(y=this->ymin[IdRegion];y<this->ymax[IdRegion];y++) for(x=this->xmin[IdRegion];x<this->xmax[IdRegion];x++) for (l=0;l<3;l++)
      this->TempVF2.Add(this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,x,y,z);
  }
  
  //3) copy the temporary smoothed VF
  DeepCopy(&this->TempVF2,VF);

}


///convolution of a 3D vector field using the predifined kernel
void MultiRegionFFTConvolver2::ConvolutionOnROI(VectorField * VF,int idRegion){
    //init temp VF
    DeepCopy(VF,&this->TempVF1);
        
    //convolve with local convolver
    this->Region_convolver[idRegion].Convolution(&this->TempVF1);
        

    // copy the temporary smoothed VF
    DeepCopy(&this->TempVF1,VF);
}



///return the number of regions considered
int MultiRegionFFTConvolver2::GetRegionNb(){
  return this->NumberOfRegions;
}




///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                  5: CLASS TO MANAGE THE MUTUAL INFORMATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///Constructor
MImanager::MImanager(){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  this->NumberOfBinsS=0;
  this->NumberOfBinsT=0;
  this->MinGreyLevelsS=0;
  this->MinGreyLevelsT=0;
  this->SizeStepsGreyLevelsS=0; 
  this->SizeStepsGreyLevelsT=0;
  this->MI=0;
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
  this->JointEntropy=0;
  this->MarginalEntropyS=0;
  this->MarginalEntropyT=0;
  this->indicatorMaskDefined=0;
  
}

///destructor
MImanager::~MImanager(void){}


//Initiate the MI manager without any mask
void MImanager::Initiate(ScalarField * SourceImage,ScalarField * TargetImage,int NbBinsSrc,int NbBinsTrg, int LocMargin){
  int i;
  int x,y,z;
  float minval,maxval;
  
  this->SrcImage=SourceImage;
  this->TrgImage=TargetImage;
  
  this->indicatorMaskDefined=0;
  
  //number of bins and margin
  this->NumberOfBinsS=NbBinsSrc+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->NumberOfBinsT=NbBinsTrg+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->Margin=LocMargin;
  
  //allocate the histograms
  this->MarginalHistogramS= new float [this->NumberOfBinsS];
  this->MarginalHistogramT= new float [this->NumberOfBinsT];
  this->JointHistogram= new float * [this->NumberOfBinsS];
  for (i=0; i<this->NumberOfBinsS; i++) 
    this->JointHistogram[i]= new float [this->NumberOfBinsT];
  
  //put to 0 the 'up-to-date indicators' (to allow the first estimations of the histograms and the MI)
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
  
  //save the size of the treated images
  this->NX=SrcImage->NX;
  this->NY=SrcImage->NY;
  this->NZ=SrcImage->NZ;
  
  //compute the minimal and maximal intensities of the source image ...
  minval=SrcImage->G(0,0,0);
  maxval=SrcImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>SrcImage->G(x,y,z)) minval=SrcImage->G(x,y,z);
      if (maxval<SrcImage->G(x,y,z)) maxval=SrcImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>SrcImage->G(x,y,0)) minval=SrcImage->G(x,y,0);
      if (maxval<SrcImage->G(x,y,0)) maxval=SrcImage->G(x,y,0);
    }
  }
  MinGreyLevelsS=minval;
  SizeStepsGreyLevelsS=(maxval*1-minval)/(this->NumberOfBinsS-4);
  
  //compute the minimal and maximal intensities of the target image
  minval=TrgImage->G(0,0,0);
  maxval=TrgImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>TrgImage->G(x,y,z)) minval=TrgImage->G(x,y,z);
      if (maxval<TrgImage->G(x,y,z)) maxval=TrgImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>TrgImage->G(x,y,0)) minval=TrgImage->G(x,y,0);
      if (maxval<TrgImage->G(x,y,0)) maxval=TrgImage->G(x,y,0);
    }
  }
  MinGreyLevelsT=minval;
  SizeStepsGreyLevelsT=(maxval*1-minval)/(this->NumberOfBinsT-4); 
}


//Initiate the MI manager with a mask. The MI and MI gradients will only be computed where the mask equals 1
void MImanager::Initiate(ScalarField * SourceImage,ScalarField * TargetImage,ScalarField * ROI_Mask,int NbBinsSrc,int NbBinsTrg, int LocMargin){
  int i;
  int x,y,z;
  float minval,maxval;
  
  //load the images and the mask
  this->SrcImage=SourceImage;
  this->TrgImage=TargetImage;
  this->Mask=ROI_Mask;
  
  //number of bins and margin
  this->NumberOfBinsS=NbBinsSrc+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->NumberOfBinsT=NbBinsTrg+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->Margin=LocMargin;
  
  //save the size of the treated images
  this->NX=SrcImage->NX;
  this->NY=SrcImage->NY;
  this->NZ=SrcImage->NZ;
  
  //check that the mask has at least 30 values equals to 1. If not, initiate the MImanager without a mask
  this->indicatorMaskDefined=0;
  for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++)
    if (fabs(Mask->G(x,y,z)-1)<0.0001) 
      this->indicatorMaskDefined++;
  
  if (this->indicatorMaskDefined<30){
    cout << "The ROI in the mask (where intensities equal 1) is too small to be considered -> MImanager initiated without a mask" << endl;
    this->Initiate(SourceImage,TargetImage,NbBinsSrc,NbBinsTrg,LocMargin);
    this->indicatorMaskDefined=0;
    return;
  }
  else{
    this->indicatorMaskDefined=1;
  }
  
  //allocate the histograms
  this->MarginalHistogramS= new float [this->NumberOfBinsS];
  this->MarginalHistogramT= new float [this->NumberOfBinsT];
  this->JointHistogram= new float * [this->NumberOfBinsS];
  for (i=0; i<this->NumberOfBinsS; i++) 
    this->JointHistogram[i]= new float [this->NumberOfBinsT];
  
  //put to 0 the 'up-to-date indicators' (to allow the first estimations of the histograms and the MI)
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
  
  
  //compute the minimal and maximal intensities of the source image ...
  minval=SrcImage->G(0,0,0);
  maxval=SrcImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,z)-1)<0.0001){
      if (minval>SrcImage->G(x,y,z)) minval=SrcImage->G(x,y,z);
      if (maxval<SrcImage->G(x,y,z)) maxval=SrcImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,0)-1)<0.0001){
      if (minval>SrcImage->G(x,y,0)) minval=SrcImage->G(x,y,0);
      if (maxval<SrcImage->G(x,y,0)) maxval=SrcImage->G(x,y,0);
    }
  }
  MinGreyLevelsS=minval;
  SizeStepsGreyLevelsS=(maxval*1-minval)/(this->NumberOfBinsS-4);
  
  //compute the minimal and maximal intensities of the target image
  minval=TrgImage->G(0,0,0);
  maxval=TrgImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,z)-1)<0.0001){
      if (minval>TrgImage->G(x,y,z)) minval=TrgImage->G(x,y,z);
      if (maxval<TrgImage->G(x,y,z)) maxval=TrgImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,0)-1)<0.0001){
      if (minval>TrgImage->G(x,y,0)) minval=TrgImage->G(x,y,0);
      if (maxval<TrgImage->G(x,y,0)) maxval=TrgImage->G(x,y,0);
    }
  }
  MinGreyLevelsT=minval;
  SizeStepsGreyLevelsT=(maxval*1-minval)/(this->NumberOfBinsT-4); 
}






//Indicate to the MI manager that the Source image has changed
void MImanager::IndicateSrcHasChanged(){
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
}

//Indicate to the MI manager that the Target image has changed
void MImanager::IndicateTrgHasChanged(){
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
}


//normalized the intensities adapted to the bins of the histograms
float MImanager::GiveFloatBinT(float intensity){
  return ((intensity-MinGreyLevelsT)/SizeStepsGreyLevelsT)+2;
}

float MImanager::GiveFloatBinS(float intensity){
  return ((intensity-MinGreyLevelsS)/SizeStepsGreyLevelsS)+2;
}


//when computing the histograms give the contribution of 'intensity' in the bin 'IdBin' for the source image
float MImanager::GiveValueParzenWindowT(float intensity,int IdBin){
  float GFBi,flId,tmp,tmp2;
  
  GFBi=GiveFloatBinT(intensity);
  flId=static_cast<float>(IdBin);
  
  return BSplineWeight(GFBi-flId);
  
  //-> Former cubic b-spline approximation:
  //if ((flId-2<GFBi)&&(GFBi<=flId-1)){
  //  return pow(GFBi-flId+2,3)/6;
  //}
  //else if ((flId-1<GFBi)&&(GFBi<=flId)){
  //  tmp=GFBi-flId+1;
  //  tmp2=tmp*tmp;
  //  return (-3*tmp2*tmp+3*tmp2+3*tmp+1)/6;
  //}
  //else if ((flId<GFBi)&&(GFBi<=flId+1)){
  //  tmp=GFBi-flId;
  //  tmp2=tmp*tmp;
  //  return (3*tmp2*tmp-6*tmp2+4)/6;
  //}
  //else if ((flId+1<GFBi)&&(GFBi<flId+2)){
  //  return pow(flId-GFBi+2,3)/6;
  //}
  //else{
  //  return 0;
  //}
  
  //-> Former linear approximation which can be inline
  //return (0.5*(2-fabs(GiveFloatBinT(intensity)-static_cast<float>(IdBin))))*static_cast<float>(fabs(GiveFloatBinT(intensity)-static_cast<float>(IdBin))<2);
}



//when computing the histograms give the contribution of 'intensity' in the bin 'IdBin' for the target image
float MImanager::GiveValueParzenWindowS(float intensity,int IdBin){
  float GFBi,flId,tmp,tmp2;
  
  GFBi=GiveFloatBinS(intensity);
  flId=static_cast<float>(IdBin);
  
  return BSplineWeight(GFBi-flId);
  
  //-> Former cubic b-spline approximation:
  //if ((flId-2<GFBi)&&(GFBi<=flId-1)){
  //  return pow(GFBi-flId+2,3)/6;
  //}
  //else if ((flId-1<GFBi)&&(GFBi<=flId)){
  //  tmp=GFBi-flId+1;
  //  tmp2=tmp*tmp;
  //  return (-3*tmp2*tmp+3*tmp2+3*tmp+1)/6;
  //}
  //else if ((flId<GFBi)&&(GFBi<=flId+1)){
  //  tmp=GFBi-flId;
  //  tmp2=tmp*tmp;
  //  return (3*tmp2*tmp-6*tmp2+4)/6;
  //}
  //else if ((flId+1<GFBi)&&(GFBi<flId+2)){
  //  return pow(flId-GFBi+2,3)/6;
  //}
  //else{
  //  return 0;
  //}
  
  //-> Former linear approximation which can be inline
  //return (0.5*(2-fabs(GiveFloatBinS(intensity)-static_cast<float>(IdBin))))*static_cast<float>(fabs(GiveFloatBinS(intensity)-static_cast<float>(IdBin))<2);
}



//Compute Joint Histogram And Entropy
void MImanager::ComputeJointHistogramAndEntropy(){
  int i,j;
  int x,y,z;
  int Total;
  int MinBinS,MinBinT;
  
  //histo
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    JointHistogram[i][j]=0.01;  //0.01 to avoid problems when computing the 
  
  
  
  if (this->indicatorMaskDefined==1){ // a mask is defined
    if (this->NZ>1){// ... 3D image
      for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,z)-1)<0.0001){
  //evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,z)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,z)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,z),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,z),j);
      }
    }
    else{// ... 2D image
      for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,0)-1)<0.0001){
        //evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,0)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,0)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,0),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,0),j);
      }
    }
  }
  else{// no mask is defined
    if (this->NZ>1){// ... 3D image
      for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
        //evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,z)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,z)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,z),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,z),j);
      }
    }
    else{// ... 2D image
      for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
        //evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,0)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,0)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,0),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,0),j);
      }
    }
  }
  
  Total=0;
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    Total+=JointHistogram[i][j];
  
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    JointHistogram[i][j]/=Total;
  
  //entropy
  this->JointEntropy=0;
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    JointEntropy+=JointHistogram[i][j]*log(JointHistogram[i][j]);
  
  //indicator
  indicatorUpdatedJointHisto=1;
}



//Compute Marginal Histogram And Entropy S
void MImanager::ComputeMarginalHistogramAndEntropyS(){
  int i,j;
  
  if (indicatorUpdatedJointHisto==0)
    this->ComputeJointHistogramAndEntropy();
  
  //histo
  for (i = 0; i < this->NumberOfBinsS; i++)
    MarginalHistogramS[i]=0;  
  
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    MarginalHistogramS[i]+=JointHistogram[i][j];
  
  //entropy
  this->MarginalEntropyS=0;
  for (i = 0; i < this->NumberOfBinsS; i++)
    MarginalEntropyS+=MarginalHistogramS[i]*log(MarginalHistogramS[i]);
  
  
  //indicator
  indicatorUpdatedSrcHisto=1;
}


//Compute Marginal Histogram And Entropy T
void MImanager::ComputeMarginalHistogramAndEntropyT(){
  int i,j;
  
  if (indicatorUpdatedJointHisto==0)
    this->ComputeJointHistogramAndEntropy();
  
  //histo
  for (i = 0; i < this->NumberOfBinsT; i++)
    MarginalHistogramT[i]=0;  
  
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    MarginalHistogramT[j]+=JointHistogram[i][j];
  
  //entropy
  this->MarginalEntropyT=0;
  for (i = 0; i < this->NumberOfBinsT; i++)
    MarginalEntropyT+=MarginalHistogramT[i]*log(MarginalHistogramT[i]);
  
  
  //indicator
  indicatorUpdatedTrgHisto=1;
}  



//returns the normalized mutual information (plus update all the histograms)
float MImanager::EvaluateMI(){

  //compute the joint histogram and joint entropy...
  if ((indicatorUpdatedMI!=1)||(indicatorUpdatedJointHisto!=1)||(indicatorUpdatedTrgHisto!=1)||(indicatorUpdatedSrcHisto!=1))
    this->ComputeJointHistogramAndEntropy();
  
  
  //compute the histgram related to the source image and related marginal entropy...
  if ((indicatorUpdatedMI!=1)||(indicatorUpdatedSrcHisto!=1))
    ComputeMarginalHistogramAndEntropyS();
  
  //compute the histgram related to the target image and related marginal entropy...
  if ((indicatorUpdatedMI!=1)||(indicatorUpdatedTrgHisto!=1))
    ComputeMarginalHistogramAndEntropyT();
  
  //compute the normalised mutual information
  if ((indicatorUpdatedMI!=1)){
    this->MI=JointEntropy-MarginalEntropyT-MarginalEntropyS;  //mutual information actually
    
    indicatorUpdatedMI=1;
  }
  
  return this->MI;
  
}



//returns the estimated gradient of normalized mutual information
//In practice: 
//  * local intensity gradients of the deformed source image are first evaluated
//  * At each point, the demons-style update is computed and return in 'Gradient'
void MImanager::EvaluateGradMI(VectorField * Gradient){
  int i, j, k;
  int x,y,z;
  float LocNorm,direcX,direcY,direcZ;
  float tmpFl,tmpX,tmpY,tmpZ;
  int SizeNgbh,SizeCheckForwardBackward;
  float S_loc,S_fwd,S_bwd,Sprime;
  float epsilon;
  float * ListIntensities_Target;
  float * ListIntensities_Centered;
  float * ListIntensities_Forward;
  float * ListIntensities_Backward;
  
  int * MinBin_Target;
  int * MinBin_Centered;
  int * MinBin_Forward;
  int * MinBin_Backward;
  
  float lambdaX;
  
  lambdaX=1;
  
  //1) Initiate and allocate
  epsilon=0.0001;
  SizeNgbh=3;
  SizeCheckForwardBackward=1;
  
  ListIntensities_Target= new float [SizeNgbh];
  ListIntensities_Centered = new float [SizeNgbh];
  ListIntensities_Forward = new float [SizeNgbh];
  ListIntensities_Backward = new float [SizeNgbh];
  
  MinBin_Target = new int [SizeNgbh];
  MinBin_Centered = new int [SizeNgbh];
  MinBin_Forward = new int [SizeNgbh];
  MinBin_Backward = new int [SizeNgbh];
  
  this->EvaluateMI();
  
  //2) Compute the gradient of the intensities in the (normaly deformed) source image   (d f(g(x,mu)) / d mu) for all x where mu is the direction in which the Gateaux derivative is the highest
  Cpt_Grad_ScalarField(this->SrcImage,Gradient);
  
  //3) Compute the amplitude of the gradient for each point in the image (minus the margins)
  if (this->NZ>1){ //3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      
      //The estimated direction of the gradients is the one of mu.
      
      //3.1) direction in which we look
      LocNorm=sqrt(Gradient->G(0,x,y,z)*Gradient->G(0,x,y,z)+Gradient->G(1,x,y,z)*Gradient->G(1,x,y,z)+Gradient->G(2,x,y,z)*Gradient->G(2,x,y,z));
      
      if (LocNorm<epsilon){ // no intensity gradient
  Gradient->P(0,0,x,y,z);
  Gradient->P(0,1,x,y,z);
  Gradient->P(0,2,x,y,z);
      }
      else{// intensity gradient OK
  direcX=Gradient->G(0,x,y,z)/LocNorm;
  direcY=Gradient->G(1,x,y,z)/LocNorm;
  direcZ=Gradient->G(2,x,y,z)/LocNorm;
  
  //3.2) intensities which will be considered to compute the gradient  ... an idea: the deformations may be weighted with the kernel size
  for (i=0;i<SizeNgbh;i++){
    tmpX=static_cast<float>(x)+(i-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
    tmpY=static_cast<float>(y)+(i-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
    tmpZ=static_cast<float>(z)+(i-(SizeNgbh/2))*direcZ*SizeCheckForwardBackward;
    ListIntensities_Target[i]=this->TrgImage->G(tmpX,tmpY,tmpZ);
    ListIntensities_Centered[i]=this->SrcImage->G(tmpX,tmpY,tmpZ);
  }
  for (i=0;i<SizeNgbh;i++){
    if (i>0) ListIntensities_Forward[i]=ListIntensities_Centered[i-1];
    if (i<SizeNgbh-1) ListIntensities_Backward[i]=ListIntensities_Centered[i+1];
  }
  
  tmpX=static_cast<float>(x)+(-1-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
  tmpY=static_cast<float>(y)+(-1-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
  tmpZ=static_cast<float>(z)+(-1-(SizeNgbh/2))*direcZ*SizeCheckForwardBackward;
  ListIntensities_Forward[0]=this->SrcImage->G(tmpX,tmpY,tmpZ);
  
  tmpX=static_cast<float>(x)+(SizeNgbh-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
  tmpY=static_cast<float>(y)+(SizeNgbh-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
  tmpZ=static_cast<float>(z)+(SizeNgbh-(SizeNgbh/2))*direcZ*SizeCheckForwardBackward;
  ListIntensities_Backward[SizeNgbh-1]=this->SrcImage->G(tmpX,tmpY,tmpZ);
  
  
  //3.3) min bins corresponding to the intensities of interest (to only check the information where it is not null)
  for (i=0;i<SizeNgbh;i++){
    MinBin_Target[i]=static_cast<int>(this->GiveFloatBinT(ListIntensities_Target[i]))-1;
    if (MinBin_Target[i]<0) MinBin_Target[i]=0;
    if (MinBin_Target[i]>=NumberOfBinsT-3) MinBin_Target[i]=NumberOfBinsT-4;
    
    MinBin_Centered[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Centered[i]))-1;
    if (MinBin_Centered[i]<0) MinBin_Centered[i]=0;
    if (MinBin_Centered[i]>=NumberOfBinsS-3) MinBin_Centered[i]=NumberOfBinsS-4;
    
    MinBin_Forward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Forward[i]))-1;
    if (MinBin_Forward[i]<0) MinBin_Forward[i]=0;
    if (MinBin_Forward[i]>=NumberOfBinsS-3) MinBin_Forward[i]=NumberOfBinsS-4;
    
    MinBin_Backward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Backward[i]))-1;
    if (MinBin_Backward[i]<0) MinBin_Backward[i]=0;
    if (MinBin_Backward[i]>=NumberOfBinsS-3) MinBin_Backward[i]=NumberOfBinsS-4;
  }
  
  
  //3.4) evaluate the derivative of the mutual information at (x,y,z) in the direction mu
  
  //3.4.1) local mutual information (in one direction and irrespective to the intensities gradient)
  S_loc=0;
  for (k=0;k<SizeNgbh;k++){
    tmpFl=0;
    for (i = MinBin_Centered[k]; i < MinBin_Centered[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
      //log_2 of the conditional probability of the bin pair [i][j]
      tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
      
      //weight of the local intensity of the target image in the parzen window
      tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
      
      //weight of the local intensity of the deformed source image in the parzen window
      tmpFl*=GiveValueParzenWindowS(ListIntensities_Centered[k],i);
      
      //update S_loc
      S_loc+=tmpFl;
    }
  }
  
  //3.4.2) local forward mutual information (in one direction and irrespective to the intensities gradient)
  S_fwd=0;
  for (k=0;k<SizeNgbh;k++){
    tmpFl=0;
    for (i = MinBin_Forward[k]; i < MinBin_Forward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
      //log_2 of the conditional probability of the bin pair [i][j]
      tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
      
      //weight of the local intensity of the target image in the parzen window
      tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
      
      //weight of the local intensity of the deformed source image in the parzen window
      tmpFl*=GiveValueParzenWindowS(ListIntensities_Forward[k],i);
      
      //update S_fwd
      S_fwd+=tmpFl;
    }
  }
  
  //3.4.3) local backward mutual information (in one direction and irrespective to the intensities gradient)
  S_bwd=0;
  for (k=0;k<SizeNgbh;k++){
    tmpFl=0;
    for (i = MinBin_Backward[k]; i < MinBin_Backward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
      //log_2 of the conditional probability of the bin pair [i][j]
      tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
      
      //weight of the local intensity of the target image in the parzen window
      tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
      
      //weight of the local intensity of the deformed source image in the parzen window
      tmpFl*=GiveValueParzenWindowS(ListIntensities_Backward[k],i);
      
      //update S_bwd
      S_bwd+=tmpFl;
    }
  }
  
  
  //3.5) define the local gradient of mutual information   (the S_... values are supposed of the same sign / normaly < 0)
  if ((S_fwd-S_loc)<(S_bwd-S_loc)){ //best score in the forward direction
    Sprime=S_fwd-S_loc;
  }
  else{
    Sprime=S_loc-S_bwd;
  }
  
  Gradient->P(Gradient->G(0,x,y,z)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),0,x,y,z);
  Gradient->P(Gradient->G(1,x,y,z)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),1,x,y,z);
  Gradient->P(Gradient->G(2,x,y,z)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),2,x,y,z);
      }
    }
  }
  else{ //3.bis: 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      
      //The estimated direction of the gradients is the one of mu.
      
      //3.bis.1) direction in which we look
      LocNorm=sqrt(Gradient->G(0,x,y,0)*Gradient->G(0,x,y,0)+Gradient->G(1,x,y,0)*Gradient->G(1,x,y,0));
      
      if (LocNorm<epsilon){ // no intensity gradient
  Gradient->P(0,0,x,y,0);
  Gradient->P(0,1,x,y,0);
  Gradient->P(0,2,x,y,0);
      }
      else{// intensity gradient OK
  direcX=Gradient->G(0,x,y,0)/LocNorm;
  direcY=Gradient->G(1,x,y,0)/LocNorm;
  direcZ=0;
  
  //3.bis.2) intensities which will be considered to compute the gradient  ... an idea: the deformations may be weighted with the kernel size
  for (i=0;i<SizeNgbh;i++){
    tmpX=static_cast<float>(x)+(i-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
    tmpY=static_cast<float>(y)+(i-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
    tmpZ=0;
    ListIntensities_Target[i]=this->TrgImage->G(tmpX,tmpY,tmpZ);
    ListIntensities_Centered[i]=this->SrcImage->G(tmpX,tmpY,tmpZ);
  }
  for (i=0;i<SizeNgbh;i++){
    if (i>0) ListIntensities_Forward[i]=ListIntensities_Centered[i-1];
    if (i<SizeNgbh-1) ListIntensities_Backward[i]=ListIntensities_Centered[i+1];
  }
  
  tmpX=static_cast<float>(x)+(-1-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
  tmpY=static_cast<float>(y)+(-1-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
  tmpZ=0;
  ListIntensities_Forward[0]=this->SrcImage->G(tmpX,tmpY,tmpZ);
  
  tmpX=static_cast<float>(x)+(SizeNgbh-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
  tmpY=static_cast<float>(y)+(SizeNgbh-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
  tmpZ=0;
  ListIntensities_Backward[SizeNgbh-1]=this->SrcImage->G(tmpX,tmpY,tmpZ);
  
  
  //3.bis.3) min bins corresponding to the intensities of interest (to only check the information where it is not null)
  for (i=0;i<SizeNgbh;i++){
    MinBin_Target[i]=static_cast<int>(this->GiveFloatBinT(ListIntensities_Target[i]))-1;
    if (MinBin_Target[i]<0) MinBin_Target[i]=0;
    if (MinBin_Target[i]>=NumberOfBinsT-3) MinBin_Target[i]=NumberOfBinsT-4;
    
    MinBin_Centered[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Centered[i]))-1;
    if (MinBin_Centered[i]<0) MinBin_Centered[i]=0;
    if (MinBin_Centered[i]>=NumberOfBinsS-3) MinBin_Centered[i]=NumberOfBinsS-4;
    
    MinBin_Forward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Forward[i]))-1;
    if (MinBin_Forward[i]<0) MinBin_Forward[i]=0;
    if (MinBin_Forward[i]>=NumberOfBinsS-3) MinBin_Forward[i]=NumberOfBinsS-4;
    
    MinBin_Backward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Backward[i]))-1;
    if (MinBin_Backward[i]<0) MinBin_Backward[i]=0;
    if (MinBin_Backward[i]>=NumberOfBinsS-3) MinBin_Backward[i]=NumberOfBinsS-4;
  }
  
  
  //3.bis.4) evaluate the derivative of the mutual information at (x,y,0) in the direction mu
  
  //3.bis.4.1) local mutual information (in one direction and irrespective to the intensities gradient)
  S_loc=0;
  for (k=0;k<SizeNgbh;k++){
    tmpFl=0;
    for (i = MinBin_Centered[k]; i < MinBin_Centered[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
      //log_2 of the conditional probability of the bin pair [i][j]
      tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
      
      //weight of the local intensity of the target image in the parzen window
      tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
      
      //weight of the local intensity of the deformed source image in the parzen window
      tmpFl*=GiveValueParzenWindowS(ListIntensities_Centered[k],i);
      
      //update S_loc
      S_loc+=tmpFl;
    }
  }
  
  //3.bis.4.2) local forward mutual information (in one direction and irrespective to the intensities gradient)
  S_fwd=0;
  for (k=0;k<SizeNgbh;k++){
    tmpFl=0;
    for (i = MinBin_Forward[k]; i < MinBin_Forward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
      //log_2 of the conditional probability of the bin pair [i][j]
      tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
      
      //weight of the local intensity of the target image in the parzen window
      tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
      
      //weight of the local intensity of the deformed source image in the parzen window
      tmpFl*=GiveValueParzenWindowS(ListIntensities_Forward[k],i);
      
      //update S_fwd
      S_fwd+=tmpFl;
    }
  }
  
  //3.bis.4.3) local backward mutual information (in one direction and irrespective to the intensities gradient)
  S_bwd=0;
  for (k=0;k<SizeNgbh;k++){
    tmpFl=0;
    for (i = MinBin_Backward[k]; i < MinBin_Backward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
      //log_2 of the conditional probability of the bin pair [i][j]
      tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
      
      //weight of the local intensity of the target image in the parzen window
      tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
      
      //weight of the local intensity of the deformed source image in the parzen window
      tmpFl*=GiveValueParzenWindowS(ListIntensities_Backward[k],i);
      
      //update S_bwd
      S_bwd+=tmpFl;
    }
  }
  
  
  //3.bis.5) define the local gradient of mutual information   (the S_... values are supposed of the same sign / normaly < 0)
  if ((S_fwd-S_loc)<(S_bwd-S_loc)){ //best score in the forward direction
    Sprime=S_fwd-S_loc;
  }
  else{
    Sprime=S_loc-S_bwd;
  }
  
  Gradient->P(Gradient->G(0,x,y,0)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),0,x,y,0);
  Gradient->P(Gradient->G(1,x,y,0)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),1,x,y,0);
      }
    }
  }
  
  
  //4) put to zero the gradients that are not taken into account...
  //4.1) ... in the margin
  if (this->NZ>1){//... 3D
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      if ((z<this->Margin)||(z>=this->NZ-this->Margin)||(y<this->Margin)||(y>=this->NY-this->Margin)||(x<this->Margin)||(x>=this->NX-this->Margin)){
        Gradient->P(0,0,x,y,z);
        Gradient->P(0,1,x,y,z);
        Gradient->P(0,2,x,y,z);
      }
    }
  }
  else{//... 2D
    for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      if ((y<this->Margin)||(y>=this->NY-this->Margin)||(x<this->Margin)||(x>=this->NX-this->Margin)){
        Gradient->P(0,0,x,y,0);
        Gradient->P(0,1,x,y,0);
        Gradient->P(0,2,x,y,0);
      }
    }
  }
  //4.2) ... outside of the ROI
  if (indicatorMaskDefined==1){
    if (this->NZ>1){//... 3D
      for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (fabs(Mask->G(x,y,z)-1)>=0.0001){
        Gradient->P(0,0,x,y,z);
        Gradient->P(0,1,x,y,z);
        Gradient->P(0,2,x,y,z);
      }
    }
    else{//... 2D
      for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (fabs(Mask->G(x,y,0)-1)>=0.0001){
        Gradient->P(0,0,x,y,0);
        Gradient->P(0,1,x,y,0);
        Gradient->P(0,2,x,y,0);
      }
    }
  }
  
  //5) delete allocated variables
  delete ListIntensities_Target;
  delete ListIntensities_Centered;
  delete ListIntensities_Forward;
  delete ListIntensities_Backward;
  
  delete MinBin_Target;
  delete MinBin_Centered;
  delete MinBin_Forward;
  delete MinBin_Backward;
}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           6: LOW LEVEL FUNCTIONS MAKING USE OF THE CLASSES ScalarField AND VectorField 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///Isotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
void Diffusion_3D(ScalarField * SField, float alpha, float dTau,int ITERATIONS_NB, float dx,float dy, float dz){
  int x, y, z;
  int NBX,NBY,NBZ;
  float *Va;
  float *Vb; 
  float *Vc;
  float *Vd;
  float *Vx;
  int n;
  int iteration;
  float Var1x,Var2x,Var3x;
  float Var1y,Var2y,Var3y;
  float Var1z,Var2z,Var3z;
  
  
  //IF alpha == 1, it can be interesting to use the following strategy (will however diffuse the image boundaries)
  //FFTconvolver3D convolver;
  //if (alpha==1){
  //  convolver.InitiateConvolver(SField->NX,SField->NY,SField->NZ,1,(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB));
  //  convolver.Convolution(SField);
  //}
  //return;
  
  
  //1) INITIALISATION
  
  //variables definition
  NBX=SField->NX;
  NBY=SField->NY;
  NBZ=SField->NZ;
  
  //precomputed values
  Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
  Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
  Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
  
  //temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
  n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
  Va= new float [n];
  Vb= new float [n];
  Vc= new float [n];
  Vd= new float [n];
  Vx= new float [n];
  
  //2) ISOTROPIC DIFFUSION
  
  for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
    //cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
    
    //2.1) diffusion - x direction
    for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
      for (x = 0; x < NBX; x++){
        Va[x+1]=Var2x;
        Vb[x+1]=Var1x;
        Vc[x+1]=Var2x;
        Vd[x+1]=SField->G(x,y,z);
      }
      Va[0]=Va[1]; Va[NBX+1]=Va[NBX]; //to avoid boundary effects
      Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX]; //to avoid boundary effects
      Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX]; //to avoid boundary effects
      Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
      for (x = 0; x < NBX; x++) SField->P(Vx[x+1],x,y,z);
    }
    
    //2.2) diffusion - y direction
    for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
      for (y = 0; y < NBY; y++){
        Va[y+1]=Var2y;
        Vb[y+1]=Var1y;
        Vc[y+1]=Var2y;
        Vd[y+1]=SField->G(x,y,z);
      }
      Va[0]=Va[1]; Va[NBY+1]=Va[NBY]; //to avoid boundary effects
      Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY]; //to avoid boundary effects
      Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY]; //to avoid boundary effects
      Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
      for (y = 0; y < NBY; y++) SField->P(Vx[y+1],x,y,z);
    }
    
    //2.3) diffusion - z direction
    if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
      for (z = 0; z < NBZ; z++){
        Va[z+1]=Var2z;
        Vb[z+1]=Var1z;
        Vc[z+1]=Var2z;
        Vd[z+1]=SField->G(x,y,z);
      }
      Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ]; //to avoid boundary effects
      Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ]; //to avoid boundary effects
      Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ]; //to avoid boundary effects
      Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
      for (z = 0; z < NBZ; z++) SField->P(Vx[z+1],x,y,z);
    }
    
  }
  
  
}


///Isotropic diffusion of 'SField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
void Diffusion_3D(ScalarField * SField,ScalarField * Mask,int MaskId, float alpha, float dTau,int ITERATIONS_NB,int optNoMaskNoDef, float dx,float dy, float dz){
  int x, y, z;
  int NBX,NBY,NBZ;
  float *Va;
  float *Vb; 
  float *Vc;
  float *Vd;
  float *Vx;
  int n;
  int iteration;
  float Var1x,Var2x,Var3x;
  float Var1y,Var2y,Var3y;
  float Var1z,Var2z,Var3z;
  float epsilon;
  float FlMaskId;
  
  //1) INITIALISATION
  
  //variables definition
  NBX=SField->NX;
  NBY=SField->NY;
  NBZ=SField->NZ;
  epsilon=0.00001;  //intensities below this value are considered as null in the mask
  FlMaskId=static_cast<float>(MaskId);
  
  //precomputed values
  Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
  Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
  Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
  
  //temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
  n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
  Va= new float [n];
  Vb= new float [n];
  Vc= new float [n];
  Vd= new float [n];
  Vx= new float [n];
  
  //2) ISOTROPIC DIFFUSION
  
  for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
    //cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
    
    //2.1) diffusion - x direction
    for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
      //main values of the diagonal matrix
      for (x = 0; x < NBX; x++){
        if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
          Va[x+1]=0; 
          Vb[x+1]=1; 
          Vc[x+1]=0;
        }
        else {
          Va[x+1]=Var2x; 
          Vb[x+1]=Var1x; 
          Vc[x+1]=Var2x;
        }
        Vd[x+1]=SField->G(x,y,z);
        
      }
      //boundaries
      Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
      Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
      Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
      Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
      //matrix inversion
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
      //store result
      for (x = 0; x < NBX; x++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) SField->P(Vx[x+1],x,y,z);
    }
    
    //2.2) diffusion - y direction
    for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
      //main values of the diagonal matrix
      for (y = 0; y < NBY; y++){
        if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
          Va[y+1]=0; 
          Vb[y+1]=1; 
          Vc[y+1]=0;
        }
        else {
          Va[y+1]=Var2y; 
          Vb[y+1]=Var1y; 
          Vc[y+1]=Var2y;
        }
        Vd[y+1]=SField->G(x,y,z);
        
      }
      //boundaries
      Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
      Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
      Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
      Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
      //matrix inversion
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
      //store result
      for (y = 0; y < NBY; y++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) SField->P(Vx[y+1],x,y,z);
    }
    
    //2.3) diffusion - z direction
    if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
      //main values of the diagonal matrix
      for (z = 0; z < NBZ; z++){
        if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
          Va[z+1]=0; 
          Vb[z+1]=1; 
          Vc[z+1]=0;
        }
        else {
          Va[z+1]=Var2z; 
          Vb[z+1]=Var1z; 
          Vc[z+1]=Var2z;
        }
        Vd[z+1]=SField->G(x,y,z);
      }
      //boundaries
      Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
      Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
      Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
      Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
      //matrix inversion
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
      //store result
      for (z = 0; z < NBZ; z++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) SField->P(Vx[z+1],x,y,z);
    }
  }
  
  
  //3) Manage the option 'optNoMaskNoDef'
  if (optNoMaskNoDef==1)
    for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
      if(Mask->G(x,y,z)<epsilon)
        SField->P(0,x,y,z);
}







///Isotropic diffusion of 'VField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
void Diffusion_3D(VectorField * VField,float alpha, float dTau,int ITERATIONS_NB,float dx,float dy, float dz){
  int x, y, z;
  int NBX,NBY,NBZ;
  float *Va;
  float *Vb; 
  float *Vc;
  float *Vd;
  float *Vx;
  int n;
  int iteration;
  float Var1x,Var2x,Var3x;
  float Var1y,Var2y,Var3y;
  float Var1z,Var2z,Var3z;
  float epsilon;
  int IdDirec;
  
  
  
  
  //IF alpha == 1, it can be interesting to use the following strategy (will however diffuse the image boundaries)
  //  FFTconvolver3D convolver;
  //  if (alpha==1){
  //  cout << "fft" << endl;
  //    convolver.InitiateConvolver(VField->NX,VField->NY,VField->NZ,1,(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB));
  //
  //    for (IdDirec=0;IdDirec<3;IdDirec++){
  //      
  //     for (z = 0; z < VField->NZ; z++) for (y = 0; y < VField->NY; y++) for (x = 0; x < VField->NX; x++) 
  //       convolver.P(VField->G(IdDirec,x,y,z), x, y, z);
  //       
  //     convolver.Convolution();
  //     
  //      for (z = 0; z < VField->NZ; z++) for (y = 0; y < VField->NY; y++) for (x = 0; x < VField->NX; x++) 
  ////        VField->P(convolver.G(x,y,z),IdDirec, x, y,  z);
  //    }
  //  }
  // return;
  
  
  //1) INITIALISATION
  
  //variables definition
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  epsilon=0.00001;  //intensities below this value are considered as null in the mask
  
  //precomputed values
  Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
  Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
  Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
  
  //temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
  n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
  Va= new float [n];
  Vb= new float [n];
  Vc= new float [n];
  Vd= new float [n];
  Vx= new float [n];
   
  //2) ISOTROPIC DIFFUSION
  
  for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
    //cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
    for (IdDirec=0;IdDirec<3;IdDirec++){
      //2.1) diffusion - x direction
      for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
        //main values of the diagonal matrix
        for (x = 0; x < NBX; x++){
          Va[x+1]=Var2x; 
          Vb[x+1]=Var1x; 
           Vc[x+1]=Var2x;
          Vd[x+1]=VField->G(IdDirec,x,y,z);
        }
        //boundaries
        Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
        Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
        Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
        Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
        //store result
        for (x = 0; x < NBX; x++) VField->P(Vx[x+1],IdDirec,x,y,z); //only store the result in the current ROI
      }
      //2.2) diffusion - y direction
      for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
        //main values of the diagonal matrix
        for (y = 0; y < NBY; y++){
          Va[y+1]=Var2y; 
           Vb[y+1]=Var1y; 
           Vc[y+1]=Var2y;
          Vd[y+1]=VField->G(IdDirec,x,y,z);
        }
        //boundaries
        Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
        Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
        Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
        Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
        //store result
        for (y = 0; y < NBY; y++)  VField->P(Vx[y+1],IdDirec,x,y,z);
      }
      
      //2.3) diffusion - z direction
      if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
        //main values of the diagonal matrix
        for (z = 0; z < NBZ; z++){
           Va[z+1]=Var2z; 
          Vb[z+1]=Var1z; 
          Vc[z+1]=Var2z;
          Vd[z+1]=VField->G(IdDirec,x,y,z);
        }
        //boundaries
        Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
        Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
        Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
        Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
        //store result
        for (z = 0; z < NBZ; z++) VField->P(Vx[z+1],IdDirec,x,y,z);
      }
    }
  }
  
}





///Isotropic diffusion of 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
///If direction > 0 then the smoothing is only performed in the direction {1=X, 2=Y or 3=Z}
void Diffusion_3D(VectorField * VField,ScalarField * Mask,int MaskId, float alpha, float dTau,int ITERATIONS_NB,  int optNoMaskNoDef,int direction, float dx,float dy, float dz){
  int x, y, z;
  int NBX,NBY,NBZ;
  float *Va;
  float *Vb; 
  float *Vc;
  float *Vd;
  float *Vx;
  int n;
  int iteration;
  float Var1x,Var2x,Var3x;
  float Var1y,Var2y,Var3y;
  float Var1z,Var2z,Var3z;
  float epsilon;
  float FlMaskId;
  int IdDirec;
  
  
  

  
  //1) INITIALISATION
  
  //variables definition
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  epsilon=0.00001;  //intensities below this value are considered as null in the mask
  FlMaskId=static_cast<float>(MaskId);
  
  //precomputed values
  Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
  Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
  Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
  
  //temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
  n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
  Va= new float [n];
  Vb= new float [n];
  Vc= new float [n];
  Vd= new float [n];
  Vx= new float [n];
   
  //2) ISOTROPIC DIFFUSION
  for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
    //cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
    for (IdDirec=0;IdDirec<3;IdDirec++) if ((direction<0)||(direction==IdDirec)) {
      //2.1) diffusion - x direction
      for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
        //main values of the diagonal matrix
        for (x = 0; x < NBX; x++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[x+1]=0; 
            Vb[x+1]=1; 
            Vc[x+1]=0;
            Vd[x+1]=VField->G(0,x,y,z);
          }
          else {
            Va[x+1]=Var2x; 
            Vb[x+1]=Var1x; 
            Vc[x+1]=Var2x;
            Vd[x+1]=VField->G(IdDirec,x,y,z);
          }
          //region boundaries
          if ((fabs(Mask->G(x,y,z)-FlMaskId)>epsilon)&&(x<NBX-3)&&(x>2)){
            if      (fabs(Mask->G(x-1,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x-1,y,z);
            else if (fabs(Mask->G(x+1,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x+1,y,z);
            else if (fabs(Mask->G(x-2,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x-2,y,z);
            else if (fabs(Mask->G(x+2,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x+2,y,z);
            else if (fabs(Mask->G(x-3,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x-3,y,z);
            else if (fabs(Mask->G(x+3,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x+3,y,z);
          }
        }
        //boundaries
        Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
        Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
        Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
        Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
        //store result
        for (x = 0; x < NBX; x++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[x+1],IdDirec,x,y,z); //only store the result in the current ROI
      }
      //2.2) diffusion - y direction
      for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
        //main values of the diagonal matrix
        for (y = 0; y < NBY; y++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[y+1]=0; 
            Vb[y+1]=1; 
            Vc[y+1]=0;
            Vd[y+1]=VField->G(0,x,y,z);
          }
          else {
            Va[y+1]=Var2y; 
            Vb[y+1]=Var1y; 
            Vc[y+1]=Var2y;
            Vd[y+1]=VField->G(IdDirec,x,y,z);
          }
          //region boundaries
          if ((fabs(Mask->G(x,y,z)-FlMaskId)>epsilon)&&(y<NBY-3)&&(y>2)){
            if      (fabs(Mask->G(x,y-1,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y-1,z);
            else if (fabs(Mask->G(x,y+1,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y+1,z);
            else if (fabs(Mask->G(x,y-2,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y-2,z);
            else if (fabs(Mask->G(x,y+2,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y+2,z);
            else if (fabs(Mask->G(x,y-3,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y-3,z);
            else if (fabs(Mask->G(x,y+3,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y+3,z);
          }
        }
        //boundaries
        Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
        Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
        Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
        Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
        //store result
        for (y = 0; y < NBY; y++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon)  VField->P(Vx[y+1],IdDirec,x,y,z);
      }
      
      //2.3) diffusion - z direction
      if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
        //main values of the diagonal matrix
        for (z = 0; z < NBZ; z++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[z+1]=0; 
            Vb[z+1]=1; 
            Vc[z+1]=0;
            Vd[z+1]=VField->G(0,x,y,z);
          }
          else {
            Va[z+1]=Var2z; 
            Vb[z+1]=Var1z; 
            Vc[z+1]=Var2z;
            Vd[z+1]=VField->G(IdDirec,x,y,z);
          }
          //region boundaries
          if ((fabs(Mask->G(x,y,z)-FlMaskId)>epsilon)&&(z<NBZ-3)&&(z>2)){
            if      (fabs(Mask->G(x,y,z-1)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z-1);
            else if (fabs(Mask->G(x,y,z+1)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z+1);
            else if (fabs(Mask->G(x,y,z-2)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z-2);
            else if (fabs(Mask->G(x,y,z+2)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z+2);
            else if (fabs(Mask->G(x,y,z-3)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z-3);
            else if (fabs(Mask->G(x,y,z+3)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z+3);
          }
        }

        //boundaries
        Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
        Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
        Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
        Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
        //store result
        for (z = 0; z < NBZ; z++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon)  VField->P(Vx[z+1],IdDirec,x,y,z);
      }
    }
  }

  //3) Manage the option 'optNoMaskNoDef'
  if (optNoMaskNoDef==1)
    for (IdDirec=0;IdDirec<3;IdDirec++) for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
      if(Mask->G(x,y,z)<epsilon)
        VField->P(0,IdDirec,x,y,z);
}









///Isotropic diffusion of 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///!!! -> Dirichlet conditions are considered at the domain boundaries !!!
void Diffusion_3Dbis(VectorField * VField,ScalarField * Mask,int MaskId,float alpha, float dTau,int ITERATIONS_NB, int optNoMaskNoDef, float dx,float dy, float dz){
  int x, y, z;
  int NBX,NBY,NBZ;
  float *Va;
  float *Vb; 
  float *Vc;
  float *Vd;
  float *Vx;
  int n;
  int iteration;
  float Var1x,Var2x,Var3x;
  float Var1y,Var2y,Var3y;
  float Var1z,Var2z,Var3z;
  float epsilon;
  float FlMaskId;
  int IdDirec;
  
  
  
  
  //1) INITIALISATION
  
  //variables definition
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  epsilon=0.00001;  //intensities below this value are considered as null in the mask
  FlMaskId=static_cast<float>(MaskId);
  
  //precomputed values
  Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
  Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
  Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
  
  //temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
  n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
  Va= new float [n];
  Vb= new float [n];
  Vc= new float [n];
  Vd= new float [n];
  Vx= new float [n];
   
  //2) ISOTROPIC DIFFUSION
  
  for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
    //cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
    for (IdDirec=0;IdDirec<3;IdDirec++){
      //2.1) diffusion - x direction
      for (z = 2; z < NBZ; z++) for (y = 0; y < NBY; y++) {
        //main values of the diagonal matrix
        for (x = 0; x < NBX; x++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[x+1]=0; 
            Vb[x+1]=1; 
            Vc[x+1]=0;
            if (optNoMaskNoDef==1) Vd[x+1]=0;
            else Vd[x+1]=VField->G(IdDirec,x,y,z);
          }
          else {
            Va[x+1]=Var2x; 
            Vb[x+1]=Var1x; 
            Vc[x+1]=Var2x;
            Vd[x+1]=VField->G(IdDirec,x,y,z);
          }
          
        }
        //boundaries
        Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
        Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
        Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
        Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
        //store result
        for (x = 0; x < NBX; x++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[x+1],IdDirec,x,y,z);
      }
      //2.2) diffusion - y direction
      for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
        //main values of the diagonal matrix
        for (y = 0; y < NBY; y++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[y+1]=0; 
            Vb[y+1]=1; 
            Vc[y+1]=0;
            if (optNoMaskNoDef==1) Vd[y+1]=0;
            else Vd[y+1]=VField->G(IdDirec,x,y,z);
          }
          else {
            Va[y+1]=Var2y; 
            Vb[y+1]=Var1y; 
            Vc[y+1]=Var2y;
            Vd[y+1]=VField->G(IdDirec,x,y,z);
          }
        }
        //boundaries
        Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
        Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
        Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
        Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
        //store result
        for (y = 0; y < NBY; y++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[y+1],IdDirec,x,y,z);
      }
      
      //2.3) diffusion - z direction
      if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
        //main values of the diagonal matrix
        for (z = 0; z < NBZ; z++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[z+1]=0; 
            Vb[z+1]=1;
            Vc[z+1]=0;
            if (optNoMaskNoDef==1) Vd[z+1]=0;
            else Vd[z+1]=VField->G(IdDirec,x,y,z);
          }
          else {
            Va[z+1]=Var2z; 
            Vb[z+1]=Var1z; 
            Vc[z+1]=Var2z;
            Vd[z+1]=VField->G(IdDirec,x,y,z);
          }
        }
        //boundaries
        Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
        Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
        Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
        Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
        //matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
        //store result
        for (z = 0; z < NBZ; z++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[z+1],IdDirec,x,y,z);
      }
    }
  }
  
  //3) Manage the option 'optNoMaskNoDef'
  if (optNoMaskNoDef==1)
    for (IdDirec=0;IdDirec<3;IdDirec++) for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
      if(Mask->G(x,y,z)<epsilon)
        VField->P(0,IdDirec,x,y,z);
  
}


///Anisotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Directions and intensities of diffusion are defined using ax, ay, az.
///The size of each voxel is defined by dx, dy, dz.
///Remark: very close to the classic paper of Perona-Malik
void anisoDiff_3D(ScalarField * SField,float ax, float ay, float az, float dTau,int ITERATIONS_NB, float dx,float dy, float dz){
  int x, y, z;
  ScalarField imageE;
  ScalarField imageO;
  int NBX,NBY,NBZ;
  double dIdx,dIdy,dIdz;
  float DivDgradI;
  float *Va;
  float *Vb; 
  float *Vc;
  float *Vd;
  float *Vx;
  int n;
  float Dxx_div_dxSq,Dyy_div_dySq,Dzz_div_dzSq;
  int iteration;
  float DivPowDxSqu,DivPowDySqu,DivPowDzSqu;
  
  //1) INITIALISATION
  
  //variables definition
  NBX=SField->NX+2;  //for boundary effects
  NBY=SField->NY+2;  //for boundary effects
  NBZ=SField->NZ+2;  //for boundary effects
  
  //temporary input and output images
  imageE.CreateVoidField(NBX,NBY,NBZ);
  imageO.CreateVoidField(NBX,NBY,NBZ);
  
  //precomputed values
  DivPowDxSqu=1./pow(dx,2);
  DivPowDySqu=1./pow(dy,2);
  DivPowDzSqu=1./pow(dz,2);
  
  
  //temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
  n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
  Va= new float [n];
  Vb= new float [n];
  Vc= new float [n];
  Vd= new float [n];
  Vx= new float [n];
  
  //2) ANISOTROPIC DIFFUSION
  
  //2.1) copy the values of the input image in a temporary 3D image
  for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
    imageE.P(SField->G(x, y, z),x+1,y+1,z+1);
  
  for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
    imageO.P(SField->G(x, y, z),x+1,y+1,z+1);
  
  
  //image extension to avoid boundary effects
  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE.P(imageE.G(1,y,z),0,y,z);
  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE.P(imageE.G(NBX-2,y,z),NBX-1,y,z);
  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO.P(imageO.G(1,y,z),0,y,z);
  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO.P(imageO.G(NBX-2,y,z),NBX-1,y,z);
  
  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,1,z),x,0,z);
  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,NBY-2,z),x,NBY-1,z);
  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,1,z),x,0,z);
  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,NBY-2,z),x,NBY-1,z);
  
  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,y,1),x,y,0);
  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,y,NBZ-2),x,y,NBZ-1);
  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,y,1),x,y,0);
  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,y,NBZ-2),x,y,NBZ-1);
  
  //2.2) diffusion in the temporary 3D image - ADI semi implicit scheme
  for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
    cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
    
    //2.2.2) diffusion - x implicit / y,z explicit
    //2.2.2.1) explicit part
    for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
      dIdy=(imageE.G(x,y+1,z)-imageE.G(x,y-1,z))/(2*dy);
      Dyy_div_dySq=(1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu;
      dIdz=(imageE.G(x,y,z+1)-imageE.G(x,y,z-1))/(2*dz);
      Dzz_div_dzSq=(1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu;
      //new value of the voxel
      DivDgradI=(imageE.G(x,y+1,z)-2*imageE.G(x,y,z)+imageE.G(x,y-1,z))*Dyy_div_dySq+
      (imageE.G(x,y,z+1)-2*imageE.G(x,y,z)+imageE.G(x,y,z-1))*Dzz_div_dzSq;
      
      imageO.P(imageE.G(x,y,z)+(dTau/3.)*DivDgradI,x,y,z);
    }
    
    //2.2.2.2) implicit part
    for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) {
      for (x = 1; x < NBX-1; x++){
        dIdx=(imageE.G(x+1,y,z)-imageE.G(x-1,y,z))/(2*dx);
        Dxx_div_dxSq=(1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu;
        Va[x+1]=(dTau/3.)*Dxx_div_dxSq;
        Vb[x+1]=-1-2*(dTau/3.)*Dxx_div_dxSq;
        Vc[x+1]=(dTau/3.)*Dxx_div_dxSq;
        Vd[x+1]=imageE.G(x,y,z); //why not imageO ???
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
      for (x = 1; x < NBX-1; x++) imageO.P(-Vx[x+1],x,y,z);
    }
    
    //2.2.3) diffusion - y implicit / x,z explicit
    //2.2.3.1) explicit part
    for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
      dIdx=(imageO.G(x+1,y,z)-imageO.G(x-1,y,z))/(2*dx);
      Dxx_div_dxSq=(1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu;
      dIdz=(imageO.G(x,y,z+1)-imageO.G(x,y,z-1))/(2*dz);
      Dzz_div_dzSq=(1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu;
      //new value of the voxel
      DivDgradI=(imageO.G(x+1,y,z)-2*imageO.G(x,y,z)+imageO.G(x-1,y,z))*Dxx_div_dxSq+
      (imageO.G(x,y,z+1)-2*imageO.G(x,y,z)+imageO.G(x,y,z-1))*Dzz_div_dzSq;
      
      imageE.P(imageO.G(x,y,z)+(dTau/3.)*DivDgradI,x,y,z);
    }
    
    //2.2.3.2) implicit part
    for (z = 1; z < NBZ-1; z++) for (x = 1; x < NBX-1; x++){
      for (y = 1; y < NBY-1; y++){
        dIdy=(imageO.G(x,y+1,z)-imageO.G(x,y-1,z))/(2*dy);
        Dyy_div_dySq=(1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu;
        Va[y+1]=(dTau/3.)*Dyy_div_dySq;
        Vb[y+1]=-1-2*(dTau/3.)*Dyy_div_dySq;
        Vc[y+1]=(dTau/3.)*Dyy_div_dySq;
        Vd[y+1]=imageO.G(x,y,z);
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBY]=Va[NBY-2]; Va[NBY+1]=Va[NBY-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBY]=Vb[NBY-2]; Vb[NBY+1]=Vb[NBY-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBY]=Vc[NBY-2]; Vc[NBY+1]=Vc[NBY-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBY]=Vd[NBY-2]; Vd[NBY+1]=Vd[NBY-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
      for (y = 1; y < NBY-1; y++) imageE.P(-Vx[y+1],x,y,z);
    }
    
    //2.2.4) diffusion - z implicit / x,y explicit
    //2.2.4.1) explicit part
    for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
      dIdx=(imageE.G(x+1,y,z)-imageE.G(x-1,y,z))/(2*dx);
      Dxx_div_dxSq=(1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu;
      dIdy=(imageE.G(x,y+1,z)-imageE.G(x,y-1,z))/(2*dy);
      Dyy_div_dySq=(1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu;
      //new value of the voxel
      DivDgradI=(imageE.G(x+1,y,z)-2*imageE.G(x,y,z)+imageE.G(x-1,y,z))*Dxx_div_dxSq+
      (imageE.G(x,y+1,z)-2*imageE.G(x,y,z)+imageE.G(x,y-1,z))*Dyy_div_dySq;
      
      imageO.P(imageE.G(x,y,z)+(dTau/3.)*DivDgradI,x,y,z);
    }
    
    //2.2.4.2) implicit part
    for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
      for (z = 1; z < NBZ-1; z++){
        dIdz=(imageE.G(x,y,z+1)-imageE.G(x,y,z-1))/(2*dz);
        Dzz_div_dzSq=(1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu;
        Va[z+1]=(dTau/3.)*Dzz_div_dzSq;
        Vb[z+1]=-1-2*(dTau/3.)*Dzz_div_dzSq;
        Vc[z+1]=(dTau/3.)*Dzz_div_dzSq;
        Vd[z+1]=imageE.G(x,y,z);
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBZ]=Va[NBZ-2]; Va[NBZ+1]=Va[NBZ-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBZ]=Vb[NBZ-2]; Vb[NBZ+1]=Vb[NBZ-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBZ]=Vc[NBZ-2]; Vc[NBZ+1]=Vc[NBZ-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBZ]=Vd[NBZ-2]; Vd[NBZ+1]=Vd[NBZ-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
      for (z = 1; z < NBZ-1; z++) imageO.P(-Vx[z+1],x,y,z);
    }
    
    //2.2.5) temporary output image is reinjected in temporary input image
    for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
      imageE.P(imageO.G(x,y,z),x,y,z);
    
  }
  //2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
  for (z = 0; z < NBZ-2; z++) for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
    SField->P(imageE.G(x+1,y+1,z+1),x,y,z);
}


//Shrink regions contained in an image/mask to single points which can be considered as the regions center. These points are estimated using iterative erosions.
void ShrinkRegionsToSinglePoints(ScalarField * Mask,float MaskId,ScalarField * CenterPointsImage){
  int x,y,z,i,j,k;
  int OnlySinglePoints,SinglePoint,BoundaryPoint;
  float epsilon;
  int TmpInt;
  
  epsilon=0.00001;
    
    
  cout << endl << "Points to add:";
    
  //1) copy Mask into CenterPointsImage and treat CenterPointsImage  (shape = 1 / background = 0 / image boundaries = 0)
  DeepCopy(Mask,CenterPointsImage);

  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
    if (fabs(CenterPointsImage->G(x,y,z)-MaskId)<epsilon)
      CenterPointsImage->P(1,x,y,z);
    else
      CenterPointsImage->P(0,x,y,z);
    }
    
   for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++) CenterPointsImage->P(0,0,y,z);
   for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++) CenterPointsImage->P(0,Mask->NX-1,y,z);
   for (y = 0; y < Mask->NY; y++) for (x = 0; x < Mask->NX; x++) CenterPointsImage->P(0,x,y,0);
   for (y = 0; y < Mask->NY; y++) for (x = 0; x < Mask->NX; x++) CenterPointsImage->P(0,x,y,Mask->NZ-1);
   for (z = 0; z < Mask->NZ; z++) for (x = 0; x < Mask->NX; x++) CenterPointsImage->P(0,x,0,z);
   for (z = 0; z < Mask->NZ; z++) for (x = 0; x < Mask->NX; x++) CenterPointsImage->P(0,x,Mask->NY-1,z);

    //CenterPointsImage->Write("MaskToShrink.nii");

    
   //2) erode the shape until we have single points only  (algorithm could probably be optimized)
   //2.1) initiate 
   OnlySinglePoints=0;
   
   while (OnlySinglePoints==0){
     OnlySinglePoints=1;
     
       //2.2) detect shape  boundaries
     for (z = 1; z < Mask->NZ-1; z++) for (y = 1; y < Mask->NY-1; y++)  for (x = 1; x < Mask->NX-1; x++) if (CenterPointsImage->G(x,y,z)>0.5){
       BoundaryPoint=0;
       if (CenterPointsImage->G(x+1,y,z)<0.5) BoundaryPoint=1;
       if (CenterPointsImage->G(x-1,y,z)<0.5) BoundaryPoint=1;
       if (CenterPointsImage->G(x,y+1,z)<0.5) BoundaryPoint=1;
       if (CenterPointsImage->G(x,y-1,z)<0.5) BoundaryPoint=1;
       if (CenterPointsImage->G(x,y,z+1)<0.5) BoundaryPoint=1;
       if (CenterPointsImage->G(x,y,z-1)<0.5) BoundaryPoint=1;
       if (BoundaryPoint==1) CenterPointsImage->P(2,x,y,z);
     }
       
      //2.3) detects points which are not isolated and only remove them
     for (z = 1; z < Mask->NZ-1; z++) for (y = 1; y < Mask->NY-1; y++)  for (x = 1; x < Mask->NX-1; x++) if (CenterPointsImage->G(x,y,z)>1.5){
         SinglePoint=1;
         if (CenterPointsImage->G(x+1,y+1,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x+1,y-1,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x+1,y,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x+1,y,z+1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x+1,y,z-1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x-1,y+1,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x-1,y-1,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x-1,y,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x-1,y,z+1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x-1,y,z-1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y+1,z-1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y+1,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y+1,z+1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y-1,z-1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y-1,z)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y-1,z+1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y,z+1)>0.5) SinglePoint=0;
         if (CenterPointsImage->G(x,y,z-1)>0.5) SinglePoint=0;
         
         if (SinglePoint==0) {
             CenterPointsImage->P(0,x,y,z);
             OnlySinglePoints=0;
         }
         else
             CenterPointsImage->P(1,x,y,z);
       }
       
   }
  
  //3) give the number of points to add
  TmpInt=0;
  for (z = 1; z < Mask->NZ-1; z++) for (y = 1; y < Mask->NY-1; y++)  for (x = 1; x < Mask->NX-1; x++) if (CenterPointsImage->G(x,y,z)>0.5)
      TmpInt++;
    cout << TmpInt << endl;
}


//compute the distance map of the edges in a 3D Mask. Greedy algorithm updating the distances in a band.
//-> 'NeaBoun' contains a vector field pointing to the nearest boundaries
//-> 'Distance' contains the distance to the nearest boundary
void Cpt_NearestBoundary(ScalarField * Mask,VectorField * NeaBoun,ScalarField * Distance){
  int x,y,z,direc,i,j,k;
  float epsilon;
  int *ListBx;
  int *ListBy;
  int *ListBz;
  int *ListBx2;
  int *ListBy2;
  int *ListBz2;
  int NbB;
  int NbB2;
  float TmpDist;
  int ix,iy,iz;
  int changes;
  int NbMaxBoundaryPts;
  float mindxdydz;
  float MaxDist;
  float dx,dy,dz;
  float ThreshDistanceNotAnyMoreTreated;
  int MaxJ;
  
  epsilon=0.00001;
  NbMaxBoundaryPts=30000000;
  
  ListBx = new int [NbMaxBoundaryPts];
  ListBy = new int [NbMaxBoundaryPts];
  ListBz = new int [NbMaxBoundaryPts];
  
  ListBx2 = new int [NbMaxBoundaryPts];
  ListBy2 = new int [NbMaxBoundaryPts];
  ListBz2 = new int [NbMaxBoundaryPts];
  
  //initate dx,dy,dz, MaxDist, NeaBoun and Distance
  dx=sqrt(Mask->Image2World[0][0]*Mask->Image2World[0][0]+Mask->Image2World[0][1]*Mask->Image2World[0][1]+Mask->Image2World[0][2]*Mask->Image2World[0][2]);
  dy=sqrt(Mask->Image2World[1][0]*Mask->Image2World[1][0]+Mask->Image2World[1][1]*Mask->Image2World[1][1]+Mask->Image2World[1][2]*Mask->Image2World[1][2]);
  dz=sqrt(Mask->Image2World[2][0]*Mask->Image2World[2][0]+Mask->Image2World[2][1]*Mask->Image2World[2][1]+Mask->Image2World[2][2]*Mask->Image2World[2][2]);
  
  MaxDist=sqrt((Mask->NZ*Mask->NZ*dz*dz)+(Mask->NY*Mask->NY*dy*dy)+(Mask->NX*Mask->NX*dx*dx));
  
  mindxdydz=dx;
  if (dy<mindxdydz) mindxdydz=dy;
  if (dz<mindxdydz) mindxdydz=dz;
  
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    Distance->P(MaxDist*MaxDist,x,y,z);
  
  for (direc = 0; direc < 3; direc++) for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    NeaBoun->P(0,direc,x,y,z);

  
  //1) FIRST BOUNDARIES
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX-1; x++){
    if(fabs(Mask->G(x,y,z)-Mask->G(x+1,y,z))>epsilon){
      Distance->P(1,x,y,z);      
      NeaBoun->Add(0.5*dx,0,x,y,z);   
      Distance->P(1,x+1,y,z);    
      NeaBoun->Add(-0.5*dx,0,x+1,y,z);
    }
  }
    
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY-1; y++)  for (x = 0; x < Mask->NX; x++){
    if(fabs(Mask->G(x,y,z)-Mask->G(x,y+1,z))>epsilon){
      Distance->P(1,x,y,z);       
      NeaBoun->Add(0.5*dy,1,x,y,z);   
      Distance->P(1,x,y+1,z);     
      NeaBoun->Add(-0.5*dy,1,x,y+1,z);
    }
  }
  
  if (Mask->NZ>2){  //2D case
    for (z = 0; z < Mask->NZ-1; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
      if (fabs(Mask->G(x,y,z)-Mask->G(x,y,z+1))>epsilon){
        Distance->P(1,x,y,z);     
        NeaBoun->Add(0.5*dz,2,x,y,z);   
        Distance->P(1,x,y,z+1);   
        NeaBoun->Add(-0.5*dz,2,x,y,z+1);
      }
    }
  }
  
  NbB=0;
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
    if (Distance->G(x,y,z)<2){
      TmpDist=sqrt(NeaBoun->G(0,x,y,z)*NeaBoun->G(0,x,y,z)+NeaBoun->G(1,x,y,z)*NeaBoun->G(1,x,y,z)+NeaBoun->G(2,x,y,z)*NeaBoun->G(2,x,y,z));
      Distance->P(TmpDist,x,y,z);
      ListBx[NbB]=x; 
      ListBy[NbB]=y; 
      ListBz[NbB]=z;
      NbB++;
    }
    else{
      NeaBoun->P(MaxDist*MaxDist,0,x,y,z);
      NeaBoun->P(MaxDist*MaxDist,1,x,y,z);
      NeaBoun->P(MaxDist*MaxDist,2,x,y,z);
    }
  }
  
  //2) DISTANCE PROPAGATION
  changes=1;
  
  ThreshDistanceNotAnyMoreTreated=-3*mindxdydz;
  while (changes==1){
    changes=0;
    NbB2=0;
    
    
    //one iteration of the propagation
    for (i=0; i<NbB; i++) {
      //init
      x=ListBx[i]; y=ListBy[i]; z=ListBz[i];
      
      //update the list with the current point for the next iteration
      if (Distance->G(x,y,z)>ThreshDistanceNotAnyMoreTreated){
        ListBx2[NbB2]=x;
        ListBy2[NbB2]=y;
        ListBz2[NbB2]=z;
        NbB2++;
      }
      
      //update the 6 ngbh of the points in the list
      ix=0;iy=0;iz=0;
      if (Mask->NZ>2) MaxJ=4;  //2D case
      else MaxJ=6;             //3D case

      for (j=0;j<MaxJ;j++){
        if (j==0) {ix=1;iy=0;iz=0;}
        if (j==1) {ix=-1;}
        if (j==2) {ix=0;iy=1;}
        if (j==3) {iy=-1;}
        if (j==4) {iy=0;iz=1;}
        if (j==5) {iz=-1;};
        
        TmpDist=sqrt((NeaBoun->G(0,x,y,z)-ix*dx)*(NeaBoun->G(0,x,y,z)-ix*dx)+(NeaBoun->G(1,x,y,z)-iy*dy)*(NeaBoun->G(1,x,y,z)-iy*dy)+(NeaBoun->G(2,x,y,z)-iz*dz)*(NeaBoun->G(2,x,y,z)-iz*dz));
        
        if ((x+ix>=0)&&(x+ix<NeaBoun->NX)&&(y+iy>=0)&&(y+iy<NeaBoun->NY)&&(z+iz>=0)&&(z+iz<NeaBoun->NZ)) if ((TmpDist<Distance->G(x+ix,y+iy,z+iz))&&(TmpDist<MaxDist)){
          if (Distance->G(x+ix,y+iy,z+iz)>MaxDist*MaxDist-2){
            ListBx2[NbB2]=x+ix;
            ListBy2[NbB2]=y+iy;
            ListBz2[NbB2]=z+iz;
            NbB2++;
            if (NbB2>NbMaxBoundaryPts-10){
              cout << "Too much points at the boundary -> limit extended" << endl;
              delete ListBx;
              delete ListBy;
              delete ListBz;
              ListBx = new int [Mask->NX*Mask->NY*Mask->NZ];
              ListBy = new int [Mask->NX*Mask->NY*Mask->NZ];
              ListBz = new int [Mask->NX*Mask->NY*Mask->NZ];
              
              for (k=0;k<NbB2;k++) ListBx[k]=ListBx2[k];
              for (k=0;k<NbB2;k++) ListBy[k]=ListBy2[k];
              for (k=0;k<NbB2;k++) ListBz[k]=ListBz2[k];
              
              delete ListBx2;
              delete ListBy2;
              delete ListBz2;
              ListBx2 = new int [Mask->NX*Mask->NY*Mask->NZ];
              ListBy2 = new int [Mask->NX*Mask->NY*Mask->NZ];
              ListBz2 = new int [Mask->NX*Mask->NY*Mask->NZ];
              
              for (k=0;k<NbB2;k++) ListBx2[k]=ListBx[k];
              for (k=0;k<NbB2;k++) ListBy2[k]=ListBy[k];
              for (k=0;k<NbB2;k++) ListBz2[k]=ListBz[k];
              
              NbMaxBoundaryPts=Mask->NX*Mask->NY*Mask->NZ;
            } 
            changes=1;
          }
          NeaBoun->P(NeaBoun->G(0,x,y,z)-ix*dx,0,x+ix,y+iy,z+iz);
          NeaBoun->P(NeaBoun->G(1,x,y,z)-iy*dy,1,x+ix,y+iy,z+iz);
          NeaBoun->P(NeaBoun->G(2,x,y,z)-iz*dz,2,x+ix,y+iy,z+iz);
          Distance->P(TmpDist,x+ix,y+iy,z+iz);
        }
      }
      
    }
    
    //prepare the next iteration of the propagation
    for(i=0;i<NbB2;i++){
      ListBx[i]=ListBx2[i];
      ListBy[i]=ListBy2[i];
      ListBz[i]=ListBz2[i];
    }
    NbB=NbB2;
    
    ThreshDistanceNotAnyMoreTreated+=mindxdydz;
  }
  
  
  //put 1 in the points of Distance in which the computations have been done and 0 elsewhere
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++)if (Distance->G(x,y,z)>MaxDist*MaxDist-2){
      Distance->P(0,x,y,z);
      NeaBoun->P(0,0,x,y,z);
      NeaBoun->P(0,1,x,y,z);
      NeaBoun->P(0,2,x,y,z);
  }

  delete ListBx;
  delete ListBy;
  delete ListBz;
  delete ListBx2;
  delete ListBy2;
  delete ListBz2;
  
}




//Compute the map of the nearsest boundaries in the 3D mask 'Mask' (propagate the nearest boundaries identification until at a maximum distance of MaxDist)
//Outputs:
//  * NeaBoun: vector field pointing the nearest boundary
//  * TempSF: In the end the points for which the nearest boundary is identified have a value equal to 1
//Option:
//If NbRegions>0 then masks with partial volume effects are treated. In that case NbRegions is the nb of regions 
//considered in the mask and IdRegions are the corresponding regions 
//If 'LowCostPropagation==1' a wise algorithm propagates the distance map at a low algorithmic cost with a small approximation (otherwise no approximation but much slower algorithm)
void Cpt_NearestBoundaryOLD(ScalarField * Mask,VectorField * NeaBoun,ScalarField * TempSF,float MaxDist,int NbRegions, float * IdRegions,float dx,float dy, float dz,int LowCostPropagation){
  int x,y,z,direc,i,j,k;
  float flx,fly,flz;
  float flxp,flyp,flzp;
  float flxm,flym,flzm;
  float flxt,flyt,flzt;
  float epsilon,epsilon2;
  int *ListBx;
  int *ListBy;
  int *ListBz;
  int *ListBx2;
  int *ListBy2;
  int *ListBz2;
  int NbB;
  int NbB2;
  float TmpDist,TmpMaxDist;
  float direcX,direcY,direcZ;
  float IdReg1,IdReg2,IdRegRef;
  float bestDifference;
  int ix,iy,iz;
  int changes;
  int NbMaxBoundaryPts;
  int NbPastIterations;
  float maxdxdydz;
  

  epsilon=0.00001;
  epsilon2=0.01;
  NbMaxBoundaryPts=30000000;
  
  ListBx = new int [NbMaxBoundaryPts];
  ListBy = new int [NbMaxBoundaryPts];
  ListBz = new int [NbMaxBoundaryPts];
  
  ListBx2 = new int [NbMaxBoundaryPts];
  ListBy2 = new int [NbMaxBoundaryPts];
  ListBz2 = new int [NbMaxBoundaryPts];
  
  
  
  
  //initate NeaBoun and TempSF
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    TempSF->P(MaxDist*MaxDist,x,y,z);
  
  for (direc = 0; direc < 3; direc++) for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    NeaBoun->P(MaxDist*MaxDist,direc,x,y,z);
  
  
  //1) FIRST BOUNDARIES
  //normalize the dx, dy, dz  (this will be undone in the 2nd step)
  maxdxdydz=dx;
  if (dy>maxdxdydz) maxdxdydz=dy;
  if (dz>maxdxdydz) maxdxdydz=dz;
  
  dx=dx/maxdxdydz;
  dy=dy/maxdxdydz;
  dz=dz/maxdxdydz;
  
  
  NbB=0;
  if (NbRegions<0){//no partial volume effect in the mask
    for (z = 1; z < Mask->NZ-2; z++) for (y = 1; y < Mask->NY-2; y++)  for (x = 1; x < Mask->NX-2; x++){
      if(fabs(Mask->G(x,y,z)-Mask->G(x+1,y,z))>epsilon){
        TempSF->P(1,x,y,z);      NeaBoun->Add(0.5*dx,0,x,y,z);   
        TempSF->P(1,x+1,y,z);    NeaBoun->Add(-0.5*dx,0,x+1,y,z);
      }
      
      if(fabs(Mask->G(x,y,z)-Mask->G(x,y+1,z))>epsilon){
        TempSF->P(1,x,y,z);       NeaBoun->Add(0.5*dy,1,x,y,z);  
        TempSF->P(1,x,y+1,z);     NeaBoun->Add(-0.5*dy,1,x,y+1,z);  
      }
      
      if(fabs(Mask->G(x,y,z)-Mask->G(x,y,z+1))>epsilon){
        TempSF->P(1,x,y,z);     NeaBoun->Add(0.5*dz,2,x,y,z);
        TempSF->P(1,x,y,z+1);   NeaBoun->Add(-0.5*dz,2,x,y,z+1);
      }
    }
    
    for (z = 1; z < Mask->NZ-2; z++) for (y = 1; y < Mask->NY-2; y++)  for (x = 1; x < Mask->NX-2; x++) if (TempSF->G(x,y,z)>2){
      TmpDist=sqrt(NeaBoun->G(0,x,y,z)*NeaBoun->G(0,x,y,z)+NeaBoun->G(1,x,y,z)*NeaBoun->G(1,x,y,z)+NeaBoun->G(2,x,y,z)*NeaBoun->G(2,x,y,z));
      if (TmpDist>epsilon){
        NeaBoun->P(0.5*NeaBoun->G(0,x,y,z)/TmpDist,0,x,y,z); 
        NeaBoun->P(0.5*NeaBoun->G(1,x,y,z)/TmpDist,1,x,y,z); 
        NeaBoun->P(0.5*NeaBoun->G(2,x,y,z)/TmpDist,2,x,y,z); 
        
        TempSF->P(TmpDist,x,y,z);
        ListBx[NbB]=x; 
        ListBy[NbB]=y; 
        ListBz[NbB]=z;
        NbB++;
      }
    }
  }
  else{//partial volume effects in the mask
    for (z = 4; z < Mask->NZ-4; z++) for (y = 4; y < Mask->NY-4; y++)  for (x = 4; x < Mask->NX-4; x++){
      flx=static_cast<float>(x);
      fly=static_cast<float>(y);
      flz=static_cast<float>(z);
      flxp=static_cast<float>(x+0.1);
      flyp=static_cast<float>(y+0.1);
      flzp=static_cast<float>(z+0.1);
      flxm=static_cast<float>(x-0.1);
      flym=static_cast<float>(y-0.1);
      flzm=static_cast<float>(z-0.1);
      
      
      if ((fabs(Mask->G(x+1,y,z)-Mask->G(x-1,y,z))>epsilon)||(fabs(Mask->G(x,y+1,z)-Mask->G(x,y-1,z))>epsilon)||(fabs(Mask->G(x,y,z+1)-Mask->G(x,y,z-1))>epsilon)){
        //identify the two regions' ID
        IdReg1=-1;
        IdReg2=-1;
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x+4,y,z)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x+4,y,z)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x-4,y,z)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x-4,y,z)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x,y+4,z)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x,y+4,z)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x,y-4,z)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x,y-4,z)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x,y,z+4)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x,y,z+4)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x,y,z-4)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x,y,z-4)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        
        if (IdReg2>0){ //otherwise the identifiers around a boundary are not identified
          //identify the grey level of the mask at the boundary
          IdRegRef=(IdReg1+IdReg2)/2;
          
          //identify the direction
          direcX=Mask->G(flxp,fly,flz)-Mask->G(flxm,fly,flz);
          direcY=Mask->G(flx,flyp,flz)-Mask->G(flx,flym,flz);
          direcZ=Mask->G(flx,fly,flzp)-Mask->G(flx,fly,flzm);
          
          
          if ((fabs(direcX)>epsilon2)||(fabs(direcY)>epsilon2)||(fabs(direcZ)>epsilon2)){ // otherwise one direction nas not been identified
            
            TmpDist=sqrt(direcX*direcX+direcY*direcY+direcZ*direcZ);
            direcX/=TmpDist;
            direcY/=TmpDist;
            direcZ/=TmpDist;
            if (Mask->G(x,y,z)>IdRegRef){
              direcX*=-1;
              direcY*=-1;
              direcZ*=-1;
            }
            
            //find the boundary
            bestDifference=fabs(Mask->G(flx,fly,flz)-IdRegRef);
            NbB2=0;
            flxt=static_cast<float>(flx+0.1*direcX);
            flyt=static_cast<float>(fly+0.1*direcY);
            flzt=static_cast<float>(flz+0.1*direcZ);
            while ((bestDifference>=fabs(Mask->G(flxt,flyt,flzt)-IdRegRef)-epsilon)&&(NbB2<30)) {
              bestDifference=fabs(Mask->G(flxt,flyt,flzt)-IdRegRef);
              NbB2++;
              flxt+=direcX;
              flyt+=direcY;
              flzt+=direcZ;
            }
            flxt-=direcX;
            flyt-=direcY;
            flzt-=direcZ;
            
            NeaBoun->P(flxt-flx,0,x,y,z); 
            NeaBoun->P(flyt-fly,1,x,y,z); 
            NeaBoun->P(flzt-flz,2,x,y,z);
            
            if ((fabs(flxt-flx)<0.05)||(fabs(flyt-fly)<0.05)||(fabs(flzt-flz)<0.05)){
              NeaBoun->P(direcX*0.05,0,x,y,z); 
              NeaBoun->P(direcY*0.05,1,x,y,z); 
              NeaBoun->P(direcZ*0.05,2,x,y,z);
            }
            
            TmpDist=((flxt-flx)*(flxt-flx)*dx*dx+(flyt-fly)*(flyt-fly)*dy*dy+(flzt-flzt)*(flzt-flz)*dz*dz);
            TempSF->P(sqrt(TmpDist),x,y,z);
            ListBx[NbB]=x; 
            ListBy[NbB]=y; 
            ListBz[NbB]=z;
            NbB++;
          }
        }
      }
    }
  }
  
  //2) DISTANCE PROPAGATION
  
  //give again the real resolution to 'dx, dy, dz' when propagating the distance
  dx=dx*maxdxdydz;
  dy=dy*maxdxdydz;
  dz=dz*maxdxdydz;
  
  //2.A) Low cost propagation
  if (LowCostPropagation==1){
    changes=1;
    NbPastIterations=0;
    while (changes==1){
      changes=0;
      NbB2=0;
      
      
      //one iteration of the propagation
      for (i=0; i<NbB; i++) {
        //init
        x=ListBx[i];
        y=ListBy[i];
        z=ListBz[i];
        
        //update the list with the current point for the next iteration
        if (TempSF->G(x,y,z)>static_cast<float>(NbPastIterations-3)){
          ListBx2[NbB2]=x;
          ListBy2[NbB2]=y;
          ListBz2[NbB2]=z;
          NbB2++;
        }
        
        //update the 6 ngbh of the points in the list
        ix=0;iy=0;iz=0;
        for (j=0;j<6;j++){
          if (j==0) {ix=1;iy=0;iz=0;}
          if (j==1) {ix=-1;}
          if (j==2) {ix=0;iy=1;}
          if (j==3) {iy=-1;}
          if (j==4) {iy=0;iz=1;}
          if (j==5) {iz=-1;};
          
          TmpDist=sqrt((NeaBoun->G(0,x,y,z)-ix*dx)*(NeaBoun->G(0,x,y,z)-ix*dx)+(NeaBoun->G(1,x,y,z)-iy*dy)*(NeaBoun->G(1,x,y,z)-iy*dy)+(NeaBoun->G(2,x,y,z)-iz*dz)*(NeaBoun->G(2,x,y,z)-iz*dz));
          
          if ((TmpDist<TempSF->G(x+ix,y+iy,z+iz))&&(TmpDist<MaxDist)) if ((x>2)&&(x<NeaBoun->NX-3)&&(y>2)&&(y<NeaBoun->NY-3)&&(z>2)&&(z<NeaBoun->NZ-3)){
            if (TempSF->G(x+ix,y+iy,z+iz)>MaxDist*MaxDist-2){
              ListBx2[NbB2]=x+ix;
              ListBy2[NbB2]=y+iy;
              ListBz2[NbB2]=z+iz;
              NbB2++;
              if (NbB2>NbMaxBoundaryPts-10){
                cout << "Too much points at the boundary -> limit extended" << endl;
                delete ListBx;
                delete ListBy;
                delete ListBz;
                ListBx = new int [Mask->NX*Mask->NY*Mask->NZ];
                ListBy = new int [Mask->NX*Mask->NY*Mask->NZ];
                ListBz = new int [Mask->NX*Mask->NY*Mask->NZ];
                
                for (k=0;k<NbB2;k++) ListBx[k]=ListBx2[k];
                for (k=0;k<NbB2;k++) ListBy[k]=ListBy2[k];
                for (k=0;k<NbB2;k++) ListBz[k]=ListBz2[k];
                
                delete ListBx2;
                delete ListBy2;
                delete ListBz2;
                ListBx2 = new int [Mask->NX*Mask->NY*Mask->NZ];
                ListBy2 = new int [Mask->NX*Mask->NY*Mask->NZ];
                ListBz2 = new int [Mask->NX*Mask->NY*Mask->NZ];
                
                for (k=0;k<NbB2;k++) ListBx2[k]=ListBx[k];
                for (k=0;k<NbB2;k++) ListBy2[k]=ListBy[k];
                for (k=0;k<NbB2;k++) ListBz2[k]=ListBz[k];
                
                NbMaxBoundaryPts=Mask->NX*Mask->NY*Mask->NZ;
              } 
              changes=1;
            }
            NeaBoun->P(NeaBoun->G(0,x,y,z)-ix*dx,0,x+ix,y+iy,z+iz);
            NeaBoun->P(NeaBoun->G(1,x,y,z)-iy*dy,1,x+ix,y+iy,z+iz);
            NeaBoun->P(NeaBoun->G(2,x,y,z)-iz*dz,2,x+ix,y+iy,z+iz);
            TempSF->P(TmpDist,x+ix,y+iy,z+iz);
          }
          
        }
        
      }
      
      //prepare the next iteration of the propagation
      for(i=0;i<NbB2;i++){
        ListBx[i]=ListBx2[i];
        ListBy[i]=ListBy2[i];
        ListBz[i]=ListBz2[i];
      }
      NbB=NbB2;
      
      NbPastIterations++;
    }
    
    
    //put 1 in the points of TempSF in which the computations have been done and 0 elsewhere
    for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
      if (TempSF->G(x,y,z)<MaxDist*MaxDist-2) TempSF->P(1,x,y,z);
      else{
        TempSF->P(0,x,y,z);
        NeaBoun->P(0,0,x,y,z);
        NeaBoun->P(0,1,x,y,z);
        NeaBoun->P(0,2,x,y,z);
      }
    }
    
    delete ListBx;
    delete ListBy;
    delete ListBz;
    delete ListBx2;
    delete ListBy2;
    delete ListBz2;
  }
  if (LowCostPropagation!=1){ //2.B) Full propagation
    
    //2.B.1) intiate the computations
    TmpMaxDist=sqrt(FLT_MAX/100);
    
    for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
      if (TempSF->G(x,y,z)==0){
        TempSF->P(TmpMaxDist,x,y,z);
        NeaBoun->P(TmpMaxDist,0,x,y,z);
        NeaBoun->P(TmpMaxDist,1,x,y,z);
        NeaBoun->P(TmpMaxDist,2,x,y,z);
        
      }
      else{
        flx=NeaBoun->G(0,x,y,z);
        fly=NeaBoun->G(1,x,y,z);
        flz=NeaBoun->G(2,x,y,z);
        
        TmpDist=sqrt(flx*flx*dx*dx+fly*fly*dy*dy+flz*flz*dz*dz);
        
        TempSF->P(TmpDist,x,y,z);
      }
  }
    
    
    //2.B.2) propagate
    changes=1;
    while (changes>0){
      //cout << changes << endl;
      changes=0;
      for (z = 1; z < Mask->NZ-1; z++) for (y = 1; y < Mask->NY-1; y++)  for (x = 1; x < Mask->NX-1; x++) for (j=0;j<6;j++){
        if (j==0) {ix=1;iy=0;iz=0;}
        if (j==1) {ix=-1;}
        if (j==2) {ix=0;iy=1;}
        if (j==3) {iy=-1;}
        if (j==4) {iy=0;iz=1;}
        if (j==5) {iz=-1;};
        
        
        TmpDist=sqrt((NeaBoun->G(0,x,y,z)-ix*dx)*(NeaBoun->G(0,x,y,z)-ix*dx)+(NeaBoun->G(1,x,y,z)-iy*dy)*(NeaBoun->G(1,x,y,z)-iy*dy)+(NeaBoun->G(2,x,y,z)-iz*dz)*(NeaBoun->G(2,x,y,z)-iz*dz));
        
        if (TmpDist<TempSF->G(x+ix,y+iy,z+iz)){
          NeaBoun->P(NeaBoun->G(0,x,y,z)-ix*dx,0,x+ix,y+iy,z+iz);
          NeaBoun->P(NeaBoun->G(1,x,y,z)-iy*dy,1,x+ix,y+iy,z+iz);
          NeaBoun->P(NeaBoun->G(2,x,y,z)-iz*dz,2,x+ix,y+iy,z+iz);
          TempSF->P(TmpDist,x+ix,y+iy,z+iz);
          changes++;
        }
      }
    }
    
    //2.B.3) post-treatment (put 1 in the points of TempSF in which the computations have been done and 0 elsewhere)
    
    for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
      if ((z==0)||(z==Mask->NZ-1)||(y==0)||(y==Mask->NY-1)||(x==0)||(x==Mask->NX-1)){
        TempSF->P(0,x,y,z);
        NeaBoun->P(0,0,x,y,z);
        NeaBoun->P(0,1,x,y,z);
        NeaBoun->P(0,2,x,y,z);
        
      }
      else{
        TempSF->P(1,x,y,z);
      }
    }
    
  }  
  
}

//Stochastically propagate the non-null values of ImageIn in a ROI.
// -> The ROI is defined where we have the non-null values of Mask.
// -> Intensity propagation takes into account the distance between neighbors in world coordinates
void PropagateNonNullValues_Sto(ScalarField * ImageIn,ScalarField * Mask){
  ScalarField ImageTmp;
  int changes;
  int x,y,z;
  int x2,y2,z2;
  float dx,dy,dz,dist;
  float Ngbh_Intensity[26];
  float Ngbh_SumDists[26];
  float Ngbh_SumAllDists;
  float ChosenCumulatedDist;
  int ChosenNgbh,Ngbh_Nb;
  
  srand(time(0));


  ImageTmp.CreateVoidField(ImageIn->NX,ImageIn->NY,ImageIn->NZ);
  
  if (ImageIn->NZ<=1) cout << "Grey level propagation was only coded for 3D images " << endl;
  
  changes=1;
  
  while (changes>0){
    changes=0;
    ImageTmp.PutToAllVoxels(0);
    
    for (z = 1; z < ImageIn->NZ-1; z++) for (y = 1; y < ImageIn->NY-1; y++)  for (x = 1; x < ImageIn->NX-1; x++) if (fabs(Mask->G(x,y,z)-1)<0.00001) if (fabs(ImageIn->G(x,y,z))<0.00001) {
      Ngbh_SumAllDists=0;
      Ngbh_Nb=0;
      for (z2 = z-1; z2 <= z+1; z2++) for (y2 = y-1; y2 <= y+1; y2++)  for (x2 = x-1; x2 <= x+1; x2++) if (fabs(Mask->G(x2,y2,z2)-1)<0.00001) if (fabs(ImageIn->G(x2,y2,z2))>0.00001){
        dx=static_cast<float>(x2-x)*ImageIn->Image2World[0][0]+static_cast<float>(y2-y)*ImageIn->Image2World[0][1]+static_cast<float>(z2-z)*ImageIn->Image2World[0][2];
        dy=static_cast<float>(x2-x)*ImageIn->Image2World[1][0]+static_cast<float>(y2-y)*ImageIn->Image2World[1][1]+static_cast<float>(z2-z)*ImageIn->Image2World[1][2];
        dz=static_cast<float>(x2-x)*ImageIn->Image2World[2][0]+static_cast<float>(y2-y)*ImageIn->Image2World[2][1]+static_cast<float>(z2-z)*ImageIn->Image2World[2][2];

        dist=sqrt((dx*dx)+(dy*dy)+(dz*dz));
        Ngbh_SumAllDists+=1/dist;
        Ngbh_Intensity[Ngbh_Nb]=ImageIn->G(x2,y2,z2);
        Ngbh_SumDists[Ngbh_Nb]=Ngbh_SumAllDists;
        
        Ngbh_Nb++;
      }
      
      if (Ngbh_Nb>0){
        ChosenCumulatedDist=static_cast<float>(Ngbh_SumAllDists-0.000000001)*static_cast<float>(rand())/static_cast<float>(RAND_MAX);
        
        ChosenNgbh=0;
        while (Ngbh_SumDists[ChosenNgbh]<ChosenCumulatedDist) ChosenNgbh++;
        
        ImageTmp.P(Ngbh_Intensity[ChosenNgbh],x,y,z);
        changes++;
        }
    }
    
    if (changes>0){
      for (z = 1; z < ImageIn->NZ-1; z++) for (y = 1; y < ImageIn->NY-1; y++)  for (x = 1; x < ImageIn->NX-1; x++) if (fabs(ImageTmp.G(x,y,z))>0.00001)
        ImageIn->P(ImageTmp.G(x,y,z),x,y,z);
      }
  
    cout << changes << " changes" << endl;
  }
  
}


//Propagate the non-null values of ImageIn in a ROI. The value given to each null value of ImageIn in the ROI is the nearest non-null value of ImageIn in world coordinates
// -> The ROI is defined where we have the non-null values of Mask.
// -> Intensity propagation takes into account the distance between neighbors in world coordinates
void PropagateNonNullValues_NN(ScalarField * ImageIn,ScalarField * Mask){
  ScalarField ImageTmp;
  int changes;
  int x,y,z;
  int x2,y2,z2;
  float dx,dy,dz,dist,minDist,bestIntensity;
  float Ngbh_Intensity[26];
  float Ngbh_SumDists[26];
  float Ngbh_SumAllDists;
  float ChosenCumulatedDist;
  int ChosenNgbh,Ngbh_Nb;
  
  
  ImageTmp.CreateVoidField(ImageIn->NX,ImageIn->NY,ImageIn->NZ);
  
  ImageTmp.PutToAllVoxels(0);
    
  for (z = 0; z < ImageIn->NZ; z++) for (y = 0; y < ImageIn->NY; y++)  for (x = 0; x < ImageIn->NX; x++) if (fabs(Mask->G(x,y,z)-1)<0.00001) if (fabs(ImageIn->G(x,y,z))<0.00001) {
    dist=-1;
    for (z2 = 0; z2 < ImageIn->NZ; z2++) for (y2 = 0; y2 < ImageIn->NY; y2++)  for (x2 = 0; x2 < ImageIn->NX; x2++) if (fabs(Mask->G(x2,y2,z2)-1)<0.00001) if (fabs(ImageIn->G(x2,y2,z2))>0.00001){
      
      dx=static_cast<float>(x2-x)*ImageIn->Image2World[0][0]+static_cast<float>(y2-y)*ImageIn->Image2World[0][1]+static_cast<float>(z2-z)*ImageIn->Image2World[0][2];
      dy=static_cast<float>(x2-x)*ImageIn->Image2World[1][0]+static_cast<float>(y2-y)*ImageIn->Image2World[1][1]+static_cast<float>(z2-z)*ImageIn->Image2World[1][2];
      dz=static_cast<float>(x2-x)*ImageIn->Image2World[2][0]+static_cast<float>(y2-y)*ImageIn->Image2World[2][1]+static_cast<float>(z2-z)*ImageIn->Image2World[2][2];
      
      dist=sqrt((dx*dx)+(dy*dy)+(dz*dz));
      
      if (dist<0) {minDist=dist;  bestIntensity=ImageIn->G(x2,y2,z2);}
      if (minDist>dist) {minDist=dist;  bestIntensity=ImageIn->G(x2,y2,z2);}
      
    }
    
    ImageTmp.P(bestIntensity,x,y,z);
  }
    
  for (z = 0; z < ImageIn->NZ; z++) for (y = 0; y < ImageIn->NY; y++)  for (x = 0; x < ImageIn->NX; x++) if (fabs(ImageTmp.G(x,y,z))>0.00001)
        ImageIn->P(ImageTmp.G(x,y,z),x,y,z);
}


//Compute the map of the nearsest boundaries in the 3D mask 'Mask' (propagate the nearest boundaries identification until at a maximum distance of MaxDist)
//Outputs:
//  * NeaBoun: vector field pointing the nearest boundary
//  * TempSF: In the end the points for which the nearest boundary is identified have a value equal to 1
//Option:
//If NbRegions>0 then masks with partial volume effects are treated. In that case NbRegions is the nb of regions 
//considered in the mask and IdRegions are the corresponding regions 
void Cpt_NearestBoundary2D(ScalarField * Mask,VectorField * NeaBoun,ScalarField * TempSF,int NbRegions, float * IdRegions,float dx,float dy){
  int x,y,direc,i,j;
  float flx,fly;
  float flxp,flyp;
  float flxm,flym;
  float flxt,flyt;
  float epsilon,epsilon2;
  int *ListBx;
  int *ListBy;
  int *ListBx2;
  int *ListBy2;
  int NbB;
  int NbB2;
  float TmpDist,TmpMaxDist;
  float direcX,direcY;
  float IdReg1,IdReg2,IdRegRef;
  float bestDifference;
  int ix,iy;
  int changes;
  int NbMaxBoundaryPts;
  float maxdxdydz;
  float MaxDist;
  float zero;
  
  zero=0;
  MaxDist=static_cast<float>(Mask->NY);
  if (static_cast<float>(Mask->NX)>MaxDist)  MaxDist=static_cast<float>(Mask->NX);
  
  epsilon=0.00001;
  epsilon2=0.01;
  NbMaxBoundaryPts=30000000;
  
  ListBx = new int [NbMaxBoundaryPts];
  ListBy = new int [NbMaxBoundaryPts];
  
  ListBx2 = new int [NbMaxBoundaryPts];
  ListBy2 = new int [NbMaxBoundaryPts];
  
  
  //normalize the dx, dy
  maxdxdydz=dx;
  if (dy>maxdxdydz) maxdxdydz=dy;
  
  dx=dx/maxdxdydz;
  dy=dy/maxdxdydz;
  
  
  //initate NeaBoun and TempSF
  for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    TempSF->P(MaxDist*MaxDist,x,y,0);
  
  for (direc = 0; direc < 3; direc++)  for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    NeaBoun->P(MaxDist*MaxDist,direc,x,y,0);
  
  
  //1) FIRST BOUNDARIES
  NbB=0;
  if (NbRegions<0){//no partial volume effect in the mask
    for (y = 1; y < Mask->NY-2; y++)  for (x = 1; x < Mask->NX-2; x++){
      if(fabs(Mask->G(x,y,0)-Mask->G(x+1,y,0))>epsilon){
        TempSF->P(1,x,y,0);      NeaBoun->Add(0.5*dx,0,x,y,0);   
        TempSF->P(1,x+1,y,0);    NeaBoun->Add(-0.5*dx,0,x+1,y,0);
      }
      
      if(fabs(Mask->G(x,y,0)-Mask->G(x,y+1,0))>epsilon){
        TempSF->P(1,x,y,0);       NeaBoun->Add(0.5*dy,1,x,y,0);  
        TempSF->P(1,x,y+1,0);     NeaBoun->Add(-0.5*dy,1,x,y+1,0);  
      }
      
    }
    
    for (y = 1; y < Mask->NY-2; y++)  for (x = 1; x < Mask->NX-2; x++) if (TempSF->G(x,y,0)>2){
      TmpDist=sqrt(NeaBoun->G(0,x,y,0)*NeaBoun->G(0,x,y,0)+NeaBoun->G(1,x,y,0)*NeaBoun->G(1,x,y,0));
      if (TmpDist>epsilon){
        NeaBoun->P(0.5*NeaBoun->G(0,x,y,0)/TmpDist,0,x,y,0); 
        NeaBoun->P(0.5*NeaBoun->G(1,x,y,0)/TmpDist,1,x,y,0); 
        
        TempSF->P(TmpDist,x,y,0);
        ListBx[NbB]=x; 
        ListBy[NbB]=y; 
        NbB++;
      }
    }
  }
  else{//partial volume effects in the mask
    for (y = 4; y < Mask->NY-4; y++)  for (x = 4; x < Mask->NX-4; x++){
      flx=static_cast<float>(x);
      fly=static_cast<float>(y);
      flxp=static_cast<float>(x+0.1);
      flyp=static_cast<float>(y+0.1);
      flxm=static_cast<float>(x-0.1);
      flym=static_cast<float>(y-0.1);
      
      
      if ((fabs(Mask->G(x+1,y,0)-Mask->G(x-1,y,0))>epsilon)||(fabs(Mask->G(x,y+1,0)-Mask->G(x,y-1,0))>epsilon)){
        //identify the two regions' ID
        IdReg1=-1;
        IdReg2=-1;
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x+4,y,0)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x+4,y,0)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x-4,y,0)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x-4,y,0)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x,y+4,0)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x,y+4,0)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        for (i=0;i<NbRegions;i++) if (fabs(Mask->G(x,y-4,0)-IdRegions[i])<epsilon2){
          if(IdReg1<0) IdReg1=IdRegions[i];
          else if ((IdReg2<0)&&(fabs(Mask->G(x,y-4,0)-IdReg1)>epsilon2)) IdReg2=IdRegions[i];
        }
        
        if (IdReg2>0){ //otherwise the identifiers around a boundary are not identified
          //identify the grey level of the mask at the boundary
          IdRegRef=(IdReg1+IdReg2)/2;
          
          //identify the direction
          direcX=Mask->G(flxp,fly,zero)-Mask->G(flxm,fly,zero);
          direcY=Mask->G(flx,flyp,zero)-Mask->G(flx,flym,zero);
          
          
          if ((fabs(direcX)>epsilon2)||(fabs(direcY)>epsilon2)){ // otherwise one direction nas not been identified
            
            TmpDist=sqrt(direcX*direcX+direcY*direcY);
            direcX/=TmpDist;
            direcY/=TmpDist;
            if (Mask->G(x,y,0)>IdRegRef){
              direcX*=-1;
              direcY*=-1;
            }
            
            //find the boundary
            bestDifference=fabs(Mask->G(flx,fly,zero)-IdRegRef);
            NbB2=0;
            flxt=static_cast<float>(flx+0.1*direcX);
            flyt=static_cast<float>(fly+0.1*direcY);
            while ((bestDifference>=fabs(Mask->G(flxt,flyt,zero)-IdRegRef)-epsilon)&&(NbB2<30)) {
              bestDifference=fabs(Mask->G(flxt,flyt,zero)-IdRegRef);
              NbB2++;
              flxt+=direcX;
              flyt+=direcY;
            }
            flxt-=direcX;
            flyt-=direcY;
            
            NeaBoun->P(flxt-flx,0,x,y,0); 
            NeaBoun->P(flyt-fly,1,x,y,0); 
            
            if ((fabs(flxt-flx)<0.05)||(fabs(flyt-fly)<0.05)){
              NeaBoun->P(direcX*0.05,0,x,y,0); 
              NeaBoun->P(direcY*0.05,1,x,y,0); 
            }
            
            TmpDist=((flxt-flx)*(flxt-flx)*dx*dx+(flyt-fly)*(flyt-fly)*dy*dy);
            TempSF->P(sqrt(TmpDist),x,y,0);
            ListBx[NbB]=x; 
            ListBy[NbB]=y; 
            NbB++;
          }
        }
      }
    }
  }
  
  //2) DISTANCE PROPAGATION
  
  //give again the real value to dx, dy
  dx=dx*maxdxdydz;
  dy=dy*maxdxdydz;
  
  
  //2.B.1) intiate the computations
  TmpMaxDist=sqrt(FLT_MAX/100);
  
  for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
    if (TempSF->G(x,y,0)==0){
      TempSF->P(TmpMaxDist,x,y,0);
      NeaBoun->P(TmpMaxDist,0,x,y,0);
      NeaBoun->P(TmpMaxDist,1,x,y,0);
      
    }
    else{
      flx=NeaBoun->G(0,x,y,0);
      fly=NeaBoun->G(1,x,y,0);
      
      TmpDist=(flx*flx*dx*dx+fly*fly*dy*dy);
      
      TempSF->P(TmpDist,x,y,0);
    }
  }
  
  
  //2.B.2) propagate
  changes=1;
  while (changes>0){
    //cout << changes << endl;
    changes=0;
    for (y = 1; y < Mask->NY-1; y++)  for (x = 1; x < Mask->NX-1; x++) for (j=0;j<4;j++){
      if (j==0) {ix=1;iy=0;}
      if (j==1) {ix=-1;}
      if (j==2) {ix=0;iy=1;}
      if (j==3) {iy=-1;}
      
      
      TmpDist=sqrt((NeaBoun->G(0,x,y,0)-ix*dx)*(NeaBoun->G(0,x,y,0)-ix*dx)+(NeaBoun->G(1,x,y,0)-iy*dy)*(NeaBoun->G(1,x,y,0)-iy*dy));
      
      if (TmpDist<TempSF->G(x+ix,y+iy,0)){
        NeaBoun->P(NeaBoun->G(0,x,y,0)-ix*dx,0,x+ix,y+iy,0);
        NeaBoun->P(NeaBoun->G(1,x,y,0)-iy*dy,1,x+ix,y+iy,0);
        TempSF->P(TmpDist,x+ix,y+iy,0);
        changes++;
      }
    }
  }
  
  //2.B.3) post-treatment (put 1 in the points of TempSF in which the computations have been done and 0 elsewhere)
  
  for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
    if ((y==0)||(y==Mask->NY-1)||(x==0)||(x==Mask->NX-1)){
      TempSF->P(0,x,y,0);
      NeaBoun->P(0,0,x,y,0);
      NeaBoun->P(0,1,x,y,0);
      NeaBoun->P(0,2,x,y,0);
      
    }
    else{
      TempSF->P(1,x,y,0);
    }
  }
  
                   
}




void Cpt_DistMap(ScalarField * SegImage,float IDregion,ScalarField * DistMap){
  int x,y,z,direc,i,j,k;
  float flx,fly,flz;
  float flxp,flyp,flzp;
  float flxm,flym,flzm;
  float flxt,flyt,flzt;
  float epsilon;
  float TmpDist,TmpMaxDist;
  float direcX,direcY,direcZ;
  float IdReg1,IdReg2,IdRegRef;
  float bestDifference;
  int ix,iy,iz;
  int changes;
  int NbPastIterations;
  float maxdxdydz;
  float MaxDist;
  VectorField NeaBoun;
  float dx,dy,dz;
  
  
  //0) INIT
  NeaBoun.CreateVoidField(SegImage->NX,SegImage->NY,SegImage->NZ);
  
  epsilon=0.01;
  
  //set the voxel sizes
  dx=sqrt(SegImage->Image2World[0][0]*SegImage->Image2World[0][0]+SegImage->Image2World[0][1]*SegImage->Image2World[0][1]+SegImage->Image2World[0][2]*SegImage->Image2World[0][2]);
  dy=sqrt(SegImage->Image2World[1][0]*SegImage->Image2World[1][0]+SegImage->Image2World[1][1]*SegImage->Image2World[1][1]+SegImage->Image2World[1][2]*SegImage->Image2World[1][2]);
  dz=sqrt(SegImage->Image2World[2][0]*SegImage->Image2World[2][0]+SegImage->Image2World[2][1]*SegImage->Image2World[2][1]+SegImage->Image2World[2][2]*SegImage->Image2World[2][2]);
  
  cout << "(dx,dy,dz) = (" << dx << "," << dy << "," << dz << ")" << endl;
  
  MaxDist=sqrt(SegImage->NX*SegImage->NX*dx*dx+SegImage->NY*SegImage->NY*dy*dy+SegImage->NZ*SegImage->NZ*dz*dz)+dx;
  
  for (z = 0; z < SegImage->NZ; z++) for (y = 0; y < SegImage->NY; y++)  for (x = 0; x < SegImage->NX; x++) DistMap->P(0,x,y,z);

  
  //1) FIRST BOUNDARIES
  
  for (z = 1; z < SegImage->NZ-1; z++) for (y = 1; y < SegImage->NY-1; y++)  for (x = 1; x < SegImage->NX-1; x++) if (fabs(SegImage->G(x,y,z)-IDregion)<epsilon){
    if(fabs(SegImage->G(x,y,z)-SegImage->G(x+1,y,z))>epsilon){
      DistMap->P(1,x,y,z);      NeaBoun.Add(0.5*dx,0,x,y,z);   
    }
    
    if(fabs(SegImage->G(x,y,z)-SegImage->G(x,y+1,z))>epsilon){
      DistMap->P(1,x,y,z);      NeaBoun.Add(0.5*dy,1,x,y,z);  
    }
    
    if(fabs(SegImage->G(x,y,z)-SegImage->G(x,y,z+1))>epsilon){
      DistMap->P(1,x,y,z);      NeaBoun.Add(0.5*dz,2,x,y,z);
    }
    
    if(fabs(SegImage->G(x,y,z)-SegImage->G(x-1,y,z))>epsilon){
      DistMap->P(1,x,y,z);      NeaBoun.Add(-0.5*dx,0,x,y,z);   
    }
    
    if(fabs(SegImage->G(x,y,z)-SegImage->G(x,y-1,z))>epsilon){
      DistMap->P(1,x,y,z);      NeaBoun.Add(-0.5*dy,1,x,y,z);  
    }
    
    if(fabs(SegImage->G(x,y,z)-SegImage->G(x,y,z-1))>epsilon){
      DistMap->P(1,x,y,z);      NeaBoun.Add(-0.5*dz,2,x,y,z);
    }
  
  
    if (fabs(DistMap->G(x,y,z)-1)<epsilon){  // in the domain and at the boundary
      TmpDist=sqrt(NeaBoun.G(0,x,y,z)*NeaBoun.G(0,x,y,z)+NeaBoun.G(1,x,y,z)*NeaBoun.G(1,x,y,z)+NeaBoun.G(2,x,y,z)*NeaBoun.G(2,x,y,z));
      
      if (TmpDist<0.00000001) { NeaBoun.Add(0.1*dx,0,x,y,z); TmpDist=0.1*dx;}
        
      DistMap->P(TmpDist,x,y,z);
    }
    else{   // in the domain but not at the boundary
      DistMap->P(MaxDist,x,y,z);
      NeaBoun.P(MaxDist,0,x,y,z);   NeaBoun.P(MaxDist,1,x,y,z);  NeaBoun.P(MaxDist,2,x,y,z);
    }
  }
  
  
  //2) DISTANCE PROPAGATION
  //At this point:
  //  * DistMap <- 0 outside of the domain / distance to the boundary at the boudary / MaxDist within the rest of the domain
  //  * NeaBoun <- 0 outside of the domain / nearest point to the boundary at the boudary / MaxDist^3 within the rest of the domain
  
  changes=1;
  while (changes>0){
    //cout << changes << endl;
    changes=0;
    for (z = 1; z < SegImage->NZ-1; z++) for (y = 1; y < SegImage->NY-1; y++)  for (x = 1; x < SegImage->NX-1; x++) if (fabs(SegImage->G(x,y,z)-IDregion)<epsilon) for (j=0;j<6;j++){
      if (j==0) {ix=1;iy=0;iz=0;}
      if (j==1) {ix=-1;}
      if (j==2) {ix=0;iy=1;}
      if (j==3) {iy=-1;}
      if (j==4) {iy=0;iz=1;}
      if (j==5) {iz=-1;};
      
      TmpDist=sqrt((NeaBoun.G(0,x,y,z)-ix*dx)*(NeaBoun.G(0,x,y,z)-ix*dx)+(NeaBoun.G(1,x,y,z)-iy*dy)*(NeaBoun.G(1,x,y,z)-iy*dy)+(NeaBoun.G(2,x,y,z)-iz*dz)*(NeaBoun.G(2,x,y,z)-iz*dz));
      
      if (TmpDist<DistMap->G(x+ix,y+iy,z+iz)){
        NeaBoun.P(NeaBoun.G(0,x,y,z)-ix*dx,0,x+ix,y+iy,z+iz);
        NeaBoun.P(NeaBoun.G(1,x,y,z)-iy*dy,1,x+ix,y+iy,z+iz);
        NeaBoun.P(NeaBoun.G(2,x,y,z)-iz*dz,2,x+ix,y+iy,z+iz);
        DistMap->P(TmpDist,x+ix,y+iy,z+iz);
        changes++;
      }
    }
  }
    
  
}



//Dilation of a scalar field
void ImageDilation(ScalarField * TreatedImage,int iterationNb){
  int it;
  int x,y,z,t;
  ScalarField TempField;
  float tmpfl;
  
  TempField.CreateVoidField(TreatedImage->NX,TreatedImage->NY,TreatedImage->NZ);
  
  
  for (it=0;it<iterationNb;it++){
    for (t=0;t<TreatedImage->NT;t++){
      //x direction
      for (z=0;z<TreatedImage->NZ;z++) for (y=0;y<TreatedImage->NY;y++){
        for (x=0;x<TreatedImage->NX;x++) {
          tmpfl=TreatedImage->G(x,y,z,t);
          if (x!=TreatedImage->NX-1) if (tmpfl<TreatedImage->G(x+1,y,z,t)) tmpfl=TreatedImage->G(x+1,y,z,t);
          if (x!=0)           if (tmpfl<TreatedImage->G(x-1,y,z,t)) tmpfl=TreatedImage->G(x-1,y,z,t);
          TempField.P(tmpfl,x,y,z);
        }
      } 
      
      //y direction
      for (z=0;z<TreatedImage->NZ;z++) for (x=0;x<TreatedImage->NX;x++){
        for (y=0;y<TreatedImage->NY;y++){
          tmpfl=TempField.G(x,y,z);
          if (y!=TreatedImage->NY-1) if (tmpfl<TempField.G(x,y+1,z)) tmpfl=TempField.G(x,y+1,z);
          if (y!=0)           if (tmpfl<TempField.G(x,y-1,z)) tmpfl=TempField.G(x,y-1,z);
          TreatedImage->P(tmpfl,x,y,z,t);
        }
      } 
      
      //z direction
      for (y=0;y<TreatedImage->NY;y++) for (x=0;x<TreatedImage->NX;x++){
        for (z=0;z<TreatedImage->NZ;z++){
          tmpfl=TreatedImage->G(x,y,z,t);
          if (z!=TreatedImage->NZ-1) if (tmpfl<TreatedImage->G(x,y,z+1,t)) tmpfl=TreatedImage->G(x,y,z+1,t);
          if (z!=0)           if (tmpfl<TreatedImage->G(x,y,z-1,t)) tmpfl=TreatedImage->G(x,y,z-1,t);
          TempField.P(tmpfl,x,y,z);
        }
      }
      
      //copy result
      for (z=0;z<TreatedImage->NZ;z++) for (y=0;y<TreatedImage->NY;y++) for (x=0;x<TreatedImage->NX;x++) TreatedImage->P(TempField.G(x,y,z),x,y,z,t);
    }
  }
}



//Erosion of a scalar field
void ImageErosion(ScalarField * TreatedImage,int iterationNb){
  int it;
  int x,y,z,t;
  ScalarField TempField;
  float tmpfl;
  
  TempField.CreateVoidField(TreatedImage->NX,TreatedImage->NY,TreatedImage->NZ);
  
  
  for (it=0;it<iterationNb;it++){
    for (t=0;t<TreatedImage->NT;t++){
      //x direction
      for (z=0;z<TreatedImage->NZ;z++) for (y=0;y<TreatedImage->NY;y++){
        for (x=0;x<TreatedImage->NX;x++) {
          tmpfl=TreatedImage->G(x,y,z,t);
          if (x!=TreatedImage->NX-1) if (tmpfl>TreatedImage->G(x+1,y,z,t)) tmpfl=TreatedImage->G(x+1,y,z,t);
          if (x!=0)           if (tmpfl>TreatedImage->G(x-1,y,z,t)) tmpfl=TreatedImage->G(x-1,y,z,t);
          TempField.P(tmpfl,x,y,z);
        }
      } 
      
      //y direction
      for (z=0;z<TreatedImage->NZ;z++) for (x=0;x<TreatedImage->NX;x++){
        for (y=0;y<TreatedImage->NY;y++){
          tmpfl=TempField.G(x,y,z);
          if (y!=TreatedImage->NY-1) if (tmpfl>TempField.G(x,y+1,z)) tmpfl=TempField.G(x,y+1,z);
          if (y!=0)           if (tmpfl>TempField.G(x,y-1,z)) tmpfl=TempField.G(x,y-1,z);
          TreatedImage->P(tmpfl,x,y,z,t);
        }
      } 
      
      //z direction
      for (y=0;y<TreatedImage->NY;y++) for (x=0;x<TreatedImage->NX;x++){
        for (z=0;z<TreatedImage->NZ;z++){
          tmpfl=TreatedImage->G(x,y,z,t);
          if (z!=TreatedImage->NZ-1) if (tmpfl>TreatedImage->G(x,y,z+1,t)) tmpfl=TreatedImage->G(x,y,z+1,t);
          if (z!=0)           if (tmpfl>TreatedImage->G(x,y,z-1,t)) tmpfl=TreatedImage->G(x,y,z-1,t);
          TempField.P(tmpfl,x,y,z);
        }
      }
      
      //copy result
      for (z=0;z<TreatedImage->NZ;z++) for (y=0;y<TreatedImage->NY;y++) for (x=0;x<TreatedImage->NX;x++) TreatedImage->P(TempField.G(x,y,z),x,y,z,t);
    }
  }
}


///Smooth the field 'SmoothedField' for fluid like regularisation with sliding conditions at th boundaries of 'LocMask'.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> 'NeaBoun' must contain a vector field of the same size as 'SmoothedField' (treated as a temporary field)
///-> 'NormalCmp'  must contain a scalar field of the same size as 'SmoothedField' (treated as a temporary field)
///-> 'TempSF'  must contain a scalar field of the same size as 'SmoothedField' (treated as a temporary field)
void SmoothNormalAndTangentContributions(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB,int NbRegions, float * IdRegions,VectorField * NeaBoun,ScalarField * NormalCmp,ScalarField * TempSF, float dx,float dy, float dz){
  int i;
  int x,y,z;
  int x2,y2,z2;
  int direction;
  float normLoc;
  float ScalProd;
  float MaxDist_DistMap;
  float DeltaTime;
  float SecureDistanceTangentContrib;
  
  SecureDistanceTangentContrib=2;
  
  //1) compute the velocity field of the nearest boundaries to the regions in LocMask
  MaxDist_DistMap=4;
  if (sqrt(2*TimeSmooth)*3>MaxDist_DistMap)
    MaxDist_DistMap=sqrt(2*TimeSmooth)*3;
  
  Cpt_NearestBoundaryOLD(LocMask,NeaBoun,TempSF,MaxDist_DistMap,NbRegions,IdRegions,dx,dy,dz); //if TempSF==1 the map is computed (0 elsewhere)
  Diffusion_3D(NeaBoun,1, 1.5,1,dx,dy,dz);
  
  //NeaBoun->Write("NeaBounX.nii","NeaBounY.nii","NeaBounZ.nii");
  
  
  //2) put to 0 the contributions at the very boundary of the regions (for numerical reasons)
  for (z = 3; z < LocMask->NZ-3; z++)  for (y = 3; y < LocMask->NY-3; y++) for (x = 3; x < LocMask->NX-3; x++){
    if ((fabs(LocMask->G(x,y,z)-LocMask->G(x-1,y,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y-1,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y,z-1))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x+1,y,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y+1,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y,z+1))>0.1)){
      for (z2 = -3; z2 <= 3; z2++)  for (y2 = -3; y2 <=3 ; y2++) for (x2 = -3; x2<=3; x2++){
        SmoothedField->P(0,0,x+x2,y+y2,z+z2);
        SmoothedField->P(0,1,x+x2,y+y2,z+z2);
        SmoothedField->P(0,2,x+x2,y+y2,z+z2);
      }
    }
  }
  
  
  
  
  //2) smooth the information in the 3 directions X, Y, Z
  for (direction=0;direction<3;direction++){
    //init the normal component
    for (z = 0; z < SmoothedField->NZ; z++) for (y = 0; y < SmoothedField->NY; y++)  for (x = 0; x < SmoothedField->NX; x++)
      NormalCmp->P(0,x,y,z);
    
    //2.1) extract the normal and tangential componant of SmoothedField (in NormalCmp and SmoothedField respectively)
    for (z = 2; z < SmoothedField->NZ-2; z++) for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++) if (TempSF->G(x,y,z)>0.5){
      
        normLoc=sqrt(NeaBoun->G(0,x,y,z)*NeaBoun->G(0,x,y,z)+NeaBoun->G(1,x,y,z)*NeaBoun->G(1,x,y,z)+NeaBoun->G(2,x,y,z)*NeaBoun->G(2,x,y,z));
        
        //2.1.1) detect the region boundaries
        TempSF->P(1+normLoc,x,y,z);
        
        //2.1.2) separate the normal and tengent contributions
        if (normLoc>0.01){
          //compute the scalar product
          ScalProd=SmoothedField->G(0,x,y,z)*NeaBoun->G(0,x,y,z)+SmoothedField->G(1,x,y,z)*NeaBoun->G(1,x,y,z)+SmoothedField->G(2,x,y,z)*NeaBoun->G(2,x,y,z);
          
          //compute the normal component
          NormalCmp->P(NeaBoun->G(direction,x,y,z)*ScalProd/(normLoc*normLoc),x,y,z);
          
          //compute the tangent component
          SmoothedField->P(SmoothedField->G(direction,x,y,z)-NormalCmp->G(x,y,z),direction,x,y,z);
        }
        if (normLoc<0.01) NormalCmp->P(0,x,y,z);
        if (normLoc<SecureDistanceTangentContrib) SmoothedField->P(0,direction,x,y,z);
    }
    
    //TempSF->Write("tempSF.nii");
    
    //2.2) diffuse the normal component in the different regions  
    //Diffusion_3D(NormalCmp,1, TimeSmooth/2,2,dx,dy,dz);
    Diffusion_3D(NormalCmp,1, TimeSmooth,1,dx,dy,dz);
    
    
    //if (direction==0) SmoothedField->Write("totoX1.nii","toto.nii","toto.nii");
    //if (direction==1) SmoothedField->Write("toto.nii","totoY2.nii","toto.nii");
    //if (direction==2) SmoothedField->Write("toto.nii","toto.nii","totoZ3.nii");
    
    //2.3) diffuse the tangential component in the different regions  

    DeltaTime=TimeSmooth/ITERATIONS_NB;
    
    for (i=0;i<NbRegions;i++){
      //Diffusion_3D(SmoothedField,LocMask,IdRegions[i], 1, DeltaTime,static_cast<int>(0.0001+TimeSmooth/DeltaTime),0,direction,dx,dy,dz);
      Diffusion_3D(SmoothedField,LocMask,IdRegions[i], 1, DeltaTime,ITERATIONS_NB,0,direction,dx,dy,dz);
    }
    
    //2.4) recompose the smoothed components
    for (z = 0; z < SmoothedField->NZ; z++) for (y = 0; y < SmoothedField->NY; y++)  for (x = 0; x < SmoothedField->NX; x++) if (TempSF->G(x,y,z)>0.5){
      SmoothedField->P(SmoothedField->G(direction,x,y,z)+NormalCmp->G(x,y,z),direction,x,y,z);
    }

  
  }
  
}






///Smooth the field 'SmoothedField' for diffusion like regularisation with sliding conditions at the boundaries of 'LocMask'.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> If 'NoDistMapEstim' == 1 the distance map (NeaBoun) is not estimated here. In that case 'TempSF' is also not reestimated.
void SmoothWithinRegions(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB, int NbRegions,float * IdRegions, float dx,float dy, float dz,int DirchletNeumann){
  int x,y,z;
  int i;
  float epsilon;
  float normLoc;
  float DeltaTime;
  float MaxDist_DistMap;
  int direction;
  int x2,y2,z2;
  
  epsilon=0.000001;
  DeltaTime=TimeSmooth/ITERATIONS_NB;
  
  //A) Dirichlet conditions (default)
  if (DirchletNeumann>=0){
    for (i=0;i<NbRegions;i++) Diffusion_3Dbis(SmoothedField,LocMask,IdRegions[i], 1, DeltaTime,ITERATIONS_NB,1,dx,dy,dz);
    
  }
  else{//B) Neumann conditions
    for (direction=0;direction<3;direction++) for (i=0;i<NbRegions;i++){
      Diffusion_3D(SmoothedField,LocMask,IdRegions[i], 1, DeltaTime,ITERATIONS_NB,0,direction,dx,dy,dz);
    }
  }
}


///Smooth the field 'SmoothedField' for diffusion like regularisation with sliding conditions at the boundaries of 'LocMask'. 
///Contrary to what is done in SmoothWithinRegions, we put to 0 the contributions too close to the bounaries.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary in 'LocMask'
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which the SmoothedField is set to zero before smoothing
void SmoothWithinRegions2(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB, int NbRegions,float * IdRegions,VectorField * NeaBoun,ScalarField * TempSF, float dx,float dy, float dz,int DistMapToEstim){
  int x,y,z;
  int x2,y2,z2;
  int i,j;
  float DeltaTime;
  int direction;
  float TmpFl;
  
  
  DeltaTime=TimeSmooth/ITERATIONS_NB;
  
  //1) compute the local distance map
  
  if (DistMapToEstim==1){
    Cpt_NearestBoundaryOLD(LocMask,NeaBoun,TempSF,100,NbRegions,IdRegions,dx,dy,dz,0);
    Diffusion_3D(NeaBoun,1, 2,2,dx,dy,dz);
  }
  
  //2) Set to zero all the contributions too close to the boundary
  if (LocMask->NZ>1){ //3D image
    for (z = 2; z < LocMask->NZ-2; z++)  for (y = 2; y < LocMask->NY-2; y++) for (x = 2; x < LocMask->NX-2; x++){
      
      TmpFl=0;
      for (z2 = -2; z2 < 3; z2++)  for (y2 = -2; y2 < 3; y2++) for (x2 = -2; x2 < 3; x2++) TmpFl+=fabs(LocMask->G(x,y,z)-LocMask->G(x+x2,y+y2,z+z2));
                 
      if (TmpFl>0.1){
        SmoothedField->P(0,0,x,y,z);
        SmoothedField->P(0,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
      }
    }
  }
  else{ //2D image
    for (y = 2; y < LocMask->NY-2; y++) for (x = 2; x < LocMask->NX-2; x++){
      if ((fabs(LocMask->G(x,y,0)-LocMask->G(x-2,y,0))>0.1)||(fabs(LocMask->G(x,y,0)-LocMask->G(x,y-2,0))>0.1)||(fabs(LocMask->G(x,y,0)-LocMask->G(x+2,y,0))>0.1)||(fabs(LocMask->G(x,y,0)-LocMask->G(x,y+2,0))>0.1)){
        SmoothedField->P(0,0,x,y,0);
        SmoothedField->P(0,1,x,y,0);
      }
    }
  }
  
  
  //3) smooth the vector field within each region
  for (direction=0;direction<3;direction++) for (i=0;i<NbRegions;i++){
    Diffusion_3D(SmoothedField,LocMask,IdRegions[i], 1, DeltaTime,ITERATIONS_NB,0,direction,dx,dy,dz);
  }
}



///Remove the normal contributions of a velocity field that are too close to a boundary
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which SmoothedField has reduced normal contirbutions
void RemoveNormalContributions(VectorField * SmoothedField,VectorField * NeaBoun,ScalarField * TempSF,float BoundaMargin, float dx,float dy, float dz,int SetBoundaryToZero){
  int x,y,z;
  int i,j;
  float epsilon;
  float normLoc;
  float normLoc3,ScalProd;
  float DeltaTime;
  float MaxDist_DistMap;
  int direction;
  int x2,y2,z2;
  float NCx,NCy,NCz;
  float NBlocX,NBlocY,NBlocZ;
  int dxyz;
  
  
  //cout << "Remove normal contributions" << endl;
  
  if (SmoothedField->NZ>1){ //3D image
    dxyz=dx;
    if (dy<dxyz) dxyz=dy;
    if (dz<dxyz) dxyz=dz;
    
    for (z = 2; z < SmoothedField->NZ-2; z++) for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++) if (TempSF->G(x,y,z)>0.5){
      //1) INIT
      //put the vector representing the nearest boundary in NCx,NCy,NCz
      NCx=NeaBoun->G(0,x,y,z);
      NCy=NeaBoun->G(1,x,y,z);
      NCz=NeaBoun->G(2,x,y,z);
      
      normLoc=sqrt(NCx*NCx+NCy*NCy+NCz*NCz); //distance to the boundary in millimiters
      
      //2) REMOVE THE NORMAL CONTRIBUTION
      if ((normLoc<=0.5*dxyz)&&(SetBoundaryToZero==1)){
        NCx=0; NCy=0; NCz=0; normLoc=0;
        if (TempSF->G(x+1,y,z)>0.5) NCx=NeaBoun->G(0,x+1,y,z); NCy=NeaBoun->G(1,x+1,y,z); NCz=NeaBoun->G(2,x+1,y,z); normLoc++;
        if (TempSF->G(x-1,y,z)>0.5) NCx=NeaBoun->G(0,x-1,y,z); NCy=NeaBoun->G(1,x-1,y,z); NCz=NeaBoun->G(2,x-1,y,z); normLoc++;
        if (TempSF->G(x,y+1,z)>0.5) NCx=NeaBoun->G(0,x,y+1,z); NCy=NeaBoun->G(1,x,y+1,z); NCz=NeaBoun->G(2,x,y+1,z); normLoc++;
        if (TempSF->G(x,y-1,z)>0.5) NCx=NeaBoun->G(0,x,y-1,z); NCy=NeaBoun->G(1,x,y-1,z); NCz=NeaBoun->G(2,x,y-1,z); normLoc++;
        if (TempSF->G(x,y,z+1)>0.5) NCx=NeaBoun->G(0,x,y,z+1); NCy=NeaBoun->G(1,x,y,z+1); NCz=NeaBoun->G(2,x,y,z+1); normLoc++;
        if (TempSF->G(x,y,z-1)>0.5) NCx=NeaBoun->G(0,x,y,z-1); NCy=NeaBoun->G(1,x,y,z-1); NCz=NeaBoun->G(2,x,y,z-1); normLoc++;
        
        if (normLoc>0.5){
          NCx/=normLoc;
          NCy/=normLoc;
          NCz/=normLoc;
          normLoc=sqrt(NCx*NCx+NCy*NCy+NCz*NCz);
        }
      }
      
      if ((normLoc<BoundaMargin)&&(normLoc>0.5*dxyz)){
        //normalise NeaBoun   (in mm)
        NCx/=normLoc;
        NCy/=normLoc;
        NCz/=normLoc;
        
        
        //compute the scalar product with SmoothedField    (in mm)
        ScalProd=SmoothedField->G(0,x,y,z)*dx*NCx+SmoothedField->G(1,x,y,z)*dy*NCy+SmoothedField->G(2,x,y,z)*dz*NCz;
        
        //add a weight to ScalProd to take into account the distance to the bounadry
        ScalProd*=((normLoc/BoundaMargin)-1)*((normLoc/BoundaMargin)-1);
        
        //compute the normal contribution of SmoothedField * its weight and transform it in voxels
        NCx*=(ScalProd)/dx;
        NCy*=(ScalProd)/dy;
        NCz*=(ScalProd)/dz;
        
        //remove the normal contribution
        SmoothedField->P(SmoothedField->G(0,x,y,z)-NCx,0,x,y,z);
        SmoothedField->P(SmoothedField->G(1,x,y,z)-NCy,1,x,y,z);
        SmoothedField->P(SmoothedField->G(2,x,y,z)-NCz,2,x,y,z);
        
      }
      
      if ((normLoc<=0.5*dxyz)&&(SetBoundaryToZero==1)){
        SmoothedField->P(0,0,x,y,z);
        SmoothedField->P(0,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
      }

    }
  }
  else{ //2D image
    dxyz=dx;
    if (dy<dxyz) dxyz=dy;
    z=0;
    for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++) if (TempSF->G(x,y,z)>0.5){
      //1) INIT
      //put the vector representing the nearest boundary in NCx,NCy,NCz
      NCx=NeaBoun->G(0,x,y,z);
      NCy=NeaBoun->G(1,x,y,z);
      
      normLoc=sqrt(NCx*NCx+NCy*NCy); //distance to the boundary in millimiters
      
      //2) REMOVE THE NORMAL CONTRIBUTION
      if ((normLoc<=0.5*dxyz)&&(SetBoundaryToZero==1)){
        NCx=0; NCy=0; normLoc=0;
        if (TempSF->G(x+1,y,z)>0.5) NCx=NeaBoun->G(0,x+1,y,z); NCy=NeaBoun->G(1,x+1,y,z); normLoc++;
        if (TempSF->G(x-1,y,z)>0.5) NCx=NeaBoun->G(0,x-1,y,z); NCy=NeaBoun->G(1,x-1,y,z); normLoc++;
        if (TempSF->G(x,y+1,z)>0.5) NCx=NeaBoun->G(0,x,y+1,z); NCy=NeaBoun->G(1,x,y+1,z); normLoc++;
        if (TempSF->G(x,y-1,z)>0.5) NCx=NeaBoun->G(0,x,y-1,z); NCy=NeaBoun->G(1,x,y-1,z); normLoc++;
        
        if (normLoc>0.5){
          NCx/=normLoc;
          NCy/=normLoc;
          normLoc=sqrt(NCx*NCx+NCy*NCy);
        }
      }
      
      if ((normLoc<BoundaMargin)&&(normLoc>0.5*dxyz)){
        //normalise NeaBoun   (in mm)
        NCx/=normLoc;
        NCy/=normLoc;
        
        //compute the scalar product with SmoothedField    (in mm)
        ScalProd=SmoothedField->G(0,x,y,z)*dx*NCx+SmoothedField->G(1,x,y,z)*dy*NCy;
        
        //add a weight to ScalProd to take into account the distance to the bounadry
        ScalProd*=((normLoc/BoundaMargin)-1)*((normLoc/BoundaMargin)-1);
        
        //compute the normal contribution of SmoothedField * its weight and transform it in voxels
        NCx*=(ScalProd)/dx;
        NCy*=(ScalProd)/dy;
        
        //remove the normal contribution
        SmoothedField->P(SmoothedField->G(0,x,y,z)-NCx,0,x,y,z);
        SmoothedField->P(SmoothedField->G(1,x,y,z)-NCy,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
      }
      
      if ((normLoc<=0.5*dxyz)&&(SetBoundaryToZero==1)){
        SmoothedField->P(0,0,x,y,z);
        SmoothedField->P(0,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
      }

    }
  }
}



///Remove the normal contributions of a velocity field that are too close to a boundary
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which SmoothedField has reduced normal contirbutions
void RemoveNormalContributionsOLD(VectorField * SmoothedField,VectorField * NeaBoun,ScalarField * TempSF,float BoundaMargin, float dx,float dy, float dz){
  int x,y,z;
  int i,j;
  float epsilon;
  float normLoc;
  float normLoc3,ScalProd;
  float DeltaTime;
  float MaxDist_DistMap;
  int direction;
  int x2,y2,z2;
  float NCx,NCy,NCz;
  float NBlocX,NBlocY,NBlocZ;
  
  
  //cout << "Remove normal contributions" << endl;
  
  //4) remove the normal contributions in the ngbh of the boundary
  if (SmoothedField->NZ>1){ //3D image
    for (z = 2; z < SmoothedField->NZ-2; z++) for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++) if (TempSF->G(x,y,z)>0.5){
      
      
      //1) INIT
      //put the vector representing the nearest boundary in NCx,NCy,NCz
      NCx=NeaBoun->G(0,x,y,z);
      NCy=NeaBoun->G(1,x,y,z);
      NCz=NeaBoun->G(2,x,y,z);
      
      normLoc=sqrt(NCx*NCx+NCy*NCy+NCz*NCz); //distance to the boundary in millimiters
      
      //2) REMOVE THE NORMAL CONTRIBUTION
      if (normLoc<0.2){
	NCx=0; NCy=0; NCz=0; normLoc=0;
	if (TempSF->G(x+1,y,z)>0.5) NCx=NeaBoun->G(0,x+1,y,z); NCy=NeaBoun->G(1,x+1,y,z); NCz=NeaBoun->G(2,x+1,y,z); normLoc++;
	if (TempSF->G(x-1,y,z)>0.5) NCx=NeaBoun->G(0,x-1,y,z); NCy=NeaBoun->G(1,x-1,y,z); NCz=NeaBoun->G(2,x-1,y,z); normLoc++;
	if (TempSF->G(x,y+1,z)>0.5) NCx=NeaBoun->G(0,x,y+1,z); NCy=NeaBoun->G(1,x,y+1,z); NCz=NeaBoun->G(2,x,y+1,z); normLoc++;
	if (TempSF->G(x,y-1,z)>0.5) NCx=NeaBoun->G(0,x,y-1,z); NCy=NeaBoun->G(1,x,y-1,z); NCz=NeaBoun->G(2,x,y-1,z); normLoc++;
	if (TempSF->G(x,y,z+1)>0.5) NCx=NeaBoun->G(0,x,y,z+1); NCy=NeaBoun->G(1,x,y,z+1); NCz=NeaBoun->G(2,x,y,z+1); normLoc++;
	if (TempSF->G(x,y,z-1)>0.5) NCx=NeaBoun->G(0,x,y,z-1); NCy=NeaBoun->G(1,x,y,z-1); NCz=NeaBoun->G(2,x,y,z-1); normLoc++;
	
	if (normLoc>0.5){
	  NCx/=normLoc;
	  NCy/=normLoc;
	  NCz/=normLoc;
	  normLoc=sqrt(NCx*NCx+NCy*NCy+NCz*NCz);
	}
	else{
        SmoothedField->P(0,0,x,y,z);
        SmoothedField->P(0,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
	}
      }
      
      if ((normLoc<BoundaMargin)&&(normLoc>0.2)){
        //normalise NeaBoun   (in mm)
        NCx/=normLoc;
        NCy/=normLoc;
        NCz/=normLoc;
        
        
        //compute the scalar product with SmoothedField    (in mm)
        ScalProd=SmoothedField->G(0,x,y,z)*dx*NCx+SmoothedField->G(1,x,y,z)*dy*NCy+SmoothedField->G(2,x,y,z)*dz*NCz;
        
        //add a weight to ScalProd to take into account the distance to the bounadry
        ScalProd*=((normLoc/BoundaMargin)-1)*((normLoc/BoundaMargin)-1);
        
        //compute the normal contribution of SmoothedField * its weight and transform it in voxels
        NCx*=(ScalProd)/dx;
        NCy*=(ScalProd)/dy;
        NCz*=(ScalProd)/dz;
        
        //remove the normal contribution
        SmoothedField->P(SmoothedField->G(0,x,y,z)-NCx,0,x,y,z);
        SmoothedField->P(SmoothedField->G(1,x,y,z)-NCy,1,x,y,z);
        SmoothedField->P(SmoothedField->G(2,x,y,z)-NCz,2,x,y,z);
        
      }
    }
  }
  else{ //2D image
    cout << "Treatment of normal contributions in the 2D case -> TO DO" << endl;
  }
}




///
///Smooth the field 'SmoothedField' for diffusion like regularisation with sliding conditions at the boundaries of 'LocMask'. 
///Contrary to what is done in SmoothWithinRegions, we reduce the contribution of the normals close to the bounary.
///-> The field is diffused during 'TimeSmooth' time units for time steps of 'DeltaTime' time units
///-> 'LocMask' contains 'NbRegions' regions. Their identifiers are listed in 'IdRegions'
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary in 'LocMask'
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadrie (in voxels) in which the SmoothedField is set to zero before smoothing
///-> If 'NoDistMapEstim' == 1 the distance map (NeaBoun) is not estimated here. In that case 'TempSF' is also not reestimated.
///DEPRECATED FUNCTION - SHOULD NOT BE USED
void SmoothWithinRegionsWithSliding(VectorField * SmoothedField,ScalarField * LocMask,float TimeSmooth,int ITERATIONS_NB, int NbRegions,float * IdRegions,VectorField * NeaBoun,ScalarField * TempSF, float dx,float dy, float dz,float BoundaMargin,int NoDistMapEstim,int BoundaryCompensation){
  int x,y,z;
  int i,j;
  float epsilon;
  float normLoc;
  float normLoc3,ScalProd;
  float DeltaTime;
  float MaxDist_DistMap;
  int direction;
  int x2,y2,z2;
  float NCx,NCy,NCz;
  float NBlocX,NBlocY,NBlocZ;
  
  
  cout << "USE OF A DEPRECATED FUNCTION TO SMOOTH A VECTOR FIELD";
  
  epsilon=0.000001;
  DeltaTime=TimeSmooth/ITERATIONS_NB;
  
  //1) compute the local distance map
  MaxDist_DistMap=BoundaMargin*1.5;
  
  if (NoDistMapEstim!=1){
    Cpt_NearestBoundaryOLD(LocMask,NeaBoun,TempSF,MaxDist_DistMap,NbRegions,IdRegions,dx,dy,dz,0);
    Diffusion_3D(NeaBoun,1, 2,2,dx,dy,dz);
  }
  
  //2) Set to zero all the contributions too close to the boundary
  if (LocMask->NZ>1){ //3D image
    for (z = 1; z < LocMask->NZ-1; z++)  for (y = 1; y < LocMask->NY-1; y++) for (x = 1; x < LocMask->NX-1; x++){
      if ((fabs(LocMask->G(x,y,z)-LocMask->G(x-1,y,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y-1,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y,z-1))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x+1,y,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y+1,z))>0.1)||(fabs(LocMask->G(x,y,z)-LocMask->G(x,y,z+1))>0.1)){
        SmoothedField->P(0,0,x,y,z);
        SmoothedField->P(0,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
      }
    }
  }
  else{ //2D image
    for (y = 1; y < LocMask->NY-1; y++) for (x = 1; x < LocMask->NX-1; x++){
      if ((fabs(LocMask->G(x,y,0)-LocMask->G(x-1,y,0))>0.1)||(fabs(LocMask->G(x,y,0)-LocMask->G(x,y-1,0))>0.1)||(fabs(LocMask->G(x,y,0)-LocMask->G(x+1,y,0))>0.1)||(fabs(LocMask->G(x,y,0)-LocMask->G(x,y+1,0))>0.1)){
        SmoothedField->P(0,0,x,y,0);
        SmoothedField->P(0,1,x,y,0);
      }
    }
  }
  
  
  //3) smooth the vector field within each region
  for (direction=0;direction<3;direction++) for (i=0;i<NbRegions;i++){
    Diffusion_3D(SmoothedField,LocMask,IdRegions[i], 1, DeltaTime,ITERATIONS_NB,0,direction,dx,dy,dz);
  }
  
  //4) remove the normal contributions in the ngbh of the boundary
  if (LocMask->NZ>1){ //3D image
    for (z = 2; z < SmoothedField->NZ-2; z++) for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++) if (TempSF->G(x,y,z)>0.5){
      
      normLoc=sqrt(NeaBoun->G(0,x,y,z)*NeaBoun->G(0,x,y,z)+NeaBoun->G(1,x,y,z)*NeaBoun->G(1,x,y,z)+NeaBoun->G(2,x,y,z)*NeaBoun->G(2,x,y,z)); //distance to the boundary in millimiters
      
      if (normLoc<0.5){
        SmoothedField->P(0,0,x,y,z);
        SmoothedField->P(0,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
      }
      else if ((normLoc<BoundaMargin)&&(BoundaryCompensation==1)){
        //normalise NeaBoun
        NCx=NeaBoun->G(0,x,y,z)/normLoc;
        NCy=NeaBoun->G(1,x,y,z)/normLoc;
        NCz=NeaBoun->G(2,x,y,z)/normLoc;
        
        //compute the scalar product with SmoothedField
        ScalProd=SmoothedField->G(0,x,y,z)*dx*NCx+SmoothedField->G(1,x,y,z)*dy*NCy+SmoothedField->G(2,x,y,z)*dz*NCz;
        
        //compute the normal contribution of SmoothedField * its weight
        if (normLoc<BoundaMargin/2){
          NCx*=(ScalProd/dx);
          NCy*=(ScalProd/dy);
          NCz*=(ScalProd/dz);
        }
        else{
          NCx*=(ScalProd/dx)*(2-2*normLoc/BoundaMargin);
          NCy*=(ScalProd/dy)*(2-2*normLoc/BoundaMargin);
          NCz*=(ScalProd/dz)*(2-2*normLoc/BoundaMargin);
        }
        
        NCx/=2;
        NCy/=2;
        NCz/=2;
        
        //remove the normal contribution
        SmoothedField->P(SmoothedField->G(0,x,y,z)-NCx,0,x,y,z);
        SmoothedField->P(SmoothedField->G(1,x,y,z)-NCy,1,x,y,z);
        SmoothedField->P(SmoothedField->G(2,x,y,z)-NCz,2,x,y,z);
      }
    }
  }
  else{ //2D image
    for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++) if (TempSF->G(x,y,0)>0.5){
      
      normLoc=sqrt(NeaBoun->G(0,x,y,0)*NeaBoun->G(0,x,y,0)+NeaBoun->G(1,x,y,0)*NeaBoun->G(1,x,y,0)); //distance to the boundary in millimiters
      
      if ((normLoc<BoundaMargin)&&(BoundaryCompensation==1)){
        SmoothedField->P(SmoothedField->G(0,x,y,0)*normLoc/BoundaMargin,0,x,y,0);
        SmoothedField->P(SmoothedField->G(1,x,y,0)*normLoc/BoundaMargin,1,x,y,0);
      }
    }
  }

}




//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
void Cpt_Grad_ScalarField(ScalarField * SField,VectorField * Gradient,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  NBX=SField->NX;
  NBY=SField->NY;
  NBZ=SField->NZ;
  NBT=SField->NT;
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
    
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)){
      Gradient->CreateVoidField(NBX,NBY,NBZ);
      cout << "Gradient added in Cpt_Grad_ScalarField\n";
    }
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) gradient in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z);
      Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z);
      Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z);
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0);
      Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)||(NBT!=Gradient->NT))
      Gradient->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //gradient in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z,t);
        Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z,t);
        Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z,t);
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0,t);
        Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0,t);
      }
    }
  }
}



//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
void Cpt_Grad_MaskedScalarField(ScalarField * SField,VectorField * Gradient,ScalarField * Mask,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  NBX=SField->NX;
  NBY=SField->NY;
  NBZ=SField->NZ;
  NBT=SField->NT;
  float epsilon;
  
  epsilon=0.0001;
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
    
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)){
      Gradient->CreateVoidField(NBX,NBY,NBZ);
      cout << "Gradient added in Cpt_Grad_ScalarField\n";
    }
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) gradient in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z);
      Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z);
      Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z);
      
      if ((fabs(Mask->G(x,y,z,t)-Mask->G(x+1,y,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x-1,y,z,t))>epsilon)) 
        Gradient->P(0,0,x,y,z);
      
      if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y+1,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y-1,z,t))>epsilon)) 
        Gradient->P(0,1,x,y,z);
      
      if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z+1,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z-1,t))>epsilon)) 
        Gradient->P(0,2,x,y,z);
      
      if (Mask->G(x,y,z,t)<epsilon){
        Gradient->P(0,0,x,y,z);
        Gradient->P(0,1,x,y,z);
        Gradient->P(0,2,x,y,z);
      } 
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0);
      Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)||(NBT!=Gradient->NT))
      Gradient->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //gradient in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z,t);
        Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z,t);
        Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z,t);
        
        if ((fabs(Mask->G(x,y,z,t)-Mask->G(x+1,y,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x-1,y,z,t))>epsilon)) 
          Gradient->P(0,0,x,y,z);
        
        if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y+1,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y-1,z,t))>epsilon)) 
          Gradient->P(0,1,x,y,z);
        
        if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z+1,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z-1,t))>epsilon)) 
          Gradient->P(0,2,x,y,z);
        
        if (Mask->G(x,y,z,t)<epsilon){
          Gradient->P(0,0,x,y,z);
          Gradient->P(0,1,x,y,z);
          Gradient->P(0,2,x,y,z);
        } 
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0,t);
        Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0,t);
      }
    }
  }
}






//Compute (d VField(X) / d x) + (d VField(Y) / d y) + (d VField(Z) / d z) and put the result in 'GradScalVF'
//where 'VField' is a vector field and 'GradScalVF' a scalar field
void Cpt_Grad_Scal_VectorField(VectorField * VField,ScalarField * GradScalVF,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float GradX,GradY,GradZ;
  
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  NBT=VField->NT;
  
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
    
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=GradScalVF->NX)||(NBY!=GradScalVF->NY)||(NBZ!=GradScalVF->NZ))
      GradScalVF->CreateVoidField(NBX,NBY,NBZ);
    
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) sum of gradients in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      GradX=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
      GradY=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
      GradZ=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
      GradScalVF->P(GradX+GradY+GradZ,x,y,z);
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      GradX=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
      GradY=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
      GradScalVF->P(GradX+GradY,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=GradScalVF->NX)||(NBY!=GradScalVF->NY)||(NBZ!=GradScalVF->NZ)||(NBT!=GradScalVF->NT))
      GradScalVF->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //sum of gradients in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        GradX=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
        GradY=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
        GradZ=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
        GradScalVF->P(GradX+GradY+GradZ,x,y,z,t);
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        GradX=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
        GradY=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
        GradScalVF->P(GradX+GradY,x,y,0,t);
      }
    }
  }
}


#ifdef COMPILE_WITH_OPENMP

//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;
  
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  NBT=VField->NT;
  
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
    
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ))
      DetJ->CreateVoidField(NBX,NBY,NBZ);
    
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //BEGIN FORK FOR THREADS
    #pragma omp parallel default(shared) private(d11,d12,d13,d21,d22,d23,d31,d32,d33,x,y,z) 
    {
      //1.2.1) sum of gradients in direction x, y, z
      #pragma omp for
      for (y=1;y<NBY-1;y++){ 
         for (z=1;z<NBZ-1;z++) for (x=1;x<NBX-1;x++){
          d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
          d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
          d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
          d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
          d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
          d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
          d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
          d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
          d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
          DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z);
        }
      }
    //END FORK FOR THREADS
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
      DetJ->P(d11*d22-d21*d12,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ)||(NBT!=DetJ->NT))
      DetJ->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //BEGIN FORK FOR THREADS
    #pragma omp parallel default(shared) private(d11,d12,d13,d21,d22,d23,d31,d32,d33,x,y,z,t) 
    {
      //1.2) Calculations
      #pragma omp for
      for (t=0;t<NBT;t++){
        //sum of gradients in direction x, y, z
        for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
          d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
          d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
          d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
          d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
          d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
          d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
          d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
          d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
          d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
          DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z,t);
        }
        
        //boundaries at 0.
        z=0;
        for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        z=NBZ-1;
        for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        y=0;
        for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        y=NBY-1;
        for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        x=0;
        for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
        x=NBX-1;
        for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
        
        //2D image case
        if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
          d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
          d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
          d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
          d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
          DetJ->P(d11*d22-d21*d12,x,y,0,t);
        }
      }
    //END FORK FOR THREADS
    }
  }
}

#else

//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;
  
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  NBT=VField->NT;
  
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
    
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ))
      DetJ->CreateVoidField(NBX,NBY,NBZ);
    
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) sum of gradients in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
      d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
      d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
      d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
      d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
      d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
      DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z);
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
      DetJ->P(d11*d22-d21*d12,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ)||(NBT!=DetJ->NT))
      DetJ->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //sum of gradients in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
        d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
        d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
        d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
        d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
        d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
        d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
        d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
        d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
        DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z,t);
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
        d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
        d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
        d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
        DetJ->P(d11*d22-d21*d12,x,y,0,t);
      }
    }
  }
}

#endif




//VF = SF1 \diamond_{Diff} SF2
//which means: VF = - SF2 \nabla SF1                    
void Diamond_Diff(ScalarField * SF1,ScalarField * SF2,VectorField * VF){
  int x,y,z;
  
  //compute the gradient of SF1
  Cpt_Grad_ScalarField(SF1,VF);
  
  //multiply the gradient by "- SF2 "
  for (z=1;z<SF1->NZ-1;z++) for (y=1;y<SF1->NY-1;y++) for (x=1;x<SF1->NX-1;x++){
    VF->P(-VF->G(x,y,z,0)*SF2->G(x,y,z),x,y,z,0);
    VF->P(-VF->G(x,y,z,1)*SF2->G(x,y,z),x,y,z,1);
    VF->P(-VF->G(x,y,z,2)*SF2->G(x,y,z),x,y,z,2);
  }
  

  //multiply the gradient by "- SF2 "   (2D VERSION)
  if (SF1->NZ==0) for (y=1;y<SF1->NY-1;y++) for (x=1;x<SF1->NX-1;x++){
    VF->P(-VF->G(x,y,0,0)*SF2->G(x,y,0),x,y,0,0);
    VF->P(-VF->G(x,y,0,1)*SF2->G(x,y,0),x,y,0,1);
  }
}



//V  = SF1 \diamond_{Trans} SF2
//which means: V =  - \int_{R^3} ( SF2 \nabla SF1  )
void Diamond_Trans(ScalarField * SF1,ScalarField * SF2,float V[3]){
  VectorField VF;
  int x,y,z;
  int NbPoints;
  
  //init
  VF.CreateVoidField(SF1->NX,SF1->NY,SF1->NZ);
  
  // compute the vector field (SF2 \nabla SF1) to average (not the fastest solution but the simple to write)
  Diamond_Diff(SF1,SF2,&VF);
  
  //compute the average of (SF2 \nabla SF1)...
  
  //... init
  NbPoints=0;
  V[0]=0;
  V[1]=0;
  V[2]=0;
  
  //...   --> 3D field
  for (z=1;z<VF.NZ-1;z++) for (y=1;y<VF.NY-1;y++) for (x=1;x<VF.NX-1;x++){
    V[0]+=VF.G(x,y,z,0);
    V[1]+=VF.G(x,y,z,1);
    V[2]+=VF.G(x,y,z,2);
    NbPoints++;
  }
  
  //...   --> 2D field
  if (VF.NZ==0) for (y=1;y<VF.NY-1;y++) for (x=1;x<VF.NX-1;x++){
    V[0]+=VF.G(x,y,0,0);
    V[1]+=VF.G(x,y,0,1);
    NbPoints++;
  }
  
  //... finalise
  V[0]/=static_cast<float>(NbPoints);
  V[1]/=static_cast<float>(NbPoints);
  V[2]/=static_cast<float>(NbPoints);

}


//T = V1 \diamond_{R^3} V2 
//which means:  T =   -0.5  (V1 \otimes V2 - V2 \otimes V1 )  
void Diamond_R3(float V1[3],float V2[3],float T[3][3]){
  T[0][0]=0;                                 T[0][1]=-0.5*(V1[0]*V2[1]-V2[0]*V1[1]);  T[0][2]=-0.5*(V1[0]*V2[2]-V2[0]*V1[2]);  
  T[1][0]=-0.5*(V1[1]*V2[0]-V2[1]*V1[0]);    T[1][1]=0;                               T[1][2]=-0.5*(V1[1]*V2[2]-V2[1]*V1[2]);  
  T[2][0]=-0.5*(V1[2]*V2[0]-V2[2]*V1[0]);    T[2][1]=-0.5*(V1[2]*V2[1]-V2[2]*V1[1]);  T[2][2]=0;  
}


//T = SF1 \diamond_{SO(3)} SF2 
//which means:  T =   -0.5  \int_{R^3} \left(    SF2(x) (  (\nabla SF1)(x) \otimes x - x \otimes (\nabla SF1)(x) )  \right) dx
//Remark: 
//  x represents where we are in the image. 
//  Here it can have two interpretations: (1) image coordinates and (2) where we are depending on the center of rotation
//  If you don't want the center of rotation to be the image coordinate (0,0,0), fill the image coordinate you want in {Cx,Cy,Cz}
void Diamond_SO3(ScalarField * SF1,ScalarField * SF2,float T[3][3],float Cx,float Cy,float Cz){
  VectorField VF;
  int x,y,z;
  int NbPoints;
  float T_temp[3][3];
  float X_temp1[3];
  float X_temp2[3];
  
  
  //compute the gradient of SF1
  Cpt_Grad_ScalarField(SF1,&VF);

  //compute the rotation matrix...
  
  //... init
  T[0][0]=0;  T[0][1]=0;  T[0][2]=0;  
  T[1][0]=0;  T[1][1]=0;  T[1][2]=0;  
  T[2][0]=0;  T[2][1]=0;  T[2][2]=0;  
  
  
  //...   --> 3D field
  for (z=1;z<VF.NZ-1;z++) for (y=1;y<VF.NY-1;y++) for (x=1;x<VF.NX-1;x++){
    X_temp1[0]=VF.G(x,y,z,0);
    X_temp1[1]=VF.G(x,y,z,1);
    X_temp1[2]=VF.G(x,y,z,2);
    
    X_temp2[0]=static_cast<float>(x)-Cx;
    X_temp2[1]=static_cast<float>(y)-Cy;
    X_temp2[2]=static_cast<float>(z)-Cz;
    
    Diamond_R3(X_temp1,X_temp2,T_temp);
    
                                          T[0][1]+=T_temp[0][1]*SF2->G(x,y,z);  T[0][2]+=T_temp[0][2]*SF2->G(x,y,z);
    T[1][0]+=T_temp[1][0]*SF2->G(x,y,z);                                        T[1][2]+=T_temp[1][2]*SF2->G(x,y,z);  
    T[2][0]+=T_temp[2][0]*SF2->G(x,y,z);  T[2][1]+=T_temp[2][1]*SF2->G(x,y,z);    
    
    NbPoints++;
  }
  
  //...   --> 2D field
  if (VF.NZ==0) for (y=1;y<VF.NY-1;y++) for (x=1;x<VF.NX-1;x++){
    X_temp1[0]=VF.G(x,y,z,0);
    X_temp1[0]=VF.G(x,y,z,1);
    X_temp1[0]=0;
    
    X_temp2[0]=static_cast<float>(x)-Cx;
    X_temp2[1]=static_cast<float>(y)-Cy;
    X_temp2[2]=0;
    
    Diamond_R3(X_temp1,X_temp2,T_temp);
    
                                         T[0][1]+=T_temp[0][1]*SF2->G(x,y,z);
    T[1][0]+=T_temp[1][0]*SF2->G(x,y,z);                                       
    
    NbPoints++;
  }
  
  //finalise
  T[0][0]/=static_cast<float>(NbPoints);  T[0][1]/=static_cast<float>(NbPoints);  T[0][2]/=static_cast<float>(NbPoints);  
  T[1][0]/=static_cast<float>(NbPoints);  T[1][1]/=static_cast<float>(NbPoints);  T[1][2]/=static_cast<float>(NbPoints);  
  T[2][0]/=static_cast<float>(NbPoints);  T[2][1]/=static_cast<float>(NbPoints);  T[2][2]/=static_cast<float>(NbPoints);  
  
  
}





#ifdef COMPILE_WITH_OPENMP

///Integrate a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptDefFromSteadyVeloField(VectorField * VeloField,VectorField * DeformationField,int log2TimeStepNb,float MultFactor){
  int x,y,z;
  int i;
  float VecTempX,VecTempY,VecTempZ;
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=MultFactor/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,VecTempX,VecTempY,VecTempZ,i) 
  {
    //initiate DeformationField
    #pragma omp for
    for (y=0;y<VeloField->NY;y++) {
      for (z=0;z<VeloField->NZ;z++) for (x=0;x<VeloField->NX;x++){
        DeformationField->P(VeloField->G(0,x,y,z)*coef,0,x,y,z);
        DeformationField->P(VeloField->G(1,x,y,z)*coef,1,x,y,z);
        DeformationField->P(VeloField->G(2,x,y,z)*coef,2,x,y,z);
      }
    }
    
    //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
    for (i=0;i<log2TimeStepNb;i++){
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          VecTempX=x+DeformationField->G(0,x,y,z);
          VecTempY=y+DeformationField->G(1,x,y,z);
          VecTempZ=z+DeformationField->G(2,x,y,z);
          
          TempVeloField1.P(DeformationField->G(0,x,y,z)+DeformationField->G(0,VecTempX,VecTempY,VecTempZ),0,x,y,z);
          TempVeloField1.P(DeformationField->G(1,x,y,z)+DeformationField->G(1,VecTempX,VecTempY,VecTempZ),1,x,y,z);
          TempVeloField1.P(DeformationField->G(2,x,y,z)+DeformationField->G(2,VecTempX,VecTempY,VecTempZ),2,x,y,z);
        }
      }
      
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          DeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
          DeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
          DeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///Integrate a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptDefFromSteadyVeloField(VectorField * VeloField,VectorField * DeformationField,int log2TimeStepNb,float MultFactor){
  int x,y,z;
  int i;
  float VecTemp[3];
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=MultFactor/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  
  //initiate DeformationField
  for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
    DeformationField->P(VeloField->G(0,x,y,z)*coef,0,x,y,z);
    DeformationField->P(VeloField->G(1,x,y,z)*coef,1,x,y,z);
    DeformationField->P(VeloField->G(2,x,y,z)*coef,2,x,y,z);
  }
  
  //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
  for (i=0;i<log2TimeStepNb;i++){
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      VecTemp[0]=x+DeformationField->G(0,x,y,z);
      VecTemp[1]=y+DeformationField->G(1,x,y,z);
      VecTemp[2]=z+DeformationField->G(2,x,y,z);
      
      TempVeloField1.P(DeformationField->G(0,x,y,z)+DeformationField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
      TempVeloField1.P(DeformationField->G(1,x,y,z)+DeformationField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
      TempVeloField1.P(DeformationField->G(2,x,y,z)+DeformationField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
    }
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      DeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
      DeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
      DeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
    }
  }
}

#endif

#ifdef COMPILE_WITH_OPENMP

///Integrate the inverse of a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptInvDefFromSteadyVeloField(VectorField * VeloField,VectorField * InvDeformationField,int log2TimeStepNb,float MultFactor){
  int x,y,z;
  int i;
  float VecTempX,VecTempY,VecTempZ;
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=MultFactor/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,VecTempX,VecTempY,VecTempZ,i) 
  {
    //initiate InvDeformationField
    #pragma omp for
    for (y=0;y<VeloField->NY;y++){
      for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
        InvDeformationField->P(-VeloField->G(0,x,y,z)*coef,0,x,y,z);
        InvDeformationField->P(-VeloField->G(1,x,y,z)*coef,1,x,y,z);
        InvDeformationField->P(-VeloField->G(2,x,y,z)*coef,2,x,y,z);
      }
    }
    
    //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
    for (i=0;i<log2TimeStepNb;i++){
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          VecTempX=x+InvDeformationField->G(0,x,y,z);
          VecTempY=y+InvDeformationField->G(1,x,y,z);
          VecTempZ=z+InvDeformationField->G(2,x,y,z);
          
          TempVeloField1.P(InvDeformationField->G(0,x,y,z)+InvDeformationField->G(0,VecTempX,VecTempY,VecTempZ),0,x,y,z);
          TempVeloField1.P(InvDeformationField->G(1,x,y,z)+InvDeformationField->G(1,VecTempX,VecTempY,VecTempZ),1,x,y,z);
          TempVeloField1.P(InvDeformationField->G(2,x,y,z)+InvDeformationField->G(2,VecTempX,VecTempY,VecTempZ),2,x,y,z);
        }
      }
      
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          InvDeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
          InvDeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
          InvDeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else


///Integrate the inverse of a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptInvDefFromSteadyVeloField(VectorField * VeloField,VectorField * InvDeformationField,int log2TimeStepNb,float MultFactor){
  int x,y,z;
  int i;
  float VecTemp[3];
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=MultFactor/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  
  //initiate InvDeformationField
  for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
    InvDeformationField->P(-VeloField->G(0,x,y,z)*coef,0,x,y,z);
    InvDeformationField->P(-VeloField->G(1,x,y,z)*coef,1,x,y,z);
    InvDeformationField->P(-VeloField->G(2,x,y,z)*coef,2,x,y,z);
  }
  
  //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
  for (i=0;i<log2TimeStepNb;i++){
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      VecTemp[0]=x+InvDeformationField->G(0,x,y,z);
      VecTemp[1]=y+InvDeformationField->G(1,x,y,z);
      VecTemp[2]=z+InvDeformationField->G(2,x,y,z);
      
      TempVeloField1.P(InvDeformationField->G(0,x,y,z)+InvDeformationField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
      TempVeloField1.P(InvDeformationField->G(1,x,y,z)+InvDeformationField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
      TempVeloField1.P(InvDeformationField->G(2,x,y,z)+InvDeformationField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
    }
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      InvDeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
      InvDeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
      InvDeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
    }
  }
}

#endif


///compute the Lie Braket of VF1 and VF2. Put the result in VF3.   (subfunction of ComposeTwoLogFieldsUsingBCH)
void LieBracket(VectorField * VF1,VectorField * VF2,VectorField * VF3){
  int x,y,z;
  //put to 0 the boundaries of VF3
  x=0;
  for (z=0;z<VF1->NZ;z++) for (y=0;y<VF1->NY;y++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  x=VF1->NX-1;
  for (z=0;z<VF1->NZ;z++) for (y=0;y<VF1->NY;y++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  y=0;
  for (z=0;z<VF1->NZ;z++)  for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  y=VF1->NY-1;
  for (z=0;z<VF1->NZ;z++)  for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  z=0;
  for (y=0;y<VF1->NY;y++) for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  z=VF1->NZ-1;
  for (y=0;y<VF1->NY;y++) for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  
  //treat the rest of the data
  for (z=1;z<VF1->NZ-1;z++) for (y=1;y<VF1->NY-1;y++) for (x=1;x<VF1->NX-1;x++){
    VF3->P(-((VF2->G(0,x+1,y,z)-VF2->G(0,x-1,y,z))*VF1->G(0,x,y,z))-(VF2->G(0,x,y+1,z)-(VF2->G(0,x,y-1,z))*VF1->G(1,x,y,z))-((VF2->G(0,x,y,z+1)-VF2->G(0,x,y,z-1))*VF1->G(2,x,y,z)),0,x,y,z);
    VF3->P(-((VF2->G(1,x+1,y,z)-VF2->G(1,x-1,y,z))*VF1->G(0,x,y,z))-(VF2->G(1,x,y+1,z)-(VF2->G(1,x,y-1,z))*VF1->G(1,x,y,z))-((VF2->G(1,x,y,z+1)-VF2->G(1,x,y,z-1))*VF1->G(2,x,y,z)),1,x,y,z);
    VF3->P(-((VF2->G(2,x+1,y,z)-VF2->G(2,x-1,y,z))*VF1->G(0,x,y,z))-(VF2->G(2,x,y+1,z)-(VF2->G(2,x,y-1,z))*VF1->G(1,x,y,z))-((VF2->G(2,x,y,z+1)-VF2->G(2,x,y,z-1))*VF1->G(2,x,y,z)),2,x,y,z);
    
    VF3->Add(((VF1->G(0,x+1,y,z)-VF1->G(0,x-1,y,z))*VF2->G(0,x,y,z))+((VF1->G(0,x,y+1,z)-VF1->G(0,x,y-1,z))*VF2->G(1,x,y,z))+((VF1->G(0,x,y,z+1)-VF1->G(0,x,y,z-1))*VF2->G(2,x,y,z)),0,x,y,z);
    VF3->Add(((VF1->G(1,x+1,y,z)-VF1->G(1,x-1,y,z))*VF2->G(0,x,y,z))+((VF1->G(1,x,y+1,z)-VF1->G(1,x,y-1,z))*VF2->G(1,x,y,z))+((VF1->G(1,x,y,z+1)-VF1->G(1,x,y,z-1))*VF2->G(2,x,y,z)),1,x,y,z);
    VF3->Add(((VF1->G(2,x+1,y,z)-VF1->G(2,x-1,y,z))*VF2->G(0,x,y,z))+((VF1->G(2,x,y+1,z)-VF1->G(2,x,y-1,z))*VF2->G(1,x,y,z))+((VF1->G(2,x,y,z+1)-VF1->G(2,x,y,z-1))*VF2->G(2,x,y,z)),2,x,y,z);
    
    VF3->P(VF3->G(0,x,y,z)/2,0,x,y,z);
    VF3->P(VF3->G(1,x,y,z)/2,1,x,y,z);
    VF3->P(VF3->G(2,x,y,z)/2,2,x,y,z);
  }
}






///RefVeloField and UpdateVeloField are two steady velocity fields. Their exponentials are respectively the current deformation
///and the update of the deformation. This function approximates the velocity field which is the log of the composition 
///between the current deformation and the update deformation (Baker-Campbell-Hausdorff formula) (cf Vercauteren MICCAI 2008)
///In output, RefVeloField is the updated velocity field. UpdateVeloField is also modified for algorithmic reasons but it
///represents nothing pertinent as an output.
void ComposeTwoLogFieldsUsingBCH(VectorField * RefVeloField,VectorField * UpdateVeloField){
  int x,y,z;
  VectorField TempVeloField1;
  VectorField TempVeloField2;
  
  //init
  TempVeloField1.CreateVoidField(RefVeloField->NX,RefVeloField->NY,RefVeloField->NZ);
  TempVeloField2.CreateVoidField(RefVeloField->NX,RefVeloField->NY,RefVeloField->NZ);
  
  
  //compute LieBracket(RefVeloField,UpdateVeloField)
  LieBracket(RefVeloField,UpdateVeloField,&TempVeloField1);
  
  //compute LieBracket(RefVeloField,LieBracket(RefVeloField,UpdateVeloField))
  LieBracket(RefVeloField,&TempVeloField1,&TempVeloField2);
  
  
  
  // estimate the new RefVeloField
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(0,x,y,z)+(TempVeloField1.G(0,x,y,z)/2)+(TempVeloField2.G(0,x,y,z)/12),0,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(1,x,y,z)+(TempVeloField1.G(1,x,y,z)/2)+(TempVeloField2.G(1,x,y,z)/12),1,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(2,x,y,z)+(TempVeloField1.G(2,x,y,z)/2)+(TempVeloField2.G(2,x,y,z)/12),2,x,y,z);
}


///RefVeloField and UpdateVeloField are two steady velocity fields. Their exponentials are respectively the current deformation
///and the update of the deformation. This function approximates the velocity field which is the log of the composition 
///between the current deformation and the update deformation
///In output, RefVeloField is the updated velocity field. UpdateVeloField is also modified for algorithmic reasons but it
///represents nothing pertinent as an output.
void ComposeTwoLogFieldsUsingSum(VectorField * RefVeloField,VectorField * UpdateVeloField){
  int x,y,z;
  
  
  // estimate the new RefVeloField
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(0,x,y,z),0,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(1,x,y,z),1,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(2,x,y,z),2,x,y,z);
}



///Compose RefField with UpdateField. The result is saved in RefField
void DisplacementFieldCompose(VectorField * RefField,VectorField * UpdateField){
  int x,y,z;
  float VecTemp[3];
  
  for (z=0;z<RefField->NZ;z++) for (y=0;y<RefField->NY;y++) for (x=0;x<RefField->NX;x++){
    VecTemp[0]=x+RefField->G(0,x,y,z);
    VecTemp[1]=y+RefField->G(1,x,y,z);
    VecTemp[2]=z+RefField->G(2,x,y,z);
    
    RefField->Add(UpdateField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
    RefField->Add(UpdateField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
    RefField->Add(UpdateField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
  }
}


///Compose InvUpdateField with InvRefField. The result is saved in InvRefField (and InvRefField)
///(we consider that an inv. disp. f. points from the target source to the source)
void InvDisplacementFieldCompose(VectorField * InvUpdateField,VectorField * InvRefField){
  int x,y,z;
  float VecTemp[3];
  
  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++){
    VecTemp[0]=x+InvUpdateField->G(0,x,y,z);
    VecTemp[1]=y+InvUpdateField->G(1,x,y,z);
    VecTemp[2]=z+InvUpdateField->G(2,x,y,z);
    
    InvUpdateField->Add(InvRefField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
    InvUpdateField->Add(InvRefField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
    InvUpdateField->Add(InvRefField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
  }

  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++) InvRefField->P(InvUpdateField->G(0,x,y,z),0,x,y,z);
  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++) InvRefField->P(InvUpdateField->G(1,x,y,z),1,x,y,z);
  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++) InvRefField->P(InvUpdateField->G(2,x,y,z),2,x,y,z);

}




///We consider here that the vectors of DispField point from 'StaticImage' to 'DeformedImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
///remark: 'DefoField' has to be sufficiently smooth to allow an accurate inversion
void ProjectImageUsingDispField(VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage, int NearestNgbh){
  int x,y,z;
  float VecTemp[3];
  float VecTemp2[3];
  int NbIt;
  int i;
  
  NbIt=3;
  
  for (z=0;z<DispField->NZ;z++) for (y=0;y<DispField->NY;y++) for (x=0;x<DispField->NX;x++){
    VecTemp[0]=x-DispField->G(0,x,y,z);
    VecTemp[1]=y-DispField->G(1,x,y,z);
    VecTemp[2]=z-DispField->G(2,x,y,z);
    
    for (i=0; i<NbIt; i++) {
      VecTemp2[0]=x-DispField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[1]=y-DispField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[2]=z-DispField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]);
      
      VecTemp[0]=VecTemp2[0];
      VecTemp[1]=VecTemp2[1];
      VecTemp[2]=VecTemp2[2];
    }
    
    if (NearestNgbh!=1) DeformedImage->P(StaticImage->G(VecTemp[0],VecTemp[1],VecTemp[2]),x,y,z);
    else DeformedImage->P(StaticImage->G(static_cast<int>(VecTemp[0]+0.5),static_cast<int>(VecTemp[1]+0.5),static_cast<int>(VecTemp[2]+0.5)),x,y,z);
  }
}


///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndDispField(float ProjectCS_2_OriginCS[4][4],VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage){
  int x,y,z;
  float VecTemp[3];
  float VecTemp2[3];
  int NbIt;
  int i;
  
  NbIt=3;
  
  for (z=0;z<DispField->NZ;z++) for (y=0;y<DispField->NY;y++) for (x=0;x<DispField->NX;x++){
    VecTemp[0]=x-DispField->G(0,x,y,z);
    VecTemp[1]=y-DispField->G(1,x,y,z);
    VecTemp[2]=z-DispField->G(2,x,y,z);
    
    for (i=0; i<NbIt; i++) {
      VecTemp2[0]=x-DispField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[1]=y-DispField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[2]=z-DispField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]);
      
      VecTemp[0]=VecTemp2[0];
      VecTemp[1]=VecTemp2[1];
      VecTemp[2]=VecTemp2[2];
    }
    
    DeformedImage->P(StaticImage->G(ProjectCS_2_OriginCS,VecTemp[0],VecTemp[1],VecTemp[2]),x,y,z);
  }
}




///We consider here that the vectors of InvDispField point from 'DeformedImage' to 'StaticImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
void ProjectImageUsingInvDispField(VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh){
  int x,y,z;
  float VecTemp[3];
  
  
  //project the image
  for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
    VecTemp[0]=x+InvDispField->G(0,x,y,z);
    VecTemp[1]=y+InvDispField->G(1,x,y,z);
    VecTemp[2]=z+InvDispField->G(2,x,y,z);
    
    
    if (NearestNgbh!=1) DeformedImage->P(StaticImage->G(VecTemp[0],VecTemp[1],VecTemp[2]),x,y,z);
    else DeformedImage->P(StaticImage->G(static_cast<int>(VecTemp[0]+0.5),static_cast<int>(VecTemp[1]+0.5),static_cast<int>(VecTemp[2]+0.5)),x,y,z);
  }
  
}

#ifdef COMPILE_WITH_OPENMP

///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndInvDispField(float ProjectCS_2_OriginCS[4][4],VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  int x,y,z;
  float x2,y2,z2;
  int x3,y3,z3;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,x2,y2,z2) 
  {
    #pragma omp for
    for (y=0;y<DeformedImage->NY;y++){
      for (z=0;z<DeformedImage->NZ;z++)  for (x=0;x<DeformedImage->NX;x++){
        x2=x+InvDispField->G(0,x,y,z);
        y2=y+InvDispField->G(1,x,y,z);
        z2=z+InvDispField->G(2,x,y,z);
        
        DeformedImage->P(StaticImage->G(ProjectCS_2_OriginCS,x2,y2,z2,0,NN),x,y,z);
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndInvDispField(float ProjectCS_2_OriginCS[4][4],VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  int x,y,z;
  float x2,y2,z2;
  
  
  //project the image
  for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
    x2=x+InvDispField->G(0,x,y,z);
    y2=y+InvDispField->G(1,x,y,z);
    z2=z+InvDispField->G(2,x,y,z);
    
    DeformedImage->P(StaticImage->G(ProjectCS_2_OriginCS,x2,y2,z2,0,NN),x,y,z);
  }
}

#endif



///... explicit name ... here DispField represents the mapping from the target c.s to the source c.s. in voxels
void ProjectImageUsingDispFieldAndInvDispField(VectorField * DispField,VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  int x,y,z;
  float x2,y2,z2;
  int x3,y3,z3;
  float VecTemp[3];
  
  //project the image using ...
  if (NN==0){  //... trilinear interpolation
    for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
      x2=x+InvDispField->G(0,x,y,z);
      y2=y+InvDispField->G(1,x,y,z);
      z2=z+InvDispField->G(2,x,y,z);
      
      DeformedImage->P(StaticImage->G(DispField->G(0,x2,y2,z2),DispField->G(1,x2,y2,z2),DispField->G(2,x2,y2,z2)),x,y,z);
    }
  }
  else{ //... nearest neighbor interpolation
    for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
      x2=x+InvDispField->G(0,x,y,z);
      y2=y+InvDispField->G(1,x,y,z);
      z2=z+InvDispField->G(2,x,y,z);
      
      x3=static_cast<int>(DispField->G(0,x2,y2,z2)+0.5);
      y3=static_cast<int>(DispField->G(1,x2,y2,z2)+0.5);
      z3=static_cast<int>(DispField->G(2,x2,y2,z2)+0.5);
      
      if (x3<0) x3=0; if (x3>=StaticImage->NX) x3=StaticImage->NX-1; 
      if (y3<0) y3=0; if (y3>=StaticImage->NY) y3=StaticImage->NY-1; 
      if (z3<0) z3=0; if (z3>=StaticImage->NZ) z3=StaticImage->NZ-1; 
      
      DeformedImage->P(StaticImage->G(x3,y3,z3),x,y,z);
    }
  }
}






///... explicit name ....
void ProjectImageUsingSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage, int NearestNgbh,float factor){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptInvDefFromSteadyVeloField(VeloField,&TempField,5,factor);
  
  //project the image
  ProjectImageUsingInvDispField(&TempField,StaticImage,DeformedImage,NearestNgbh);
  
}


///... explicit name ....
void ProjectImageUsingAffineTransfoAndSteadyVeloField(float ProjectCS_2_OriginCS[4][4],VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN,float factor){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptInvDefFromSteadyVeloField(VeloField,&TempField,5,factor);
  
  //project the image
  ProjectImageUsingAffineTransfoAndInvDispField(ProjectCS_2_OriginCS,&TempField,StaticImage,DeformedImage,NN);
}

///... explicit name .... here DispField represents the mapping from the target c.s to the source c.s. in voxels
void ProjectImageUsingDispFieldAndSteadyVeloField(VectorField * DispField,VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN,float factor){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptInvDefFromSteadyVeloField(VeloField,&TempField,5,factor);
  
  //project the image
  ProjectImageUsingDispFieldAndInvDispField(DispField,&TempField,StaticImage,DeformedImage,NN);
  
}





///... explicit name ....
void ProjectImageUsingInverseSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage, int NearestNgbh,float factor){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptDefFromSteadyVeloField(VeloField,&TempField,5,factor);
  
  //project the image
  ProjectImageUsingInvDispField(&TempField,StaticImage,DeformedImage,NearestNgbh);

}


///Image transformation using the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfo(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,ScalarField * TransformedImage){
  int x,y,z,t;
  
  for (t = 0; t < TransformedImage->NT; t++) for (z = 0; z < TransformedImage->NZ; z++) for (y = 0; y < TransformedImage->NY; y++) for (x = 0; x < TransformedImage->NX; x++){
    TransformedImage->P(ImagToPropag->G(ProjectCS_2_OriginCS,x,y,z,t),x,y,z,t);
  }
}


//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'
void CptMappingFromVeloField(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps,float DeltaX){
  float VecTemp[3];
  float VecTemp2[3];
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  
  //0) INITIALISATION
  
  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //allocate memory in GradScalVF if not done
  if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ)||(NBT!=Map->NT))
    Map->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  
  //1) MAPPING AT THE REFERENCE TIME SUBDIVISION
  if ((refTimeStep<0)||(refTimeStep>NBT-1)) refTimeStep=0;
  
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    Map->P(MappingAtRefTimeStep->G(0,x,y,z),0,x,y,z,refTimeStep);
    Map->P(MappingAtRefTimeStep->G(1,x,y,z),1,x,y,z,refTimeStep);
    Map->P(MappingAtRefTimeStep->G(2,x,y,z),2,x,y,z,refTimeStep);
  }
  
  
  
  //2) FORWARD MAPPING for the time subdivisions > refTimeStep
  for (t=refTimeStep+1;t<NBT;t++){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        //find the original coordinates
        Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[1]=VeloField->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[2]=VeloField->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
  }
  
  //3) BACKWARD MAPPING for the time subdivisions < refTimeStep
  for (t=refTimeStep-1;t>=0;t--){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[1]=VeloField->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[2]=VeloField->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
  }
}


#ifdef COMPILE_WITH_OPENMP

//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'
void CptMappingFromVeloField_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps,float DeltaX){
  float VecTempX,VecTempY,VecTempZ;
  int NBX,NBY,NBZ,NBT;
  int x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  
  //0) INITIALISATION
  
  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  
  //allocate memory in GradScalVF if not done
  if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ)||(NBT!=Map->NT))
    Map->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,VecTempX,VecTempY,VecTempZ,t) 
  {
    //1) MAPPING AT THE REFERENCE TIME SUBDIVISION
    if ((refTimeStep<0)||(refTimeStep>NBT-1)) refTimeStep=0;
    
    #pragma omp for
    for (y=0;y<NBY;y++){
      for (z=0;z<NBZ;z++)  for (x=0;x<NBX;x++){
        Map->P(static_cast<float>(x),0,x,y,z,refTimeStep);
        Map->P(static_cast<float>(y),1,x,y,z,refTimeStep);
        Map->P(static_cast<float>(z),2,x,y,z,refTimeStep);
      }
    }
    
  
    //2) FORWARD MAPPING for the time subdivisions > refTimeStep
    for (t=refTimeStep+1;t<NBT;t++){
      #pragma omp for
      for (y=0;y<NBY;y++){
        for (z=0;z<NBZ;z++)  for (x=0;x<NBX;x++){
          VecTempX=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempY=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempZ=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
          
          //find the original coordinates
          Map->P(Map->G(0,x-VecTempX,y-VecTempY,z-VecTempZ,t-1),0,x,y,z,t);
          Map->P(Map->G(1,x-VecTempX,y-VecTempY,z-VecTempZ,t-1),1,x,y,z,t);
          Map->P(Map->G(2,x-VecTempX,y-VecTempY,z-VecTempZ,t-1),2,x,y,z,t);
        }
      }
    }
    
  
    //3) BACKWARD MAPPING for the time subdivisions < refTimeStep
    for (t=refTimeStep-1;t>=0;t--){
      #pragma omp for
      for (y=0;y<NBY;y++){
        for (z=0;z<NBZ;z++)  for (x=0;x<NBX;x++){
          VecTempX=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempY=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempZ=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
          
          Map->P(Map->G(0,x+VecTempX,y+VecTempY,z+VecTempZ,t+1),0,x,y,z,t);
          Map->P(Map->G(1,x+VecTempX,y+VecTempY,z+VecTempZ,t+1),1,x,y,z,t);
          Map->P(Map->G(2,x+VecTempX,y+VecTempY,z+VecTempZ,t+1),2,x,y,z,t);
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'
void CptMappingFromVeloField_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps,float DeltaX){
  float VecTemp[3];
  float VecTemp2[3];
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  
  //0) INITIALISATION
  
  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //allocate memory in GradScalVF if not done
  if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ)||(NBT!=Map->NT))
    Map->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  
  //1) MAPPING AT THE REFERENCE TIME SUBDIVISION
  if ((refTimeStep<0)||(refTimeStep>NBT-1)) refTimeStep=0;
  
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    Map->P(static_cast<float>(x),0,x,y,z,refTimeStep);
    Map->P(static_cast<float>(y),1,x,y,z,refTimeStep);
    Map->P(static_cast<float>(z),2,x,y,z,refTimeStep);
  }
  
  
  
  //2) FORWARD MAPPING for the time subdivisions > refTimeStep
  for (t=refTimeStep+1;t<NBT;t++){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        //find the original coordinates
        Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[1]=VeloField->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[2]=VeloField->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
  }
  
  //3) BACKWARD MAPPING for the time subdivisions < refTimeStep
  for (t=refTimeStep-1;t>=0;t--){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[1]=VeloField->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[2]=VeloField->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
  }
}

#endif


//Do the same thing as the initialisation of CptMappingFromVeloField -> load the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
//The function 'CptMappingFromVeloField2_Increment' then allows to increment the field backward or forward according to 'VeloField'
void CptMappingFromVeloField2_Init(VectorField * MappingAtRefTimeStep,VectorField * Map){
  int NBX,NBY,NBZ;
  int x,y,z;
  
  //0) INITIALISATION
  
  //initialisation
  NBX=MappingAtRefTimeStep->NX;
  NBY=MappingAtRefTimeStep->NY;
  NBZ=MappingAtRefTimeStep->NZ;
  
  //allocate memory in GradScalVF if not done
  if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ))
    Map->CreateVoidField(NBX,NBY,NBZ);
  
  
  //1) MAPPING AT THE REFERENCE TIME SUBDIVISION
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    Map->P(MappingAtRefTimeStep->G(0,x,y,z),0,x,y,z);
    Map->P(MappingAtRefTimeStep->G(1,x,y,z),1,x,y,z);
    Map->P(MappingAtRefTimeStep->G(2,x,y,z),2,x,y,z);
  }
  
}



//same as CptMappingFromVeloField2_Init with an identity mapping
void CptMappingFromVeloField2_Init_IniIdMap(VectorField * Map){
  int x,y,z;
  
  for (z=0;z<Map->NZ;z++) for (y=0;y<Map->NY;y++) for (x=0;x<Map->NX;x++){
    Map->P(static_cast<float>(x),0,x,y,z);
    Map->P(static_cast<float>(y),1,x,y,z);
    Map->P(static_cast<float>(z),2,x,y,z);
  }
}




//Do the same thing as an incrementation of CptMappingFromVeloField from 'CurrentTimeStep' and backward (BackwardOrForward==-1) or forward (BackwardOrForward==1) 
//The function 'CptMappingFromVeloField2_Init' is then supposed to have loaded the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
void CptMappingFromVeloField2_Increment(VectorField * VeloField,VectorField * Map,int CurrentTimeStep,int BackwardOrForward,int ConvergenceSteps,float DeltaX){
  float VecTemp[3];
  float VecTemp2[3];
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  VectorField NewMap;
  
  //1) INITIALISATION
  
  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  NewMap.CreateVoidField(NBX,NBY,NBZ);
  
  t=CurrentTimeStep;
  
  
  if (BackwardOrForward==1){  //2.1) FORWARD MAPPING for the time subdivisions > refTimeStep
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        //find the original coordinates
        NewMap.P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),0,x,y,z);
        NewMap.P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),1,x,y,z);
        NewMap.P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),2,x,y,z);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[1]=VeloField->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[2]=VeloField->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        NewMap.P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),0,x,y,z);
        NewMap.P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),1,x,y,z);
        NewMap.P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),2,x,y,z);
      }
    }
  }
  else{  //2.2) BACKWARD MAPPING for the time subdivisions < refTimeStep
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        NewMap.P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),0,x,y,z);
        NewMap.P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),1,x,y,z);
        NewMap.P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),2,x,y,z);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[1]=VeloField->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[2]=VeloField->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        NewMap.P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),0,x,y,z);
        NewMap.P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),1,x,y,z);
        NewMap.P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),2,x,y,z);
      }
    }
  }
  
  
  //3) Compute the new map
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    Map->P(NewMap.G(0,x,y,z),0,x,y,z);
    Map->P(NewMap.G(1,x,y,z),1,x,y,z);
    Map->P(NewMap.G(2,x,y,z),2,x,y,z);
  }
  
}





//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialMapping' which is the partial mapping of 'MappingAtRefTimeStep' from the time 
//subdivision 'refTimeStep' due to the contribution of 'PartialVeloField'. Note, that an Identity mapping 'MappingId' is //also defined in the inputs (to avoid defining it each time the function is used)
void CptPartialMappingFromVeloFields(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * MappingId,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  float x1,y1,z1;
  float x2,y2,z2;
  float x3,y3,z3;
  float x4,y4,z4;
  float xS,yS,zS;
  int t;
  float DeltaT_div_DeltaX;
  VectorField TotalBmap;  //total Backward mapping from the [new time_sub-1] to 0
  VectorField TotalBmap2;  //total Backward mapping from the [new time_sub] to 0
  VectorField TotalFmap;  //total forward mapping from the [new time_sub-1] to 0
  VectorField TotalFmap2;  //total forwrd mapping from the [new time_sub] to 0
  
  //1) INITIALISATION
  //1.1) constants
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in GradScalVF if not done
  if ((NBX!=PartialMapping->NX)||(NBY!=PartialMapping->NY)||(NBZ!=PartialMapping->NZ)||(NBT!=PartialMapping->NT))
    PartialMapping->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //1.3) allocate memory for TotalBmapand TotalFmap
  TotalBmap.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalBmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalFmap.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalFmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
  
  
  //1.4) PartialMapping at the reference time subdivision
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    PartialMapping->P(MappingAtRefTimeStep->G(0,x,y,z),0,x,y,z,refTimeStep);
    PartialMapping->P(MappingAtRefTimeStep->G(1,x,y,z),1,x,y,z,refTimeStep);
    PartialMapping->P(MappingAtRefTimeStep->G(2,x,y,z),2,x,y,z,refTimeStep);
  }
  
  //2) PARTIAL FORWARD MAPPING FOR t>refTimeStep
  for (t=refTimeStep+1;t<NBT;t++){
    
    //2.1) compute the total backward mapping from t-1 to 0
    CptMappingFromVeloField(t-1,MappingId,VeloField,&TotalBmap,ConvergenceSteps);
    
    //2.2) compute the total backward mapping from t to 0
    CptMappingFromVeloField(t,MappingId,VeloField,&TotalBmap2,ConvergenceSteps);
    
    //2.3) compute partial forward map at t from the one at t-1
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //2.3.1) first estimation
      //2.3.1.a) first guess of where the information comes from at t-1
      x1=PartialMapping->G(0,x,y,z,t-1); y1=PartialMapping->G(1,x,y,z,t-1); z1=PartialMapping->G(2,x,y,z,t-1);
      x2=TotalBmap.G(0,x1,y1,z1,0);   y2=TotalBmap.G(1,x1,y1,z1,0);   z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
      
      xS=x-PartialVeloField->G(0,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      yS=y-PartialVeloField->G(1,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      zS=z-PartialVeloField->G(2,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      
      //2.3.1.b) first transport of the information
      PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
      PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
      PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
      
      
      //2.3.1) leap frog style improvement of the estimation
      for (i=0;i<ConvergenceSteps*2;i++){
        //2.3.2.a) where the information comes from at t-1
        x1=PartialMapping->G(0,xS,yS,zS,t-1); y1=PartialMapping->G(1,xS,yS,zS,t-1); z1=PartialMapping->G(2,xS,yS,zS,t-1);
        x2=TotalBmap.G(0,x1,y1,z1,0);      y2=TotalBmap.G(1,x1,y1,z1,0);      z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
        
        x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
        x4=TotalBmap2.G(0,x3,y3,z3,0);  y4=TotalBmap2.G(1,x3,y3,z3,0);  z4=TotalBmap2.G(2,x3,y3,z3,0); //TotalBmap2 -> t
        
        xS=x-(PartialVeloField->G(0,x2,y2,z2,t-1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        yS=y-(PartialVeloField->G(1,x2,y2,z2,t-1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        zS=z-(PartialVeloField->G(2,x2,y2,z2,t-1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        
        //2.3.1.b) update the transport of the information
        PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
        PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
        PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
      }
    }
  }
  
  //3) PARTIAL BACWARD MAPPING FOR t<refTimeStep
  for (t=refTimeStep-1;t>=0;t--){      //not tested
    //3.1) compute the total forward mapping from 0 to t+1
    CptMappingFromVeloField(t+1,MappingId,VeloField,&TotalFmap,ConvergenceSteps);
    
    //3.2) compute the total forward mapping from 0 to t
    CptMappingFromVeloField(t,MappingId,VeloField,&TotalFmap2,ConvergenceSteps);
    
    //3.3) compute partial forward map at t from the one at t+1
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //3.3.1) first estimation
      //3.3.1.a) first guess of where the information comes from at t-1
      
      x1=PartialMapping->G(0,x,y,z,t+1); y1=PartialMapping->G(1,x,y,z,t+1); z1=PartialMapping->G(2,x,y,z,t+1);
      x2=TotalFmap.G(0,x1,y1,z1,NBT-1);   y2=TotalFmap.G(1,x1,y1,z1,NBT-1);   z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t+1
      
      xS=x+PartialVeloField->G(0,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      yS=y+PartialVeloField->G(1,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      zS=z+PartialVeloField->G(2,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      
      //3.3.1.b) first transport of the information
      PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
      PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
      PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
      
      //3.3.1) leap frog style improvement of the estimation
      for (i=0;i<ConvergenceSteps*2;i++){
        //3.3.2.a) where the information comes from at t-1
        x1=PartialMapping->G(0,xS,yS,zS,t+1); y1=PartialMapping->G(1,xS,yS,zS,t+1); z1=PartialMapping->G(2,xS,yS,zS,t+1);
        x2=TotalFmap.G(0,x1,y1,z1,NBT-1);      y2=TotalFmap.G(1,x1,y1,z1,NBT-1);      z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t-1
        
        x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
        x4=TotalFmap2.G(0,x3,y3,z3,NBT-1);  y4=TotalFmap2.G(1,x3,y3,z3,NBT-1);  z4=TotalFmap2.G(2,x3,y3,z3,NBT-1); //TotalFmap2 -> t
        
        xS=x+(PartialVeloField->G(0,x2,y2,z2,t+1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        yS=y+(PartialVeloField->G(1,x2,y2,z2,t+1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        zS=z+(PartialVeloField->G(2,x2,y2,z2,t+1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        
        //3.3.1.b) update the transport of the information
        PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
        PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
        PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
      }
    }
  }
}



//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialMapping' which is the partial mapping of 'MappingAtRefTimeStep' from the time 
//subdivision 'refTimeStep' due to the contribution of 'PartialVeloField'.
void CptPartialMappingFromVeloFields_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  float x1,y1,z1;
  float x2,y2,z2;
  float x3,y3,z3;
  float x4,y4,z4;
  float xS,yS,zS;
  int t;
  float DeltaT_div_DeltaX;
  VectorField TotalBmap;  //total Backward mapping from the [new time_sub-1] to 0
  VectorField TotalBmap2;  //total Backward mapping from the [new time_sub] to 0
  VectorField TotalFmap;  //total forward mapping from the [new time_sub-1] to 0
  VectorField TotalFmap2;  //total forwrd mapping from the [new time_sub] to 0
  
  //1) INITIALISATION
  //1.1) constants
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in GradScalVF if not done
  if ((NBX!=PartialMapping->NX)||(NBY!=PartialMapping->NY)||(NBZ!=PartialMapping->NZ)||(NBT!=PartialMapping->NT))
    PartialMapping->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //1.3) allocate memory for TotalBmapand TotalFmap
  TotalBmap.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalBmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalFmap.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalFmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
  
  
  //1.4) PartialMapping at the reference time subdivision
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    PartialMapping->P(static_cast<float>(x),0,x,y,z,refTimeStep);
    PartialMapping->P(static_cast<float>(y),1,x,y,z,refTimeStep);
    PartialMapping->P(static_cast<float>(z),2,x,y,z,refTimeStep);
  }
  
  //2) PARTIAL FORWARD MAPPING FOR t>refTimeStep
  for (t=refTimeStep+1;t<NBT;t++){
    
    //2.1) compute the total backward mapping from t-1 to 0
    CptMappingFromVeloField_IniIdMap(t-1,VeloField,&TotalBmap,ConvergenceSteps);
    
    //2.2) compute the total backward mapping from t to 0
    CptMappingFromVeloField_IniIdMap(t,VeloField,&TotalBmap2,ConvergenceSteps);
    
    //2.3) compute partial forward map at t from the one at t-1
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //2.3.1) first estimation
      //2.3.1.a) first guess of where the information comes from at t-1
      x1=PartialMapping->G(0,x,y,z,t-1); y1=PartialMapping->G(1,x,y,z,t-1); z1=PartialMapping->G(2,x,y,z,t-1);
      x2=TotalBmap.G(0,x1,y1,z1,0);   y2=TotalBmap.G(1,x1,y1,z1,0);   z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
      
      xS=x-PartialVeloField->G(0,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      yS=y-PartialVeloField->G(1,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      zS=z-PartialVeloField->G(2,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      
      //2.3.1.b) first transport of the information
      PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
      PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
      PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
      
      
      //2.3.1) leap frog style improvement of the estimation
      for (i=0;i<ConvergenceSteps*2;i++){
        //2.3.2.a) where the information comes from at t-1
        x1=PartialMapping->G(0,xS,yS,zS,t-1); y1=PartialMapping->G(1,xS,yS,zS,t-1); z1=PartialMapping->G(2,xS,yS,zS,t-1);
        x2=TotalBmap.G(0,x1,y1,z1,0);      y2=TotalBmap.G(1,x1,y1,z1,0);      z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
        
        x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
        x4=TotalBmap2.G(0,x3,y3,z3,0);  y4=TotalBmap2.G(1,x3,y3,z3,0);  z4=TotalBmap2.G(2,x3,y3,z3,0); //TotalBmap2 -> t
        
        xS=x-(PartialVeloField->G(0,x2,y2,z2,t-1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        yS=y-(PartialVeloField->G(1,x2,y2,z2,t-1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        zS=z-(PartialVeloField->G(2,x2,y2,z2,t-1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        
        //2.3.1.b) update the transport of the information
        PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
        PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
        PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
      }
    }
  }
  
  //3) PARTIAL BACWARD MAPPING FOR t<refTimeStep
  for (t=refTimeStep-1;t>=0;t--){      //not tested
    //3.1) compute the total forward mapping from 0 to t+1
    CptMappingFromVeloField_IniIdMap(t+1,VeloField,&TotalFmap,ConvergenceSteps);
    
    //3.2) compute the total forward mapping from 0 to t
    CptMappingFromVeloField_IniIdMap(t,VeloField,&TotalFmap2,ConvergenceSteps);
    
    //3.3) compute partial forward map at t from the one at t+1
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //3.3.1) first estimation
      //3.3.1.a) first guess of where the information comes from at t-1
      
      x1=PartialMapping->G(0,x,y,z,t+1); y1=PartialMapping->G(1,x,y,z,t+1); z1=PartialMapping->G(2,x,y,z,t+1);
      x2=TotalFmap.G(0,x1,y1,z1,NBT-1);   y2=TotalFmap.G(1,x1,y1,z1,NBT-1);   z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t+1
      
      xS=x+PartialVeloField->G(0,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      yS=y+PartialVeloField->G(1,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      zS=z+PartialVeloField->G(2,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      
      //3.3.1.b) first transport of the information
      PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
      PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
      PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
      
      //3.3.1) leap frog style improvement of the estimation
      for (i=0;i<ConvergenceSteps*2;i++){
        //3.3.2.a) where the information comes from at t-1
        x1=PartialMapping->G(0,xS,yS,zS,t+1); y1=PartialMapping->G(1,xS,yS,zS,t+1); z1=PartialMapping->G(2,xS,yS,zS,t+1);
        x2=TotalFmap.G(0,x1,y1,z1,NBT-1);      y2=TotalFmap.G(1,x1,y1,z1,NBT-1);      z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t-1
        
        x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
        x4=TotalFmap2.G(0,x3,y3,z3,NBT-1);  y4=TotalFmap2.G(1,x3,y3,z3,NBT-1);  z4=TotalFmap2.G(2,x3,y3,z3,NBT-1); //TotalFmap2 -> t
        
        xS=x+(PartialVeloField->G(0,x2,y2,z2,t+1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        yS=y+(PartialVeloField->G(1,x2,y2,z2,t+1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        zS=z+(PartialVeloField->G(2,x2,y2,z2,t+1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        
        //3.3.1.b) update the transport of the information
        PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
        PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
        PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
      }
    }
  }
}




//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialLocMap' which is the partial mapping ONLY AT 'TargetSubdiv' FROM 'SourceSubdiv' due to the contribution of PartialVeloField.
//-> PartialLocMap therefore represents where are the coordinates of the points of time subdivision 'SourceSubdiv' when transported on time subdivision 'TargetSubdiv'. Note that we consider here an identity mapping at 'SourceSubdiv'.
void ComputeLagrangianPartialMapping(int SourceSubdiv,int TargetSubdiv,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialLocMap,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float x1,y1,z1;
  float x2,y2,z2;
  float x3,y3,z3;
  float DeltaT_div_DeltaX;
  float DX,DY,DZ;
  float DX2,DY2,DZ2;
  
  //1) initialisation
  //1.1) constants
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in GradScalVF if not done
  if ((NBX!=PartialLocMap->NX)||(NBY!=PartialLocMap->NY)||(NBZ!=PartialLocMap->NZ)||(NBT!=1))
    PartialLocMap->CreateVoidField(NBX,NBY,NBZ,1);
  
  //2) Compute the transportation of the points of time SourceSubdiv to TargetSubdiv
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    
    //initial coordinates
    x1=x*1.;  y1=y*1.;  z1=z*1.;  //for the transportation in the complete velocity field
    x2=x*1.;  y2=y*1.;  z2=z*1.;  //for the transportation in the partial velocity field
    
    //transportation...
    //...forward
    if (SourceSubdiv<TargetSubdiv) for (t=SourceSubdiv;t<TargetSubdiv;t++){
      x3=x1; y3=y1; z3=z1;
      
      DX=VeloField->G(0,x3,y3,z3,t)*DeltaT_div_DeltaX;
      DY=VeloField->G(1,x3,y3,z3,t)*DeltaT_div_DeltaX;
      DZ=VeloField->G(2,x3,y3,z3,t)*DeltaT_div_DeltaX;
      
      DX2=(VeloField->G(0,x3,y3,z3,t)+VeloField->G(0,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
      DY2=(VeloField->G(1,x3,y3,z3,t)+VeloField->G(1,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
      DZ2=(VeloField->G(2,x3,y3,z3,t)+VeloField->G(2,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
      DX=DX2;
      DY=DY2;
      DZ=DZ2;
      
      x1=x1+DX;
      y1=y1+DY;
      z1=z1+DZ;
      
      x2=x2+(PartialVeloField->G(0,x3,y3,z3,t)+PartialVeloField->G(0,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
      y2=y2+(PartialVeloField->G(1,x3,y3,z3,t)+PartialVeloField->G(1,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
      z2=z2+(PartialVeloField->G(2,x3,y3,z3,t)+PartialVeloField->G(2,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
    }
    
    //...backward
    if (SourceSubdiv<TargetSubdiv) for (t=SourceSubdiv;t>TargetSubdiv;t--){ //not tested
      x3=x1; y3=y1; z3=z1;
      
      DX=VeloField->G(0,x3,y3,z3,t)*DeltaT_div_DeltaX;
      DY=VeloField->G(1,x3,y3,z3,t)*DeltaT_div_DeltaX;
      DZ=VeloField->G(2,x3,y3,z3,t)*DeltaT_div_DeltaX;
      
      DX2=(VeloField->G(0,x3,y3,z3,t)+VeloField->G(0,x3-DX,y3-DY,z3-DZ,t-1))*DeltaT_div_DeltaX/2.;
      DY2=(VeloField->G(1,x3,y3,z3,t)+VeloField->G(1,x3-DX,y3-DY,z3-DZ,t-1))*DeltaT_div_DeltaX/2.;
      DZ2=(VeloField->G(2,x3,y3,z3,t)+VeloField->G(2,x3-DX,y3-DY,z3-DZ,t-1))*DeltaT_div_DeltaX/2.;
      DX=DX2;
      DY=DY2;
      DZ=DZ2;
      
      x1=x1-DX;
      y1=y1-DY;
      z1=z1-DZ;
      
      x2=x2-(PartialVeloField->G(0,x3,y3,z3,t)+PartialVeloField->G(0,x1,y1,z1,t-1))*DeltaT_div_DeltaX/2.;
      y2=y2-(PartialVeloField->G(1,x3,y3,z3,t)+PartialVeloField->G(1,x1,y1,z1,t-1))*DeltaT_div_DeltaX/2.;
      z2=z2-(PartialVeloField->G(2,x3,y3,z3,t)+PartialVeloField->G(2,x1,y1,z1,t-1))*DeltaT_div_DeltaX/2.;
    }
    
    
    //save where the point x,y,z is transported
    PartialLocMap->P(x2,0,x,y,z);
    PartialLocMap->P(y2,1,x,y,z);
    PartialLocMap->P(z2,2,x,y,z);
  }
}



void EulerScheme(int refTimeStep,VectorField * VeloField,VectorField * Map,float DeltaX){
	float VecTemp[3];
	int NBX,NBY,NBZ,NBT;
	int x,y,z;
	int t;
	float DeltaT_div_DeltaX;
	
	//0) INITIALISATION
	
	//initialisation
	NBX=VeloField->NX;
	NBY=VeloField->NY;
	NBZ=VeloField->NZ;
	NBT=VeloField->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
	//allocate memory in GradScalVF if not done
	if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ)||(NBT!=Map->NT))
		Map->CreateVoidField(NBX,NBY,NBZ,NBT);
	
	
	//1) MAPPING AT THE REFERENCE TIME SUBDIVISION
	if ((refTimeStep<0)||(refTimeStep>NBT-1)) refTimeStep=0;
	
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
		Map->P(static_cast<float>(x),0,x,y,z,refTimeStep);
		Map->P(static_cast<float>(y),1,x,y,z,refTimeStep);
		Map->P(static_cast<float>(z),2,x,y,z,refTimeStep);
	}
	
	//2) Euler Scheme
	for (t=refTimeStep+1;t<NBT;t++)
    {
            // simple integration scheme (centered in time)
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
            {
				VecTemp[0]=VeloField->G(0,x,y,z,t-1)*DeltaT_div_DeltaX;
				VecTemp[1]=VeloField->G(1,x,y,z,t-1)*DeltaT_div_DeltaX;
				VecTemp[2]=VeloField->G(2,x,y,z,t-1)*DeltaT_div_DeltaX;
				
				//find the original coordinates
				Map->P(Map->G(0,x,y,z,t-1)+VecTemp[0],0,x,y,z,t);
				Map->P(Map->G(1,x,y,z,t-1)+VecTemp[1],1,x,y,z,t);
				Map->P(Map->G(2,x,y,z,t-1)+VecTemp[2],2,x,y,z,t);
            }
	}
}




#ifdef COMPILE_WITH_OPENMP

//Compute the projection of a 3D image 'ImagToPropag' using Mapping 'Map'.
//The image is projected at the time step 'TimeStepProj' of 'Map' and stored in 'ImageTimeT'.
//
//Importantly, the Mapping 'Map' should be an identity transformation at the time step 't' where 'ImagToPropag' is.
//It should also represent a forward mapping after 't' and a backward mapping before 't'.
void Project3Dimage(ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,t) 
  {
    #pragma omp for
    for (y = 0; y < ImageTimeT->NY; y++){
      for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++)  for (x = 0; x < ImageTimeT->NX; x++){
        ImageTimeT->P(ImagToPropag->G(Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
      }
    }
  //END FORK FOR THREADS
  }
}

#else

//Compute the projection of a 3D image 'ImagToPropag' using Mapping 'Map'.
//The image is projected at the time step 'TimeStepProj' of 'Map' and stored in 'ImageTimeT'.
//
//Importantly, the Mapping 'Map' should be an identity transformation at the time step 't' where 'ImagToPropag' is.
//It should also represent a forward mapping after 't' and a backward mapping before 't'.
void Project3Dimage(ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
  for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    ImageTimeT->P(ImagToPropag->G(Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
  }
}

#endif


#ifdef COMPILE_WITH_OPENMP

//same as above but the coordinates are transformed with the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfoAndTimeDepVF(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
 //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,t) 
  {
    #pragma omp for
    for (y = 0; y < ImageTimeT->NY; y++){
      for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++)  for (x = 0; x < ImageTimeT->NX; x++){
        ImageTimeT->P(ImagToPropag->G(ProjectCS_2_OriginCS,Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
      }
    }
  //END FORK FOR THREADS
  }
}

#else

//same as above but the coordinates are transformed with the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfoAndTimeDepVF(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
  for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    ImageTimeT->P(ImagToPropag->G(ProjectCS_2_OriginCS,Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
  }
}

#endif



//same as above but the coordinates are transformed with a displacement field in voxels (from 'Map' c. s. to 'ImagToPropag' c. s.)
void Project3DImageUsingDispFieldAndTimeDepVF(VectorField * DispField,ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  float x2,y2,z2;
  
  for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    x2=Map->G(0,x,y,z,TimeStepProj);
    y2=Map->G(1,x,y,z,TimeStepProj);
    z2=Map->G(2,x,y,z,TimeStepProj);
    
    ImageTimeT->P(ImagToPropag->G(DispField->G(0,x2,y2,z2),DispField->G(1,x2,y2,z2),DispField->G(2,x2,y2,z2),t),x,y,z,t);
  }
}




//same as above but with vector fields (NOT REORIENTED!!!)
void Project3Dimage(VectorField * VFToPropag,VectorField * Map,VectorField * VFTimeT,int TimeStepProj){
  int x,y,z,t,i;
  
  for (i = 0; i < 3; i++) for (t = 0; t < VFTimeT->NT; t++) for (z = 0; z < VFTimeT->NZ; z++) for (y = 0; y < VFTimeT->NY; y++) for (x = 0; x < VFTimeT->NX; x++){
    VFTimeT->P(VFToPropag->G(i,Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),i,x,y,z,t);
  }
}





///By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
///'VeloField4Measure' in the length of the flow from each point of the field. The length of flow
///is projected AT T=0 and returned in the 3D scalar field 'LengthOfFlow'
/// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
/// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
///   is computed.
void CptLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps,float DeltaX){
  float VecTemp[3];
  float VecTemp2[3];
  float NormVecTemp;
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  ScalarField PrevLengthOfFlow;
  
  //1) INITIALISATION
  //1.1) field size
  NBX=VeloField4Flow->NX;
  NBY=VeloField4Flow->NY;
  NBZ=VeloField4Flow->NZ;
  NBT=VeloField4Flow->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in LengthOfFlow at time t if not done
  if ((NBX!=LengthOfFlow->NX)||(NBY!=LengthOfFlow->NY)||(NBZ!=LengthOfFlow->NZ))
    LengthOfFlow->CreateVoidField(NBX,NBY,NBZ);
  
  //1.3) allocate memory of the PrevLengthOfFlow at time t+1
  PrevLengthOfFlow.CreateVoidField(NBX,NBY,NBZ);
  
  //1.4) PrevLengthOfFlow at the last time subdivision
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) PrevLengthOfFlow.P(0.,x,y,z);
  
  //JO at the other time subdivisions
  for (t=NBT-2;t>=0;t--){
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //init
      VecTemp[0]=0.; 
      VecTemp[1]=0.;
      VecTemp[2]=0.;
      
      //convergence
      for (i=0;i<ConvergenceSteps;i++){
        VecTemp2[0]=VeloField4Flow->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
        VecTemp2[1]=VeloField4Flow->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
        VecTemp2[2]=VeloField4Flow->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
        
        VecTemp[0]=(VecTemp2[0]+VeloField4Flow->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VecTemp2[1]+VeloField4Flow->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VecTemp2[2]+VeloField4Flow->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
      }
      
      //compute the lenght
      NormVecTemp =(float)pow((double)VeloField4Measure->G(0,x,y,z,t)*DeltaT_div_DeltaX,2.);
      NormVecTemp+=(float)pow((double)VeloField4Measure->G(1,x,y,z,t)*DeltaT_div_DeltaX,2.);
      NormVecTemp+=(float)pow((double)VeloField4Measure->G(2,x,y,z,t)*DeltaT_div_DeltaX,2.);
      NormVecTemp=sqrt(NormVecTemp);
      
      LengthOfFlow->P(NormVecTemp+PrevLengthOfFlow.G(x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),x,y,z);
    }
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
      PrevLengthOfFlow.P(LengthOfFlow->G(x,y,z),x,y,z);
  }
}











///By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
///'VeloField4Measure' AT THE CURRENT TIME in the length of the flow from each point of the field. The length of flow
///is returned in the 3D+t scalar field 'LengthOfFlow'
/// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
/// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
///   is computed.
void CptEvoLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps,float DeltaX){
  float VecTemp[3];
  float VecTemp2[3];
  float VecTemp3[3];
  float TmpFl;
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  
  //1) INITIALISATION
  //1.1) field size
  NBX=VeloField4Flow->NX;
  NBY=VeloField4Flow->NY;
  NBZ=VeloField4Flow->NZ;
  NBT=VeloField4Flow->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in LengthOfFlow at time t if not done
  if ((NBX!=LengthOfFlow->NX)||(NBY!=LengthOfFlow->NY)||(NBZ!=LengthOfFlow->NZ)||(NBT!=LengthOfFlow->NT))
    LengthOfFlow->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //1.3) JO at the first time subdivision
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
    LengthOfFlow->P(0.,x,y,z,0);
  
  //   ScalarField toto;
  //   toto.Read("AOD_Deformation1to2.nii");
  //   for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
  //         LengthOfFlow->P(toto.G(x,y,z,14),x,y,z,0);
  
  
  //2) Computation of the amplitude of the deformations
  for (t=1;t<NBT;t++){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField4Flow->G(0,x,y,z,t-1)+VeloField4Flow->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField4Flow->G(1,x,y,z,t-1)+VeloField4Flow->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField4Flow->G(2,x,y,z,t-1)+VeloField4Flow->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        //propagate the lenght of flow
        VecTemp3[0]=(VeloField4Measure->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp3[1]=(VeloField4Measure->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp3[2]=(VeloField4Measure->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        TmpFl =(float)pow(VecTemp3[0],2);
        TmpFl+=(float)pow(VecTemp3[1],2);
        TmpFl+=(float)pow(VecTemp3[2],2);
        TmpFl=sqrt(TmpFl);
        TmpFl+=LengthOfFlow->G(x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
        LengthOfFlow->P(TmpFl,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField4Flow->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[1]=VeloField4Flow->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[2]=VeloField4Flow->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField4Flow->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField4Flow->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField4Flow->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //propagate the lenght of flow
        VecTemp3[0]=(VeloField4Measure->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp3[1]=(VeloField4Measure->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp3[2]=(VeloField4Measure->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        
        TmpFl =(float)pow(VecTemp3[0],2);
        TmpFl+=(float)pow(VecTemp3[1],2);
        TmpFl+=(float)pow(VecTemp3[2],2);
        TmpFl=sqrt(TmpFl);
        TmpFl+=LengthOfFlow->G(x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
        LengthOfFlow->P(TmpFl,x,y,z,t);
      }
    }
  }
}


void SimpleEulerStep(VectorField *VectorField1,VectorField *VectorField2, int time, float dt){
  cout << "SimpleEulerStep <- TO DO" << endl;
}




///compute the L_2 norm of the difference between two scalar fields
float CalcSqrtSumOfSquaredDif(ScalarField * I1,ScalarField * I2){
  int x,y,z;
  float L2_norm,tmp;
  
  L2_norm=0.;
  
  for (z=0;z<I1->NZ;z++) for (y=0;y<I1->NY;y++) for (x=0;x<I1->NX;x++){
    tmp=(I1->G(x,y,z)-I2->G(x,y,z));
    L2_norm+=tmp*tmp;
  }
  
  L2_norm=sqrt(L2_norm);
  
  return L2_norm;
}




///...
void TransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,int t,float DeltaX)
{	int x,y,z;
	int NBX,NBY,NBZ;
	float d11,d12,d13,d21,d22,d23,d31,d32,d33;
	float temp;
	NBX=TempInvDiffeo->NX;
	NBY=TempInvDiffeo->NY;
	NBZ=TempInvDiffeo->NZ;
	for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
		d11=(TempInvDiffeo->G(0,x+1,y,z,t)-TempInvDiffeo->G(0,x-1,y,z,t))/(2.*DeltaX);
		d12=(TempInvDiffeo->G(0,x,y+1,z,t)-TempInvDiffeo->G(0,x,y-1,z,t))/(2.*DeltaX);
		d13=(TempInvDiffeo->G(0,x,y,z+1,t)-TempInvDiffeo->G(0,x,y,z-1,t))/(2.*DeltaX);
		d21=(TempInvDiffeo->G(1,x+1,y,z,t)-TempInvDiffeo->G(1,x-1,y,z,t))/(2.*DeltaX);
		d22=(TempInvDiffeo->G(1,x,y+1,z,t)-TempInvDiffeo->G(1,x,y-1,z,t))/(2.*DeltaX);
		d23=(TempInvDiffeo->G(1,x,y,z+1,t)-TempInvDiffeo->G(1,x,y,z-1,t))/(2.*DeltaX);
		d31=(TempInvDiffeo->G(2,x+1,y,z,t)-TempInvDiffeo->G(2,x-1,y,z,t))/(2.*DeltaX);
		d32=(TempInvDiffeo->G(2,x,y+1,z,t)-TempInvDiffeo->G(2,x,y-1,z,t))/(2.*DeltaX);
		d33=(TempInvDiffeo->G(2,x,y,z+1,t)-TempInvDiffeo->G(2,x,y,z-1,t))/(2.*DeltaX);
		temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0);
		Momentum->P( temp* (d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13)),x,y,z);
	}
	//1.2.2) boundaries at 0.
	z=0;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	z=NBZ-1;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	y=0;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	y=NBY-1;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	x=0;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	x=NBX-1;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	
	//1.2.3) 2D image case
	//float max=0.0;
	if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
		d11=(TempInvDiffeo->G(0,x+1,y,0,t)-TempInvDiffeo->G(0,x-1,y,0,t))/(2.*DeltaX);
		d12=(TempInvDiffeo->G(0,x,y+1,0,t)-TempInvDiffeo->G(0,x,y-1,0,t))/(2.*DeltaX);
		d21=(TempInvDiffeo->G(1,x+1,y,0,t)-TempInvDiffeo->G(1,x-1,y,0,t))/(2.*DeltaX);
		d22=(TempInvDiffeo->G(1,x,y+1,0,t)-TempInvDiffeo->G(1,x,y-1,0,t))/(2.*DeltaX);
		temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,0,t),TempInvDiffeo->G(1,x,y,0,t),TempInvDiffeo->G(2,x,y,0,t),0);
		//if (max<abs(temp*(d11*d22-d21*d12))){max=abs(temp*(d11*d22-d21*d12));}
		Momentum->P(temp*(d11*d22-d21*d12),x,y,0);
	}
	//cout << "c'est penible  "<< Momentum->GetMaxAbsVal() <<"\n";
}

/// Bug - Not VaLidated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void AddTransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,float cste, int t)
{  int x,y,z;
  int NBX,NBY,NBZ;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;
  float temp;
  NBX=TempInvDiffeo->NX;
  NBY=TempInvDiffeo->NY;
  NBZ=TempInvDiffeo->NZ;
  for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
  {
    d11=(TempInvDiffeo->G(0,x+1,y,z,t)-TempInvDiffeo->G(0,x-1,y,z,t))/(2.*DeltaX);
    d12=(TempInvDiffeo->G(0,x,y+1,z,t)-TempInvDiffeo->G(0,x,y-1,z,t))/(2.*DeltaX);
    d13=(TempInvDiffeo->G(0,x,y,z+1,t)-TempInvDiffeo->G(0,x,y,z-1,t))/(2.*DeltaX);
    d21=(TempInvDiffeo->G(1,x+1,y,z,t)-TempInvDiffeo->G(1,x-1,y,z,t))/(2.*DeltaX);
    d22=(TempInvDiffeo->G(1,x,y+1,z,t)-TempInvDiffeo->G(1,x,y-1,z,t))/(2.*DeltaX);
    d23=(TempInvDiffeo->G(1,x,y,z+1,t)-TempInvDiffeo->G(1,x,y,z-1,t))/(2.*DeltaX);
    d31=(TempInvDiffeo->G(2,x+1,y,z,t)-TempInvDiffeo->G(2,x-1,y,z,t))/(2.*DeltaX);
    d32=(TempInvDiffeo->G(2,x,y+1,z,t)-TempInvDiffeo->G(2,x,y-1,z,t))/(2.*DeltaX);
    d33=(TempInvDiffeo->G(2,x,y,z+1,t)-TempInvDiffeo->G(2,x,y,z-1,t))/(2.*DeltaX);
    temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0);
    Momentum->Add(cste* temp* (d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13)),x,y,z);
  }
  //1.2.2) boundaries at 0.
  z=0;
  for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
  z=NBZ-1;
  for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
  y=0;
  for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
  y=NBY-1;
  for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
  x=0;
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
  x=NBX-1;
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
  
  //1.2.3) 2D image case
  if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
  {
    d11=(TempInvDiffeo->G(0,x+1,y,0,t)-TempInvDiffeo->G(0,x-1,y,0,t))/(2.*DeltaX);
    d12=(TempInvDiffeo->G(0,x,y+1,0,t)-TempInvDiffeo->G(0,x,y-1,0,t))/(2.*DeltaX);
    d21=(TempInvDiffeo->G(1,x+1,y,0,t)-TempInvDiffeo->G(1,x-1,y,0,t))/(2.*DeltaX);
    d22=(TempInvDiffeo->G(1,x,y+1,0,t)-TempInvDiffeo->G(1,x,y-1,0,t))/(2.*DeltaX);
    temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,0,t),TempInvDiffeo->G(1,x,y,0,t));
    Momentum->Add(cste*temp*(d11*d22-d21*d12),x,y,0);
  }
}

///...
void TransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image, int t)
{
  int x,y,z;
  for (z = 0; z < InitialImage->NZ; z++) for (y = 0; y < InitialImage->NY; y++) for (x = 0; x < InitialImage->NX; x++)
  {
    Image->P(InitialImage->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0),x,y,z,0);
  }
}

///...
void AddTransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image, float cste, int t)
{
  int x,y,z;
  for (z = 0; z < InitialImage->NZ; z++) for (y = 0; y < InitialImage->NY; y++) for (x = 0; x < InitialImage->NX; x++)
  {
    Image->Add(cste * InitialImage->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0),x,y,z,0);
  }
}

///...
void DeepCopy(VectorField *VectorField1,VectorField *VectorField2,int t)
{
  int i,x,y,z;
  for (i=0;i<3;i++)
  {
    for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
    {
      VectorField2->P(VectorField1->G(i,x,y,z),i,x,y,z,t);
    }
  }
}
void DeepCopy(VectorField *VectorField,ScalarField* ScalarField,int direc,int t)
{
  int i,x,y,z;
  for (i=0;i<3;i++)
  {
    for (z = 0; z < ScalarField->NZ; z++) for (y = 0; y < ScalarField->NY; y++) for (x = 0; x < ScalarField->NX; x++)
    {
      ScalarField->P(VectorField->G(direc,x,y,z,t),x,y,z,0);
    }
  }
}


///...
void DeepCopy(ScalarField *ScalarField1,ScalarField *ScalarField2,int t)
{
  int x,y,z;
  for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
  {
    ScalarField2->P(ScalarField1->G(x,y,z),x,y,z,t);
  }  
}

///...
void DeepCopyMultiplyVectorField(VectorField *VectorField1,VectorField *VectorField2, float cste)
{
  int i,x,y,z,t;
  for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++) for (t=0; t<VectorField1->NT; t++)
  {
    VectorField2->P(cste*VectorField1->G(i,x,y,z,t),i,x,y,z,t);
  }  
}


/// Compute the scalar product between the vector fields and put it in ScalarField0 (for which NT=1)
void ScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t,float cste)
{
  int i,x,y,z;
  float temp;
  for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
  {
    temp=0.0;
    for (i=0;i<3;i++){temp+=VectorField1->G(i,x,y,z,t)*VectorField2->G(i,x,y,z,t);}
    ScalarField0->P(cste*temp,x,y,z);
  }
}


/// Compute the scalar product between the scalar fields and put it in ScalarField0 (for which NT=1)
void ScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t,float cste)
{
  int x,y,z;
  for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
  {
    ScalarField0->P(cste * ScalarField1->G(x,y,z,t)*ScalarField2->G(x,y,z,t),x,y,z);
  }  
}


/// Add the scalar product between the vector fields to ScalarField0 (for which NT=1)
void AddScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t)
{
  int i,x,y,z;
  for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
  {
    for (i=0;i<3;i++){ScalarField0->Add(VectorField1->G(i,x,y,z,t)*VectorField2->G(i,x,y,z,t),x,y,z);}
  }  
}


/// Add the scalar product between the scalar fields to ScalarField0 (for which NT=1)
void AddScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t)
{
  int x,y,z;
  for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
  {
    ScalarField0->Add(ScalarField1->G(x,y,z,t)*ScalarField2->G(x,y,z,t),x,y,z);
  }  
}


/// Add  ScalarField1 at time t1 to ScalarField2 at time t2
void AddScalarField(ScalarField *ScalarField1, ScalarField *ScalarField2,float cste, int t1,int t2)
{
  int x,y,z;
  for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
  {
    ScalarField2->Add(cste*ScalarField1->G(x,y,z,t1),x,y,z,t2);
  }  
}
/// Add  ScalarField1 at time t1 to ScalarField2 at time t2
void AddVectorField(VectorField *VectorField1, VectorField *VectorField2,float cste, int t1,int t2)
{
  int i,x,y,z;
  for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
  {
    VectorField2->Add(cste*VectorField1->G(i,x,y,z,t1),i,x,y,z,t2);
  }  
}
/// Multiply a vector field by the cste
void MultiplyVectorField(VectorField *VectorField1, float cste,int t)
{
  int i,x,y,z;
  for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
  {
    VectorField1->P(cste*VectorField1->G(i,x,y,z,t),i,x,y,z,t);
  }  
}

/// Sum two vector fields and put it in Output.
void SumVectorField(VectorField *VectorField1, VectorField *VectorField2, VectorField *Output, int t1, int t2, int t3, float cste1 ,float cste2)
{
  int i,x,y,z;
  for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
  {
    Output->P(cste1 * VectorField1->G(i,x,y,z,t1) + cste2 * VectorField2->G(i,x,y,z,t2),i,x,y,z,t3);
  }
}
/// compute the product element by element along each dimension
void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2)
{
  int i,x,y,z;
  for (i=0;i<3;i++) for (z = 0; z < ScalarField->NZ; z++) for (y = 0; y < ScalarField->NY; y++) for (x = 0; x < ScalarField->NX; x++)
  {
    VectorField2->P(VectorField1->G(i,x,y,z)*ScalarField->G(x,y,z),i,x,y,z);
  }
}

void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2, int time)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < ScalarField->NZ; z++) for (y = 0; y < ScalarField->NY; y++) for (x = 0; x < ScalarField->NX; x++)
	{
		VectorField2->P(VectorField1->G(i,x,y,z)*ScalarField->G(x,y,z),i,x,y,z,time);
	}
}

void AddProduct(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < ScalarField->NZ; z++) for (y = 0; y < ScalarField->NY; y++) for (x = 0; x < ScalarField->NX; x++)
	{
		VectorField2->Add(VectorField1->G(i,x,y,z)*ScalarField->G(x,y,z),i,x,y,z);
	}
}


///...
float DotProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, int t1,int t2)
{
  float result;
  int x,y,z;
  
  result=0.0;
  for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
    result += ScalarField1->G(x,y,z,t1) * ScalarField2->G(x,y,z,t2);
  
  return result;
}

float DotProduct(VectorField *ScalarField1, VectorField *ScalarField2, int t1,int t2)
{
	float result;
	int i,x,y,z;
  
	result=0.0;
  for (i=0;i<3;i++)	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
		result += ScalarField1->G(i,x,y,z,t1) * ScalarField2->G(i,x,y,z,t2);
	
  return result;
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                   7: OTHER FUNCTIONS OF SCIENTIFIC COMPUTATION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// 7.1) ++++++++++++++++++ small matrices ++++++++++++++++++ 


/// Solve the problem: MX=D where D is a known vector, M a tridiagonal matrix and X the unknown vector.
/// Inputs are a,b,c,d,n where M(i,i)=b(i), M(i,i-1)=a(i), M(i,i+1)=c(i), D(i)=d(i), D in R^n and M in R^n*R^n.
/// Output is X where X in R^n.  Warning: will modify c and d! */
void TridiagonalSolveFloat(const float *a, const float *b, float *c, float *d, float *x, int n){
  int i;
  double id;
  
  // Modify the coefficients.
  c[0] /= b[0];                       // Division by zero risk. 
  d[0] /= b[0];                       // Division by zero would imply a singular matrix. 
  for(i = 1; i < n; i++){
    id = (b[i] - c[i-1] * a[i]);      // Division by zero risk. 
    c[i] /= id;                       // Last value calculated is redundant. 
    d[i] = (d[i] - d[i-1] * a[i])/id;
  }
  
  // Now back substitute. 
  x[n - 1] = d[n - 1];
  for(i = n - 2; i >= 0; i--)
    x[i] = d[i] - c[i] * x[i + 1];
}


//Perform the eigenvalue decomposition of a 3*3 matrix
//Adapted from the algorithm having the same name in 'numerical recipes'.
//Input:  The 3*3 matrix 'MatIni' that has to be symetric.
//Ouputs: 'ValP' is a vector of size 3 which contains the eigen values (in decreasing order). 
//        'VecP' is a 3*3 matrix containg the eigen vectors in columns.
#define ROTATE3(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
void jacobi3(float **MatIni,float *ValP, float **VecP){
  int j,iq,ip,i; 
  float tresh,theta,tau,t,sm,s,h,g,c;
  float b[4];
  float z[4];
  float a[4][4];   //correspond a MatIni
  float d[4];    //correspond a ValP
  float v[4][4];   //correspond a VecP
  int vTri1,vTri2;
  float TempF;
  int n;
  
  
  n=3;
  for(i=0;i<n;i++) for(j=0;j<n;j++) a[i+1][j+1]=MatIni[i][j];
  
  
  //algo de numerical recipes
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++)
      for (iq=ip+1;iq<=n;iq++)
        sm += fabs(a[ip][iq]);
    
    if (sm == 0.0) {
      //adaptation des valeurs de l'algo de numerical recipes aux valeurs de sortie
      for(i=0;i<n;i++) ValP[i]=d[i+1];
      for(i=0;i<n;i++) for(j=0;j<n;j++) MatIni[i][j]=a[i+1][j+1];
      for(i=0;i<n;i++) for(j=0;j<n;j++) VecP[i][j]=v[i+1][j+1];
      
      //tri des donnees
      for(vTri1=0;vTri1<n-1;vTri1++) for(vTri2=vTri1+1;vTri2<n;vTri2++) if (ValP[vTri1]<ValP[vTri2]){
        TempF=ValP[vTri1]; ValP[vTri1]=ValP[vTri2]; ValP[vTri2]=TempF;
        for(i=0;i<n;i++) { TempF=VecP[i][vTri1]; VecP[i][vTri1]=VecP[i][vTri2]; VecP[i][vTri2]=TempF;}
      }
      
      return;
    }
    if (i < 4) tresh=0.2*sm/(n*n);
    else tresh=0.0;
    
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
        g=100.0*fabs(a[ip][iq]);
        
        if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((float)(fabs(h)+g) == (float)fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a[ip][iq];
          z[ip] -= h; z[iq] += h; d[ip] -= h; d[iq] += h; a[ip][iq]=0.0;
          for (j=1;j<=ip-1;j++) { ROTATE3(a,j,ip,j,iq) }
          for (j=ip+1;j<=iq-1;j++) { ROTATE3(a,ip,j,j,iq) }
          for (j=iq+1;j<=n;j++) { ROTATE3(a,ip,j,iq,j) }
          for (j=1;j<=n;j++) { ROTATE3(v,j,ip,j,iq)}
        }
      }
    }
    for (ip=1;ip<=n;ip++) { b[ip] += z[ip]; d[ip]=b[ip]; z[ip]=0.0; }
  }
  printf("Too many iterations in the routine jacobi\n");
  
  //adaptation des valeurs de l'algo de numerical recipes aux valeurs de sortie
  for(i=0;i<n;i++) ValP[i]=d[i+1];
  for(i=0;i<n;i++) for(j=0;j<n;j++) MatIni[i][j]=a[i+1][j+1];
  for(i=0;i<n;i++) for(j=0;j<n;j++) VecP[i][j]=v[i+1][j+1];
  
  //tri des donnees
  for(vTri1=0;vTri1<n-1;vTri1++) for(vTri2=vTri1+1;vTri2<n;vTri2++) if (ValP[vTri1]<ValP[vTri2]){
    TempF=ValP[vTri1]; ValP[vTri1]=ValP[vTri2]; ValP[vTri2]=TempF;
    for(i=0;i<n;i++) { TempF=VecP[i][vTri1]; VecP[i][vTri1]=VecP[i][vTri2]; VecP[i][vTri2]=TempF;}
  }
  
}



///compute two orthogonal vectors tvec1 and tvec2 in R^3 which are orthogonal to nvec
///the norm of tvec1 and tvec2 is defined as equal to the one of nvec
void CptVecsTangentPlane(float nvec[3],float tvec1[3],float tvec2[3]){
  float epsilon,dist;
  
  dist=sqrt(nvec[0]*nvec[0]+nvec[1]*nvec[1]+nvec[2]*nvec[2]);
  epsilon=dist/100;
  
  //define two orthogonal directions in the plan where transformations are allowed...
  //... vec1
  if (fabs(nvec[0])<epsilon){
    tvec1[0]=1; tvec1[1]=0; tvec1[2]=0;
  }
  else{
    tvec1[0]=0;
    
    if (fabs(nvec[1])<epsilon){
      tvec1[1]=1; tvec1[2]=0; 
    }
    else{
      tvec1[1]=-nvec[2]/nvec[1];
      tvec1[2]=1;
    }
  }
  
  VecNormalize(tvec1,dist);
  
  //... vec2
  tvec2[0]=(nvec[1]*tvec1[2])-(nvec[2]*tvec1[1]);
  tvec2[1]=-(nvec[0]*tvec1[2])+(nvec[2]*tvec1[0]);
  tvec2[2]=(nvec[0]*tvec1[1])-(nvec[1]*tvec1[0]);
  
  VecNormalize(tvec2,dist);
}



///normalize a vector
void VecNormalize(float vec[3],float norm){
  float tmpfl;
        
  tmpfl=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  vec[0]*=norm/tmpfl;
  vec[1]*=norm/tmpfl;
  vec[2]*=norm/tmpfl;
}


///compute the determinant of a 3*3 matrix
float determinant_3t3matrix(float m[3][3]){
  float deter;
  deter=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]);
  deter-=m[1][0]*(m[0][1]*m[2][2]-m[2][1]*m[0][2]);
  deter+=m[2][0]*(m[0][1]*m[1][2]-m[1][1]*m[0][2]);
  return deter;
}


///compute the comatrix of a 3*3 matrix
void comatrix_3t3matrix(float m1[3][3],float m2[3][3]){
  m2[0][0]=(m1[1][1]*m1[2][2])-(m1[2][1]*m1[1][2]);
  m2[1][0]=(m1[2][1]*m1[0][2])-(m1[0][1]*m1[2][2]);
  m2[2][0]=(m1[0][1]*m1[1][2])-(m1[1][1]*m1[0][2]);
  m2[0][1]=(m1[2][0]*m1[1][2])-(m1[1][0]*m1[2][2]);
  m2[1][1]=(m1[0][0]*m1[2][2])-(m1[2][0]*m1[0][2]);
  m2[2][1]=(m1[1][0]*m1[0][2])-(m1[0][0]*m1[1][2]);
  m2[0][2]=(m1[1][0]*m1[2][1])-(m1[2][0]*m1[1][1]);
  m2[1][2]=(m1[2][0]*m1[0][1])-(m1[0][0]*m1[2][1]);
  m2[2][2]=(m1[0][0]*m1[1][1])-(m1[1][0]*m1[0][1]);
  
}

///Estimate the exponential of a 3*3 matrix   (Checked with matlab and OK)
void Exponential_3t3matrix(float m1[3][3],float m2[3][3]){
  float tempMat[3][3];
  float tempMat2[3][3];
  int k,facto;
  
  //init  (k=0)
  m2[0][0]=1; m2[0][1]=0; m2[0][2]=0; 
  m2[1][0]=0; m2[1][1]=1; m2[1][2]=0; 
  m2[2][0]=0; m2[2][1]=0; m2[2][2]=1; 
  
  tempMat[0][0]=m1[0][0]; tempMat[0][1]=m1[0][1]; tempMat[0][2]=m1[0][2];
  tempMat[1][0]=m1[1][0]; tempMat[1][1]=m1[1][1]; tempMat[1][2]=m1[1][2];
  tempMat[2][0]=m1[2][0]; tempMat[2][1]=m1[2][1]; tempMat[2][2]=m1[2][2];
  
  //estimation...
  facto=1;
  for (k=1;k<15;k++){
    //...update m2
    m2[0][0]+=tempMat[0][0]/facto; m2[0][1]+=tempMat[0][1]/facto; m2[0][2]+=tempMat[0][2]/facto;
    m2[1][0]+=tempMat[1][0]/facto; m2[1][1]+=tempMat[1][1]/facto; m2[1][2]+=tempMat[1][2]/facto;
    m2[2][0]+=tempMat[2][0]/facto; m2[2][1]+=tempMat[2][1]/facto; m2[2][2]+=tempMat[2][2]/facto;
    
    //...update facto
    facto*=k+1;
    
    //...update tempMat
    mult_3t3mat_3t3mat(tempMat,m1,tempMat2);
    
    tempMat[0][0]=tempMat2[0][0]; tempMat[0][1]=tempMat2[0][1]; tempMat[0][2]=tempMat2[0][2];
    tempMat[1][0]=tempMat2[1][0]; tempMat[1][1]=tempMat2[1][1]; tempMat[1][2]=tempMat2[1][2];
    tempMat[2][0]=tempMat2[2][0]; tempMat[2][1]=tempMat2[2][1]; tempMat[2][2]=tempMat2[2][2];
  }
}


///transpose a 3*3 matrix
void transpose_3t3matrix(float m1[3][3],float m2[3][3]){
  int i,j;
  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      m2[i][j]=m1[j][i];
}


///multiply two 3*3 matrices
void  mult_3t3mat_3t3mat(float m1[3][3], float m2[3][3], float MatRes[3][3]){
  int i,j,k;
  
  for (i=0;i<3;i++)
    for(j=0;j<3;j++){
      MatRes[i][j]=0;
      for(k=0;k<3;k++) MatRes[i][j]+=m1[i][k]*m2[k][j];
    }
}



///inverse of a 3*3 matrix
void invert_3t3matrix(float m1[3][3],float m2[3][3]){
  float det;
  int i,j;
  
  det=determinant_3t3matrix(m1);
  comatrix_3t3matrix(m1,m2);
  transpose_3t3matrix(m2,m1);
  
  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      m2[i][j]=m1[i][j]/det;
  
}

///multiply a vector of size 3 and a 3*3 matrix
void mult_3t3mat_3vec(float mat[3][3], float vectIni[3], float vectRes[3]){
  int lig,col;
  
  for (lig=0;lig<3;lig++){
    vectRes[lig]=0;
    for(col=0;col<3;col++)
      vectRes[lig]+=mat[lig][col]*vectIni[col];
  }
}


///In a 3D domain, project a point "Pt2proj" to a line defined by a point "LinePt" and a vector "LineVec". Result is saved in "ProjPt"
float Project_point_to_a_line(float Pt2proj[3],float LinePt[3],float LineVec[3],float ProjPt[3]){
  float scalarProd; 
  
  //make sure LineVec is normalized
  VecNormalize(LineVec,1);
  
  //scalar product computation 
  scalarProd=((Pt2proj[0]-LinePt[0])*LineVec[0])+((Pt2proj[1]-LinePt[1])*LineVec[1])+((Pt2proj[2]-LinePt[2])*LineVec[2]);
  
  //estimate the projected point
  ProjPt[0]=LinePt[0]+scalarProd*LineVec[0];
  ProjPt[1]=LinePt[1]+scalarProd*LineVec[1];
  ProjPt[2]=LinePt[2]+scalarProd*LineVec[2];
  
  return 0;
}




///scalar product between two vectors of size 3
float scalarProd_3vec(float v1[3],float v2[3]){
  float scalarProd;
  scalarProd=(v1[0]*v2[0])+(v1[1]*v2[1])+(v1[2]*v2[2]);
  return scalarProd;
}

///cross product between two vectors of size 3
void crossProd_3vec(float v1[3],float v2[3],float OutV[3]){
  OutV[0]=(v1[1]*v2[2])-(v1[2]*v2[1]);
  OutV[1]=-(v1[0]*v2[2])+(v1[2]*v2[0]);
  OutV[2]=(v1[0]*v2[1])-(v1[1]*v2[0]);
}

///project a vector 'vec' to a plan defined by the vectors (Pvec1,Pvec2)
void Project_3vec_plan(float Pvec1[3],float Pvec2[3],float vec[3]){
  float ps1,ps2;
  
  ps1=scalarProd_3vec(vec,Pvec1);
  ps2=scalarProd_3vec(vec,Pvec2);
  
  vec[0]=ps1*Pvec1[0]+ps2*Pvec2[0];
  vec[1]=ps1*Pvec1[1]+ps2*Pvec2[1];
  vec[2]=ps1*Pvec1[2]+ps2*Pvec2[2];
}

///inverse of a quaternion
void invert_4t4quaternion(float q1[4][4],float q2[4][4]){
  float r11,r12,r13,r21,r22,r23,r31,r32,r33,v1,v2,v3 , deti ;
  
  //algorithm inspired from the one of nifti_io.c
  
  //  INPUT MATRIX IS:  
  r11 = q1[0][0]; r12 = q1[0][1]; r13 = q1[0][2];  // [ r11 r12 r13 v1 ] 
  r21 = q1[1][0]; r22 = q1[1][1]; r23 = q1[1][2];  // [ r21 r22 r23 v2 ] 
  r31 = q1[2][0]; r32 = q1[2][1]; r33 = q1[2][2];  // [ r31 r32 r33 v3 ] 
  v1  = q1[0][3]; v2  = q1[1][3]; v3  = q1[2][3];  // [  0   0   0   1 ] 
  
  deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
  +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;
  
  if( deti != 0.0 ) deti = 1.0 / deti ;
  
  q2[0][0] = deti*( r22*r33-r32*r23) ;
  q2[0][1] = deti*(-r12*r33+r32*r13) ;
  q2[0][2] = deti*( r12*r23-r22*r13) ;
  q2[0][3] = deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3
                    -r22*v1*r33-r32*r13*v2+r32*v1*r23) ;
  
  q2[1][0] = deti*(-r21*r33+r31*r23) ;
  q2[1][1] = deti*( r11*r33-r31*r13) ;
  q2[1][2] = deti*(-r11*r23+r21*r13) ;
  q2[1][3] = deti*( r11*r23*v3-r11*v2*r33-r21*r13*v3
                    +r21*v1*r33+r31*r13*v2-r31*v1*r23) ;
  
  q2[2][0] = deti*( r21*r32-r31*r22) ;
  q2[2][1] = deti*(-r11*r32+r31*r12) ;
  q2[2][2] = deti*( r11*r22-r21*r12) ;
  q2[2][3] = deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3
                    -r21*r32*v1-r31*r12*v2+r31*r22*v1) ;
  
  q2[3][0] = q2[3][1] = q2[3][2] = 0.0 ;
  q2[3][3] = (deti == 0.0) ? 0.0 : 1.0 ; // failure flag if deti == 0 
  
}


///multiply a vector of size 4 and a 4*4 matrix
void mult_4t4mat_4vec(float mat[4][4], float vectIni[4], float vectRes[4]){
  int lig,col;
  
  for (lig=0;lig<4;lig++){
    vectRes[lig]=0;
    for(col=0;col<4;col++)
      vectRes[lig]+=mat[lig][col]*vectIni[col];
  }
}

///multiply two 4*4 matrix: mat_i1 * mat_i2 -> mat_o
void mult_quat4t4mat_quat4t4mat(float mat_i1[4][4], float mat_i2[4][4], float mat_o[4][4]){
  int o1,o2,i;
  
  for (o1=0;o1<4;o1++) for (o2=0;o2<4;o2++){
    mat_o[o1][o2]=0;
    for (i=0;i<4;i++) mat_o[o1][o2]+=mat_i1[o1][i]*mat_i2[i][o2];
  }
}


///read a 4*4 matrix in a text file
void Read_quat4t4mat(char * FileName,float locmat[4][4]){
  int i,j;
  FILE *DataFile;
  
  //open file
  DataFile=fopen(FileName,"r");
  
  //go to the file beginning
  fseek(DataFile,0,SEEK_SET);
  
  //read the 4*4 matrix
  fscanf(DataFile,"%f  %f  %f  %f",&locmat[0][0],&locmat[0][1],&locmat[0][2],&locmat[0][3]);
  fscanf(DataFile,"%f  %f  %f  %f",&locmat[1][0],&locmat[1][1],&locmat[1][2],&locmat[1][3]);
  fscanf(DataFile,"%f  %f  %f  %f",&locmat[2][0],&locmat[2][1],&locmat[2][2],&locmat[2][3]);
  fscanf(DataFile,"%f  %f  %f  %f",&locmat[3][0],&locmat[3][1],&locmat[3][2],&locmat[3][3]);
  
  //close file
  fclose(DataFile);
}

///write a 4*4 matrix in a text file
void Write_quat4t4mat(char * FileName,float locmat[4][4]){
  int i,j;
  FILE *DataFile;
  
  //open file
  DataFile=fopen(FileName,"w");
  
  //read the 4*4 matrix
  fprintf(DataFile,"%f  %f  %f  %f\n",locmat[0][0],locmat[0][1],locmat[0][2],locmat[0][3]);
  fprintf(DataFile,"%f  %f  %f  %f\n",locmat[1][0],locmat[1][1],locmat[1][2],locmat[1][3]);
  fprintf(DataFile,"%f  %f  %f  %f\n",locmat[2][0],locmat[2][1],locmat[2][2],locmat[2][3]);
  fprintf(DataFile,"%f  %f  %f  %f\n",locmat[3][0],locmat[3][1],locmat[3][2],locmat[3][3]);
  
  //close file
  fclose(DataFile);
}

/// 7.2) ++++++++++++++++++ Square matrices ++++++++++++++++++ 


///constructors
SquareMatrix::SquareMatrix(void){
  this->NI=0;
  this->NJ=0;
  this->CopiedByRef=0;
}


///destructor
SquareMatrix::~SquareMatrix(void){
  if (this->CopiedByRef==0) 
    if ((this->LocMatrix!=NULL)&&(this->NI>0)){
      //cout << "delete a SquareMatrix" << endl;
      delete this->LocMatrix;
      }
  
  this->NI=0;
  this->NJ=0;
}


///read the matrix in CSV files with space delimiters
void SquareMatrix::Read(char * FileName){
  FILE *DataFile;
  float LocValue;
  int NbPts;
  
  //open file
  DataFile=fopen(FileName,"r");
  
  //count the number of floats in the file
  fseek(DataFile,0,SEEK_SET);
  NbPts=0;
  while(!feof(DataFile)){
    fscanf(DataFile,"%f ",&LocValue);
    NbPts++;
  }
  
  //compute the matrix size and allocate space for the data
  this->NI=static_cast<int>(floor(sqrt(static_cast<float>(NbPts))));
  this->NJ=static_cast<int>(floor(sqrt(static_cast<float>(NbPts))));
  cout << FileName << " is read as a " << this->NI << "*" << this->NJ << " matrix" << endl;
  
  this->LocMatrix = new float [NbPts];
  
  //read the data
  fseek(DataFile,0,SEEK_SET);
  NbPts=0;
  while(!feof(DataFile)){
    fscanf(DataFile,"%f ",&LocValue);
    LocMatrix[NbPts]=LocValue;
    NbPts++;
  }
  
  //close DataFile
  fclose(DataFile);
}

///allocate memory for a n*n matrix
void SquareMatrix::CreateVoidMatrix(int n){
  
  this->NI=n;
  this->NJ=n;
  this->CopiedByRef=0;
  
  this->LocMatrix = new float [n*n];
  
}


///put all values to zero
void SquareMatrix::SetToZero(){
  int i;
  for (i=0;i<this->NI*this->NJ;i++) this->LocMatrix[i]=0;
  }

///make an eye matrix
void SquareMatrix::Eye(){
  int i;
  
  this->SetToZero();
  
  for (i=0;i<this->NI;i++) this->P(1,i,i);
  }


///write the matrix in CSV files with space delimiters
void SquareMatrix::Write(char * FileName){
  FILE *DataFile;
  float LocValue;
  int NbPts;
  int i,j;
  
  //open file
  DataFile=fopen(FileName,"w");
  
  //write the data
  for (i=0;i<this->NI;i++){
    for (j=0;j<this->NJ;j++){
      fprintf(DataFile,"%.3e ",this->G(i,j));
      //fprintf(DataFile,"%5f ",this->G(i,j));
    }
    fprintf(DataFile,"\n");
  }

  //close DataFile
  fclose(DataFile);
}


///write the matrix in CSV files with space delimiters
void SquareMatrix::Show(){
  int i,j;
  
  cout << endl << this->NI << "*" << this->NJ << " matrix:" << endl;
  
  for (i=0;i<this->NI;i++){
    for (j=0;j<this->NJ;j++){
      cout <<  fixed << setprecision( 2 ) << this->G(i,j) << " ";
    }
    cout << endl;
  }
  cout << endl;
}



///copy a matrix by reference
void SquareMatrix::CopyByRef(SquareMatrix * copiedMatrix){
    
    if (this->NI>0) this->~SquareMatrix();
    
    this->NI=copiedMatrix->NI;
    this->NJ=copiedMatrix->NJ;
    this->CopiedByRef=1;
    this->LocMatrix=copiedMatrix->LocMatrix;
    //cout << "Matrix copied by reference" << endl;
}


//copy another square matrix.
void SquareMatrix::Copy(SquareMatrix * copiedMatrix){
  int i;
  //make sure that memory is well allocated to copy the matrix
  if (this->CopiedByRef==1) this->~SquareMatrix();
  
  if (this->NI==0) this->CreateVoidMatrix(copiedMatrix->NI);
  
  if (this->NI!=copiedMatrix->NI){
      this->~SquareMatrix();
      this->CreateVoidMatrix(copiedMatrix->NI);
    }
    
  //copy matrix
  for (i=0;i<this->NI*this->NJ;i++) this->LocMatrix[i]=copiedMatrix->LocMatrix[i];
  }


///transpose a matrix
void SquareMatrix::Transpose(){
  int i,j;
  float tmpFl;
  
  for (i=0;i<this->NI;i++) for (j=i+1;j<this->NJ;j++){
    tmpFl=this->G(i,j);
    this->P(this->G(j,i),i,j);
    this->P(tmpFl,j,i);
  }
}

///symmetrize a matrix -> Mat <- (Mat + Mat')/2
void SquareMatrix::symmetrize(){
  int i,j;
  float tmpFl;
  
  for (i=0;i<this->NI;i++) for (j=i+1;j<this->NJ;j++){
    tmpFl=(this->G(i,j)+this->G(j,i))/2;
    this->P(tmpFl,i,j);
    this->P(tmpFl,j,i);
  }
}

  


#ifdef COMPILE_WITH_OPENMP

///Matrix/Matrix multiplication: (this M2)
///The result is stored in this
void SquareMatrix::Mult(SquareMatrix * M2){
  float **tmpVec;
  int i,j,k;
  int MaxThreadNB,tid;

  //get the maximum number of threads allowed here
  MaxThreadNB=omp_get_max_threads();
  
  //allocate memory for each temporary row of the new matrix
  tmpVec = new float* [MaxThreadNB];
  for (i=0;i<MaxThreadNB;i++) 
    tmpVec[i] = new float [this->NJ];
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j,k,tid) num_threads(MaxThreadNB)
  {  
    //matrix multiplication
    #pragma omp for
    for (i=0;i<this->NI;i++){
      
      tid = omp_get_thread_num();
      //cout << "Row " << i << " treated in thread " << tid << endl;
      
      for (j=0;j<M2->NJ;j++){
        //init row treatment
        tmpVec[tid][j]=0;
        
        //define multiplied row
        for (k=0;k<this->NJ;k++) tmpVec[tid][j]+=this->G(i,k)*M2->G(k,j);
      }
      
      //copy modified row to row i in this
      for (j=0;j<this->NJ;j++) this->P(tmpVec[tid][j],i,j);
    }
    
  //END FORK FOR THREADS
  }
  
  //delete allocated memory
  for (i=0;i<MaxThreadNB;i++) delete tmpVec[i];
  delete tmpVec;
}




///Put the outer-product of two vectors in the square matrix
///The size of the vectors must be coherent with the matrix
void SquareMatrix::PutOuterProduct(float * v1,float * v2){
  int i,j;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j)
  {
    #pragma omp for 
    for (i=0;i<this->NI;i++){ 
      for (j=0;j<this->NJ;j++) this->P(v1[i]*v2[j],i,j);
    }
  //END FORK FOR THREADS
  }
}

///Add the outer-product of two vectors in the square matrix
///The size of the vectors must be coherent with the matrix
void SquareMatrix::AddOuterProduct(float * v1,float * v2){
  int i,j;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j)
  {
    #pragma omp for 
    for (i=0;i<this->NI;i++){
      for (j=0;j<this->NJ;j++) this->PlusEq(i,j,v1[i]*v2[j]);
    }
  //END FORK FOR THREADS
  }
}



///Addition with another square matrix
///The result is stored in this
void SquareMatrix::Add(SquareMatrix * M2){
  int i,j;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j)
  {
    #pragma omp for 
    for (i=0;i<this->NI;i++){
      for (j=0;j<this->NJ;j++) this->PlusEq(i,j,M2->G(i,j));
    }
  //END FORK FOR THREADS
  }
}
 
///Multiplication by a scalar
void SquareMatrix::Mult(float coefMult){
  int i,j;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j)
  {
    #pragma omp for 
    for (i=0;i<this->NI;i++){
      for (j=0;j<this->NJ;j++) this->TimesEq(i,j,coefMult);
    }
  //END FORK FOR THREADS
  }
}

///weighted sum of two matrices: this <- (weight_this * this) + (weight_M2 * M2)
void SquareMatrix::WeightedSum(float weight_this,SquareMatrix * M2,float weight_M2){
  int i,j;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j)
  {
    #pragma omp for 
    for (i=0;i<this->NI;i++){
      for (j=0;j<this->NJ;j++) this->LocMatrix[j+this->NI*i]=  (weight_this*this->LocMatrix[j+this->NI*i])  +  (weight_M2*M2->LocMatrix[j+this->NI*i]);
    }
  //END FORK FOR THREADS
  }
}


///Addition with a scalar
void SquareMatrix::Add(float b){
  int i,j;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j)
  {
    #pragma omp for 
    for (i=0;i<this->NI;i++){
      for (j=0;j<this->NJ;j++) this->PlusEq(i,j,b);
    }
  //END FORK FOR THREADS
  }
}



///Matrix/Vector multiplication. Return the product in V
void SquareMatrix::Mult(float * V){
  float *tmpVec;
  int i,k;

  tmpVec = new float [this->NJ];
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,k)
  {  
    #pragma omp for
    for (i=0;i<this->NI;i++){
      tmpVec[i]=0;
      for (k=0;k<this->NI;k++) tmpVec[i]+=this->G(i,k)*V[k];
    }
  
  //END FORK FOR THREADS
  }
  
  for (i=0;i<this->NI;i++) V[i]=tmpVec[i];
  
  delete tmpVec;
}



#else



///Matrix/Matrix multiplication: (this M2)
///The result is stored in this
void SquareMatrix::Mult(SquareMatrix * M2){
  float *tmpVec;
  int i,j,k;
  
  tmpVec = new float [this->NJ];
  
  for (i=0;i<this->NI;i++){
    for (j=0;j<M2->NJ;j++){
      //init row treatment
      tmpVec[j]=0;
      
      //define multiplied row
      for (k=0;k<this->NJ;k++) tmpVec[j]+=this->G(i,k)*M2->G(k,j);
    }
    
    //copy modified row to row i in this
    for (j=0;j<this->NJ;j++) this->P(tmpVec[j],i,j);
  }
  
  delete tmpVec;
}

///Put the outer-product of two vectors in the square matrix
///The size of the vectors must be coherent with the matrix
void SquareMatrix::PutOuterProduct(float * v1,float * v2){
  int i,j;
  for (i=0;i<this->NI;i++) for (j=0;j<this->NJ;j++) this->P(v1[i]*v2[j],i,j);
  }

///Add the outer-product of two vectors in the square matrix
///The size of the vectors must be coherent with the matrix
void SquareMatrix::AddOuterProduct(float * v1,float * v2){
  int i,j;
  for (i=0;i<this->NI;i++) for (j=0;j<this->NJ;j++) this->PlusEq(i,j,v1[i]*v2[j]);
  }



///Addition with another square matrix
///The result is stored in this
void SquareMatrix::Add(SquareMatrix * M2){
  int i,j;
  
  for (i=0;i<this->NI;i++) for (j=0;j<this->NJ;j++) this->PlusEq(i,j,M2->G(i,j));
  }
 
///Multiplication by a scalar
void SquareMatrix::Mult(float coefMult){
  int i,j;
  
  for (i=0;i<this->NI;i++) for (j=0;j<this->NJ;j++) this->TimesEq(i,j,coefMult);
}

///weighted sum of two matrices: this <- (weight_this * this) + (weight_M2 * M2)
void SquareMatrix::WeightedSum(float weight_this,SquareMatrix * M2,float weight_M2){
  int i,j;
  
  for (i=0;i<this->NI;i++) for (j=0;j<this->NJ;j++) 
    this->LocMatrix[j+this->NI*i]=  (weight_this*this->LocMatrix[j+this->NI*i])  +  (weight_M2*M2->LocMatrix[j+this->NI*i]);
}



///Addition with a scalar
void SquareMatrix::Add(float b){
  int i,j;
  
  for (i=0;i<this->NI;i++) for (j=0;j<this->NJ;j++) this->PlusEq(i,j,b);
}

///Matrix/Vector multiplication. Return the product in V
void SquareMatrix::Mult(float * V){
  float *tmpVec;
  int i,k;
  
  tmpVec = new float [this->NJ];
  
  for (i=0;i<this->NI;i++){
    tmpVec[i]=0;
    for (k=0;k<this->NI;k++) tmpVec[i]+=this->G(i,k)*V[k];
  }
  
  for (i=0;i<this->NI;i++) V[i]=tmpVec[i];
  
  delete tmpVec;
  }


#endif








///weighted sum of this with identity matrix: this <- (weight_this * this) + (weight_M2 * Id)
void SquareMatrix::WeightedSumWithId(float weight_this,float weight_M2){
  int i,j;
  
  if (weight_this!=static_cast<float>(1)) this->Mult(weight_this);
  
  for (i=0;i<this->NI;i++) 
    this->Add(weight_M2,i,i);
}




///sum and max of the values or absolute values
float SquareMatrix::SumValues(){
  float tmpFl;   int i;
  tmpFl=0; for (i=0;i<this->NI*this->NJ;i++) tmpFl+=this->LocMatrix[i];
  return tmpFl;
  }
  
float SquareMatrix::MaxValue(){
  float tmpFl;  int i;
  tmpFl=this->LocMatrix[0];  for (i=1;i<this->NI*this->NJ;i++) if (tmpFl<this->LocMatrix[i]) tmpFl=this->LocMatrix[i];
  return tmpFl;
  }
  
float SquareMatrix::SumAbsValues(){
  float tmpFl;  int i;
  tmpFl=0; for (i=0;i<this->NI*this->NJ;i++) tmpFl+=fabs(this->LocMatrix[i]);
  return tmpFl;
  }
  
float SquareMatrix::MaxAbsValue(){
  float tmpFl;  int i;
  tmpFl=fabs(this->LocMatrix[0]);  for (i=1;i<this->NI*this->NJ;i++) if (tmpFl<fabs(this->LocMatrix[i])) tmpFl=fabs(this->LocMatrix[i]);
  return tmpFl;
  }
  
  
///normalize the matrix so that the average sum of the columns is 1
void SquareMatrix::MatrixNormalization(){
    this->Mult(static_cast<float>(this->NJ)/this->SumAbsValues());
    }
  
  
///Perform a Householder reduction of a real symmetric matrix z[][]
///Outputs:
///  -> z[][] is replaced by the orthogonal matrix effecting the transformation. 
///  -> d[] returns the diagonal elements of the tri-diagonal matrix
///  -> e[] the off-diagonal elements with e[0] = 0.
///
///-> Adapted from www.physics.sdsu.edu/~johnson/phys580/tred2.c which was very slightly adapted from tred2() in numerical recipies in C++ 3rd edition (differences where checked)
///-> Validated on small matrices - results differ from those obtained using the hess function in matlab but make sense when reconstructing the original matrix
void tred2(SquareMatrix * z, float *d, float *e){
  int l,k,j,i;
  float scale,hh,h,g,f;
  int n;
  
  n=z->NI;
  
  cout << "Householder transformation" << endl;
  //cout << "Householder reduction - STEP 1" << endl;
  for(i=n-1;i>0;i--){
    l = i-1;
    h = 0.0;
    scale = 0.0;
    if(l > 0){
      for(k=0;k<=l;k++) scale += fabs(z->G(i,k));
      if(scale == 0.0)               //Skip transformation
        e[i] = z->G(i,l);
      else{
        for(k=0;k<=l;k++){
          z->DivEq(i,k,scale);           //Use scaled z's for transformation
          h+=z->G(i,k)*z->G(i,k);      //Form \sigma in h
        }
        f=z->G(i,l);
        g=(f>=0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h-=f*g;                    //Now h is equation (11.3.4)
        z->P(f-g,i,l);               //Store u in row i of z
        f=0.0;
        
        for(j=0;j<=l;j++){
          z->P(z->G(i,j)/h,j,i);        //Store u/H in column i of z
          g=0.0;                    //Form an element A \dot u in g
          for(k=0;k<=j;k++) g+=z->G(j,k)*z->G(i,k);
          for(k=j+1;k<=l;k++) g+=z->G(k,j)*z->G(i,k);
          e[j]=g/h;                 //Form element of p in temporarily unused element of e
          f+=e[j]*z->G(i,j);
        }
        hh=f/(h+h);                //Form K, equation (11.3.11)
        for(j=0;j<=l;j++){          //Form q and store e overwritting p
          f=z->G(i,j);
          e[j]=g=e[j]-hh*f;
          for(k=0;k<=j;k++)         //Reduce z, equation (11.3.13)
            z->MinusEq(j,k,f*e[k]+g*z->G(i,k));
        }
      }
    }
    else{
      e[i]=z->G(i,l);
    }
    d[i]=h;
  }  // end i-loop
  
  d[0]=0.0;
  e[0]=0.0;
  
  //cout << "Householder reduction - STEP 2" << endl;
  for(i=0;i<n;i++){       //Begin accumulation of transformation matrices
    l=i-1;
    if(d[i]){         //This block is skipped when d[i]==0
      for(j=0;j<=l;j++) {
        g= 0.0;
        for(k=0;k<=l;k++) g+=z->G(i,k)*z->G(k,j);  //Use u and u/H stored in z to form P \dot Q
        for (k=0;k<=l;k++) z->MinusEq(k,j,g*z->G(k,i));
      }
    }
    d[i]=z->G(i,i);
    z->P(1.0,i,i);            //reset row and column of z to identity matrix for next iteration
    for(j=0;j<=l;j++) { z->P(0.0,j,i);  z->P(0.0,i,j); }
  }
  //cout << "Householder reduction - done" << endl;
}



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
void tqli(float *d, float *e, SquareMatrix * z){
  int m,l,iter,i,k;
  float s,r,p,g,f,dd,c,b;
  int n;
  
  cout << "QL  algorithm" << endl;
  
  n=z->NI;
  
  //convenient to renumber the elements of e
  for(i=1;i<n;i++) e[i-1] = e[i];
  e[n]=0.0;
   
  //main loop
  //cout << "QL  algorithm with implicit shifts" << endl;
  for(l=0;l<n;l++){
    iter=0;
    do{
      //look for a single small subdiagonal element to split the matrix
      for(m=l;m<n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if((double)(fabs(e[m])+dd) == dd) break;  //Slightly different to NR3-C++
      }
      if(m!=l){
        //check whether the maximal number of iterations is reached
        //if(iter++==30){
        //  cout << "\n\nToo many iterations in tqli.\n" ;
        //  exit(1);
        //}
        //form shift
        g = (d[l+1] - d[l])/(2.0 * e[l]);
        r = pythag(g,1.0);
        //This is d_m-k_s
        g = d[m]-d[l]+e[l]/(g+nr_SIGN(r,g));
        s = c = 1.0;
        p = 0.0;
        //A plane rotation as in the original QL, folled by Givens rotations to restore tridiagonal form
        for (i=m-1;i>=l;i--){
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          //recover from underflow
          if(r==0.0){
            d[i+1]-=p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
          //Form eigenvectors
          for(k=0;k<n;k++){
            f=z->G(k,i+1);
            z->P(s*z->G(k,i)+c*f,k,i+1);
            z->P(c*z->G(k,i)-s*f,k,i);
          }
        }
        if(r == 0.0 && i >= l) continue;
        d[l]-=p;
        e[l]=g;
        e[m]=0.0;
      }
    } while(m != l);
  }
  
  //cout << "QL  algorithm with implicit shifts - done" << endl;
}

//rank eigenvalues and eigenvetors with decreasing eigenvalues 
void rankEigenvalues(float *eigenvalues, SquareMatrix * eigenvectors){
  float topEigenValue;
  int topEigenValueID;
  float tmpFl;
  int k,k2,i,n;
  int GoSwap;
  
  n=eigenvectors->NI;
  
  for(k=0;k<n-1;k++){
    //get the location of the highest eigenvalue after k (->topEigenValueID)
    topEigenValue=fabs(eigenvalues[k]);
    topEigenValueID=k;
    GoSwap=0;
    
    for(k2=k+1;k2<n;k2++) if (topEigenValue<fabs(eigenvalues[k2])){
      topEigenValue=fabs(eigenvalues[k2]);
      topEigenValueID=k2;
      GoSwap=1;
    }
    
    if (GoSwap==1){
      //swap the eigenvalues   k <-> topEigenValueID
      tmpFl=eigenvalues[k];
      eigenvalues[k]=eigenvalues[topEigenValueID];
      eigenvalues[topEigenValueID]=tmpFl;
      
      //swap corresponding eigenvectors
      for(i=0;i<n;i++){
        tmpFl=eigenvectors->G(i,k);
        eigenvectors->P(eigenvectors->G(i,topEigenValueID),i,k);
        eigenvectors->P(tmpFl,i,topEigenValueID);
      }
    }
  }
}

//show recomposed matrix after a PCA
void ShowRecomposedMatrixAfterPCA(float *eigenvalues, SquareMatrix * eigenvectors){
  int i,j,k,n;
  float tmpFl;
  
  n=eigenvectors->NI;
  
  cout << "Recomposed matrix:" << endl;
  for (i=0;i<n;i++){
     for (j=0;j<n;j++){
       tmpFl=0;
       for (k=0;k<n;k++) tmpFl+=eigenvectors->G(i,k)*eigenvectors->G(j,k)*eigenvalues[k];
       
       cout << tmpFl << " ";
       }
       cout << endl;
     }
}


//Construct the QR decomposition of a. Outputs are the orthogonal matrix Q and the upper triangluar matrix R.
//1 is returned if a singularity is encountered during the decomposition, but the decomposition is still completed in this case; 0 is returned otherwise
//Adapted from the algorithm QRdcmp of numerical recipes in c++ - third edition
//Works very well on small examples
int QRdecomp(SquareMatrix * a,SquareMatrix * q, SquareMatrix * r){
  int i,j,k;
  float scale, sigma, sum, tau;
  int sing;
  int n;
  float * d;
  float * c;
  
  //1) init
  sing=0;
  n=a->NI;
  q->SetToZero();
  r->Copy(a);
  
  d = new float [a->NI];
  c = new float [a->NI];
  
  //2) QR algorithm   -  note that what is called q here is actually q transpose
  for (k=0;k<n-1;k++){
    scale=0.0;
    for (i=k;i<n;i++) scale = nr_MAX(scale,fabs(r->G(i,k)));
    //singular case
    if (scale == 0.0){
      sing=1;
      c[k]=d[k]=0.0;
    }
    else{ //form Q_k and Q_k \dot A
      for (i=k;i<n;i++) r->DivEq(i,k,scale);
      for (sum=0.0,i=k;i<n;i++) sum+=r->G(i,k)*r->G(i,k);
      sigma=nr_SIGN(sqrt(sum),r->G(k,k));
      r->PlusEq(k,k,sigma);
      c[k]=sigma*r->G(k,k);
      d[k] = -scale*sigma;
      for (j=k+1;j<n;j++){
        for (sum=0.0,i=k;i<n;i++) sum+= r->G(i,k)*r->G(i,j);
        tau=sum/c[k];
        for (i=k;i<n;i++) r->MinusEq(i,j,tau*r->G(i,k));
      }
    }
  }
  d[n-1]=r->G(n-1,n-1);
  if (d[n-1]==0.0) sing=1;
  for (i=0;i<n;i++){                    //Form Q' explicitly
    for (j=0;j<n;j++) q->P(0.0,i,j);
    q->P(1.0,i,i);
  }
  for (k=0;k<n-1;k++){
    if (c[k]!=0.0){
      for (j=0;j<n;j++){
        sum=0.0;
        for(i=k;i<n;i++)
          sum+= r->G(i,k)*q->G(i,j);
        sum/=c[k];
        for (i=k;i<n;i++)
          q->MinusEq(i,j,sum*r->G(i,k));
      }
    }
  }
  for (i=0;i<n;i++){                    //Form R explicitly
    r->P(d[i],i,i);
    for (j=0;j<i;j++) r->P(0.0,i,j);
  }
  
  //3) transpose q' to output q
  q->Transpose();
  
  return sing;
}


//perform the PCA of a reasonably large SYMMETRIC matrix
void LargeMatrixPCA(SquareMatrix * eigenvectors, float *eigenvalues){
  float *e;
  e = new float [eigenvectors->NI];
  
  tred2(eigenvectors,eigenvalues,e);
  tqli(eigenvalues,e,eigenvectors);
  rankEigenvalues(eigenvalues,eigenvectors);
  }


/// 7.3) ++++++++++++++++++ Other scientific computation stuffs ++++++++++++++++++


///minmod function
float minmodfunc(float a,float b){
    float sa,sb;
    float minabsab;
    
    sa=1-2*(static_cast<float>(a<0));
    sb=1-2*(static_cast<float>(b<0));
    
    minabsab=fabs(a);
    if (fabs(b)<minabsab) minabsab=fabs(b);
    
    return minabsab*(sa+sb)/2;
}

///computes sqrt(a^2+b^2) without destructive underflow or overflow
///From numerical recipies (2nd edition)
float pythag(float a, float b){
  float absa,absb;
  
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

///Generate a random number following a normal distribution
float sampleNormal(){
  double u,v,r,c;
  u = 2.0 * static_cast<double>(rand()) / (RAND_MAX) - 1;
  v = 2.0 * static_cast<double>(rand()) / (RAND_MAX) - 1;
  r = u * u + v * v;
  if (r <0.01 || r > 0.99) return sampleNormal();
  c = sqrt(-2 * log(r) / r);
  return static_cast<float>(u * c);
  }

///Generate a random number following an uniform distribution  (in [0,1])
float sampleUniform(){
  return static_cast<float>(rand())/static_cast<float>(RAND_MAX);
  }

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                            8: LANDMARKS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// 8.1) ++++++++++++++++++ point landmarks ++++++++++++++++++ 

///constructor
LDMK_Points::LDMK_Points(void){
  this->LDMK_Points_Nb=0;
}

///destructor
LDMK_Points::~LDMK_Points(void){
  if ((this->Lx!=NULL)&&(this->LDMK_Points_Nb>0)) delete this->Lx;
  if ((this->Ly!=NULL)&&(this->LDMK_Points_Nb>0)) delete this->Ly;
  if ((this->Lz!=NULL)&&(this->LDMK_Points_Nb>0)) delete this->Lz;
  if ((this->Val!=NULL)&&(this->LDMK_Points_Nb>0)) delete this->Val;
  this->LDMK_Points_Nb=0;
}


///return the X coordinate of a Landmark
float LDMK_Points::GetX(int Id){
  
  if (Id<0) Id=0;
  if (Id>=this->LDMK_Points_Nb) Id=this->LDMK_Points_Nb-1;
  
  if (this->LDMK_Points_Nb==0) return 0;
  else return Lx[Id];
}  

///return the Y coordinate of a Landmark
float LDMK_Points::GetY(int Id){
  
  if (Id<0) Id=0;
  if (Id>=this->LDMK_Points_Nb) Id=this->LDMK_Points_Nb-1;
  
  if (this->LDMK_Points_Nb==0) return 0;
  else return Ly[Id];
}  

///return the Z coordinate of a Landmark
float LDMK_Points::GetZ(int Id){
  
  if (Id<0) Id=0;
  if (Id>=this->LDMK_Points_Nb) Id=this->LDMK_Points_Nb-1;
  
  if (this->LDMK_Points_Nb==0) return 0;
  else return Lz[Id];
}  

///return the value of a Landmark
float LDMK_Points::GetVal(int Id){
  
  if (Id<0) Id=0;
  if (Id>=this->LDMK_Points_Nb) Id=this->LDMK_Points_Nb-1;
  
  if (this->LDMK_Points_Nb==0) return 0;
  else return this->Val[Id];
}  


///Get the number of LDMK_Points
int LDMK_Points::Get_LDMK_PointsNumber(void){
  return this->LDMK_Points_Nb;
}  



///read LDMK_Points in a CSV file
///If (withValues!=1) the file should have this format:
///[Point 1: x]\t[Point 1: y]\t[Point 1:  z]
///[Point 2: x]\t[Point 2: y]\t[Point 2:  z]
///...
///or, if (withValues==1):
///[Point 1: x]\t[Point 1: y]\t[Point 1:  z]\t[Point 1:  Val]
///[Point 2: x]\t[Point 2: y]\t[Point 2:  z]\t[Point 2:  Val]
///...
void LDMK_Points::Read(char * FileName,int withValues){
  FILE *DataFile;
  float LocValue;
  int Nb_LDMK_Points;
  float fl1,fl2,fl3,fl4;
  
  //open file
  DataFile=fopen(FileName,"r");
  
  //1) Check the number of points and allocate memory
  //1.1) go to the file beginning
  fseek(DataFile,0,SEEK_SET);
  
  //1.2) read the file
  Nb_LDMK_Points=0;
  if (withValues!=1){  //1.2.1) no values
    while(!feof(DataFile)){
      fscanf(DataFile,"%f  %f  %f",&fl1,&fl2,&fl3);
      Nb_LDMK_Points++;
    }
  }
  else{  //1.2.2) with values
    while(!feof(DataFile)){
      fscanf(DataFile,"%f  %f  %f  %f",&fl1,&fl2,&fl3,&fl4);
      Nb_LDMK_Points++;
    }
  }
  
  this->LDMK_Points_Nb=Nb_LDMK_Points-1;
  
  //1.3) allocate the memory for the LDMK_Points
  this->Lx = new float [this->LDMK_Points_Nb];
  this->Ly = new float [this->LDMK_Points_Nb];
  this->Lz = new float [this->LDMK_Points_Nb];
  this->Val = new float [this->LDMK_Points_Nb];

  
  
  //2) fill the class
  //2.1) go to the file beginning
  fseek(DataFile,0,SEEK_SET);
  
  //2.2) read the file
  Nb_LDMK_Points=0;
  if (withValues!=1){  //2.2.1) no values
    while(!feof(DataFile)){
      fscanf(DataFile,"%f  %f  %f",&this->Lx[Nb_LDMK_Points],&this->Ly[Nb_LDMK_Points],&this->Lz[Nb_LDMK_Points]);
      this->Val[Nb_LDMK_Points]=1;
      Nb_LDMK_Points++;
    }
  }
  else{  //2.2.2) with values
    while(!feof(DataFile)){
      fscanf(DataFile,"%f  %f  %f  %f",&this->Lx[Nb_LDMK_Points],&this->Ly[Nb_LDMK_Points],&this->Lz[Nb_LDMK_Points],&this->Val[Nb_LDMK_Points]);
      Nb_LDMK_Points++;
    }
  }
  
  //close DataFile
  fclose(DataFile);
}


///Write LDMK_Points in a CSV file
void LDMK_Points::Write(char * FileName){
  int i;
  ofstream myfile;
  
  myfile.open(FileName);
  for (i=0;i<this->LDMK_Points_Nb;i++){
    myfile << static_cast<float>(this->Lx[i]) << " "  << static_cast<float>(this->Ly[i]) << " "  << static_cast<float>(this->Lz[i]) << " "  << static_cast<float>(this->Val[i]) << endl;
   }
  myfile.close();
}


///Show LDMK_Points
void LDMK_Points::Show(){
  int i;
  for (i=0;i<this->LDMK_Points_Nb;i++)
    cout << this->Lx[i] << " "  << this->Ly[i] << " "  << this->Lz[i] << " "  << this->Val[i] << endl;
  cout << "-> " << this->LDMK_Points_Nb << " points" << endl;
}




///translate the points
void LDMK_Points::Translate(float DecX,float DecY,float DecZ){
  int i;
  
  for (i=0;i<this->LDMK_Points_Nb;i++){
    this->Lx[i]+=DecX;
    this->Ly[i]+=DecY;
    this->Lz[i]+=DecZ;
    }
}

///read the LDMK_Points as the non-null values of a scalar field
void LDMK_Points::ReadInScalarField(ScalarField * img3d){
  int i,j,k;
  float epsilon;
  int PointsNb;
  FILE *dataFile;
  float i_fl,j_fl,k_fl;
  double x_loc,y_loc,z_loc,gl_loc;
  
  //1) count the number of non-null points in the image
  epsilon=0.0001;
  
  PointsNb=0;
  for(k=0;k<img3d->NZ;k++)
      for(j=0;j<img3d->NY;j++)
          for(i=0;i<img3d->NX;i++)
              if (fabs(img3d->G(i,j,k))>epsilon)
                  PointsNb++;
  
  
  //2) allocate the memory for the LDMK_Points
  this->LDMK_Points_Nb=PointsNb;
  
  this->Lx = new float [this->LDMK_Points_Nb];
  this->Ly = new float [this->LDMK_Points_Nb];
  this->Lz = new float [this->LDMK_Points_Nb];
  this->Val = new float [this->LDMK_Points_Nb];
  
  //3) fill the class
  PointsNb=0;
  for(k=0;k<img3d->NZ;k++)
      for(j=0;j<img3d->NY;j++)
          for(i=0;i<img3d->NX;i++)
              if (fabs(img3d->G(i,j,k))>epsilon){
                    i_fl=static_cast<float>(i);
                    j_fl=static_cast<float>(j);
                    k_fl=static_cast<float>(k);
                    this->Lx[PointsNb]=i_fl*img3d->Image2World[0][0]+j_fl*img3d->Image2World[0][1]+k_fl*img3d->Image2World[0][2]+img3d->Image2World[0][3];
                    this->Ly[PointsNb]=i_fl*img3d->Image2World[1][0]+j_fl*img3d->Image2World[1][1]+k_fl*img3d->Image2World[1][2]+img3d->Image2World[1][3];
                    this->Lz[PointsNb]=i_fl*img3d->Image2World[2][0]+j_fl*img3d->Image2World[2][1]+k_fl*img3d->Image2World[2][2]+img3d->Image2World[2][3];
                    this->Val[PointsNb]=img3d->G(i,j,k);
                    PointsNb++;
                }
    
}



///transform the LDMK_Points coordinates and diameters from voxels to mm according to the image 2 world properties of RefSF 
void LDMK_Points::VoxelsToMillimeters(ScalarField * RefSF){
  int i;
  float x_new,y_new,z_new;
  
  for (i=0;i<this->LDMK_Points_Nb;i++){
    x_new=this->Lx[i]*RefSF->Image2World[0][0]+this->Ly[i]*RefSF->Image2World[0][1]+this->Lz[i]*RefSF->Image2World[0][2]+RefSF->Image2World[0][3];
    y_new=this->Lx[i]*RefSF->Image2World[1][0]+this->Ly[i]*RefSF->Image2World[1][1]+this->Lz[i]*RefSF->Image2World[1][2]+RefSF->Image2World[1][3];
    z_new=this->Lx[i]*RefSF->Image2World[2][0]+this->Ly[i]*RefSF->Image2World[2][1]+this->Lz[i]*RefSF->Image2World[2][2]+RefSF->Image2World[2][3];

    this->Lx[i]=x_new;
    this->Ly[i]=y_new;
    this->Lz[i]=z_new;
  }

}

///transform the LDMK_Points coordinates and diameters from mm to voxels according to the image 2 world properties of RefSF 
void LDMK_Points::MillimetersToVoxels(ScalarField * RefSF){
  int i;
  float x_new,y_new,z_new;
  
  for (i=0;i<this->LDMK_Points_Nb;i++){
    x_new=this->Lx[i]*RefSF->World2Image[0][0]+this->Ly[i]*RefSF->World2Image[0][1]+this->Lz[i]*RefSF->World2Image[0][2]+RefSF->World2Image[0][3];
    y_new=this->Lx[i]*RefSF->World2Image[1][0]+this->Ly[i]*RefSF->World2Image[1][1]+this->Lz[i]*RefSF->World2Image[1][2]+RefSF->World2Image[1][3];
    z_new=this->Lx[i]*RefSF->World2Image[2][0]+this->Ly[i]*RefSF->World2Image[2][1]+this->Lz[i]*RefSF->World2Image[2][2]+RefSF->World2Image[2][3];

    this->Lx[i]=x_new;
    this->Ly[i]=y_new;
    this->Lz[i]=z_new;
  }
}


/// 8.2) ++++++++++++++++++ curve landmarks ++++++++++++++++++ 



///constructor
LDMK_Curves::LDMK_Curves(void){
  this->NbSeg=0;
}

///alternative constructor
LDMK_Curves::LDMK_Curves(int SegNb, int ElNb){
  int i,j;
  this->NbSeg=SegNb;
  
  //allocate the segments
  this->NbEl = new int [2*SegNb];
  this->x = new float* [2*SegNb];
  this->y = new float* [2*SegNb];
  this->z = new float* [2*SegNb];
  this->d = new float* [2*SegNb];
  
  //allocate the elements
  for (i=0;i<SegNb;i++){
    this->NbEl[i]=ElNb;
    this->x[i] = new float [ElNb];
    this->y[i] = new float [ElNb];
    this->z[i] = new float [ElNb];
    this->d[i] = new float [ElNb];
  }
  
  //give default values to the elements
  for (i=0;i<SegNb;i++) for (j=0;j<this->NbEl[i];j++){
    this->x[i][j] = 0;
    this->y[i][j] = 0;
    this->z[i][j] = 0;
    this->d[i][j] = 1;
  }
}

///destructor
LDMK_Curves::~LDMK_Curves(void){
  int i,j;
  
  if (this->NbSeg>0){
    for (i=0;i<this->NbSeg;i++) if (this->NbEl[i]>0){
        delete this->x[i];
        delete this->y[i];
        delete this->z[i];
        delete this->d[i];
    }
    
    delete this->NbEl;
  }
  this->NbSeg=0;
}

///read a LDMK_Curves in a mv3d file
void LDMK_Curves::Read(char * Mv3dName){
  FILE *DataFile;
  char CT;
  char CharTestPrec;
  int NoLoc,i,j,k,m,n;
  double xLoc,yLoc,zLoc,dLoc;

  //open the file containing the curves (at the mv3d format)
  DataFile=fopen(Mv3dName,"rb");

  //1) count the number of segments
  cout << "Load the curves of "<< Mv3dName << endl;

  fseek(DataFile,0,SEEK_SET);  //start at the beginning of the file
  CT=fgetc(DataFile);
    
  while(!feof(DataFile)){
    //read the current character and store the previous one
    CharTestPrec=CT;
    CT=fgetc(DataFile);   //move to next character in DataFile
    
    //test
    if (CharTestPrec=='\n')
      if ((CT=='0')||(CT=='1')||(CT=='2')||(CT=='3')||(CT=='4')||(CT=='5')||(CT=='6')||(CT=='7')||(CT=='8')||(CT=='9')){
        fseek(DataFile,-1,SEEK_CUR);
        fscanf(DataFile,"%d  %lf  %lf  %lf  %lf",&NoLoc,&xLoc,&yLoc,&zLoc,&dLoc);  //move to next raw
        fseek(DataFile,-2,SEEK_CUR);   //move backward of 2 characters to reinitiate properly CT and CharTestPrec
        }
    }

  NoLoc++; //the last NoLoc is the number of segments minus one
  this->NbSeg=NoLoc;   //contains the number of segments

  //2) memory allocation (1/2)
  this->NbEl = new int [2*this->NbSeg];     //we allocate two times the amount of required memory, just in case...
  this->x = new float* [2*this->NbSeg];
  this->y = new float* [2*this->NbSeg];
  this->z = new float* [2*this->NbSeg];
  this->d = new float* [2*this->NbSeg];
  
  for (i=0;i<2*this->NbSeg;i++) this->NbEl[i]=0;
    
  //3) read the number of elements in each segment
  fseek(DataFile,0,SEEK_SET);  //go back to the beginning of the file
  CT=fgetc(DataFile);
    
  while(!feof(DataFile)){
    //read the current character and store the previous one
    CharTestPrec=CT;
    CT=fgetc(DataFile);   //move to next character in DataFile
    
    //test
    if (CharTestPrec=='\n')
      if ((CT=='0')||(CT=='1')||(CT=='2')||(CT=='3')||(CT=='4')||(CT=='5')||(CT=='6')||(CT=='7')||(CT=='8')||(CT=='9')){
        fseek(DataFile,-1,SEEK_CUR);
        fscanf(DataFile,"%d  %lf  %lf  %lf  %lf",&NoLoc,&xLoc,&yLoc,&zLoc,&dLoc);  //move to next raw
        fseek(DataFile,-2,SEEK_CUR);   //move backward of 2 characters to reinitiate properly CT and CharTestPrec
        this->NbEl[NoLoc]++;  //we add an element to the current segment
        }
    }


  //4) memory allocation (2/2)

  for(i=0;i<this->NbSeg;i++){
      this->x[i] = new float [this->NbEl[i]];
      this->y[i] = new float [this->NbEl[i]];
      this->z[i] = new float [this->NbEl[i]];
      this->d[i] = new float [this->NbEl[i]];
    }
    
  for (i=0;i<this->NbSeg;i++)
    this->NbEl[i]=0;
    
  //5) load the values of each element

  fseek(DataFile,0,SEEK_SET);  //go back to the beginning of the file
  CT=fgetc(DataFile);
  fseek(DataFile,1,SEEK_CUR);
    
  while(!feof(DataFile)){
    //read the current character and store the previous one
    CharTestPrec=CT;
    CT=fgetc(DataFile);   //move to next character in DataFile
    
    //test
    if (CharTestPrec=='\n')
      if ((CT=='0')||(CT=='1')||(CT=='2')||(CT=='3')||(CT=='4')||(CT=='5')||(CT=='6')||(CT=='7')||(CT=='8')||(CT=='9')){
        fseek(DataFile,-1,SEEK_CUR);
        fscanf(DataFile,"%d  %lf  %lf  %lf  %lf",&NoLoc,&xLoc,&yLoc,&zLoc,&dLoc);  //move to next raw
        fseek(DataFile,-2,SEEK_CUR);   //move backward of 2 characters to reinitiate properly CT and CharTestPrec
        this->x[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(xLoc);
        this->y[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(yLoc);
        this->z[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(zLoc);
        this->d[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(dLoc);
  //cout << NoLoc << " " << this->NbEl[NoLoc] << " " << this->x[NoLoc][this->NbEl[NoLoc]] << " "  << this->y[NoLoc][this->NbEl[NoLoc]] << " "  << this->z[NoLoc][this->NbEl[NoLoc]]  << " "  << this->d[NoLoc][this->NbEl[NoLoc]] << endl;
        this->NbEl[NoLoc]++;
        }
    }


  //close DataFile
  fclose(DataFile);
}



///read a LDMK_Curves in a raw ascii file -> The file only contains the coordinates
///If SubSamplingFactor is higher than 1, one point of RawFileName is loaded every 'SubSamplingFactor' points
void LDMK_Curves::ReadInRawFile(char * RawFileName,int SubSamplingFactor){
  FILE *DataFile;
  char CT;
  char CharTestPrec;
  int NoLoc,i,j,k,m,n;
  float xLoc,yLoc,zLoc;
  int count;
  
  //0) open the file containing the curves (at the mv3d format)
  DataFile=fopen(RawFileName,"ra");

  //1) count the number of segments
  cout << "Load the curves of "<< RawFileName << endl;

  fseek(DataFile,0,SEEK_SET);  //start at the beginning of the file
  
  NoLoc=0;
  while(!feof(DataFile)){
    fscanf(DataFile,"%f  %f  %f",&xLoc,&yLoc,&zLoc);
    NoLoc++;
    }

  //2) memory allocation
  this->NbSeg=1;                            //contains the number of segments
  this->NbEl = new int [2*this->NbSeg];     //we allocate two times the amount of required memory, just in case...
  this->x = new float* [2*this->NbSeg];
  this->y = new float* [2*this->NbSeg];
  this->z = new float* [2*this->NbSeg];
  this->d = new float* [2*this->NbSeg];

  this->NbEl[0]=NoLoc;
  this->x[0] = new float [this->NbEl[0]];
  this->y[0] = new float [this->NbEl[0]];
  this->z[0] = new float [this->NbEl[0]];
  this->d[0] = new float [this->NbEl[0]];
    
    
  //3) load the values of each element and eventually subsample the curve

  fseek(DataFile,0,SEEK_SET);  //go back to the beginning of the file
    
  NoLoc=0;
  count=0;
  while(!feof(DataFile)){
    count++;
    fscanf(DataFile,"%f  %f  %f",&xLoc,&yLoc,&zLoc);
    if (count==SubSamplingFactor){
      this->x[0][NoLoc]=static_cast<float>(xLoc);
      this->y[0][NoLoc]=static_cast<float>(yLoc);
      this->z[0][NoLoc]=static_cast<float>(zLoc);
      this->d[0][NoLoc]=1;
      NoLoc++;
      count=0;
    }
  }
  
  this->NbEl[0]=NoLoc;
  
  //close DataFile
  fclose(DataFile);
}




///create a LDMK_Curves containing one curve with NbElem elements
void LDMK_Curves::Create_LDMK_1_Curve(int NbElem){
  int i;

  this->NbSeg=1;

  this->NbEl = new int [2];
  this->x = new float* [2];
  this->y = new float* [2];
  this->z = new float* [2];
  this->d = new float* [2];
  
  this->NbEl[0]=NbElem;
  this->x[0] = new float [NbElem];
  this->y[0] = new float [NbElem];
  this->z[0] = new float [NbElem];
  this->d[0] = new float [NbElem];
  this->NbEl[1]=0;
  
  for (i=0;i<this->NbEl[0];i++){
    this->x[0][i] = 0;
    this->y[0][i] = 0;
    this->z[0][i] = 0;
    this->d[0][i] = 0;
  }
}



///create a LDMK_Curves containing NbCurves void curves (no elements)
void LDMK_Curves::Create_LDMK_N_VoidCurves(int N){
  int i;

  this->NbSeg=N;

  this->NbEl = new int [2*N+1];
  this->x = new float* [2*N+1];
  this->y = new float* [2*N+1];
  this->z = new float* [2*N+1];
  this->d = new float* [2*N+1];
  
  for (i=0;i<2*N+1;i++)
    this->NbEl[i]=0;
}




///write a LDMK_Curves in a mv3d file.
///Set Preserve_IDs to 1 to preserve the original segment identifiers (default). They are optimaly resampled otherwise.
void LDMK_Curves::Write(char * Mv3dName, int Preserve_IDs){
  int i,iBis,j,tempInt;
  FILE *DataFile;
  int LocNbSeg,LocNbEl;
  int temp;
  int * NbSegEnds;
  int CaseSegEndI,CaseSegEndJ;
  int SegEndI,SegEndJ;
  int Nb1,Nb2,Nb3,Nb4,Nb5,Nb6,Nb7;
  int NbIntersections;

  //1) open the file in which the data will be saved
  DataFile=fopen(Mv3dName,"w");

  //2) header

  //2.1 - total number of segments and elements
  LocNbEl=0;
  LocNbSeg=0;
  for (i=0;i<2*this->NbSeg;i++){    //the "2 *" is in case there are additional segments
    if (this->NbEl[i]!=0){
      for (j=0;j<this->NbEl[i];j++)
  LocNbEl++;
      LocNbSeg++;
      }
    }
  
  //2.2 - count the number of intersections
  NbSegEnds = new int [4*(this->NbSeg)];  // we consider a maximum of 4*[the segments number]  for the number of segment-ends

  for (i=0;i<2*this->NbSeg;i++){
    NbSegEnds[2*i]=0;
    NbSegEnds[2*i+1]=0;
    }

  for (i=0;i<2*this->NbSeg-1;i++) if (this->NbEl[i]!=0) for(CaseSegEndI=0;CaseSegEndI<2;CaseSegEndI++){
    for (j=i+1;j<2*this->NbSeg;j++) if (this->NbEl[j]!=0) for(CaseSegEndJ=0;CaseSegEndJ<2;CaseSegEndJ++){
      
      //number of the current element
      SegEndI=(this->NbEl[i]-1)*CaseSegEndI;
      SegEndJ=(this->NbEl[j]-1)*CaseSegEndJ;
      
      
      //add the number of segment-ends if an intersection is found
      if (fabs(this->x[i][SegEndI]-this->x[j][SegEndJ])<0.1)
      if (fabs(this->y[i][SegEndI]-this->y[j][SegEndJ])<0.1)
      if (fabs(this->z[i][SegEndI]-this->z[j][SegEndJ])<0.1){
    NbSegEnds[2*i+CaseSegEndI]++;
  NbSegEnds[2*j+CaseSegEndJ]++;
  }
      }
    }


  //for each segment-end, the number N of segment arriving on a node is: NbSegEnds[X]= sum_{i=1}^{i=N-1} i = N*(N-1)/2
  Nb2=0; Nb3=0; Nb4=0; Nb5=0; Nb6=0; Nb7=0;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for(CaseSegEndI=0;CaseSegEndI<2;CaseSegEndI++){
    if (NbSegEnds[2*i+CaseSegEndI]==0)  NbSegEnds[2*i+CaseSegEndI]=0;
    else if (NbSegEnds[2*i+CaseSegEndI]==1)   {NbSegEnds[2*i+CaseSegEndI]=2; Nb2++;}   //
    else if (NbSegEnds[2*i+CaseSegEndI]<=3)   {NbSegEnds[2*i+CaseSegEndI]=3; Nb3++;}   //
    else if (NbSegEnds[2*i+CaseSegEndI]<=6)   {NbSegEnds[2*i+CaseSegEndI]=4; Nb4++;}   // -> we use <= et not == for potential loops
    else if (NbSegEnds[2*i+CaseSegEndI]<=10)  {NbSegEnds[2*i+CaseSegEndI]=5; Nb5++;}  //
    else if (NbSegEnds[2*i+CaseSegEndI]<=15)  {NbSegEnds[2*i+CaseSegEndI]=6; Nb6++;}  //
    else if (NbSegEnds[2*i+CaseSegEndI]<=21)  {NbSegEnds[2*i+CaseSegEndI]=7; Nb7++;}  //
    }

  NbIntersections=Nb2/2+Nb3/3+Nb4/4+Nb5/5+Nb6/6+Nb7/7;
    

  //2.3) write header
  fprintf(DataFile,"# MicroVisu3D file\n");
  fprintf(DataFile,"# Number of lines   %d\n",LocNbSeg);
  fprintf(DataFile,"# Number of points  %d\n",LocNbEl);
  fprintf(DataFile,"# Number of inter.  %d\n",NbIntersections);
  fprintf(DataFile,"#\n");
  fprintf(DataFile,"# No    x    y    z    d\n");
  fprintf(DataFile,"#\n");

  //3) save data


  if (Preserve_IDs==1){
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
      for (j=0;j<this->NbEl[i];j++)
  fprintf(DataFile,"%d  %lf  %lf  %lf  %lf\n",i,static_cast<double>(this->x[i][j]),static_cast<double>(this->y[i][j]),static_cast<double>(this->z[i][j]),static_cast<double>(this->d[i][j]));
      fprintf(DataFile,"\n");
      }
    }
  else{
    tempInt=0;
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
      for (j=0;j<this->NbEl[i];j++)
  fprintf(DataFile,"%d  %lf  %lf  %lf  %lf\n",tempInt,static_cast<double>(this->x[i][j]),static_cast<double>(this->y[i][j]),static_cast<double>(this->z[i][j]),static_cast<double>(this->d[i][j]));
      tempInt++;
      fprintf(DataFile,"\n");
      }
    }

  fclose(DataFile);
}

///Export LDMK_Curves in a vtk file
void LDMK_Curves::ExportAsVtk(char * VtkName,int ShowVolume){
  int i,iBis,j,tempInt;
  FILE *DataFile;
  int LocNbSeg,LocNbEl,LocNbLinks;
  int temp;
  int * NbSegEnds;
  int CaseSegEndI,CaseSegEndJ;
  int SegEndI,SegEndJ;
  int Nb1,Nb2,Nb3,Nb4,Nb5,Nb6,Nb7;
  int NbIntersections;
  float xloc,yloc,zloc,rloc;
  int NbEl;
  
  //1) open the file in which the data will be saved
  DataFile=fopen(VtkName,"w");

  //2) init
  
  //2.1) total number of segments and elements
  LocNbEl=0;
  LocNbSeg=0;
  LocNbLinks=0;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){    //the "2 *" is in case there are additional segments
      for (j=0;j<this->NbEl[i];j++){
  LocNbEl++;
  if (j!=this->NbEl[i]-1) LocNbLinks++;
      }
      LocNbSeg++;
    }

  

  //2.2) write header
  fprintf(DataFile,"# vtk DataFile Version 3.0\n");
  fprintf(DataFile,"vtk output\n");
  fprintf(DataFile,"ASCII\n");
  
  //3) save data
  
  //3.1) save the points
  fprintf(DataFile,"DATASET POLYDATA\n");
  fprintf(DataFile,"POINTS %d float\n",7*LocNbEl);

  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
      xloc=this->x[i][j];
      yloc=this->y[i][j];
      zloc=this->z[i][j];
      rloc=this->d[i][j]/2;
      fprintf(DataFile,"%f  %f  %f\n",xloc+rloc,yloc,zloc);
      fprintf(DataFile,"%f  %f  %f\n",xloc-rloc,yloc,zloc);
      fprintf(DataFile,"%f  %f  %f\n",xloc,yloc+rloc,zloc);
      fprintf(DataFile,"%f  %f  %f\n",xloc,yloc-rloc,zloc);
      fprintf(DataFile,"%f  %f  %f\n",xloc,yloc,zloc+rloc);
      fprintf(DataFile,"%f  %f  %f\n",xloc,yloc,zloc-rloc);
      fprintf(DataFile,"%f  %f  %f\n",xloc,yloc,zloc);
  }
  
  fprintf(DataFile,"\n");

  
  //3.2) save the polygons
  if (ShowVolume==1){
    fprintf(DataFile,"POLYGONS %d %d\n",8*LocNbEl,32*LocNbEl);
  
    tempInt=0;
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+1,tempInt+3,tempInt+4);
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+1,tempInt+2,tempInt+4);
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+0,tempInt+3,tempInt+4);
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+0,tempInt+2,tempInt+4);
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+1,tempInt+3,tempInt+5);
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+1,tempInt+2,tempInt+5);
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+0,tempInt+3,tempInt+5);
        fprintf(DataFile,"3  %d  %d  %d\n",tempInt+0,tempInt+2,tempInt+5);
        tempInt+=7;
    }
    
    fprintf(DataFile,"\n");
  }
   
  //3.3) save the curves
  fprintf(DataFile,"LINES %d %d\n",LocNbLinks,3*LocNbLinks);

  tempInt=0;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
      if (j<this->NbEl[i]-1) fprintf(DataFile,"2  %d  %d\n",tempInt+6,tempInt+13);
      tempInt+=7;
  }
  
  //3.4) colors representing the diameter
  NbEl=0;
	for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++) NbEl+=7;
  
  fprintf(DataFile,"\n");
  fprintf(DataFile,"POINT_DATA %d\n",NbEl);
  fprintf(DataFile,"SCALARS volume float\n");
  fprintf(DataFile,"LOOKUP_TABLE default\n");

  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
      fprintf(DataFile,"%.2lf\n",this->d[i][j]);
      fprintf(DataFile,"%.2lf\n",this->d[i][j]);
      fprintf(DataFile,"%.2lf\n",this->d[i][j]);
      fprintf(DataFile,"%.2lf\n",this->d[i][j]);
      fprintf(DataFile,"%.2lf\n",this->d[i][j]);
      fprintf(DataFile,"%.2lf\n",this->d[i][j]);
      fprintf(DataFile,"%.2lf\n",this->d[i][j]);
  }


  
  fprintf(DataFile,"\n");
  
  fclose(DataFile);
}


///Delete a segment
void LDMK_Curves::DeleteSegment(int Seg){
  delete this->x[Seg];
  delete this->y[Seg];
  delete this->z[Seg];
  delete this->d[Seg];
  this->NbEl[Seg]=0;
}


///Reduce the number of elements in a segment - consider the elements between FirstEl and LastEl (included)
void LDMK_Curves::ReduceElNumber(int IdSeg, int FirstEl, int LastEl){
  int i;
  
  //tests about the maximum and minimum values
  if (FirstEl>LastEl){
    i=FirstEl;
    FirstEl=LastEl;
    LastEl=i;
  }
  if (FirstEl<0) FirstEl=0;
  if (LastEl>this->NbEl[IdSeg]-1) LastEl=this->NbEl[IdSeg]-1;
  
  //move all data if necessary
  if (FirstEl>0){
    for (i=FirstEl;i<=LastEl;i++){
      x[IdSeg][i-FirstEl]=x[IdSeg][i];
      y[IdSeg][i-FirstEl]=y[IdSeg][i];
      z[IdSeg][i-FirstEl]=z[IdSeg][i];
      d[IdSeg][i-FirstEl]=d[IdSeg][i];
    }
  }
  
  //cut the segment
  if ((FirstEl>0)||(this->NbEl[IdSeg]-1>LastEl)) this->NbEl[IdSeg]=LastEl-FirstEl+1;
}

///Merge two segments at their nearest extremity. The new segment is saved in Seg1.
void LDMK_Curves::MergeSegments(int Seg1, int Seg2){
  float * LocX;
  float * LocY;
  float * LocZ;
  float * LocD;
  int Side1,Side2,Side1_star,Side2_star;
  float tempFL,BestTempFL;
  int i;
  int NewElNb;
  
  if ((this->NbEl[Seg1]==0)||(this->NbEl[Seg2]==0)) return;

  
  //find the nearest connection
  for (Side1=0;Side1<2;Side1++) for (Side2=0;Side2<2;Side2++){
    tempFL=(this->x[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->x[Seg2][(this->NbEl[Seg2]-1)*Side2])*(this->x[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->x[Seg2][(this->NbEl[Seg2]-1)*Side2]);
    tempFL+=(this->y[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->y[Seg2][(this->NbEl[Seg2]-1)*Side2])*(this->y[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->y[Seg2][(this->NbEl[Seg2]-1)*Side2]);
    tempFL+=(this->z[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->z[Seg2][(this->NbEl[Seg2]-1)*Side2])*(this->z[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->z[Seg2][(this->NbEl[Seg2]-1)*Side2]);
  
    if ((Side1==0)&&(Side2==0)){
      Side1_star=Side1;
      Side2_star=Side2;
      BestTempFL=tempFL;
    }
    
    if (BestTempFL>tempFL){
      Side1_star=Side1;
      Side2_star=Side2;
      BestTempFL=tempFL;
    }
  }
  
  
  //allocate memory for the temporary segment
  NewElNb=this->NbEl[Seg1]+this->NbEl[Seg2]-1;
  LocX=new float [NewElNb];
  LocY=new float [NewElNb];
  LocZ=new float [NewElNb];
  LocD=new float [NewElNb];
  
  //fill the temporary segment
  if (Side1_star==1){
    for (i=0;i<this->NbEl[Seg1];i++){
      LocX[i]=this->x[Seg1][i];
      LocY[i]=this->y[Seg1][i];
      LocZ[i]=this->z[Seg1][i];
      LocD[i]=this->d[Seg1][i];
    }
  }
  else{
    for (i=0;i<this->NbEl[Seg1];i++){
      LocX[i]=this->x[Seg1][this->NbEl[Seg1]-1-i];
      LocY[i]=this->y[Seg1][this->NbEl[Seg1]-1-i];
      LocZ[i]=this->z[Seg1][this->NbEl[Seg1]-1-i];
      LocD[i]=this->d[Seg1][this->NbEl[Seg1]-1-i];
    }
  }

  if (Side2_star==0){
    for (i=1;i<this->NbEl[Seg2];i++){
      LocX[this->NbEl[Seg1]+i-1]=this->x[Seg2][i];
      LocY[this->NbEl[Seg1]+i-1]=this->y[Seg2][i];
      LocZ[this->NbEl[Seg1]+i-1]=this->z[Seg2][i];
      LocD[this->NbEl[Seg1]+i-1]=this->d[Seg2][i];
    }
  }
  else{
    for (i=1;i<this->NbEl[Seg2];i++){
      LocX[this->NbEl[Seg1]+i-1]=this->x[Seg2][this->NbEl[Seg2]-1-i];
      LocY[this->NbEl[Seg1]+i-1]=this->y[Seg2][this->NbEl[Seg2]-1-i];
      LocZ[this->NbEl[Seg1]+i-1]=this->z[Seg2][this->NbEl[Seg2]-1-i];
      LocD[this->NbEl[Seg1]+i-1]=this->d[Seg2][this->NbEl[Seg2]-1-i];
    }
  }
  
  
  //delete Seg1 and Seg2
  this->DeleteSegment(Seg1);
  this->DeleteSegment(Seg2);
  
  //rebuild Seg1
  this->NbEl[Seg1]=NewElNb;
  this->x[Seg1] = new float [NewElNb];
  this->y[Seg1] = new float [NewElNb];
  this->z[Seg1] = new float [NewElNb];
  this->d[Seg1] = new float [NewElNb];
  
  for (i=0;i<NewElNb;i++){
    this->x[Seg1][i]=LocX[i];
    this->y[Seg1][i]=LocY[i];
    this->z[Seg1][i]=LocZ[i];
    this->d[Seg1][i]=LocD[i];
  }
  
  //dealloc temporary data
  delete LocX;
  delete LocY;
  delete LocZ;
  delete LocD;

}




///Resample the number of elements in a segment. The resampled segment will be roughly homogeneously samped in space.
void LDMK_Curves::ResampleSegment(int IdSeg, int NewElNb){
  float * NewX;
  float * NewY;
  float * NewZ;
  float * NewD;
  int i;
  float curloc;
  float TotalDistToOrigin,CurrentDistToOrigin,TargetDistToOrigin;
  int   Current_Element;
  int ResamplingFactor;
  float im1_x,i_x,im1_y,i_y,im1_z,i_z,im1_d,i_d;

  
  if ((this->NbEl[IdSeg]<2)) return;
  
  ResamplingFactor=10000;
  
  
  //1) INIT  -> allocate memory for the temporary segment
  NewX=new float [NewElNb];
  NewY=new float [NewElNb];
  NewZ=new float [NewElNb];
  NewD=new float [NewElNb];
  
  
  //2) resample the curve using B-spline interpolation and homogeneous steps
  
  //2.1) Compute the total length of the original dense curve
  TotalDistToOrigin=0;
  i_x=BSplineCurveSampler(this->x[IdSeg],this->NbEl[IdSeg],0);
  i_y=BSplineCurveSampler(this->y[IdSeg],this->NbEl[IdSeg],0);
  i_z=BSplineCurveSampler(this->z[IdSeg],this->NbEl[IdSeg],0);
  i_d=BSplineCurveSampler(this->d[IdSeg],this->NbEl[IdSeg],0);
  
  
  for (i=1;i<=(this->NbEl[IdSeg]-1)*ResamplingFactor;i++){
    curloc=static_cast<float>(i)/static_cast<float>((this->NbEl[IdSeg]-1)*ResamplingFactor);
    
    im1_x=i_x;
    im1_y=i_y;
    im1_z=i_z;
    im1_d=i_d;
  
    i_x=BSplineCurveSampler(this->x[IdSeg],this->NbEl[IdSeg],curloc);
    i_y=BSplineCurveSampler(this->y[IdSeg],this->NbEl[IdSeg],curloc);
    i_z=BSplineCurveSampler(this->z[IdSeg],this->NbEl[IdSeg],curloc);
    i_d=BSplineCurveSampler(this->d[IdSeg],this->NbEl[IdSeg],curloc);
    
    TotalDistToOrigin+=sqrt( ((i_x-im1_x)*(i_x-im1_x)) + ((i_y-im1_y)*(i_y-im1_y)) +  ((i_z-im1_z)*(i_z-im1_z))   );
  }
  
  
  //2.2) Compute the total length of the original dense curve
  CurrentDistToOrigin=0;
  
  i_x=BSplineCurveSampler(this->x[IdSeg],this->NbEl[IdSeg],0);
  i_y=BSplineCurveSampler(this->y[IdSeg],this->NbEl[IdSeg],0);
  i_z=BSplineCurveSampler(this->z[IdSeg],this->NbEl[IdSeg],0);
  i_d=BSplineCurveSampler(this->d[IdSeg],this->NbEl[IdSeg],0);
  
  NewX[0]=i_x;
  NewY[0]=i_y;
  NewZ[0]=i_z;
  NewD[0]=i_d;
  
  Current_Element=1;
  TargetDistToOrigin=static_cast<float>(Current_Element)*TotalDistToOrigin/static_cast<float>(NewElNb-1);
  
  for (i=1;i<=(this->NbEl[IdSeg]-1)*ResamplingFactor;i++){
    curloc=static_cast<float>(i)/static_cast<float>((this->NbEl[IdSeg]-1)*ResamplingFactor);
    
    im1_x=i_x;
    im1_y=i_y;
    im1_z=i_z;
    im1_d=i_d;
    
    i_x=BSplineCurveSampler(this->x[IdSeg],this->NbEl[IdSeg],curloc);
    i_y=BSplineCurveSampler(this->y[IdSeg],this->NbEl[IdSeg],curloc);
    i_z=BSplineCurveSampler(this->z[IdSeg],this->NbEl[IdSeg],curloc);
    i_d=BSplineCurveSampler(this->d[IdSeg],this->NbEl[IdSeg],curloc);
    
    CurrentDistToOrigin+=sqrt( ((i_x-im1_x)*(i_x-im1_x)) + ((i_y-im1_y)*(i_y-im1_y)) +  ((i_z-im1_z)*(i_z-im1_z))   );
    
    
    if (CurrentDistToOrigin>TargetDistToOrigin){
      NewX[Current_Element]=i_x;
      NewY[Current_Element]=i_y;
      NewZ[Current_Element]=i_z;
      NewD[Current_Element]=i_d;
      
      Current_Element++;
      if (Current_Element>NewElNb-1) Current_Element=NewElNb-1;
      TargetDistToOrigin=static_cast<float>(Current_Element)*TotalDistToOrigin/static_cast<float>(NewElNb-1);;
    }
  }
  
  
  NewX[NewElNb-1]=BSplineCurveSampler(this->x[IdSeg],this->NbEl[IdSeg],1);
  NewY[NewElNb-1]=BSplineCurveSampler(this->y[IdSeg],this->NbEl[IdSeg],1);
  NewZ[NewElNb-1]=BSplineCurveSampler(this->z[IdSeg],this->NbEl[IdSeg],1);
  NewD[NewElNb-1]=BSplineCurveSampler(this->d[IdSeg],this->NbEl[IdSeg],1);

  
  //3) Copy the resampled segment and clean-up memory
  
  //3.1) delete IdSeg
  this->DeleteSegment(IdSeg);
  
  //3.2) rebuild IdSeg
  this->NbEl[IdSeg]=NewElNb;
  this->x[IdSeg] = new float [NewElNb];
  this->y[IdSeg] = new float [NewElNb];
  this->z[IdSeg] = new float [NewElNb];
  this->d[IdSeg] = new float [NewElNb];
  
  for (i=0;i<NewElNb;i++){
    this->x[IdSeg][i]=NewX[i];
    this->y[IdSeg][i]=NewY[i];
    this->z[IdSeg][i]=NewZ[i];
    this->d[IdSeg][i]=NewD[i];
  }
  
  //3.3) dealloc temporary data
  delete NewX;
  delete NewY;
  delete NewZ;
  delete NewD;
  
}


///return a weight for BSpline interpolation. x is in [-1.5,1.5] and the weights are centered in 0   
///-> it will then be necessary to use the values at pts: [pt-2] [pt-1] [pt] [pt+1] [pt+2] for interpolation
///+++ properly assessed +++
float BSplineWeight(float x){
  float u;
  
  u=(x+1.5)*0.75/3;
  
  if (u<0) return 0;
  else if (u<0.25) return (8*u*u);
  else if (u<0.5) return (-1.5+12*u-16*u*u);
  else if (u<0.75) return (4.5-12*u+8*u*u);
  else return 0;
}


///return the value of a curve 'densified' using B-spline interpolation at a given location
///location should have a value in [0 , 1]
float BSplineCurveSampler(float * Values, int NbValues,float location){
  int intLocation;
  float residualLocation;
  float interpolatedValue;
  
  location*=NbValues-1;
  
  //manage locations outside of the domain
  if (location<0) location=0;
  if (location>NbValues-1) location=NbValues-1;

  intLocation=static_cast<int>(location);
  residualLocation=location-intLocation;

  
  //interpolation depending on the location due to the treatment of boundaries
  if (intLocation==0){//location between 0 and 1
    interpolatedValue=Values[intLocation]*BSplineWeight(-1-residualLocation) +
                      Values[intLocation]*BSplineWeight(-residualLocation) +
                      Values[intLocation+1]*BSplineWeight(1-residualLocation) +
                      Values[intLocation+2]*BSplineWeight(2-residualLocation);
  }
  else if (intLocation<NbValues-2){//location between 1 and NbValues-2
    interpolatedValue=Values[intLocation-1]*BSplineWeight(-1-residualLocation) +
                      Values[intLocation]*BSplineWeight(-residualLocation) +
                      Values[intLocation+1]*BSplineWeight(1-residualLocation) +
                      Values[intLocation+2]*BSplineWeight(2-residualLocation);
  }
  else if (intLocation<NbValues-1){//location between NbValues-2 and NbValues-1
    interpolatedValue=Values[intLocation-1]*BSplineWeight(-1-residualLocation) +
                      Values[intLocation]*BSplineWeight(-residualLocation) +
                      Values[intLocation+1]*BSplineWeight(1-residualLocation) +
                      Values[intLocation+1]*BSplineWeight(2-residualLocation);
  }
  else{//location equals NbValues-1
    interpolatedValue=Values[intLocation-1]*BSplineWeight(-1-residualLocation) +
                      Values[intLocation]*BSplineWeight(-residualLocation) +
                      Values[intLocation]*BSplineWeight(1-residualLocation) +
                      Values[intLocation]*BSplineWeight(2-residualLocation);
  }
  
  return interpolatedValue;
}



///Generate a segment of size 'NbElements' with only null values
///The ID of the segment in 'this' is returned. -1 is returned if there is no more available space
int LDMK_Curves::GenerateVoidSegment(int NbElements){
  int i,CurrentSegment;
  
  //1) find the id of the segment
  i=this->GetSegNumber();
  CurrentSegment=-1;
  while((i<2*this->GetSegNumber())&&(CurrentSegment==-1)){
      if (this->GetElNumber(i)==0){
          CurrentSegment=i;
      }
  i++;
  }
  
  if (CurrentSegment==-1){
        cout << "No more free segments to fill a gap" << endl;
        return -1;
  }
  
  //2) allocate memory for the new segment
  this->NbEl[CurrentSegment]=NbElements;
  this->x[CurrentSegment] = new float [NbElements];
  this->y[CurrentSegment] = new float [NbElements];
  this->z[CurrentSegment] = new float [NbElements];
  this->d[CurrentSegment] = new float [NbElements];
  
  for (i=0;i<NbElements;i++){
    this->x[CurrentSegment][i]=0;
    this->y[CurrentSegment][i]=0;
    this->z[CurrentSegment][i]=0;
    this->d[CurrentSegment][i]=1;
  }
  
  return CurrentSegment;
}




///count the number of segments related to the segment end SegEnd (= 0 pr 1) of the segment Seg
int LDMK_Curves::CountLinkedSeg(int Seg,int SegEnd,float epsilon){
  int i,SideI;
  float LocX,LocY,LocZ;
  float LocX2,LocY2,LocZ2;
  int count;
  float tmpFl;
  
  if (SegEnd==0){ LocX=this->x[Seg][0];                 LocY=this->y[Seg][0];                 LocZ=this->z[Seg][0]; }
  if (SegEnd==1){ LocX=this->x[Seg][this->NbEl[Seg]-1]; LocY=this->y[Seg][this->NbEl[Seg]-1]; LocZ=this->z[Seg][this->NbEl[Seg]-1]; }
  
  count=0;
  for (i=0;i<2*this->NbSeg;i++) if ((this->NbEl[i]!=0)&&(i!=Seg)) for (SideI=0;SideI<2;SideI++){
    if (SideI==0){ LocX2=this->x[i][0];               LocY2=this->y[i][0];               LocZ2=this->z[i][0]; }
    if (SideI==1){ LocX2=this->x[i][this->NbEl[i]-1]; LocY2=this->y[i][this->NbEl[i]-1]; LocZ2=this->z[i][this->NbEl[i]-1]; }

    tmpFl=sqrt(((LocX-LocX2)*(LocX-LocX2))+((LocY-LocY2)*(LocY-LocY2))+((LocZ-LocZ2)*(LocZ-LocZ2)));
    
    if (tmpFl<epsilon)count++;
  }

  return count;

}

///Clean-up a moderately large LDMK_Curves structure
/// -> two segments ends are supposed linked if their distance is less than epsilon
/// -> merge the segments linked by a node with only two segments
/// -> remove the segments linked to only one node and for which the node has more than two segments related to other nodes
void LDMK_Curves::CleanUp(float epsilon){
  int changes;
  int i,j,j_star;
  int SideI, SideJ;
  float LocX,LocY,LocZ;
  float LocX2,LocY2,LocZ2;
  float LocX3,LocY3,LocZ3;
  float LocX4,LocY4,LocZ4;
  int count,count2;
  float tmpFl,tmpFl2,tmpFl3,tmpFl4;
  
  cout << "clean-up the network" << endl;
   
  //1) Remove isolated points
  //cout << "Remove isolated points" << endl;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]==1) this->DeleteSegment(i);
  
  //2) remove the segments linked to only one node and for which the node has more than two segments related to other nodes
  //cout << "Useless segments removal" << endl;
  int NgbhNbS0,NgbhNbS1,RefSide;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) if (this->NbEl[i]<10){
    NgbhNbS0=CountLinkedSeg(i,0,epsilon);
    NgbhNbS1=CountLinkedSeg(i,1,epsilon);
    
    //only one node related to something
    if ( ((NgbhNbS0==0)&&(NgbhNbS1!=0))  || ((NgbhNbS1==0)&&(NgbhNbS0!=0)) ){
      if (NgbhNbS0==0) RefSide=1;
      else RefSide=0;
      
      //find the ngbhs of 'i/RefSide' and treat them
      LocX=this->x[i][(this->NbEl[i]-1)*RefSide];
      LocY=this->y[i][(this->NbEl[i]-1)*RefSide];
      LocZ=this->z[i][(this->NbEl[i]-1)*RefSide];
      
      count=0;
      
      for (j=0;j<2*this->NbSeg;j++) if ((this->NbEl[j]!=0)&&(i!=j)) for (SideJ=0;SideJ<2;SideJ++) if (fabs(LocX-this->x[j][(this->NbEl[j]-1)*SideJ])<epsilon){
        LocX2=this->x[j][(this->NbEl[j]-1)*SideJ];
        LocY2=this->y[j][(this->NbEl[j]-1)*SideJ];
        LocZ2=this->z[j][(this->NbEl[j]-1)*SideJ];
        
        tmpFl=sqrt(((LocX-LocX2)*(LocX-LocX2))+((LocY-LocY2)*(LocY-LocY2))+((LocZ-LocZ2)*(LocZ-LocZ2)));
  
        //'i/RefSide' is linked to 'j/SideJ'
        if (tmpFl<epsilon){
          if (SideJ==0) count2=CountLinkedSeg(j,1,epsilon);
          if (SideJ==1) count2=CountLinkedSeg(j,0,epsilon);
          
          //'j/abs(1-SideJ)' is related to something
          if (count2>0) count++;
        }
      }
      
      if (count>=2){
  this->DeleteSegment(i);
  //cout << "Delete " << i << endl;
      }
    }
  }
  
  //3) remove the loops
  //cout << "Useless segments removal" << endl;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
        LocX=this->x[i][0];
        LocY=this->y[i][0];
        LocZ=this->z[i][0]; 
        LocX2=this->x[i][this->NbEl[i]-1];
        LocY2=this->y[i][this->NbEl[i]-1];
        LocZ2=this->z[i][this->NbEl[i]-1]; 
        tmpFl=sqrt(((LocX-LocX2)*(LocX-LocX2))+((LocY-LocY2)*(LocY-LocY2))+((LocZ-LocZ2)*(LocZ-LocZ2)));
        
        if (tmpFl<epsilon) this->DeleteSegment(i);
  }
  
  //4) merge the segments linked by a node with only two segments
  changes=1;
  while (changes==1){
    //cout << "Iteration of segments merging" << endl;
    changes=0;
    
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
      for (SideI=0;SideI<2;SideI++){
        //define the segment end treated
        if (SideI==0){ LocX=this->x[i][0];               LocY=this->y[i][0];               LocZ=this->z[i][0]; }
        if (SideI==1){ LocX=this->x[i][this->NbEl[i]-1]; LocY=this->y[i][this->NbEl[i]-1]; LocZ=this->z[i][this->NbEl[i]-1]; }
        
        //count the number of segment-ends linked to LocX,LocY,LocZ
        count=0;
        for (j=0;j<2*this->NbSeg;j++) if (j!=i) if (this->NbEl[j]!=0) if ((fabs(LocX-this->x[j][this->NbEl[j]-1])<epsilon)||(fabs(LocX-this->x[j][0])<epsilon)) for (SideJ=0;SideJ<2;SideJ++){
          //define the segment end treated
          if (SideJ==0){ LocX2=this->x[j][0];               LocY2=this->y[j][0];               LocZ2=this->z[j][0]; }
          if (SideJ==1){ LocX2=this->x[j][this->NbEl[j]-1]; LocY2=this->y[j][this->NbEl[j]-1]; LocZ2=this->z[j][this->NbEl[j]-1]; }
    
          //compute the distance and make the test
          tmpFl=sqrt(((LocX-LocX2)*(LocX-LocX2))+((LocY-LocY2)*(LocY-LocY2))+((LocZ-LocZ2)*(LocZ-LocZ2)));
    
          if (tmpFl<epsilon){
            count++;
            j_star=j;
          }
        }
  
        //merge the segments if only one segment is related to segment i
        if (count==1){
          this->MergeSegments(i,j_star);
          //cout << "Merge " << i << " " << j_star << endl;
          changes=1;
          SideI=2;
          SideJ=2;
        }
      }
    }
  }
  
  //5) remove segments which have the same extremities as another segment
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=i+1;j<2*this->NbSeg;j++) if (this->NbEl[j]!=0)
        if ((fabs(this->x[i][0]-this->x[j][0])<epsilon)||
            (fabs(this->x[i][0]-this->x[j][this->NbEl[j]-1])<epsilon)||
            (fabs(this->x[i][this->NbEl[i]-1]-this->x[j][0])<epsilon)||
            (fabs(this->x[i][this->NbEl[i]-1]-this->x[j][this->NbEl[j]-1])<epsilon)){
        
        LocX=this->x[i][0];
        LocY=this->y[i][0];
        LocZ=this->z[i][0]; 
        
        LocX2=this->x[i][this->NbEl[i]-1];
        LocY2=this->y[i][this->NbEl[i]-1];
        LocZ2=this->z[i][this->NbEl[i]-1];
        
        LocX3=this->x[j][0];
        LocY3=this->y[j][0];
        LocZ3=this->z[j][0]; 
        
        LocX4=this->x[j][this->NbEl[j]-1];
        LocY4=this->y[j][this->NbEl[j]-1];
        LocZ4=this->z[j][this->NbEl[j]-1];
        
        
        tmpFl=sqrt(((LocX-LocX3)*(LocX-LocX3))+((LocY-LocY3)*(LocY-LocY3))+((LocZ-LocZ3)*(LocZ-LocZ3)));
        tmpFl2=sqrt(((LocX2-LocX4)*(LocX2-LocX4))+((LocY2-LocY4)*(LocY2-LocY4))+((LocZ2-LocZ4)*(LocZ2-LocZ4)));
        
        tmpFl3=sqrt(((LocX-LocX4)*(LocX-LocX4))+((LocY-LocY4)*(LocY-LocY4))+((LocZ-LocZ4)*(LocZ-LocZ4)));
        tmpFl4=sqrt(((LocX2-LocX3)*(LocX2-LocX3))+((LocY2-LocY3)*(LocY2-LocY3))+((LocZ2-LocZ3)*(LocZ2-LocZ3)));
        
        if (((tmpFl<epsilon)&&(tmpFl2<epsilon))||((tmpFl3<epsilon)&&(tmpFl4<epsilon))) this->DeleteSegment(j);
  }

  
  
}



///Smooth the segments
void LDMK_Curves::Smooth(int ItNb){
  int it,i,j;
  float x_new,y_new,z_new,d_new;
  float x_svg,y_svg,z_svg,d_svg;
  
  
  for (it=0;it<ItNb;it++){
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=1;j<this->NbEl[i]-1;j++){
      if (j==1){
        x_svg=this->x[i][0];
        y_svg=this->y[i][0];
        z_svg=this->z[i][0];
        d_svg=this->d[i][0];
      }
      
      x_new=(x_svg+this->x[i][j]+this->x[i][j+1])/3;
      y_new=(y_svg+this->y[i][j]+this->y[i][j+1])/3;
      z_new=(z_svg+this->z[i][j]+this->z[i][j+1])/3;
      d_new=(d_svg+this->d[i][j]+this->d[i][j+1])/3;
      
      x_svg=this->x[i][j];
      y_svg=this->y[i][j];
      z_svg=this->z[i][j];
      d_svg=this->d[i][j];
      
      this->x[i][j]=x_new;
      this->y[i][j]=y_new;
      this->z[i][j]=z_new;
      this->d[i][j]=d_new;
    }
  }

}


///Smooth only one segments
void LDMK_Curves::SmoothOneSegment(int SegmentID,int ItNb){
  int it,j;
  float x_new,y_new,z_new,d_new;
  float x_svg,y_svg,z_svg,d_svg;
  
  
  for (it=0;it<ItNb;it++){
    if (this->NbEl[SegmentID]!=0) for (j=1;j<this->NbEl[SegmentID]-1;j++){
      if (j==1){
        x_svg=this->x[SegmentID][0];
        y_svg=this->y[SegmentID][0];
        z_svg=this->z[SegmentID][0];
        d_svg=this->d[SegmentID][0];
      }
      
      x_new=(x_svg+this->x[SegmentID][j]+this->x[SegmentID][j+1])/3;
      y_new=(y_svg+this->y[SegmentID][j]+this->y[SegmentID][j+1])/3;
      z_new=(z_svg+this->z[SegmentID][j]+this->z[SegmentID][j+1])/3;
      d_new=(d_svg+this->d[SegmentID][j]+this->d[SegmentID][j+1])/3;
      
      x_svg=this->x[SegmentID][j];
      y_svg=this->y[SegmentID][j];
      z_svg=this->z[SegmentID][j];
      d_svg=this->d[SegmentID][j];
      
      this->x[SegmentID][j]=x_new;
      this->y[SegmentID][j]=y_new;
      this->z[SegmentID][j]=z_new;
      this->d[SegmentID][j]=d_new;
    }
  }
}


///transform the LDMK_Curves coordinates and diameters from voxels to mm according to the image 2 world properties of RefSF 
void LDMK_Curves::VoxelsToMillimeters(ScalarField * RefSF){
  int i,j;
  float x_new,y_new,z_new;
  float x_mm,y_mm,z_mm,vox_mm,three;
  
  three=3;
  
  x_mm=sqrt(RefSF->Image2World[0][0]*RefSF->Image2World[0][0]+RefSF->Image2World[0][1]*RefSF->Image2World[0][1]+RefSF->Image2World[0][2]*RefSF->Image2World[0][2]);
  y_mm=sqrt(RefSF->Image2World[1][0]*RefSF->Image2World[1][0]+RefSF->Image2World[1][1]*RefSF->Image2World[1][1]+RefSF->Image2World[1][2]*RefSF->Image2World[1][2]);
  z_mm=sqrt(RefSF->Image2World[2][0]*RefSF->Image2World[2][0]+RefSF->Image2World[2][1]*RefSF->Image2World[2][1]+RefSF->Image2World[2][2]*RefSF->Image2World[2][2]);
  vox_mm=(x_mm+y_mm+z_mm)/three;
  
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    x_new=this->x[i][j]*RefSF->Image2World[0][0]+this->y[i][j]*RefSF->Image2World[0][1]+this->z[i][j]*RefSF->Image2World[0][2]+RefSF->Image2World[0][3];
    y_new=this->x[i][j]*RefSF->Image2World[1][0]+this->y[i][j]*RefSF->Image2World[1][1]+this->z[i][j]*RefSF->Image2World[1][2]+RefSF->Image2World[1][3];
    z_new=this->x[i][j]*RefSF->Image2World[2][0]+this->y[i][j]*RefSF->Image2World[2][1]+this->z[i][j]*RefSF->Image2World[2][2]+RefSF->Image2World[2][3];

    this->x[i][j]=x_new;
    this->y[i][j]=y_new;
    this->z[i][j]=z_new;
    this->d[i][j]=this->d[i][j]*vox_mm;
  }

}

///transform the LDMK_Curves coordinates and diameters from mm to voxels according to the image 2 world properties of RefSF 
void LDMK_Curves::MillimetersToVoxels(ScalarField * RefSF){
  int i,j;
  float x_new,y_new,z_new;
  float x_mm,y_mm,z_mm,vox_mm,three;
  
  three=3;

  x_mm=sqrt(RefSF->Image2World[0][0]*RefSF->Image2World[0][0]+RefSF->Image2World[0][1]*RefSF->Image2World[0][1]+RefSF->Image2World[0][2]*RefSF->Image2World[0][2]);
  y_mm=sqrt(RefSF->Image2World[1][0]*RefSF->Image2World[1][0]+RefSF->Image2World[1][1]*RefSF->Image2World[1][1]+RefSF->Image2World[1][2]*RefSF->Image2World[1][2]);
  z_mm=sqrt(RefSF->Image2World[2][0]*RefSF->Image2World[2][0]+RefSF->Image2World[2][1]*RefSF->Image2World[2][1]+RefSF->Image2World[2][2]*RefSF->Image2World[2][2]);
  vox_mm=(x_mm+y_mm+z_mm)/three;
  
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    x_new=this->x[i][j]*RefSF->World2Image[0][0]+this->y[i][j]*RefSF->World2Image[0][1]+this->z[i][j]*RefSF->World2Image[0][2]+RefSF->World2Image[0][3];
    y_new=this->x[i][j]*RefSF->World2Image[1][0]+this->y[i][j]*RefSF->World2Image[1][1]+this->z[i][j]*RefSF->World2Image[1][2]+RefSF->World2Image[1][3];
    z_new=this->x[i][j]*RefSF->World2Image[2][0]+this->y[i][j]*RefSF->World2Image[2][1]+this->z[i][j]*RefSF->World2Image[2][2]+RefSF->World2Image[2][3];

    this->x[i][j]=x_new;
    this->y[i][j]=y_new;
    this->z[i][j]=z_new;
    this->d[i][j]=this->d[i][j]/vox_mm;
  }
}


///Affine transformation (TransfoMat is from the source/template to the target)
void LDMK_Curves::AffineTransfo(float TransfoMat[4][4]){
  int i,j;
  float x_new,y_new,z_new;
  
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    x_new=this->x[i][j]*TransfoMat[0][0]+this->y[i][j]*TransfoMat[0][1]+this->z[i][j]*TransfoMat[0][2]+TransfoMat[0][3];
    y_new=this->x[i][j]*TransfoMat[1][0]+this->y[i][j]*TransfoMat[1][1]+this->z[i][j]*TransfoMat[1][2]+TransfoMat[1][3];
    z_new=this->x[i][j]*TransfoMat[2][0]+this->y[i][j]*TransfoMat[2][1]+this->z[i][j]*TransfoMat[2][2]+TransfoMat[2][3];

    this->x[i][j]=x_new;
    this->y[i][j]=y_new;
    this->z[i][j]=z_new;
  }
  
}

  
///translate a network
void LDMK_Curves::Translate(float DecX,float DecY,float DecZ){
  int i,j;
  
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    this->x[i][j]+=DecX;
    this->y[i][j]+=DecY;
    this->z[i][j]+=DecZ;
  }
  
}


///estimate network ROI
void LDMK_Curves::EstimateROI(float * LowerX,float * LowerY,float * LowerZ,float * WidthX,float * WidthY,float * WidthZ){
  int i,j;
  int InitMade;
  float minX,minY,minZ;
  float maxX,maxY,maxZ;
  
  
  InitMade=0;
    
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    
    if (InitMade==0){
      InitMade=1;
      minX=this->GetX(i,j);        minY=this->GetY(i,j);        minZ=this->GetZ(i,j);
      maxX=this->GetX(i,j);        maxY=this->GetY(i,j);        maxZ=this->GetZ(i,j);
      }
    
    if (this->GetX(i,j)>maxX) maxX=this->GetX(i,j);
    if (this->GetY(i,j)>maxY) maxY=this->GetY(i,j);
    if (this->GetZ(i,j)>maxZ) maxZ=this->GetZ(i,j);
    if (this->GetX(i,j)<minX) minX=this->GetX(i,j);
    if (this->GetY(i,j)<minY) minY=this->GetY(i,j);
    if (this->GetZ(i,j)<minZ) minZ=this->GetZ(i,j);
  }
  
  *LowerX=minX;
  *LowerY=minY;
  *LowerZ=minZ;
  *WidthX=maxX-minX;
  *WidthY=maxY-minY;
  *WidthZ=maxZ-minZ;
}
  
   

//Private function used for the skeletonization  -  Rotates the 3*3*3 scalar field CubeIn and save the result in CubeOut. 
//Options are:
// opt1==1 => swap z <-> z
// opt1==2 => swap z <-> -z
// opt1==3 => swap z <-> x
// opt1==4 => swap z <-> -x
// opt1==5 => swap z <-> y
// opt1==6 => swap z <-> -y
// opt2==1 => Rotation of 0 degrees around z-axis
// opt2==2 => Rotation of 90 degrees around z-axis
// opt2==3 => Rotation of 180 degrees around z-axis
// opt2==4 => Rotation of 270 degrees around z-axis
void LDMK_Curves::NgbhTransform(ScalarField * CubeIn,ScalarField * CubeTemp,ScalarField * CubeOut,int opt1,int opt2){
int i,j,k;

if (opt1<1) opt1=1;
if (opt1>6) opt1=1;
if (opt2<1) opt2=1;
if (opt2>4) opt2=1;

//rotate around z-axis
switch (opt2){
        case 1 : for (i=0;i<3;i++){
                        CubeTemp->P(CubeIn->G(0,0,i),0,0,i); CubeTemp->P(CubeIn->G(1,0,i),1,0,i); CubeTemp->P(CubeIn->G(2,0,i),2,0,i);
                        CubeTemp->P(CubeIn->G(0,1,i),0,1,i); CubeTemp->P(CubeIn->G(1,1,i),1,1,i); CubeTemp->P(CubeIn->G(2,1,i),2,1,i);
                        CubeTemp->P(CubeIn->G(0,2,i),0,2,i); CubeTemp->P(CubeIn->G(1,2,i),1,2,i); CubeTemp->P(CubeIn->G(2,2,i),2,2,i);
                        }
                break;
        case 2 : for (i=0;i<3;i++){
                        CubeTemp->P(CubeIn->G(0,2,i),0,0,i); CubeTemp->P(CubeIn->G(0,1,i),1,0,i); CubeTemp->P(CubeIn->G(0,0,i),2,0,i);
                        CubeTemp->P(CubeIn->G(1,2,i),0,1,i); CubeTemp->P(CubeIn->G(1,1,i),1,1,i); CubeTemp->P(CubeIn->G(1,0,i),2,1,i);
                        CubeTemp->P(CubeIn->G(2,2,i),0,2,i); CubeTemp->P(CubeIn->G(2,1,i),1,2,i); CubeTemp->P(CubeIn->G(2,0,i),2,2,i);
                        }
                break;
        case 3 : for (i=0;i<3;i++){
                        CubeTemp->P(CubeIn->G(2,2,i),0,0,i); CubeTemp->P(CubeIn->G(1,2,i),1,0,i); CubeTemp->P(CubeIn->G(0,2,i),2,0,i);
                        CubeTemp->P(CubeIn->G(2,1,i),0,1,i); CubeTemp->P(CubeIn->G(1,1,i),1,1,i); CubeTemp->P(CubeIn->G(0,1,i),2,1,i);
                        CubeTemp->P(CubeIn->G(2,0,i),0,2,i); CubeTemp->P(CubeIn->G(1,0,i),1,2,i); CubeTemp->P(CubeIn->G(0,0,i),2,2,i);
                        }
                break;
        case 4 : for (i=0;i<3;i++){
                        CubeTemp->P(CubeIn->G(2,0,i),0,0,i); CubeTemp->P(CubeIn->G(2,1,i),1,0,i); CubeTemp->P(CubeIn->G(2,2,i),2,0,i);
                        CubeTemp->P(CubeIn->G(1,0,i),0,1,i); CubeTemp->P(CubeIn->G(1,1,i),1,1,i); CubeTemp->P(CubeIn->G(1,2,i),2,1,i);
                        CubeTemp->P(CubeIn->G(0,0,i),0,2,i); CubeTemp->P(CubeIn->G(0,1,i),1,2,i); CubeTemp->P(CubeIn->G(0,2,i),2,2,i);
                        }
                break;
        }

//swap axes
switch (opt1){
        case 1 : for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) CubeOut->P(CubeTemp->G(k,j,i),k,j,i);
                break;
        case 2 : for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) CubeOut->P(CubeTemp->G(2-k,j,2-i),k,j,i);
                break;
        case 3 : for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) CubeOut->P(CubeTemp->G(k,i,2-j),k,j,i);
                break;
        case 4 : for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) CubeOut->P(CubeTemp->G(k,2-i,j),k,j,i);
                break;
        case 5 : for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) CubeOut->P(CubeTemp->G(i,j,2-k),k,j,i);
                break;
        case 6 : for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) CubeOut->P(CubeTemp->G(2-i,j,k),k,j,i);
                break;
        }

}


//Private function used for the skeletonization  -  To know whether a point is a centerline point, test a 3*3*3 ScalarField with the filters of 
//"A 3d 6 subiteration thinning algorithm for extracting medial lines / Kalman Palagyi / //Pattern recogntion letters 19 (1998)"
// -> 6 filters are implemented (FilterID=1,2,3,4,5,6)
// -> return 1 if we have a clue that the point is not a centerline point and 0 otherwise (different tests are performed)
int LDMK_Curves::IsNotCenterlinePt(ScalarField * CubeIn,int FilterID){
int result;

result=1;

switch (FilterID){
  case 1 : 
    //manage 'fix' points
    if (CubeIn->G(0,0,2)==1){result=0; break;} if (CubeIn->G(1,0,2)==1){result=0; break;} if (CubeIn->G(2,0,2)==1){result=0; break;}
    if (CubeIn->G(0,1,2)==1){result=0; break;} if (CubeIn->G(1,1,2)==1){result=0; break;} if (CubeIn->G(2,1,2)==1){result=0; break;}
    if (CubeIn->G(0,2,2)==1){result=0; break;} if (CubeIn->G(1,2,2)==1){result=0; break;} if (CubeIn->G(2,2,2)==1){result=0; break;}
    if (CubeIn->G(1,1,1)==0){result=0; break;} if (CubeIn->G(1,1,0)==0){result=0; break;}
    //manage the "at least 1 point"
    if (CubeIn->G(0,0,0)==1) break; if (CubeIn->G(0,1,0)==1) break; if (CubeIn->G(0,2,0)==1) break;
    if (CubeIn->G(2,0,0)==1) break; if (CubeIn->G(2,1,0)==1) break; if (CubeIn->G(2,2,0)==1) break;
    if (CubeIn->G(1,0,0)==1) break; if (CubeIn->G(1,2,0)==1) break; if (CubeIn->G(0,0,1)==1) break; 
    if (CubeIn->G(0,1,1)==1) break; if (CubeIn->G(0,2,1)==1) break; if (CubeIn->G(2,0,1)==1) break; 
    if (CubeIn->G(2,1,1)==1) break; if (CubeIn->G(2,2,1)==1) break; if (CubeIn->G(1,0,1)==1) break; 
    if (CubeIn->G(1,2,1)==1) break;
    result=0;
    break;
  case 2 : 
    //manage 'fix' points
    if (CubeIn->G(0,1,2)==1){result=0; break;} if (CubeIn->G(1,1,2)==1){result=0; break;} if (CubeIn->G(2,1,2)==1){result=0; break;}
    if (CubeIn->G(0,2,2)==1){result=0; break;} if (CubeIn->G(1,2,2)==1){result=0; break;} if (CubeIn->G(2,2,2)==1){result=0; break;}
    if (CubeIn->G(1,1,1)==0){result=0; break;} if (CubeIn->G(1,0,1)==0){result=0; break;} if (CubeIn->G(1,1,0)==0){result=0; break;}
    break;
  case 3 : 
    //manage 'fix' points
    if (CubeIn->G(0,2,2)==1){result=0; break;} if (CubeIn->G(1,2,2)==1){result=0; break;} if (CubeIn->G(0,1,2)==1){result=0; break;}
    if (CubeIn->G(1,1,2)==1){result=0; break;} if (CubeIn->G(1,1,1)==0){result=0; break;} if (CubeIn->G(2,1,1)==0){result=0; break;}
    if (CubeIn->G(1,0,1)==0){result=0; break;} if (CubeIn->G(1,1,0)==0){result=0; break;}
    break;
  case 4 : 
    //manage 'fix' points
    if (CubeIn->G(0,0,2)==1){result=0; break;} if (CubeIn->G(1,0,2)==1){result=0; break;} if (CubeIn->G(2,0,2)==0){result=0; break;}
    if (CubeIn->G(0,1,2)==1){result=0; break;} if (CubeIn->G(1,1,2)==1){result=0; break;} if (CubeIn->G(2,1,2)==1){result=0; break;}
    if (CubeIn->G(0,2,2)==1){result=0; break;} if (CubeIn->G(1,2,2)==1){result=0; break;} if (CubeIn->G(2,2,2)==1){result=0; break;}
    if (CubeIn->G(2,0,1)==0){result=0; break;} if (CubeIn->G(1,1,1)==0){result=0; break;} if (CubeIn->G(1,1,0)==0){result=0; break;}
    break;
  case 5 : 
    //manage 'fix' points
    if (CubeIn->G(0,0,2)==1){result=0; break;} if (CubeIn->G(1,0,2)==1){result=0; break;} if (CubeIn->G(2,0,2)==1){result=0; break;}
    if (CubeIn->G(0,1,2)==1){result=0; break;} if (CubeIn->G(1,1,2)==1){result=0; break;} if (CubeIn->G(2,1,2)==1){result=0; break;}
    if (CubeIn->G(0,2,2)==1){result=0; break;} if (CubeIn->G(1,2,2)==1){result=0; break;} if (CubeIn->G(2,2,2)==1){result=0; break;}
    if (CubeIn->G(0,2,1)==1){result=0; break;} if (CubeIn->G(1,2,1)==1){result=0; break;} if (CubeIn->G(2,2,1)==1){result=0; break;}
    if (CubeIn->G(0,2,0)==1){result=0; break;} if (CubeIn->G(1,2,0)==1){result=0; break;} if (CubeIn->G(2,2,0)==1){result=0; break;}
    if (CubeIn->G(1,1,1)==0){result=0; break;} if (CubeIn->G(1,0,0)==0){result=0; break;} if (CubeIn->G(1,1,0)==1){result=0; break;}
    //manage the "at least 1 point"
    if (CubeIn->G(0,1,0)==1) break; if (CubeIn->G(0,0,0)==1) break; if (CubeIn->G(2,1,0)==1) break;
    if (CubeIn->G(2,0,0)==1) break; if (CubeIn->G(2,0,1)==1) break; if (CubeIn->G(1,0,1)==1) break;
    if (CubeIn->G(0,1,1)==1) break; if (CubeIn->G(0,0,1)==1) break; if (CubeIn->G(2,1,1)==1) break; 
    result=0;
    break;
  case 6 : 
    //manage 'fix' points
    if (CubeIn->G(0,0,2)==1){result=0; break;} if (CubeIn->G(1,0,2)==1){result=0; break;} if (CubeIn->G(2,0,2)==1){result=0; break;}
    if (CubeIn->G(0,1,2)==1){result=0; break;} if (CubeIn->G(1,1,2)==1){result=0; break;} if (CubeIn->G(2,1,2)==1){result=0; break;}
    if (CubeIn->G(0,2,2)==1){result=0; break;} if (CubeIn->G(1,2,2)==1){result=0; break;} if (CubeIn->G(2,2,2)==1){result=0; break;}
    if (CubeIn->G(0,2,0)==1){result=0; break;} if (CubeIn->G(0,1,0)==1){result=0; break;} if (CubeIn->G(1,2,0)==1){result=0; break;}
    if (CubeIn->G(1,1,0)==1){result=0; break;} if (CubeIn->G(0,2,1)==1){result=0; break;} if (CubeIn->G(1,2,1)==1){result=0; break;}
    if (CubeIn->G(0,1,1)==1){result=0; break;} if (CubeIn->G(1,1,1)==0){result=0; break;} if (CubeIn->G(1,0,0)==0){result=0; break;}
    if (CubeIn->G(2,1,0)==0){result=0; break;}
    break;
  }

return result;
}


///generate the LDMK_Curves by skeletonizing the img3d (shape=1 / backgroud=0)
///The 18-neighborhood is considered for the shape and the 6-neighborhood is considered for the background
///Remark 1: Radii represent the nearest boundary in mm
///Remark 2: Coordinates are in mm
///Remark 3: Post-treatments are usually necessary to merge the segments (often subdivided into several segments) and useless segments
void LDMK_Curves::Skeletonize(ScalarField * img3d){
  int i,j,k,l;
  int i2,j2,k2;
  int i3,j3,k3;
  int i4,j4,k4;
  ScalarField TestField;
  ScalarField TestField2;
  ScalarField TestField3;
  ScalarField DistanceMap;
  ScalarField DistanceMapMM;
  int TmpCount;
  int GlobalChanges;
  int Direc;
  float minNgbh;
  float tmpFl;
  FILE *dataFile;
  int cpt1,cpt2,testBr;
  int NbNodes;
  int NbLeaves;
  int heapSize,HeapLoc;
  int Px[50000];
  int Py[50000];
  int Pz[50000];
  float epsilon=0.01;
  double NodeGC_X,NodeGC_Y,NodeGC_Z;
  double x_new,y_new,z_new;
  int MaxZ,MinZ;
  char mv3dFileName[256];
  char PrefixTmpFile[256];
  char SuffixTmpFile[256];
  char locString[10];
  
  //STEP 0 : INITS

  //0.1) create temporary 3*3*3 fields
  TestField.CreateVoidField(3,3,3);
  TestField2.CreateVoidField(3,3,3);
  TestField3.CreateVoidField(3,3,3);
  
  //0.2) inititiate the non euclidian distance map (city block distance) and the euclidian distance in mm
  DistanceMapMM.CreateVoidField(img3d->NX,img3d->NY,img3d->NZ);
  Cpt_DistMap(img3d,1,&DistanceMapMM);
  
  DistanceMap.CreateVoidField(img3d->NX,img3d->NY,img3d->NZ);
  
  for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) if (img3d->G(k,j,i)>0.5) DistanceMap.P(10000,k,j,i);
  
  for (i=0;i<img3d->NZ;i++) for (j=0;j<img3d->NY;j++) img3d->P(0,0,j,i);
  for (i=0;i<img3d->NZ;i++) for (j=0;j<img3d->NY;j++) img3d->P(0,img3d->NX-1,j,i);
  for (i=0;i<img3d->NZ;i++) for (k=0;k<img3d->NX;k++) img3d->P(0,k,0,i);
  for (i=0;i<img3d->NZ;i++) for (k=0;k<img3d->NX;k++) img3d->P(0,k,img3d->NY-1,i);
  for (j=0;j<img3d->NY;j++) for (k=0;k<img3d->NX;k++) img3d->P(0,k,j,0);
  for (j=0;j<img3d->NY;j++) for (k=0;k<img3d->NX;k++) img3d->P(0,k,j,img3d->NZ-1);
      
  //0.3) generate a string for a temporary mv3d file
  static const char alphanum[] = "0123456789";
  
  for (int i = 0; i < 10; ++i) locString[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
  
  strcpy(PrefixTmpFile,"tempFile");
  strcpy(SuffixTmpFile,".mv3d");
  strcat(PrefixTmpFile,locString);
  strcat(PrefixTmpFile,SuffixTmpFile);
  
  strcpy(mv3dFileName,PrefixTmpFile);
  
  //STEP 1 : SURFACE SKELETONIZATION
  printf("Surface Skeletonization\n");
  
  GlobalChanges=1; 
  Direc=-1;
  while (GlobalChanges!=0){
    //label voxels at the edge of the surface to skeletonize
    GlobalChanges=0;
  
    Direc++; if (Direc==6) Direc=0;    //Direc -> erosion direction
    for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) if (fabs(img3d->G(k,j,i)-1)<epsilon){
      if (Direc==0) if (fabs(img3d->G(k,j,i+1))<epsilon) img3d->P(2,k,j,i);
      if (Direc==1) if (fabs(img3d->G(k,j+1,i))<epsilon) img3d->P(2,k,j,i);
      if (Direc==2) if (fabs(img3d->G(k+1,j,i))<epsilon) img3d->P(2,k,j,i);
      if (Direc==3) if (fabs(img3d->G(k,j,i-1))<epsilon) img3d->P(2,k,j,i);
      if (Direc==4) if (fabs(img3d->G(k,j-1,i))<epsilon) img3d->P(2,k,j,i);
      if (Direc==5) if (fabs(img3d->G(k-1,j,i))<epsilon) img3d->P(2,k,j,i);
      if (fabs(img3d->G(k,j,i)-2)<epsilon) GlobalChanges=1;
      }
    
    //set the voxel to 0 if the local point is not a centerline point
    for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) if (fabs(img3d->G(k,j,i)-2)<epsilon){
      //update distance map at local point
      minNgbh=DistanceMap.G(k+1,j,i);
      if (DistanceMap.G(k-1,j,i)<minNgbh) minNgbh=DistanceMap.G(k-1,j,i);
      if (DistanceMap.G(k,j+1,i)<minNgbh) minNgbh=DistanceMap.G(k,j+1,i);
      if (DistanceMap.G(k,j-1,i)<minNgbh) minNgbh=DistanceMap.G(k,j-1,i);
      if (DistanceMap.G(k,j,i+1)<minNgbh) minNgbh=DistanceMap.G(k,j,i+1);
      if (DistanceMap.G(k,j,i-1)<minNgbh) minNgbh=DistanceMap.G(k,j,i-1);
      DistanceMap.P(minNgbh+1,k,j,i);
      
      //copy local ngbh in TestField
      for (i2=-1;i2<2;i2++)  for (j2=-1;j2<2;j2++) for (k2=-1;k2<2;k2++)
        TestField.P(static_cast<float>(img3d->G(k+k2,j+j2,i+i2)>0.5),k2+1,j2+1,i2+1);
  
  
      //test whether we have a centerline point
      if (Direc==0) for (cpt1=1;cpt1<=4;cpt1++){
          this->NgbhTransform(&TestField,&TestField3,&TestField2,1,cpt1);
          for (cpt2=1;cpt2<=6;cpt2++){testBr=this->IsNotCenterlinePt(&TestField2,cpt2);if (testBr==1) {img3d->P(0,k,j,i) ; break;}}
        }
        
      if (Direc==3) for (cpt1=1;cpt1<=4;cpt1++){
          this->NgbhTransform(&TestField,&TestField3,&TestField2,2,cpt1);
          for (cpt2=1;cpt2<=6;cpt2++){ testBr=this->IsNotCenterlinePt(&TestField2,cpt2); if (testBr==1) {img3d->P(0,k,j,i) ; break;}}
        }
      
      if (Direc==2) for (cpt1=1;cpt1<=4;cpt1++){
          this->NgbhTransform(&TestField,&TestField3,&TestField2,5,cpt1);
          for (cpt2=1;cpt2<=6;cpt2++){ testBr=this->IsNotCenterlinePt(&TestField2,cpt2); if (testBr==1) {img3d->P(0,k,j,i) ; break;}}
        }
      
      if (Direc==5) for (cpt1=1;cpt1<=4;cpt1++){
          this->NgbhTransform(&TestField,&TestField3,&TestField2,6,cpt1);
          for (cpt2=1;cpt2<=6;cpt2++){ testBr=this->IsNotCenterlinePt(&TestField2,cpt2); if (testBr==1) {img3d->P(0,k,j,i) ; break;}}
        }
  
      if (Direc==1) for (cpt1=1;cpt1<=4;cpt1++){
          this->NgbhTransform(&TestField,&TestField3,&TestField2,3,cpt1);
          for (cpt2=1;cpt2<=6;cpt2++){ testBr=this->IsNotCenterlinePt(&TestField2,cpt2); if (testBr==1) {img3d->P(0,k,j,i) ; break;}}
        }
  
      if (Direc==4) for (cpt1=1;cpt1<=4;cpt1++){
          this->NgbhTransform(&TestField,&TestField3,&TestField2,4,cpt1);
          for (cpt2=1;cpt2<=6;cpt2++){ testBr=this->IsNotCenterlinePt(&TestField2,cpt2); if (testBr==1) {img3d->P(0,k,j,i) ; break;}}
        }
      }
  
    }
  
  
  //copy the euclidian distances in DistanceMap
  for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) DistanceMap.P(DistanceMapMM.G(k,j,i),k,j,i);
  
  
  //special test to detect 2D graphs
  
  MaxZ=0;
  MinZ=img3d->NZ-1;
  
  for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) if (img3d->G(k,j,i)>0.5) {
    if (i>MaxZ) MaxZ=i;
    if (i<MinZ) MinZ=i;
  }
  
  if (MaxZ-MinZ<3){
    cout << "2D or almost 2D shape is skeletonised -> our technique to estimate skeleton diameters is only appropriate for 3D images -> all diameters are set to 3" << endl;
    for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++)
          DistanceMap.P(3,k,j,i);
  }
      
      
  //STEP 2 : GRAPH CREATION
  printf("Graph creation\n");
  
  //2.1 - open the temporary  file and write its header
  
  dataFile=fopen(mv3dFileName,"w");
  
  fprintf(dataFile,"# MicroVisu3D file\n");
  fprintf(dataFile,"# Number of lines   %d\n",1);
  fprintf(dataFile,"# Number of points  %d\n",1);
  fprintf(dataFile,"# Number of inter.  %d\n",1);
  fprintf(dataFile,"#\n");
  fprintf(dataFile,"# No                x                y                z                d\n");
  fprintf(dataFile,"#\n");
  TmpCount=0;
  
  //2.2 - clean-up the image and define the list of nodes and leaves
  
  //2.2.1 - clean-up the image
  for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++){
    if (img3d->G(k,j,i)>0.5) 
      img3d->P(1,k,j,i);
    else
      img3d->P(0,k,j,i);
  }
  
  //put to zero the points at the domain boundary
  for(i=0;i<img3d->NZ;i++) for(j=0;j<img3d->NY;j++){ img3d->P(0,0,j,i);  img3d->P(0,img3d->NX-1,j,i); }
  for(i=0;i<img3d->NZ;i++) for(k=0;k<img3d->NX;k++){ img3d->P(0,k,0,i);  img3d->P(0,k,img3d->NY-1,i); }
  for(j=0;j<img3d->NY;j++) for(k=0;k<img3d->NX;k++){ img3d->P(0,k,j,0);  img3d->P(0,k,j,img3d->NZ-1); }
  
  
  //2.2.2 Detect the nodes / leaves / points within the segments
    
  NbNodes=0;
  NbLeaves=0;
  
  
  for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) if (fabs(img3d->G(k,j,i)-1)<epsilon){
    
    TmpCount=0;
    for(i2=-1;i2<=1;i2++) for(j2=-1;j2<=1;j2++) for(k2=-1;k2<=1;k2++) if (img3d->G(k+k2,j+j2,i+i2)>0.5) TmpCount++;
  
    //delete isolated points
    if (TmpCount==1) img3d->P(0,k,j,i);
    
    //detect the leaves (TmpCount-1 est le nombre de 26-voisins)
    if (TmpCount==2){
      img3d->P(10,k,j,i);
      NbLeaves++;
      }
    
    //find the points within segments
    if (TmpCount==3){
      img3d->P(5,k,j,i);
      NbLeaves++;
      }
      
    //detect the nodes
    if (TmpCount>=4){
      img3d->P(20,k,j,i);
      NbNodes++;
      }
    }
    
  //img3d->Write("toto1.nii");
  
  //2.2.3 - Extract the nodes and connected sets of 'nodes'
  TmpCount=0;
  for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) if (fabs(img3d->G(k,j,i)-20)<epsilon){
    i4=i; j4=j; k4=k;
    
      //2.2.3.1 : in case the node is not alone, construct a "connected set of nodes"
      Px[0]=k4; Py[0]=j4; Pz[0]=i4; 
      heapSize=1;
      HeapLoc=0;
      img3d->P(0,Px[HeapLoc],Py[HeapLoc],Pz[HeapLoc]);
      while ((HeapLoc<heapSize)&&(heapSize<49999)){
        for(i2=-1;i2<=1;i2++) for(j2=-1;j2<=1;j2++) for(k2=-1;k2<=1;k2++){
          if ((fabs(img3d->G(Px[HeapLoc]+k2,Py[HeapLoc]+j2,Pz[HeapLoc]+i2)-20)<epsilon)&&(heapSize<49999)){ 
            img3d->P(0,Px[HeapLoc]+k2,Py[HeapLoc]+j2,Pz[HeapLoc]+i2);
            Px[heapSize]=Px[HeapLoc]+k2;
            Py[heapSize]=Py[HeapLoc]+j2;
            Pz[heapSize]=Pz[HeapLoc]+i2;
            heapSize++;
            }
          }
        HeapLoc++;
        }
      if (heapSize>49998) printf("A heap is too large\n");
      
      
      //2.2.3.2 : compute the center of the connected set of nodes
      
      NodeGC_X=0; NodeGC_Y=0; NodeGC_Z=0;
      
      for (HeapLoc=0;HeapLoc<heapSize;HeapLoc++){
        NodeGC_X+=static_cast<double>(Px[HeapLoc]);
        NodeGC_Y+=static_cast<double>(Py[HeapLoc]);
        NodeGC_Z+=static_cast<double>(Pz[HeapLoc]);
      }
      
      NodeGC_X/=static_cast<double>(heapSize);
      NodeGC_Y/=static_cast<double>(heapSize);
      NodeGC_Z/=static_cast<double>(heapSize);
      
      //2.2.3.3 : detect the leaves or segment points connected to the 'node'/'connected set of nodes'
      for (HeapLoc=0;HeapLoc<heapSize;HeapLoc++){
        for(i2=-1;i2<=1;i2++) for(j2=-1;j2<=1;j2++) for(k2=-1;k2<=1;k2++){
          if ((fabs(img3d->G(Px[HeapLoc]+k2,Py[HeapLoc]+j2,Pz[HeapLoc]+i2)-5)<epsilon)||   
          (fabs(img3d->G(Px[HeapLoc]+k2,Py[HeapLoc]+j2,Pz[HeapLoc]+i2)-10)<epsilon)){
            x_new=static_cast<double>(NodeGC_X)*static_cast<double>(img3d->Image2World[0][0])+static_cast<double>(NodeGC_Y)*static_cast<double>(img3d->Image2World[0][1])+static_cast<double>(NodeGC_Z)*static_cast<double>(img3d->Image2World[0][2])+static_cast<double>(img3d->Image2World[0][3]);
            y_new=static_cast<double>(NodeGC_X)*static_cast<double>(img3d->Image2World[1][0])+static_cast<double>(NodeGC_Y)*static_cast<double>(img3d->Image2World[1][1])+static_cast<double>(NodeGC_Z)*static_cast<double>(img3d->Image2World[1][2])+static_cast<double>(img3d->Image2World[1][3]);
            z_new=static_cast<double>(NodeGC_X)*static_cast<double>(img3d->Image2World[2][0])+static_cast<double>(NodeGC_Y)*static_cast<double>(img3d->Image2World[2][1])+static_cast<double>(NodeGC_Z)*static_cast<double>(img3d->Image2World[2][2])+static_cast<double>(img3d->Image2World[2][3]);
            fprintf(dataFile,"%d        %lf        %lf        %lf        %lf\n",TmpCount,x_new,y_new,z_new,static_cast<double>(DistanceMap.G(NodeGC_X,NodeGC_Y,NodeGC_Z)));
            x_new=static_cast<double>(Px[HeapLoc]+k2)*static_cast<double>(img3d->Image2World[0][0])+static_cast<double>(Py[HeapLoc]+j2)*static_cast<double>(img3d->Image2World[0][1])+static_cast<double>(Pz[HeapLoc]+i2)*static_cast<double>(img3d->Image2World[0][2])+static_cast<double>(img3d->Image2World[0][3]);
            y_new=static_cast<double>(Px[HeapLoc]+k2)*static_cast<double>(img3d->Image2World[1][0])+static_cast<double>(Py[HeapLoc]+j2)*static_cast<double>(img3d->Image2World[1][1])+static_cast<double>(Pz[HeapLoc]+i2)*static_cast<double>(img3d->Image2World[1][2])+static_cast<double>(img3d->Image2World[1][3]);
            z_new=static_cast<double>(Px[HeapLoc]+k2)*static_cast<double>(img3d->Image2World[2][0])+static_cast<double>(Py[HeapLoc]+j2)*static_cast<double>(img3d->Image2World[2][1])+static_cast<double>(Pz[HeapLoc]+i2)*static_cast<double>(img3d->Image2World[2][2])+static_cast<double>(img3d->Image2World[2][3]);
            fprintf(dataFile,"%d        %lf        %lf        %lf        %lf\n",TmpCount,x_new,y_new,z_new,static_cast<double>(DistanceMap.G(Px[HeapLoc]+k2,Py[HeapLoc]+j2,Pz[HeapLoc]+i2)));
            TmpCount++;
            fprintf(dataFile,"\n");
            if (fabs(img3d->G(Px[HeapLoc]+k2,Py[HeapLoc]+j2,Pz[HeapLoc]+i2)-5)<epsilon)
              img3d->P(10,Px[HeapLoc]+k2,Py[HeapLoc]+j2,Pz[HeapLoc]+i2);
          }
        }
      }
      
      
  }
  
  //2.2.3 - Extract the segments
  
  for(i=1;i<img3d->NZ-1;i++) for(j=1;j<img3d->NY-1;j++) for(k=1;k<img3d->NX-1;k++) if (fabs(img3d->G(k,j,i)-10)<epsilon){
      i4=i; j4=j; k4=k;
      
      //2.2.3.1 : follow the segment
      tmpFl=1;
      while ((tmpFl>0)&&(tmpFl<9)){
              
        x_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[0][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[0][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[0][2])+static_cast<double>(img3d->Image2World[0][3]);
        y_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[1][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[1][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[1][2])+static_cast<double>(img3d->Image2World[1][3]);
        z_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[2][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[2][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[2][2])+static_cast<double>(img3d->Image2World[2][3]);
        fprintf(dataFile,"%d        %lf        %lf        %lf        %lf\n",TmpCount,x_new,y_new,z_new,static_cast<double>(DistanceMap.G(k4,j4,i4)));
        img3d->P(0,k4,j4,i4);
  
        tmpFl=0;
        for(i2=-1;i2<=1;i2++) for(j2=-1;j2<=1;j2++) for(k2=-1;k2<=1;k2++) if ((img3d->G(k4+k2,j4+j2,i4+i2)>tmpFl)){
          i3=i4+i2; j3=j4+j2; k3=k4+k2;
          tmpFl=img3d->G(k4+k2,j4+j2,i4+i2);
          }
        
        i4=i3; j4=j3; k4=k3;
        }
      
      //2.2.3.2 : a leaf ends the segment
      if(fabs(tmpFl-10)<epsilon){
        x_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[0][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[0][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[0][2])+static_cast<double>(img3d->Image2World[0][3]);
        y_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[1][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[1][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[1][2])+static_cast<double>(img3d->Image2World[1][3]);
        z_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[2][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[2][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[2][2])+static_cast<double>(img3d->Image2World[2][3]);
        fprintf(dataFile,"%d        %lf        %lf        %lf        %lf\n",TmpCount,x_new,y_new,z_new,static_cast<double>(DistanceMap.G(k4,j4,i4)));
        img3d->P(0,k4,j4,i4);
        }
  
      //2.2.3.3 : a node ends the segment
      if(fabs(tmpFl-20)<epsilon){
        x_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[0][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[0][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[0][2])+static_cast<double>(img3d->Image2World[0][3]);
        y_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[1][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[1][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[1][2])+static_cast<double>(img3d->Image2World[1][3]);
        z_new=static_cast<double>(k4)*static_cast<double>(img3d->Image2World[2][0])+static_cast<double>(j4)*static_cast<double>(img3d->Image2World[2][1])+static_cast<double>(i4)*static_cast<double>(img3d->Image2World[2][2])+static_cast<double>(img3d->Image2World[2][3]);
        fprintf(dataFile,"%d        %lf        %lf        %lf        %lf\n",TmpCount,x_new,y_new,z_new,static_cast<double>(DistanceMap.G(k4,j4,i4)));
        img3d->P(0,k4,j4,i4);
        cout << "Unexpected node when extracting a segment" << endl;
      }
      
      //2.2.3.4 : finalize the segment
      TmpCount++;
      fprintf(dataFile,"\n");    
  
  }
  
  
  //3 close the temporary file
  fclose(dataFile);

  cout << mv3dFileName << " is saved" << endl;

  
  //4 load the network in the temporary file
  this->Read(mv3dFileName);
}


  
  

///Reconstruct the 3D volume form the network (shape=1 / backgroud=0)
///warning: all intensities of RefImageDomain will be modified so that they are the reconstructed volume in the end of the computations
void LDMK_Curves::Generate3DVolume(char * OutputImageName,ScalarField * RefImageDomain){
  float size_x,size_y,size_z;
  float x_vox,y_vox,z_vox;
  float d_vox_x,d_vox_y,d_vox_z;
  
  int x_vox_i,y_vox_i,z_vox_i;
  float d_vox_x_i,d_vox_y_i,d_vox_z_i;
  
  int i,j;
  int x,y,z;
  
  RefImageDomain->PutToAllVoxels(0);
  
  size_x=sqrt((RefImageDomain->World2Image[0][0]*RefImageDomain->World2Image[0][0])+(RefImageDomain->World2Image[1][0]*RefImageDomain->World2Image[1][0])+(RefImageDomain->World2Image[2][0]*RefImageDomain->World2Image[2][0]));
  size_y=sqrt((RefImageDomain->World2Image[0][1]*RefImageDomain->World2Image[0][1])+(RefImageDomain->World2Image[1][1]*RefImageDomain->World2Image[1][1])+(RefImageDomain->World2Image[2][1]*RefImageDomain->World2Image[2][1]));
  size_z=sqrt((RefImageDomain->World2Image[0][2]*RefImageDomain->World2Image[0][2])+(RefImageDomain->World2Image[1][2]*RefImageDomain->World2Image[1][2])+(RefImageDomain->World2Image[2][2]*RefImageDomain->World2Image[2][2]));
  
  //cout << size_x << " " <<  size_y << " " <<  size_z << endl;
  
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    x_vox=static_cast<float>(this->x[i][j])*RefImageDomain->World2Image[0][0]+static_cast<double>(this->y[i][j])*RefImageDomain->World2Image[0][1]+static_cast<double>(this->z[i][j])*RefImageDomain->World2Image[0][2]+RefImageDomain->World2Image[0][3];
    y_vox=static_cast<float>(this->x[i][j])*RefImageDomain->World2Image[1][0]+static_cast<double>(this->y[i][j])*RefImageDomain->World2Image[1][1]+static_cast<double>(this->z[i][j])*RefImageDomain->World2Image[1][2]+RefImageDomain->World2Image[1][3];
    z_vox=static_cast<float>(this->x[i][j])*RefImageDomain->World2Image[2][0]+static_cast<double>(this->y[i][j])*RefImageDomain->World2Image[2][1]+static_cast<double>(this->z[i][j])*RefImageDomain->World2Image[2][2]+RefImageDomain->World2Image[2][3];

    d_vox_x=size_x*this->d[i][j];
    d_vox_y=size_y*this->d[i][j];
    d_vox_z=size_z*this->d[i][j];
    
    if (d_vox_x<0.99) d_vox_x=0.99;
    if (d_vox_y<0.99) d_vox_y=0.99;
    if (d_vox_z<0.99) d_vox_z=0.99;
    
    x_vox_i=static_cast<int>(x_vox+0.5);
    y_vox_i=static_cast<int>(y_vox+0.5);
    z_vox_i=static_cast<int>(z_vox+0.5);
    
    d_vox_x_i=static_cast<int>(d_vox_x+0.5);
    d_vox_y_i=static_cast<int>(d_vox_y+0.5);
    d_vox_z_i=static_cast<int>(d_vox_z+0.5);
    
    //cout << this->d[i][j] << " -> " << d_vox_x_i << " " <<  d_vox_y_i << " " <<  d_vox_z_i << endl;
    
    for (z=z_vox_i-d_vox_z_i ; z<=z_vox_i+d_vox_z_i ; z++) 
      for (y=y_vox_i-d_vox_y_i ; y<=y_vox_i+d_vox_y_i ; y++)
        for (x=x_vox_i-d_vox_x_i ; x<=x_vox_i+d_vox_x_i ; x++)
          RefImageDomain->P(1,x,y,z);
  }

  RefImageDomain->Write(OutputImageName);
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                              9: B Spline basis related to a vector field
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///9.1: constructor and destructor

///constructor
BS_GlobalConvolver::BS_GlobalConvolver(void){
  this->NbCtPtsX=0;
  this->NbCtPtsY=0;
  this->NbCtPtsZ=0;
  this->NbCtPtsT=0;
  this->StepSize=0;
  
  this->NodesWeigther_Type=0;
  this->NodesWeigther_NbNodes=0;
}

///destructor
BS_GlobalConvolver::~BS_GlobalConvolver(void){
  this->NbCtPtsX=0;
  this->NbCtPtsY=0;
  this->NbCtPtsZ=0;
  this->NbCtPtsT=0;
  this->StepSize=0;

  this->NodesWeigther_Type=0;
  this->NodesWeigther_NbNodes=0;
}

///9.2: initialization of the class

#ifdef COMPILE_WITH_OPENMP

///put a the same value at every points of the scalar field
void BS_GlobalConvolver::LinkWithVecField(VectorField * LinkedVF,int GridStepSize)
{
  //set the grid step size
  this->StepSize=GridStepSize;
  
  //link the reference vecto field to the spline basis   -> so we can make this->RefVecField->G(0,x,y,z) and stuffs like that
  this->RefVecField=LinkedVF;
  
  //initiate the grid size
   this->NbCtPtsX=this->RefVecField->NX/this->StepSize;
   this->NbCtPtsY=this->RefVecField->NY/this->StepSize;
   if (this->RefVecField->NZ==1) 
     this->NbCtPtsZ=1;  //2D image
   else
     this->NbCtPtsZ=this->RefVecField->NZ/this->StepSize; //3D image
   this->NbCtPtsT=this->RefVecField->NT;
   
   //initiate the residual field
   this->ResidualField.CreateVoidField(this->RefVecField->NX,this->RefVecField->NY,this->RefVecField->NZ,this->RefVecField->NT);
   
   //initiate the grid of vectors
    this->VecAtControlPts.CreateVoidField(this->NbCtPtsX,this->NbCtPtsY,this->NbCtPtsZ,this->NbCtPtsT);
    this->TmpVecAtControlPts.CreateVoidField(this->NbCtPtsX,this->NbCtPtsY,this->NbCtPtsZ,this->NbCtPtsT);  //to parallelize on the t axis !!!

   //initiate the spline convolver
   this->SplineConvolver.InitiateSplineConvolver(this->RefVecField->NX,this->RefVecField->NY,this->RefVecField->NZ,this->StepSize);
}

#else

///put a the same value at every points of the scalar field
void BS_GlobalConvolver::LinkWithVecField(VectorField * LinkedVF,int GridStepSize)
{
  //set the grid step size
  this->StepSize=GridStepSize;
  
  //link the reference vecto field to the spline basis   -> so we can make this->RefVecField->G(0,x,y,z) and stuffs like that
  this->RefVecField=LinkedVF;
  
  //initiate the grid size
   this->NbCtPtsX=this->RefVecField->NX/this->StepSize;
   this->NbCtPtsY=this->RefVecField->NY/this->StepSize;
   if (this->RefVecField->NZ==1) 
     this->NbCtPtsZ=1;  //2D image
   else
     this->NbCtPtsZ=this->RefVecField->NZ/this->StepSize; //3D image
   this->NbCtPtsT=this->RefVecField->NT;
   
   //initiate the residual field
   this->ResidualField.CreateVoidField(this->RefVecField->NX,this->RefVecField->NY,this->RefVecField->NZ,this->RefVecField->NT);
   
   //initiate the grid of vectors
    this->VecAtControlPts.CreateVoidField(this->NbCtPtsX,this->NbCtPtsY,this->NbCtPtsZ,this->NbCtPtsT);
    this->TmpVecAtControlPts.CreateVoidField(this->NbCtPtsX,this->NbCtPtsY,this->NbCtPtsZ);  //to save memory  !!!

   //initiate the spline convolver
   this->SplineConvolver.InitiateSplineConvolver(this->RefVecField->NX,this->RefVecField->NY,this->RefVecField->NZ,this->StepSize);
}

#endif




///load a matrix which will weight the terms of the basis with a matrix-vector multiplication
///Weigths of the basis are considered as a vector with (first x) / (then y) / (then z) / (then t) / (then direction)
///Basis terms are weighted when using 'NodesWeight'
void  BS_GlobalConvolver::Load_MatrixNodesWeigther(char * Mfile){
  
  //1) read the matrix
  this->NodesWeigther_MatrixStyle.Read(Mfile);
  
  //2) evaluate the number of nodes from the matrix
  this->NodesWeigther_NbNodes=this->NodesWeigther_MatrixStyle.NI/3;  //each node has 3 DOF because it contains a 3D vector
  
  //3) check that the size is OK with the basis
  cout << "Nb nodes found in " << Mfile << ": " << this->NodesWeigther_NbNodes << endl;
  cout << "Nb nodes in the basis: " << this->NbCtPtsX*this->NbCtPtsY*this->NbCtPtsZ << endl;

  if (this->NodesWeigther_NbNodes!=this->NbCtPtsX*this->NbCtPtsY*this->NbCtPtsZ){
    cout << "-> Basis to project the vector field must be defined before loading the weighting matrix OR no coherency" << endl;
    cout << "         => " << this->NodesWeigther_NbNodes << " nodes vs " << this->NbCtPtsX*this->NbCtPtsY*this->NbCtPtsZ << "nodes" << endl;
    cout << "         => nothing is defined to weight the nodes (use the option -GrdStep_No_M to generate a template)" << endl;
    this->NodesWeigther_Type=0;
  }
  else{
    cout << "-> Coherency OK :-)" << endl;
    this->NodesWeigther_Type=1;
    }
    
}


///same as Load_MatrixNodesWeigther but the matrix does not come from a file but is copied by reference
void BS_GlobalConvolver::Load_MatrixNodesWeigther_usingCopyByRef(SquareMatrix * RefMatrix, int Verbose){
    
    //1) copy the address of the matrix
    this->NodesWeigther_MatrixStyle.CopyByRef(RefMatrix);
    
    //2) evaluate the number of nodes from the matrix
    this->NodesWeigther_NbNodes=this->NodesWeigther_MatrixStyle.NI/3;  //each node has 3 DOF because it contains a 3D vector
    
    //3) check that the size is OK with the basis
    if (Verbose==1) cout << "Nb nodes in M: " << this->NodesWeigther_NbNodes << endl;
    if (Verbose==1) cout << "Nb nodes in the basis: " << this->NbCtPtsX*this->NbCtPtsY*this->NbCtPtsZ << endl;
    
    if (this->NodesWeigther_NbNodes!=this->NbCtPtsX*this->NbCtPtsY*this->NbCtPtsZ){
        cout << "-> Basis to project the vector field must be defined before loading the weighting matrix OR no coherency" << endl;
        cout << "         => " << this->NodesWeigther_NbNodes << " nodes vs " << this->NbCtPtsX*this->NbCtPtsY*this->NbCtPtsZ << "nodes" << endl;
        cout << "         => nothing is defined to weight the nodes (use the option -GrdStep_No_M to generate a template)" << endl;
        this->NodesWeigther_Type=0;
    }
    else{
        if (Verbose==1) cout << "-> Coherency OK :-)" << endl;
        this->NodesWeigther_Type=1;
    }
}



///initiate the convolver
///-> link a bspline basis with the VECTOR FIELD 'LinkedVF' which will be smoothed by the convolver
///-> define the grid step size of a bspline basis on which the VF will be projected
///-> Either copy by reference or load in an ascii file the square matrix (linear operator) which smoothes 
///     the weights of the basis nodes (must be coherent with the VF and grid step size)
void BS_GlobalConvolver::InitiateConvolver(VectorField * LinkedVF,int GridStepSize,char * Mfile){
  this->LinkWithVecField(LinkedVF,GridStepSize);
  this->Load_MatrixNodesWeigther(Mfile); 
}
void BS_GlobalConvolver::InitiateConvolver(VectorField * LinkedVF,int GridStepSize,SquareMatrix * RefMatrix){
  this->LinkWithVecField(LinkedVF,GridStepSize);
  this->Load_MatrixNodesWeigther_usingCopyByRef(RefMatrix);
}
void BS_GlobalConvolver::InitiateConvolver(VectorField * LinkedVF,int GridStepSize){
  this->LinkWithVecField(LinkedVF,GridStepSize);
}

///9.3: smoothing stuffs

///compress the information of RefVectorField at the control points VecAtControlPts
void  BS_GlobalConvolver::RefVecField_2_VecAtControlPts(){
 int x,y,z,t,d;
 
 for (d=0;d<3;d++) for(t=0;t<this->NbCtPtsT;t++){
   for(z=0;z<this->RefVecField->NZ;z++) for(y=0;y<this->RefVecField->NY;y++) for(x=0;x<this->RefVecField->NX;x++)
     this->SplineConvolver.P(this->RefVecField->G(d,x,y,z,t),x,y,z);
   
   this->SplineConvolver.Convolution();
   
   for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++)
    this->VecAtControlPts.P(this->SplineConvolver.G(x*this->StepSize,y*this->StepSize,z*this->StepSize),d,x,y,z,t);
 }
}

///uncompress the information at the control points VecAtControlPts in RefVectorField 
void  BS_GlobalConvolver::VecAtControlPts_2_RefVecField(){
 int x,y,z,t,d;
 float SZ3;

 for (d=0;d<3;d++) for(t=0;t<this->NbCtPtsT;t++){
   
   for(z=0;z<this->RefVecField->NZ;z++) for(y=0;y<this->RefVecField->NY;y++) for(x=0;x<this->RefVecField->NX;x++)
     this->SplineConvolver.P(0,x,y,z);
   
   SZ3=static_cast<float>(this->StepSize*this->StepSize*this->StepSize);
   
   if (this->RefVecField->NZ==1) SZ3=static_cast<float>(this->StepSize*this->StepSize);
   
   
   for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++)
     this->SplineConvolver.P(this->VecAtControlPts.G(d,x,y,z,t)*SZ3,x*this->StepSize,y*this->StepSize,z*this->StepSize);
   
   this->SplineConvolver.Convolution();
   
   for(z=0;z<this->RefVecField->NZ;z++) for(y=0;y<this->RefVecField->NY;y++) for(x=0;x<this->RefVecField->NX;x++)
     this->RefVecField->P(this->SplineConvolver.G(x,y,z),d,x,y,z,t);
   
 }

}



///add the values of other control points to the local ones   (note that no test is done to check the size consistency between the grids)
void  BS_GlobalConvolver::Add_VecAtControlPts(BS_GlobalConvolver * OtherBasis){
  int x,y,z,t,d;
  for (d=0;d<3;d++) for(t=0;t<this->NbCtPtsT;t++) for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++)
    this->VecAtControlPts.Add(OtherBasis->G(d,x,y,z,t),d,x,y,z,t);
}

///multiply the control points by the coef MultCoef
void  BS_GlobalConvolver::MultCoef(float MultCoef){
  int x,y,z,t,d;
  for (d=0;d<3;d++) for(t=0;t<this->NbCtPtsT;t++) for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++)
    this->VecAtControlPts.P(this->VecAtControlPts.G(d,x,y,z,t)*MultCoef,d,x,y,z,t);
}



#ifdef COMPILE_WITH_OPENMP


///weight the terms of the basis with the strategy predefined in 'Load_MatrixNodesWeigther'
///Nothing is made if no strategy is defined
void  BS_GlobalConvolver::NodesWeight(){
  int x,y,z,t,d;
  int x2,y2,z2,d2;
  int i,j;
  float tmpFl;
  
    if (this->NodesWeigther_Type==1){  //1) matrix style weigth
      //BEGIN FORK FOR THREADS
      #pragma omp parallel default(shared) private(t,d,z,y,x,i,j,d2,z2,y2,x2,tmpFl) 
      {  
        #pragma omp for
        for(t=0;t<this->NbCtPtsT;t++){
            //1.1) S1
            i=0;
            for (d=0;d<3;d++)  for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++){
                tmpFl=0;
                j=0;
                for (d2=0;d2<3;d2++)  for(z2=0;z2<this->NbCtPtsZ;z2++) for(y2=0;y2<this->NbCtPtsY;y2++) for(x2=0;x2<this->NbCtPtsX;x2++){
                    tmpFl+=this->VecAtControlPts.G(d2,x2,y2,z2,t)*this->NodesWeigther_MatrixStyle.G(i,j);
                    j++;
                }
                this->TmpVecAtControlPts.P(tmpFl,d,x,y,z,t);
                i++;
            }
            
            //1.2) S2
            for (d=0;d<3;d++)  for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++)
                this->VecAtControlPts.P(this->TmpVecAtControlPts.G(d,x,y,z,t),d,x,y,z,t);
        }
      
      
      //END FORK FOR THREADS
      }
    
    }
  else{  //2) no weigth  (identity matrix)
    //cout << "Nothing was defined to weight the node vectors"  << endl;
  }
}


#else
 
///weight the terms of the basis with the strategy predefined in 'Load_MatrixNodesWeigther'
///Nothing is made if no strategy is defined
void  BS_GlobalConvolver::NodesWeight(){
  int x,y,z,t,d;
  int x2,y2,z2,d2;
  int i,j;
  float tmpFl;  
  
    if (this->NodesWeigther_Type==1){  //1) matrix style weigth
        for(t=0;t<this->NbCtPtsT;t++){
            //1.1) S1
            i=0;
            for (d=0;d<3;d++)  for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++){
                tmpFl=0;
                j=0;
                for (d2=0;d2<3;d2++)  for(z2=0;z2<this->NbCtPtsZ;z2++) for(y2=0;y2<this->NbCtPtsY;y2++) for(x2=0;x2<this->NbCtPtsX;x2++){
                    tmpFl+=this->VecAtControlPts.G(d2,x2,y2,z2,t)*this->NodesWeigther_MatrixStyle.G(i,j);
                    j++;
                }
                this->TmpVecAtControlPts.P(tmpFl,d,x,y,z);
                i++;
            }
            
            //1.2) S2
            for (d=0;d<3;d++)  for(z=0;z<this->NbCtPtsZ;z++) for(y=0;y<this->NbCtPtsY;y++) for(x=0;x<this->NbCtPtsX;x++)
                this->VecAtControlPts.P(this->TmpVecAtControlPts.G(d,x,y,z),d,x,y,z,t);
        }  
    }
  else{  //2) no weigth  (identity matrix)
    cout << "Nothing was defined to weight the node vectors"  << endl;
  }
}

#endif


///perform the global smoothing
void  BS_GlobalConvolver::SmoothLinkedVF(){
  int i,t,x,y,z;
  
  //1) Project the linked VF on the grid
  this->RefVecField_2_VecAtControlPts();

  //2) Compute the residual information lost when projecting the linked VF on the grid
  DeepCopyMultiplyVectorField(this->RefVecField,&this->ResidualField, 1);
  
  this->VecAtControlPts_2_RefVecField();  //modifies the linked VF but the former values of the linked VF are stored in ResidualField

  for (i=0;i<3;i++) for (t=0;t<this->ResidualField.NT;t++)  for (z = 0; z < this->ResidualField.NZ; z++) for (y = 0; y < this->ResidualField.NY; y++) for (x = 0; x < this->ResidualField.NX; x++){
    this->ResidualField.P(this->ResidualField.G(i,x,y,z,t) - this->RefVecField->G(i,x,y,z,t),i,x,y,z,t);
  }
  
  //3) weight the project values
  this->NodesWeight();
  
  //4) recompute GradE after the weighting process
  this->VecAtControlPts_2_RefVecField();
  
  
  //5) add the residual
  for (i=0;i<3;i++) for (t=0;t<this->ResidualField.NT;t++)  for (z = 0; z < this->ResidualField.NZ; z++) for (y = 0; y < this->ResidualField.NY; y++) for (x = 0; x < this->ResidualField.NX; x++){
    this->RefVecField->P(this->RefVecField->G(i,x,y,z,t) + this->ResidualField.G(i,x,y,z,t),i,x,y,z,t);
  }
}

///9.4: write an identity M which can be used for a given velocity field and a given grid size + other writing functions

///Write an identity template matrix which could be used with 'Load_WeightingMatrix'
void  BS_GlobalConvolver::WriteTemplateNodesWeighterIdMat(VectorField * LinkedVF,int GridStepSize,char * Mfile){
  FILE *DataFile;
  float LocValue;
  int ReqNbNodes;
  int i,j;
  
  //1) half-initiate the class with the linked velocity field and the grid step size
  this->LinkWithVecField(LinkedVF,GridStepSize);
  
  //2) compute the number of nodes required
  ReqNbNodes=this->NbCtPtsX*this->NbCtPtsY*this->NbCtPtsZ;
  
  if (ReqNbNodes==0){
    cout << "Basis to project the vector field must be defined before evaluating the weighting matrix size" << endl;
    return;
    }
  
  //3) write matrix template
  
  //3.1) open file
  DataFile=fopen(Mfile,"w");
    
   //3.2) write the data
   for (i=0;i<ReqNbNodes*3;i++){
      for (j=0;j<ReqNbNodes*3;j++){
        if (j!=0)
          fprintf(DataFile," ");
        if (i==j)
          fprintf(DataFile,"%f",1.);
        else
          fprintf(DataFile,"%f",0.);
      }
    fprintf(DataFile,"\n");
  }
  
  //3.3) close DataFile
  fclose(DataFile);
  
  
  //4) Message
  cout << "Identity matrix was returned in "<< Mfile << endl;
}



void BS_GlobalConvolver::WriteControlPointWeights(char * FileNameX,char * FileNameY,char * FileNameZ){
  this->VecAtControlPts.Write(FileNameX,FileNameY,FileNameZ);
  }

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                     10: SURFACE MESHES (with triangles)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// 10.1) ++++++++++++++++ Class nodes list  (simplified compared to UtilzGraph) ++++++++++++++++ 


/// Constructor
NodesList::NodesList(){
  this->NbNodes=0;
}
  
/// Destructor
NodesList::~NodesList(){
  if (this->NbNodes>0){
    this->NbNodes=0;
    delete this->Node;
    if (this->WithWeight==1) delete this->Weight;
  }
}


///allocate memory for the list
///size; size of the list
///with_weight: allocate memory for weights if ==1 (in addition to the node identifiers)
void NodesList::allocate(int size,int with_weight){
  int i;
  
  //if the list is already allocated, delete it
  if (this->NbNodes>0){
    cout << "An existing NodesList is redefined" << endl;
    this->NbNodes=0;
    delete this->Node;
    delete this->Weight;
  }
  
  //set the number of nodes
  this->NbNodes=size;
  
  //allocate the memory
  this->Node = new int [this->NbNodes];
  if (with_weight==1){
    this->WithWeight=1;
    this->Weight = new float [this->NbNodes];
  }
  else{
    this->WithWeight=0;
    this->Weight = new float [1];
  }
  
  //fill with default values
  for (i=0;i<this->NbNodes;i++) this->Node[i]=0;
  if (this->WithWeight==1){
    for (i=0;i<this->NbNodes;i++) this->Weight[i]=0;
  }
  else
    this->Weight[0]=0;
}


///Reduce the number of nodes in the list (additional memory is still allocated but we cannot access data which may not be defined)
inline void NodesList::ReduceNbNodes(int size){
  if (size<this->NbNodes) this->NbNodes=size; 
}

///Get and put functions
inline int NodesList::GetNbNodes(){
  return this->NbNodes;
}

inline int NodesList::GetNode(int entry){
  return this->Node[entry];
}

inline void NodesList::PutNode(int entry, int ID_node){
  this->Node[entry]=ID_node;
}

inline float NodesList::GetWeight(int entry){
  return this->WithWeight*this->Weight[this->WithWeight*entry];
}

inline void NodesList::PutWeight(int entry, float locWeight){
  this->Weight[this->WithWeight*entry]=locWeight;
}



/// 10.2) ++++++++++++++++ Class Graph (simplified compared to UtilzGraph)  +++++++++++++++++

//constructor
Graph::Graph(void){
  this->NbNodes=0;
  this->NbEdges=0;
  this->isSymmetric=0;
};

//destructor
Graph::~Graph(void){};





/// 10.3) ++++++++++++++++ Class ValuedGraph +++++++++++++++++
/// -> simplified compared to UtilzGraph but with AffineTransfo and NonRigidTransfo


/// Constructor
ValuedGraph::ValuedGraph(){
  this->KnownNodeCoordinates=0;
  this->KnownTriangularFaces=0;
  this->NbNodes=0;
  this->NbTriangularFaces=0;
  }
  
/// Destructor
ValuedGraph::~ValuedGraph(){
  int i;
  
  if (this->NbNodes>0){
    delete this->NodeValues;
    
    this->Connections.~Graph();  //not sure it will make anything
    
    if (this->KnownNodeCoordinates==1){
      delete this->NodeCoordX;
      delete this->NodeCoordY;
      delete this->NodeCoordZ;
      }
    
    if (this->KnownTriangularFaces==1){
      for (i=0;i<this->NbTriangularFaces;i++) 
        delete this->TriangularFaces[i];
      
      delete this->TriangularFaces;
      }
    }
  
  
  this->KnownNodeCoordinates=0;
  this->KnownTriangularFaces=0;
  this->NbNodes=0;
  this->NbTriangularFaces=0;
  }
  

///Read the graph in an ascci POLYDATA vtk file containing a mesh with values at each node. 
///The imported mesh is a graph with weigths which are the distances between the nodes
///Returns 1 if the vtk file seems to be read normally and 0 otherwise.
int ValuedGraph::ReadMeshInVtkFile(char MeshFileName[256]){
  FILE *fichdata;
  int tmpInt,tmpInt2,LocNbNodes,i,j;
  int found;
  int Node1,Node2,Node3;
  float tmpFl,locX,locY,locZ;
  char str [80];
  int LocNode,NgbhNode1,NgbhNode2;
  int EverythingIsFine;
  float distance;
  
  //open the file
  fichdata = fopen (MeshFileName, "r" );
  EverythingIsFine=1;
  
  //1) load the node coordinates
  fseek(fichdata,0,SEEK_SET);
  
  found=0;
  while(!feof(fichdata)){
    //1.1) find the beginning of the part
    tmpInt=fscanf (fichdata, "%s", str);
    
    if (strcmp(str, "POINTS") == 0){ //read the node coordinates once found
      found=1; 
      
      //1.2) read the local header
      tmpInt=fscanf (fichdata, "%d", &LocNbNodes);
      cerr << "Point coordinates found (" << LocNbNodes << " coordinates expected)" << endl;
      tmpInt=fscanf (fichdata, "%s", str);
      
      //1.3) allocate memory for the graph
      this->NodeCoordX = new float [LocNbNodes];
      this->NodeCoordY = new float [LocNbNodes];
      this->NodeCoordZ = new float [LocNbNodes];
      
      //1.4) read the node coordinates
      for (i=0;i<LocNbNodes;i++){
          tmpInt=fscanf(fichdata,"%f",&locX); this->NodeCoordX[i]=-locX;
          tmpInt=fscanf(fichdata,"%f",&locY); this->NodeCoordY[i]=-locY;  // WARNING: X and Y axes are INVERTED
          tmpInt=fscanf(fichdata,"%f",&locZ); this->NodeCoordZ[i]=locZ;
          //printf ("Node: %f %f %f\n",locX,locY,locZ);
      }
      
    this->KnownNodeCoordinates=1;
    this->NbNodes=LocNbNodes;
    }
  }
  
  if (found==0){
    cerr << "Point coordinates not found" << endl; 
    EverythingIsFine=0;
    }
  
  //2) load the nodes connectivity (edges)
  fseek(fichdata,0,SEEK_SET);
  
  found=0;
  while(!feof(fichdata)){
    //2.1) find the beginning of the part
    tmpInt=fscanf (fichdata, "%s", str);
    
    if (strcmp(str, "POLYGONS") == 0){ //read the node connectivity once found
      found=1;
      
      //2.2) read the local header
      tmpInt=fscanf (fichdata, "%d", &this->NbTriangularFaces);
      tmpInt=fscanf (fichdata, "%d", &tmpInt2);
      cerr << "Faces found (" << this->NbTriangularFaces << " faces expected)" << endl;
      
      //2.3) allocate memory for the graph
      this->TriangularFaces = new int * [this->NbTriangularFaces];
      for (i=0;i<this->NbTriangularFaces;i++)
        this->TriangularFaces[i] = new int [3];
      
      
      //2.4) read the node composing the faces
      for (i=0;i<this->NbTriangularFaces;i++){
          tmpInt=fscanf (fichdata, "%d", &tmpInt2);
          tmpInt=fscanf (fichdata, "%d", &Node1); this->TriangularFaces[i][0]=Node1;
          tmpInt=fscanf (fichdata, "%d", &Node2); this->TriangularFaces[i][1]=Node2;
          tmpInt=fscanf (fichdata, "%d", &Node3); this->TriangularFaces[i][2]=Node3;
          //printf ("Faces: %d %d %d\n",Node1,Node2,Node3);
      }
      
    this->KnownTriangularFaces=1;
    }
  }
  
  if (found==0){
    cerr << "Faces not found" << endl;
    EverythingIsFine=0;
    }

  
  
  //3) load the values attributed to each node
  fseek(fichdata,0,SEEK_SET);
  
  found=0;
  while(!feof(fichdata)){
    //3.1) find the beginning of the part
    tmpInt=fscanf (fichdata, "%s", str);
    
    if (strcmp(str, "POINT_DATA") == 0){ //read the values once found
      found=1;
      
      //3.2) read the local header
      tmpInt=fscanf (fichdata, "%d", &tmpInt2);
      
      if (LocNbNodes==tmpInt2){
        LocNbNodes=tmpInt2;
        cerr << "Node values found (" << LocNbNodes << " values expected)" << endl;
      }
      else{
        EverythingIsFine=0;
        cerr << "The number of node values is not the same as the number of nodes!" << endl;
      }
      
      tmpInt=fscanf (fichdata, "%s", str);  //SCALARS
      tmpInt=fscanf (fichdata, "%s", str);  //scalars 
      tmpInt=fscanf (fichdata, "%s", str);  //float
      tmpInt=fscanf (fichdata, "%s", str);  //LOOKUP_TABLE 
      tmpInt=fscanf (fichdata, "%s", str);  //default
      
      
      //3.3) allocate memory for the node values
      this->NodeValues = new float [LocNbNodes];
      
      //3.4) read the node values
      for (i=0;i<LocNbNodes;i++){
          tmpInt=fscanf(fichdata,"%f",&tmpFl);
          this->NodeValues[i]=tmpFl;
          
          //printf ("Node value: %f\n",tmpFl);
      }
      
      this->NbNodes=LocNbNodes;
    }
  }
  
  if (found==0){
    cerr << "Node values not found" << endl;
    EverythingIsFine=0;
    }


  //4) convert faces into an actual graph
  
  //4.1) fill known entries and allocate memory for the lists of neighbors
  this->Connections.NbNodes=LocNbNodes;
  this->Connections.isSymmetric=1;
  this->Connections.NgbhsOfNode = new NodesList [LocNbNodes];
  
  //4.2) count the nb of neighbors of each node and allocate corresponding memory
  int * NodeNgbhNb;
  NodeNgbhNb = new int  [LocNbNodes];
  for (i=0;i<LocNbNodes;i++) NodeNgbhNb[i]=0;
  
  for (i=0;i<this->NbTriangularFaces;i++){
    NodeNgbhNb[this->TriangularFaces[i][0]]+=2;
    NodeNgbhNb[this->TriangularFaces[i][1]]+=2;
    NodeNgbhNb[this->TriangularFaces[i][2]]+=2;
    }   //remark: too many connections will be probably estimated as a connection may be observed in 2 faces (and not only 1). This will be post-treated.
  
  for (i=0;i<LocNbNodes;i++) 
    this->Connections.NgbhsOfNode[i].allocate(NodeNgbhNb[i],1);
  
  
  //4.3) fill the connections
  this->Connections.NbEdges=0;       //number of directed edges in the graph
  
  for (i=0;i<LocNbNodes;i++) NodeNgbhNb[i]=0;
  
  for (i=0;i<this->NbTriangularFaces;i++){
    //4.3.1) connections of this->TriangularFaces[i][0]
    LocNode=this->TriangularFaces[i][0]; NgbhNode1=this->TriangularFaces[i][1]; NgbhNode2=this->TriangularFaces[i][2];
    
    found=0;
    for (j=0;j<NodeNgbhNb[LocNode];j++) if(this->Connections.NgbhsOfNode[LocNode].GetNode(j)==NgbhNode1) found=1;
    if (found==0){
      this->Connections.NgbhsOfNode[LocNode].PutNode(NodeNgbhNb[LocNode],NgbhNode1); 
      distance=(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode1])*(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode1]);
      distance+=(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode1])*(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode1]);
      distance+=(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode1])*(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode1]);
      distance=sqrt(distance);
      this->Connections.NgbhsOfNode[LocNode].PutWeight(NodeNgbhNb[LocNode],distance);
      NodeNgbhNb[LocNode]++; 
      this->Connections.NbEdges++;
      }
    
    found=0;
    for (j=0;j<NodeNgbhNb[LocNode];j++) if(this->Connections.NgbhsOfNode[LocNode].GetNode(j)==NgbhNode2) found=1;
    if (found==0){
      this->Connections.NgbhsOfNode[LocNode].PutNode(NodeNgbhNb[LocNode],NgbhNode2);
      distance=(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode2])*(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode2]);
      distance+=(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode2])*(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode2]);
      distance+=(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode2])*(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode2]);
      distance=sqrt(distance);
      this->Connections.NgbhsOfNode[LocNode].PutWeight(NodeNgbhNb[LocNode],distance);
      NodeNgbhNb[LocNode]++;
      this->Connections.NbEdges++;
      }
    
    //4.3.2) connections of this->TriangularFaces[i][1]
    LocNode=this->TriangularFaces[i][1]; NgbhNode1=this->TriangularFaces[i][0]; NgbhNode2=this->TriangularFaces[i][2];
    
    found=0;
    for (j=0;j<NodeNgbhNb[LocNode];j++) if(this->Connections.NgbhsOfNode[LocNode].GetNode(j)==NgbhNode1) found=1;
    if (found==0){
      this->Connections.NgbhsOfNode[LocNode].PutNode(NodeNgbhNb[LocNode],NgbhNode1); 
      distance=(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode1])*(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode1]);
      distance+=(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode1])*(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode1]);
      distance+=(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode1])*(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode1]);
      distance=sqrt(distance);
      this->Connections.NgbhsOfNode[LocNode].PutWeight(NodeNgbhNb[LocNode],distance);
      NodeNgbhNb[LocNode]++; 
      this->Connections.NbEdges++;
      }
    
    found=0;
    for (j=0;j<NodeNgbhNb[LocNode];j++) if(this->Connections.NgbhsOfNode[LocNode].GetNode(j)==NgbhNode2) found=1;
    if (found==0){
      this->Connections.NgbhsOfNode[LocNode].PutNode(NodeNgbhNb[LocNode],NgbhNode2);
      distance=(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode2])*(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode2]);
      distance+=(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode2])*(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode2]);
      distance+=(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode2])*(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode2]);
      distance=sqrt(distance);
      this->Connections.NgbhsOfNode[LocNode].PutWeight(NodeNgbhNb[LocNode],distance);
      NodeNgbhNb[LocNode]++;
      this->Connections.NbEdges++;
      }
    
    //4.3.3) connections of this->TriangularFaces[i][2]
    LocNode=this->TriangularFaces[i][2]; NgbhNode1=this->TriangularFaces[i][1]; NgbhNode2=this->TriangularFaces[i][0];

    found=0;
    for (j=0;j<NodeNgbhNb[LocNode];j++) if(this->Connections.NgbhsOfNode[LocNode].GetNode(j)==NgbhNode1) found=1;
    if (found==0){
      this->Connections.NgbhsOfNode[LocNode].PutNode(NodeNgbhNb[LocNode],NgbhNode1); 
      distance=(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode1])*(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode1]);
      distance+=(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode1])*(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode1]);
      distance+=(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode1])*(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode1]);
      distance=sqrt(distance);
      this->Connections.NgbhsOfNode[LocNode].PutWeight(NodeNgbhNb[LocNode],distance);
      NodeNgbhNb[LocNode]++; 
      this->Connections.NbEdges++;
      }
    
    found=0;
    for (j=0;j<NodeNgbhNb[LocNode];j++) if(this->Connections.NgbhsOfNode[LocNode].GetNode(j)==NgbhNode2) found=1;
    if (found==0){
      this->Connections.NgbhsOfNode[LocNode].PutNode(NodeNgbhNb[LocNode],NgbhNode2);
      distance=(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode2])*(this->NodeCoordX[LocNode]-this->NodeCoordX[NgbhNode2]);
      distance+=(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode2])*(this->NodeCoordY[LocNode]-this->NodeCoordY[NgbhNode2]);
      distance+=(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode2])*(this->NodeCoordZ[LocNode]-this->NodeCoordZ[NgbhNode2]);
      distance=sqrt(distance);
      this->Connections.NgbhsOfNode[LocNode].PutWeight(NodeNgbhNb[LocNode],distance);
      NodeNgbhNb[LocNode]++;
      this->Connections.NbEdges++;
      }
  }
  
  //4.4) reduce the number of neigbors to its actual number
  for (i=0;i<LocNbNodes;i++) 
    this->Connections.NgbhsOfNode[i].ReduceNbNodes(NodeNgbhNb[i]);
  
  //close the file
  fclose(fichdata);
  
  cerr << "The graph contains " <<  this->Connections.NbNodes << " nodes and " << this->Connections.NbEdges << " edges." << endl;
  
  return EverythingIsFine;
}

   
///Save the graph as a mesh in an ascci POLYDATA vtk file. 
///The graph will only be saved if 'KnownNodeCoordinates' and 'KnownTriangularFaces' are equal to 1.
///Returns 1 if the graph is saved and 0 otherwise.
int ValuedGraph::WriteMeshInVtkFile(char MeshFileName[256]){
  FILE *outFile;
  int i,j;
  
  //1) Init
  
  //check whether the node coordinates and the triangular faces are known
  if ((KnownNodeCoordinates==0)||(KnownTriangularFaces==0))
    return 0;
  
  //open file
  outFile=fopen(MeshFileName,"w");

  //save the header
  fprintf(outFile,"# vtk DataFile Version 3.0\n");
  fprintf(outFile,"vtk output\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"\n");

  //2) save the nodes
  fprintf(outFile,"DATASET POLYDATA\n");
  fprintf(outFile,"POINTS %d float\n",this->NbNodes);

  for (i=0;i<this->NbNodes;i++)
    fprintf(outFile,"%f %f %f\n",-this->NodeCoordX[i],-this->NodeCoordY[i],this->NodeCoordZ[i]);  // WARNING: X and Y axes are INVERTED
  
  
  //3) save the faces
  fprintf(outFile,"POLYGONS %d %d\n",this->NbTriangularFaces,4*this->NbTriangularFaces);

  for (i=0;i<this->NbTriangularFaces;i++)
    fprintf(outFile,"3 %d %d %d\n",this->TriangularFaces[i][0],this->TriangularFaces[i][1],this->TriangularFaces[i][2]);
  
  //4) save the intensities
  fprintf(outFile,"POINT_DATA %d\n",this->NbNodes);
  fprintf(outFile,"SCALARS volume float\n");
  fprintf(outFile,"LOOKUP_TABLE default\n");

  for (i=0;i<this->NbNodes;i++)
    fprintf(outFile,"%f ",this->NodeValues[i]);

  //5) finalize
  fclose(outFile);
  
  return 1;
}



///Affine transformation of the mesh coordinates  (if KnownNodeCoordinates==1)
///DefQuaternion is a 4*4 matrix: [R_xx R_xy R_xz T_x ; R_yx R_yy R_yz T_y ; R_zx R_zy R_zz T_z ; 0 0 0 1]
void ValuedGraph::AffineTransfo(float DefQuaternion[4][4]){
  int i;
  float tmpX,tmpY,tmpZ;
  
  for (i=0;i<this->NbNodes;i++){
    tmpX=this->NodeCoordX[i]*DefQuaternion[0][0]+this->NodeCoordY[i]*DefQuaternion[0][1]+this->NodeCoordZ[i]*DefQuaternion[0][2]+DefQuaternion[0][3];
    tmpY=this->NodeCoordX[i]*DefQuaternion[1][0]+this->NodeCoordY[i]*DefQuaternion[1][1]+this->NodeCoordZ[i]*DefQuaternion[1][2]+DefQuaternion[1][3];
    tmpZ=this->NodeCoordX[i]*DefQuaternion[2][0]+this->NodeCoordY[i]*DefQuaternion[2][1]+this->NodeCoordZ[i]*DefQuaternion[2][2]+DefQuaternion[2][3];
    NodeCoordX[i]=tmpX;
    NodeCoordY[i]=tmpY;
    NodeCoordZ[i]=tmpZ;
    }
}


///Non rigid deformation with a DisplacementField
///-> to make sense the DisplacementField should cover the domain of the ValuedGraph (according to its image to world properties).
///-> Coordinates of 'this' will be transformed according to the local displacement vector of the DisplacementField
void ValuedGraph::NonRigidTransfo(VectorField * DisplacementField){
  int i;
  float tmpX,tmpY,tmpZ;
  float tmpX2,tmpY2,tmpZ2;
  
  
  for (i=0;i<this->NbNodes;i++){
    tmpX=this->NodeCoordX[i]*DisplacementField->World2Image[0][0]+this->NodeCoordY[i]*DisplacementField->World2Image[0][1]+this->NodeCoordZ[i]*DisplacementField->World2Image[0][2]+DisplacementField->World2Image[0][3];
    tmpY=this->NodeCoordX[i]*DisplacementField->World2Image[1][0]+this->NodeCoordY[i]*DisplacementField->World2Image[1][1]+this->NodeCoordZ[i]*DisplacementField->World2Image[1][2]+DisplacementField->World2Image[1][3];
    tmpZ=this->NodeCoordX[i]*DisplacementField->World2Image[2][0]+this->NodeCoordY[i]*DisplacementField->World2Image[2][1]+this->NodeCoordZ[i]*DisplacementField->World2Image[2][2]+DisplacementField->World2Image[2][3];
    
    this->NodeCoordX[i]+=DisplacementField->G(0,tmpX,tmpY,tmpZ);
    this->NodeCoordY[i]+=DisplacementField->G(1,tmpX,tmpY,tmpZ);
    this->NodeCoordZ[i]+=DisplacementField->G(2,tmpX,tmpY,tmpZ);
    }
}




///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                     11: EMPIRICAL MODE DECOMPOSITION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///Compute the local minima and maxima from the image
void ExtractLocalExtrema2(ScalarField * InputImag,ScalarField * LocalMax, ScalarField * LocalMin){
  float epsilon,tmpfl;
  int x,y,z;
  int LocCount;
  int NbLocMax,NbLocMin;
  
  epsilon=0.000000001;
  
  //extract local max
  for (z=0;z<InputImag->NZ;z++) for (y=0;y<InputImag->NY;y++) for (x=0;x<InputImag->NX;x++) LocalMax->P(0,x,y,z);
  
  NbLocMax=0;
  for (z=1;z<InputImag->NZ-1;z++) for (y=1;y<InputImag->NY-1;y++) for (x=1;x<InputImag->NX-1;x++){
    LocCount=0;
    tmpfl=InputImag->G(x,y,z)+epsilon;
    if (InputImag->G(x,y+1,z-1)>tmpfl) LocCount++;
    if (InputImag->G(x,y+1,z)>tmpfl)   LocCount++;
    if (InputImag->G(x,y+1,z+1)>tmpfl) LocCount++;
    if (InputImag->G(x,y-1,z-1)>tmpfl) LocCount++;
    if (InputImag->G(x,y-1,z)>tmpfl)   LocCount++;
    if (InputImag->G(x,y-1,z+1)>tmpfl) LocCount++;
    if (InputImag->G(x,y,z+1)>tmpfl)   LocCount++;
    if (InputImag->G(x,y,z-1)>tmpfl)   LocCount++;
    if (InputImag->G(x-1,y,z)>tmpfl)   LocCount++;
    if (InputImag->G(x-1,y+1,z)>tmpfl) LocCount++;
    if (InputImag->G(x-1,y-1,z)>tmpfl) LocCount++;
    if (InputImag->G(x-1,y,z+1)>tmpfl) LocCount++;
    if (InputImag->G(x-1,y,z-1)>tmpfl) LocCount++;
    if (InputImag->G(x+1,y,z)>tmpfl)   LocCount++;
    if (InputImag->G(x+1,y+1,z)>tmpfl) LocCount++;
    if (InputImag->G(x+1,y-1,z)>tmpfl) LocCount++;
    if (InputImag->G(x+1,y,z+1)>tmpfl) LocCount++;
    if (InputImag->G(x+1,y,z-1)>tmpfl) LocCount++;
    
    if (LocCount>17) {LocalMax->P(tmpfl,x,y,z); NbLocMax++;}
    }
      
  //extract local min
  for (z=0;z<InputImag->NZ;z++) for (y=0;y<InputImag->NY;y++) for (x=0;x<InputImag->NX;x++) LocalMin->P(0,x,y,z);
  
  NbLocMin=0;
  for (z=1;z<InputImag->NZ-1;z++) for (y=1;y<InputImag->NY-1;y++) for (x=1;x<InputImag->NX-1;x++){
    LocCount=0;
    tmpfl=InputImag->G(x,y,z)-epsilon;
    if (InputImag->G(x,y+1,z-1)<tmpfl) LocCount++;
    if (InputImag->G(x,y+1,z)<tmpfl)   LocCount++;
    if (InputImag->G(x,y+1,z+1)<tmpfl) LocCount++;
    if (InputImag->G(x,y-1,z-1)<tmpfl) LocCount++;
    if (InputImag->G(x,y-1,z)<tmpfl)   LocCount++;
    if (InputImag->G(x,y-1,z+1)<tmpfl) LocCount++;
    if (InputImag->G(x,y,z+1)<tmpfl)   LocCount++;
    if (InputImag->G(x,y,z-1)<tmpfl)   LocCount++;
    if (InputImag->G(x-1,y,z)<tmpfl)   LocCount++;
    if (InputImag->G(x-1,y+1,z)<tmpfl) LocCount++;
    if (InputImag->G(x-1,y-1,z)<tmpfl) LocCount++;
    if (InputImag->G(x-1,y,z+1)<tmpfl) LocCount++;
    if (InputImag->G(x-1,y,z-1)<tmpfl) LocCount++;
    if (InputImag->G(x+1,y,z)<tmpfl)   LocCount++;
    if (InputImag->G(x+1,y+1,z)<tmpfl) LocCount++;
    if (InputImag->G(x+1,y-1,z)<tmpfl) LocCount++;
    if (InputImag->G(x+1,y,z+1)<tmpfl) LocCount++;
    if (InputImag->G(x+1,y,z-1)<tmpfl) LocCount++;
    
    if (LocCount>17){ LocalMin->P(tmpfl,x,y,z); NbLocMin++;}
    }
    
    cout << NbLocMax << " Loc. Max. and " << NbLocMin <<  " Loc. Min.  " << endl;
    
    //LocalMax->Write("LocMax.nii");
    //LocalMin->Write("LocMin.nii");
}



///Compute the local minima and maxima from the image
void ExtractLocalExtrema(ScalarField * InputImag,ScalarField * TmpImag,ScalarField * LocalMax, ScalarField * LocalMin){
  float epsilon,tmpfl;
  int x,y,z;
  
  epsilon=0.00001;
  
  //1) dilation of InputImag
  //1.1) x direction
  for (z=0;z<InputImag->NZ;z++) for (y=0;y<InputImag->NY;y++){
    for (x=0;x<InputImag->NX;x++) {
      tmpfl=InputImag->G(x,y,z);
      if (x!=InputImag->NX-1) if (tmpfl<InputImag->G(x+1,y,z)) tmpfl=InputImag->G(x+1,y,z);
      if (x!=0)           if (tmpfl<InputImag->G(x-1,y,z)) tmpfl=InputImag->G(x-1,y,z);
      LocalMax->P(tmpfl,x,y,z);
    }
  } 
      
  //1.2) y direction
  for (z=0;z<InputImag->NZ;z++) for (x=0;x<InputImag->NX;x++){
    for (y=0;y<InputImag->NY;y++){
      tmpfl=LocalMax->G(x,y,z);
      if (y!=LocalMax->NY-1) if (tmpfl<LocalMax->G(x,y+1,z)) tmpfl=LocalMax->G(x,y+1,z);
      if (y!=0)           if (tmpfl<LocalMax->G(x,y-1,z)) tmpfl=LocalMax->G(x,y-1,z);
      TmpImag->P(tmpfl,x,y,z);
    }
  } 
      
  //1.3) z direction
  for (y=0;y<InputImag->NY;y++) for (x=0;x<InputImag->NX;x++){
    for (z=0;z<InputImag->NZ;z++){
      tmpfl=TmpImag->G(x,y,z);
      if (z!=TmpImag->NZ-1) if (tmpfl<TmpImag->G(x,y,z+1)) tmpfl=TmpImag->G(x,y,z+1);
      if (z!=0)           if (tmpfl<TmpImag->G(x,y,z-1)) tmpfl=TmpImag->G(x,y,z-1);
      LocalMax->P(tmpfl,x,y,z);
    }
  }
  
  //2) extract local maxima
  for (z=0;z<InputImag->NZ;z++) for (y=0;y<InputImag->NY;y++) for (x=0;x<InputImag->NX;x++){
    if (fabs(LocalMax->G(x,y,z)-InputImag->G(x,y,z))<epsilon)
      LocalMax->P(InputImag->G(x,y,z),x,y,z);
    else
      LocalMax->P(0,x,y,z);
    }

  //3) dilation of InputImag
  //3.1) x direction
  for (z=0;z<InputImag->NZ;z++) for (y=0;y<InputImag->NY;y++){
    for (x=0;x<InputImag->NX;x++) {
      tmpfl=InputImag->G(x,y,z);
      if (x!=InputImag->NX-1) if (tmpfl>InputImag->G(x+1,y,z)) tmpfl=InputImag->G(x+1,y,z);
      if (x!=0)           if (tmpfl>InputImag->G(x-1,y,z)) tmpfl=InputImag->G(x-1,y,z);
      LocalMin->P(tmpfl,x,y,z);
    }
  } 
      
  //3.2) y direction
  for (z=0;z<InputImag->NZ;z++) for (x=0;x<InputImag->NX;x++){
    for (y=0;y<InputImag->NY;y++){
      tmpfl=LocalMin->G(x,y,z);
      if (y!=LocalMin->NY-1) if (tmpfl>LocalMin->G(x,y+1,z)) tmpfl=LocalMin->G(x,y+1,z);
      if (y!=0)           if (tmpfl>LocalMin->G(x,y-1,z)) tmpfl=LocalMin->G(x,y-1,z);
      TmpImag->P(tmpfl,x,y,z);
    }
  } 
      
  //3.3) z direction
  for (y=0;y<InputImag->NY;y++) for (x=0;x<InputImag->NX;x++){
    for (z=0;z<InputImag->NZ;z++){
      tmpfl=TmpImag->G(x,y,z);
      if (z!=TmpImag->NZ-1) if (tmpfl>TmpImag->G(x,y,z+1)) tmpfl=TmpImag->G(x,y,z+1);
      if (z!=0)           if (tmpfl>TmpImag->G(x,y,z-1)) tmpfl=TmpImag->G(x,y,z-1);
      LocalMin->P(tmpfl,x,y,z);
    }
  }
  
  //4) extract local maxima
  for (z=0;z<InputImag->NZ;z++) for (y=0;y<InputImag->NY;y++) for (x=0;x<InputImag->NX;x++){
    if (fabs(LocalMin->G(x,y,z)-InputImag->G(x,y,z))<epsilon)
      LocalMin->P(InputImag->G(x,y,z),x,y,z);
    else
      LocalMin->P(0,x,y,z);
    }
}



///Estimate a good standard deviation of the Jacobian to perform normalized convolution of the input local extra map
float EstimateGoodSigma(ScalarField * LocalMax, ScalarField * LocalMin){
  float epsilon,GoodSigma,tmpFl,tmpFl2,tmpFl3;
  int NonNullMin,NonNullMax,NbPoints;
  int LargestAxis1,LargestAxis2,LargestAxis3,tmpInt;
  int x,y,z;
  float OneThird;
    
  OneThird=static_cast<float>(0.333);
    
  epsilon=0.00001;
  
  //1) compute NonNullMin, NonNullMax...
  //...
  NbPoints=LocalMax->NZ*LocalMax->NY*LocalMax->NX;
  
  //...
  NonNullMax=0;
  for (z=0;z<LocalMax->NZ;z++) for (y=0;y<LocalMax->NY;y++) for (x=0;x<LocalMax->NX;x++)
    if (fabs(LocalMax->G(x,y,z))>epsilon)
      NonNullMax++;

  //...
  NonNullMin=0;
  for (z=0;z<LocalMin->NZ;z++) for (y=0;y<LocalMin->NY;y++) for (x=0;x<LocalMin->NX;x++)
    if (fabs(LocalMin->G(x,y,z))>epsilon)
      NonNullMin++;
      
  //2) compute the largest axes
  LargestAxis1=LocalMax->NX;  LargestAxis2=LocalMax->NY;  LargestAxis3=LocalMax->NZ;
  if (LargestAxis1<LargestAxis2){tmpInt=LargestAxis1; LargestAxis1=LargestAxis2; LargestAxis2=tmpInt;}
  if (LargestAxis1<LargestAxis3){tmpInt=LargestAxis1; LargestAxis1=LargestAxis3; LargestAxis3=tmpInt;}
  if (LargestAxis2<LargestAxis3){tmpInt=LargestAxis2; LargestAxis2=LargestAxis3; LargestAxis3=tmpInt;}
  
  
  //3) compute a good sigma for...  (min/max are supposed homogenously distributed in space)
  if ((LargestAxis1/5>LargestAxis2)&&(LargestAxis1/5>LargestAxis3)){ //1D or almost 1D domain
    tmpFl=static_cast<float>(LargestAxis1);
    tmpFl2=static_cast<float>(NonNullMax);
    tmpFl3=static_cast<float>(NonNullMin);
    GoodSigma=(2*tmpFl)/(tmpFl2+tmpFl3);
    }
  else if ((LargestAxis1/5>LargestAxis3)){ //2D or almost 2D domain
    tmpFl=static_cast<float>(LargestAxis1+LargestAxis2)/2;
    tmpFl2=sqrt(static_cast<float>(NonNullMax));
    tmpFl3=sqrt(static_cast<float>(NonNullMin));
    GoodSigma=(2*tmpFl)/(tmpFl2+tmpFl3);
    }
  else{ //3D domain
    tmpFl=static_cast<float>(LargestAxis1+LargestAxis2+LargestAxis3)/3;
    tmpFl2=pow(static_cast<float>(NonNullMax),OneThird);
    tmpFl3=pow(static_cast<float>(NonNullMin),OneThird);
    GoodSigma=(2*tmpFl)/(tmpFl2+tmpFl3);
    }
  
  GoodSigma/=2;
  
  //4) return the sigma
  return GoodSigma;
  
}



//Add simulated Local Extrema to a map of Local extrema to numerically stabilize the normalized convolution
void AddSimulatedLocalExtrema(ScalarField * LocalExtrema,ScalarField * TmpImag,ScalarField * PointsToAdd){
  int x,y,z;
  int x2,y2,z2;
  int OnlySinglePoints,SinglePoint,BoundaryPoint;
  float epsilon;
  int NbLocalExtrema,LocExtToDetect,DetectedLocExt;
  float tmpFl;
  int iteration;
  float Value_LocExtremum1,Value_LocExtremum2,Value_LocExtremum3;
  float Wght_LocExtremum1,Wght_LocExtremum2,Wght_LocExtremum3;
  
  epsilon=0.00001;
  
  //1) init
  
  //1.1) count the number of local extrema 
  NbLocalExtrema=0;
  for (z=0;z<LocalExtrema->NZ;z++) for (y=0;y<LocalExtrema->NY;y++) for (x=0;x<LocalExtrema->NX;x++) 
    if (fabs(LocalExtrema->G(x,y,z))>epsilon)
      NbLocalExtrema++;
  
  //cout << "Nb local extrema detected = " << NbLocalExtrema << endl;
  
  //1.2) special case where there is no local extrema
  if (NbLocalExtrema==0){
    cout << "No local extrema to extrapolate/interpolate -> all set to zero" << endl;
    for (z=0;z<LocalExtrema->NZ;z++) for (y=0;y<LocalExtrema->NY;y++) for (x=0;x<LocalExtrema->NX;x++) 
      LocalExtrema->P(0,x,y,z);
      return;
    }
    
  //1.3) special case where there is one local extrema
  if (NbLocalExtrema==1){
    for (z=0;z<LocalExtrema->NZ;z++) for (y=0;y<LocalExtrema->NY;y++) for (x=0;x<LocalExtrema->NX;x++) 
      if (fabs(LocalExtrema->G(x,y,z))>epsilon)
        tmpFl=LocalExtrema->G(x,y,z);
    
    for (z=0;z<LocalExtrema->NZ;z++) for (y=0;y<LocalExtrema->NY;y++) for (x=0;x<LocalExtrema->NX;x++)
        LocalExtrema->P(tmpFl,x,y,z);
    
    cout << "Only one local extremum -> all set to its value" << endl;
    return;
    }
    
  //1.4) define the number of local extrema to detect for each extrapolation/interpolation
  if (NbLocalExtrema==2) LocExtToDetect=2;
  if (NbLocalExtrema>=3) LocExtToDetect=3;
  
  //cout << "For each simulated local extremum, " << LocExtToDetect << " points are used for interpolation (" << NbLocalExtrema << " local extrema in total)" << endl;
  
  //2) loop on the points to add
  for (z=0;z<PointsToAdd->NZ;z++) for (y=0;y<PointsToAdd->NY;y++) for (x=0;x<PointsToAdd->NX;x++) if (PointsToAdd->G(x,y,z)>0.5){
    
    //2.1) clean-up the temporary image and initiate the neighbors search
    for (z2=0;z2<TmpImag->NZ;z2++) for (y2=0;y2<TmpImag->NY;y2++) for (x2=0;x2<TmpImag->NX;x2++) TmpImag->P(0,x2,y2,z2);
    TmpImag->P(1,x,y,z);
    
    //2.2) detect the closest local extrema
    DetectedLocExt=0;
    iteration=0;
    
    while (DetectedLocExt<LocExtToDetect){
      //cout << "iteration " << iteration << endl;
      
      //2.2.1) find the points at the external boundary of the propagated domain
      for (z2=1;z2<TmpImag->NZ-1;z2++) for (y2=1;y2<TmpImag->NY-1;y2++) for (x2=1;x2<TmpImag->NX-1;x2++) if (fabs(TmpImag->G(x2,y2,z2)-1)<epsilon){
        if (fabs(TmpImag->G(x2+1,y2,z2))<epsilon) TmpImag->P(2,x2+1,y2,z2);
        if (fabs(TmpImag->G(x2-1,y2,z2))<epsilon) TmpImag->P(2,x2-1,y2,z2);
        if (fabs(TmpImag->G(x2,y2+1,z2))<epsilon) TmpImag->P(2,x2,y2+1,z2);
        if (fabs(TmpImag->G(x2,y2-1,z2))<epsilon) TmpImag->P(2,x2,y2-1,z2);
        if (fabs(TmpImag->G(x2,y2,z2+1))<epsilon) TmpImag->P(2,x2,y2,z2+1);
        if (fabs(TmpImag->G(x2,y2,z2-1))<epsilon) TmpImag->P(2,x2,y2,z2-1);
        }
      iteration++;
      
        if (iteration==TmpImag->NZ+TmpImag->NY+TmpImag->NX){
            cout << "Not enough points for interpolation have not been found... algorithm stopped" << endl;
            exit(0);
        }
      
      //2.2.2) check whether a local extremum is found
      for (z2=0;z2<TmpImag->NZ;z2++) for (y2=0;y2<TmpImag->NY;y2++) for (x2=0;x2<TmpImag->NX;x2++) if (fabs(TmpImag->G(x2,y2,z2)-2)<epsilon){
        if (fabs(LocalExtrema->G(x2,y2,z2))>epsilon){
          if (DetectedLocExt==0) {Value_LocExtremum1=LocalExtrema->G(x2,y2,z2); Wght_LocExtremum1=1/(static_cast<float>(iteration));}
          if (DetectedLocExt==1) {Value_LocExtremum2=LocalExtrema->G(x2,y2,z2); Wght_LocExtremum2=1/(static_cast<float>(iteration));}
          if (DetectedLocExt==2) {Value_LocExtremum3=LocalExtrema->G(x2,y2,z2); Wght_LocExtremum3=1/(static_cast<float>(iteration));}
          DetectedLocExt++;
        }
          TmpImag->P(1,x2,y2,z2);
        }
      }  
      
    //2.3) add the point
    if (LocExtToDetect==2)
      tmpFl=((Value_LocExtremum1*Wght_LocExtremum1)+(Value_LocExtremum2*Wght_LocExtremum2))/(Wght_LocExtremum1+Wght_LocExtremum2);
    
    if (LocExtToDetect==3)
      tmpFl=((Value_LocExtremum1*Wght_LocExtremum1)+(Value_LocExtremum2*Wght_LocExtremum2)+(Value_LocExtremum3*Wght_LocExtremum3))/(Wght_LocExtremum1+Wght_LocExtremum2+Wght_LocExtremum3);
    
    LocalExtrema->P(tmpFl,x,y,z);
  }
  

}



///Compute the normalized convolution of local minima/maxima map
void NormalizedConvolution(ScalarField * LocalExtrema,ScalarField * InputImag,LightFFTconvolver3D * TmpConvolver,ScalarField * TmpImag,float sigma,ScalarField * OutputConvol){
  int x,y,z;
  float epsilon,threshZero;
  int OK;
  float val;
  int NbNghb;
  int NoUnstableRegion;
  epsilon=0.00001;
  threshZero=epsilon*0.001;
  
  //1) init the convolver
  TmpConvolver->InitiateConvolver(LocalExtrema->NX,LocalExtrema->NY,LocalExtrema->NZ,1,sigma,sigma,sigma);
  
  //2) check if we need to add fake LocalExtrema to make the n.c. numercially stable
  NoUnstableRegion=0;
  while (NoUnstableRegion==0){
    
    //2.1) copy LocalExtrema + treat the TmpImag
    for (z=0;z<TmpImag->NZ;z++) for (y=0;y<TmpImag->NY;y++) for (x=0;x<TmpImag->NX;x++){
      if (fabs(LocalExtrema->G(x,y,z))>epsilon)
          TmpImag->P(1,x,y,z);
      else
          TmpImag->P(0,x,y,z);
      }
    
    //2.2) smooth the TmpImag
    TmpConvolver->Convolution(TmpImag);
   
    //2.3) identify regions which would lead to unstabilities (denominator is too close to zero)
    NoUnstableRegion=1;
    for (z=0;z<OutputConvol->NZ;z++) for (y=0;y<OutputConvol->NY;y++) for (x=0;x<OutputConvol->NX;x++){
      if (fabs(TmpImag->G(x,y,z))<threshZero){
        TmpImag->P(1,x,y,z);
        if ((x>0)&&(x<OutputConvol->NX-1)&&(y>0)&&(y<OutputConvol->NY-1)&&(z>0)&&(z<OutputConvol->NZ-1)) NoUnstableRegion=0;
        }
      else
        TmpImag->P(0,x,y,z);
      }
      
    //2.4) Add fake LocalExtrema if necessary (i.e. control points)
    if (NoUnstableRegion==0){
      ShrinkRegionsToSinglePoints(TmpImag,1,OutputConvol);  // OutputConvol is used here as a temporary image which gets the center points
      AddSimulatedLocalExtrema(LocalExtrema,TmpImag,OutputConvol);
      }
    }
  
  //3) perform the actual normalized convolution
  //3.1) copy LocalExtrema to the future numerators and denominators + treat the TmpImag
  DeepCopy(LocalExtrema,OutputConvol);
  
  for (z=0;z<TmpImag->NZ;z++) for (y=0;y<TmpImag->NY;y++) for (x=0;x<TmpImag->NX;x++){
    if (fabs(OutputConvol->G(x,y,z))>epsilon)
      TmpImag->P(1,x,y,z);
    else
      TmpImag->P(0,x,y,z);
    }

  //3.2) smooth the numerator
  TmpConvolver->Convolution(OutputConvol);
  
  //3.3) smooth the denominator
  TmpConvolver->Convolution(TmpImag);
  
  //3.4) normalize the convolution where we do not have to divide by zero
  for (z=0;z<OutputConvol->NZ;z++) for (y=0;y<OutputConvol->NY;y++) for (x=0;x<OutputConvol->NX;x++){
      OutputConvol->P(OutputConvol->G(x,y,z)/TmpImag->G(x,y,z),x,y,z);
      if (fabs(TmpImag->G(x,y,z))<threshZero*0.1){
       cout << "Nan expected" << endl;
      }
      }
  }
  
  
  
  
  
///Estimate the IMF from the input image
void EstimateIMF(ScalarField * InputImag,int NbIt,ScalarField * TmpImag1,ScalarField * LocalMax,ScalarField * LocalMin,ScalarField * TmpImag4,LightFFTconvolver3D * TmpConvolver,ScalarField * OutputIMF){
  int it;
  float EstimatedSigma;
  char FileName[256];

  int x,y,z;

  
  //1) copy InputImag as the 1st estimate of the IMF
  DeepCopy(InputImag,OutputIMF);
  
  
  //2) main loop
  for (it=0;it<NbIt;it++){
    cout << "Iteration " << it+1 << " -> ";
    //2.1) extract local extrema of the current estimate of the IMF
    ExtractLocalExtrema2(OutputIMF,LocalMax,LocalMin);
    
    //strcpy(FileName,"LocalMax.nii");
    //TmpImag2->Write(FileName);
    //strcpy(FileName,"LocalMin.nii");
    //TmpImag3->Write(FileName);
  
    //2.2) estimate a good sigma to perform normalized convolution
    EstimatedSigma=EstimateGoodSigma(LocalMax,LocalMin);
    cout << "  Estimated Sigma = " << EstimatedSigma;
    
    //2.3) perform the  normalized convolution
    NormalizedConvolution(LocalMax,InputImag,TmpConvolver,TmpImag1,EstimatedSigma,TmpImag4);
    DeepCopy(TmpImag4,LocalMax);
    NormalizedConvolution(LocalMin,InputImag,TmpConvolver,TmpImag1,EstimatedSigma,TmpImag4);
    DeepCopy(TmpImag4,LocalMin);

    //strcpy(FileName,"NC_Max.nii");
    //TmpImag2->Write(FileName);
    //strcpy(FileName,"NC_Min.nii");
    //TmpImag3->Write(FileName);

    
    //2.4) remove low freq from the estimated IMF
    for (z=0;z<OutputIMF->NZ;z++) for (y=0;y<OutputIMF->NY;y++) for (x=0;x<OutputIMF->NX;x++)
      OutputIMF->Add(-(LocalMax->G(x,y,z)+LocalMin->G(x,y,z))/2,x,y,z);
    
     cout << endl;
    }
  cout << endl;
  }


  
///Perform an EMD
void PerformEMD(char InputImagFile[256],char OutputPrefix[256]){
  ScalarField InputImag;
  ScalarField TmpImag;
  ScalarField TmpImag2;
  ScalarField LocalMax;
  ScalarField LocalMin;
  ScalarField OutputIMF;
  ScalarField Residue;
  LightFFTconvolver3D TmpConvolver;
  char OutputFile[256];
  char Suffix1[256];
  char Suffix2[256];
  char Suffix3[256];
  float EstimatedSigma;
  int NbIterForSifting,NbModesTocompute;
  int i,x,y,z;
  
  NbIterForSifting=1;
  NbModesTocompute=3;
  
  //1) read input image and allocate memory for other images 
  InputImag.Read(InputImagFile);
  TmpImag.Read(InputImagFile);
  TmpImag2.Read(InputImagFile);
  LocalMax.Read(InputImagFile);
  LocalMin.Read(InputImagFile);
  OutputIMF.Read(InputImagFile);
  
  //2) perform EMD
  for (i=0;i<NbModesTocompute;i++){
    //2.1) estimate the IMF
    EstimateIMF(&InputImag,NbIterForSifting,&TmpImag,&LocalMax,&LocalMin,&TmpImag2,&TmpConvolver,&OutputIMF);
    
    //2.2) save current IMF
    strcpy(Suffix1,"_IMF");
    if (i==0) strcpy(Suffix2,"1");
    if (i==1) strcpy(Suffix2,"2");
    if (i==2) strcpy(Suffix2,"3");
    if (i==3) strcpy(Suffix2,"4");
    if (i==4) strcpy(Suffix2,"5");
    if (i==5) strcpy(Suffix2,"6");
    if (i==6) strcpy(Suffix2,"7");
    if (i==7) strcpy(Suffix2,"8");
    if (i==8) strcpy(Suffix2,"9");
    strcpy(Suffix3,".nii");
    
    strcpy(OutputFile,OutputPrefix);
    strcat(OutputFile,Suffix1);
    strcat(OutputFile,Suffix2);
    strcat(OutputFile,Suffix3);
    
    OutputIMF.Write(OutputFile);
    
    //2.3) compute the residue
    for (z=0;z<OutputIMF.NZ;z++) for (y=0;y<OutputIMF.NY;y++) for (x=0;x<OutputIMF.NX;x++)
      InputImag.Add(-OutputIMF.G(x,y,z),x,y,z);
  
    //TmpConvolver.InitiateConvolver(OutputIMF.NX,OutputIMF.NY,OutputIMF.NZ,1,0.5,0.5,0.5);  //to remove
    //TmpConvolver.Convolution(&InputImag); //to remove

  
  }
  
  //3) save the residue  
  strcpy(Suffix1,"_residue.nii");
  strcpy(OutputFile,OutputPrefix);
  strcat(OutputFile,Suffix1);
  InputImag.Write(OutputFile);
  
  }






///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                     12: metamorphoses
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///update CurrentMetamo using a pseudo-gradient method to manage the non-differentiability at 0 -- here the regularisation is everywhere the same (lambdaMetamo)
void PseudoGradMetamoUpdate(ScalarField * CurrentMetamo,ScalarField * UpdateForcesMetamo,float lambdaMetamo,float refMultFactor,float ThreshZero){
    int x,y,z;
    float signMF;
    
    for (z = 0; z < CurrentMetamo->NZ; z++)  for (y = 0; y < CurrentMetamo->NY; y++) for (x = 0; x < CurrentMetamo->NX; x++){
        if (fabs(CurrentMetamo->G(x,y,z))<ThreshZero){ //1) MetamoField is zero
            if (-UpdateForcesMetamo->G(x,y,z)>lambdaMetamo)
                CurrentMetamo->Add(refMultFactor*(-UpdateForcesMetamo->G(x,y,z)-lambdaMetamo),x,y,z);
            if (-UpdateForcesMetamo->G(x,y,z)<-lambdaMetamo)
                CurrentMetamo->Add(refMultFactor*(-UpdateForcesMetamo->G(x,y,z)+lambdaMetamo),x,y,z);
        }
        else{ //2) MetamoField is non-null
            signMF=(1-2*(static_cast<float>(CurrentMetamo->G(x,y,z)<0)));
            CurrentMetamo->Add(refMultFactor*(-UpdateForcesMetamo->G(x,y,z)-signMF*lambdaMetamo),x,y,z);
        }
    }
}

///update CurrentMetamo using a pseudo-gradient method to manage the non-differentiability at 0 -- here the regularisation depends on space
void PseudoGradMetamoUpdate(ScalarField * CurrentMetamo,ScalarField * UpdateForcesMetamo,ScalarField *  lambdaMetamoField,float refMultFactor,float ThreshZero){
    int x,y,z;
    float signMF;
    
    for (z = 0; z < CurrentMetamo->NZ; z++)  for (y = 0; y < CurrentMetamo->NY; y++) for (x = 0; x < CurrentMetamo->NX; x++){
        if (fabs(CurrentMetamo->G(x,y,z))<ThreshZero){ //1) MetamoField is zero
            if (-UpdateForcesMetamo->G(x,y,z)>lambdaMetamoField->G(x,y,z))
                CurrentMetamo->Add(refMultFactor*(-UpdateForcesMetamo->G(x,y,z)-lambdaMetamoField->G(x,y,z)),x,y,z);
            if (-UpdateForcesMetamo->G(x,y,z)<-lambdaMetamoField->G(x,y,z))
                CurrentMetamo->Add(refMultFactor*(-UpdateForcesMetamo->G(x,y,z)+lambdaMetamoField->G(x,y,z)),x,y,z);
        }
        else{ //2) MetamoField is non-null
            signMF=(1-2*(static_cast<float>(CurrentMetamo->G(x,y,z)<0)));
            CurrentMetamo->Add(refMultFactor*(-UpdateForcesMetamo->G(x,y,z)-signMF*lambdaMetamoField->G(x,y,z)),x,y,z);
        }
    }
}


/// For metamorphoses, compute regularisation forces using the 1st order derivatives of CurrentMetamo
void CptMetamoReg1stOrderDerivatives(ScalarField * CurrentMetamo,ScalarField *  lambdaMetamoField,ScalarField * tempFl,ScalarField * tempFl2,float locLambda){
    int x,y,z;
    float locSign;
    float Dxp,Dxm,Dyp,Dym,Dzp,Dzm;
    float mmDxpDxm,mmDypDym,mmDzpDzm;
    float nDx,nDy,nDz;
    float epsilon;
    
    epsilon=0.000001;
    
    lambdaMetamoField->PutToAllVoxels(0);
    
    //1) x-direction
    tempFl->PutToAllVoxels(0);
    tempFl2->PutToAllVoxels(0);
     
    for (z = 1; z < CurrentMetamo->NZ-1; z++) for (y = 1; y < CurrentMetamo->NY-1; y++) for (x = 0; x < CurrentMetamo->NX-1; x++){
                Dxp=CurrentMetamo->G(x+1,y,z)-CurrentMetamo->G(x,y,z);
                Dyp=CurrentMetamo->G(x,y+1,z)-CurrentMetamo->G(x,y,z);
                Dym=CurrentMetamo->G(x,y,z)-CurrentMetamo->G(x,y-1,z);
                Dzp=CurrentMetamo->G(x,y,z+1)-CurrentMetamo->G(x,y,z);
                Dzm=CurrentMetamo->G(x,y,z)-CurrentMetamo->G(x,y,z-1);
                
                mmDypDym=minmodfunc(Dyp,Dym);
                mmDzpDzm=minmodfunc(Dzp,Dzm);
                
                nDx=sqrt((Dxp*Dxp)+(mmDypDym*mmDypDym)+(mmDzpDzm*mmDzpDzm)+epsilon);
                
                tempFl->P(locLambda*(Dxp/nDx),x,y,z);
                tempFl2->P(Dxp,x,y,z);
            }
    
    for (z = 1; z < CurrentMetamo->NZ-1; z++) for (y = 1; y < CurrentMetamo->NY-1; y++)  for (x = 0; x < CurrentMetamo->NX-1; x++){
      if (fabs(tempFl->G(x,y,z)-tempFl->G(x-1,y,z))>fabs(tempFl2->G(x,y,z)-tempFl2->G(x-1,y,z))/4)
        lambdaMetamoField->Add(-((tempFl2->G(x,y,z)-tempFl2->G(x-1,y,z))/4),x,y,z);                  
      else                                                                                      
        lambdaMetamoField->Add(-(tempFl->G(x,y,z)-tempFl->G(x-1,y,z)),x,y,z);                      
    }                                                                                           
    
    
    //2) y-direction
    tempFl->PutToAllVoxels(0);
    tempFl2->PutToAllVoxels(0);

    for (z = 1; z < CurrentMetamo->NZ-1; z++) for (y = 0; y < CurrentMetamo->NY-1; y++) for (x = 1; x < CurrentMetamo->NX-1; x++){
                Dxp=CurrentMetamo->G(x+1,y,z)-CurrentMetamo->G(x,y,z);
                Dxm=CurrentMetamo->G(x,y,z)-CurrentMetamo->G(x-1,y,z);
                Dyp=CurrentMetamo->G(x,y+1,z)-CurrentMetamo->G(x,y,z);
                Dzp=CurrentMetamo->G(x,y,z+1)-CurrentMetamo->G(x,y,z);
                Dzm=CurrentMetamo->G(x,y,z)-CurrentMetamo->G(x,y,z-1);
                
                mmDxpDxm=minmodfunc(Dxp,Dxm);
                mmDzpDzm=minmodfunc(Dzp,Dzm);
                
                nDy=sqrt((Dyp*Dyp)+(mmDxpDxm*mmDxpDxm)+(mmDzpDzm*mmDzpDzm)+epsilon);
                
                tempFl->P(locLambda*(Dyp/nDy),x,y,z);
                tempFl2->P(Dyp,x,y,z);
            }
    
    for (z = 1; z < CurrentMetamo->NZ-1; z++) for (y = 1; y < CurrentMetamo->NY-1; y++)  for (x = 0; x < CurrentMetamo->NX-1; x++){
      if (fabs(tempFl->G(x,y,z)-tempFl->G(x,y-1,z))>fabs(tempFl2->G(x,y,z)-tempFl2->G(x,y-1,z))/4)
        lambdaMetamoField->Add(-((tempFl2->G(x,y,z)-tempFl2->G(x,y-1,z))/4),x,y,z);                  
      else                                                                                      
        lambdaMetamoField->Add(-(tempFl->G(x,y,z)-tempFl->G(x,y-1,z)),x,y,z);                      
    }                                                                                           
    
    
    //3) z-direction
    tempFl->PutToAllVoxels(0);
    tempFl2->PutToAllVoxels(0);

    for (z = 0; z < CurrentMetamo->NZ-1; z++) for (y = 1; y < CurrentMetamo->NY-1; y++) for (x = 1; x < CurrentMetamo->NX-1; x++){
                Dxp=CurrentMetamo->G(x+1,y,z)-CurrentMetamo->G(x,y,z);
                Dxm=CurrentMetamo->G(x,y,z)-CurrentMetamo->G(x-1,y,z);
                Dyp=CurrentMetamo->G(x,y+1,z)-CurrentMetamo->G(x,y,z);
                Dym=CurrentMetamo->G(x,y,z)-CurrentMetamo->G(x,y-1,z);
                Dzp=CurrentMetamo->G(x,y,z+1)-CurrentMetamo->G(x,y,z);
                
                mmDxpDxm=minmodfunc(Dxp,Dxm);
                mmDypDym=minmodfunc(Dyp,Dym);
                
                nDz=sqrt((Dzp*Dzp)+(mmDypDym*mmDypDym)+(mmDxpDxm*mmDxpDxm)+epsilon);
                
                tempFl->P(locLambda*(Dzp/nDz),x,y,z);
                tempFl2->P(Dzp,x,y,z);
            }
    
    for (z = 1; z < CurrentMetamo->NZ-1; z++) for (y = 1; y < CurrentMetamo->NY-1; y++)  for (x = 0; x < CurrentMetamo->NX-1; x++){
      if (fabs(tempFl->G(x,y,z)-tempFl->G(x,y,z-1))>fabs(tempFl2->G(x,y,z)-tempFl2->G(x,y,z-1))/4)
        lambdaMetamoField->Add(-((tempFl2->G(x,y,z)-tempFl2->G(x,y,z-1))/4),x,y,z);                  
      else                                                                                      
        lambdaMetamoField->Add(-(tempFl->G(x,y,z)-tempFl->G(x,y,z-1)),x,y,z);                      
    }                                                                                           
    
    
}



