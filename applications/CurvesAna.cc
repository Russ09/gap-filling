/*=========================================================================

Author: Laurent Risser

Disclaimer: This software has been developed for research purposes only, and hence should 
not be used as a diagnostic tool. In no event shall the authors or distributors
be liable to any direct, indirect, special, incidental, or consequential 
damages arising of the use of this software, its documentation, or any 
derivatives thereof, even if the authors have been advised of the possibility 
of such damage. 

=========================================================================*/


#include <SciCalcPack.h>


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 NEW CLASS FOR THE TENSOR VOTING (TV functions outside of the class are related to the 2009 algorithm)
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


//class definition
class TensorVoting2{
private:
  
protected:
  /// Protected functions 
  virtual void GetNetworkInClass(LDMK_Curves *, LDMK_Points *);
  virtual void InitFields(int,int,int);
  virtual void InitLists();
  virtual void FillStaticLists();
  virtual void GenerateEnhancementField(float,int);
  virtual void InsertAllStickFields(float,float,int);
  virtual void InsertStickField(int,float,float,int);
  virtual void InsertPlateChainAnchors(float);
  virtual void CptSaliencyMap();
  virtual void FillTheGaps();
  virtual void TestPath(int IdSegEnd,int *TargetType, int *ID_target, int *PathLength,float LengthPathSteps=1, float ConsideredGLThresh=0.05, float NoiseLevel=0);
  virtual void TestNoisyPaths(int,int *, int *,int *,float);
  virtual void IdentifyAndValidateSegmentAndElement(int,int,int *,int *,int);
  virtual void IdentifyAndValidateSegmentEnd(int,int,int,int *);
  virtual int  CreatePath(int,int,int,int,int,int);
  virtual void CalcParamStickField(float,float,float *,float *,int *);
  
  /// Protected parameters
  TensorField TensField;
  ScalarField Functional;
  ScalarField EnhancementField;
  ScalarField denominatorField;
  ScalarField NominatorField;
  IntScalarField NetworkLabels;
  
  LDMK_Curves * Network;
  LDMK_Points * NonNullGLs;
  
  int MaxPathLength;
  int SizeSmallSegment;
  int SizeSmallExtremity;
  
  int SizeListIslands;
  int SizeListSegmentEnds;
  
  int * ListFilledGaps;
  unsigned char * NbConnectedSegEnds;
  int ** ListSegmentEnds;
  int * ListIslands;
  float ** Path;
  float ** TmpPath;
  
  
public:
  
  LDMK_Curves DeltaNetwork;  // all new segments  (can be usefull to visualize what has been done)
  
  /// Constructor
  TensorVoting2();
  
  /// Destructor
  ~TensorVoting2();
  
  /// Run the tensor voting
  virtual void PerformTensorVoting2(LDMK_Curves * InputNetwork, LDMK_Points *  Input_NonNullGLs,float dista,float angl);
};



//constructor and destructor
TensorVoting2::TensorVoting2(void){}

TensorVoting2::~TensorVoting2(void){}


//copy by address the networks in the class
void TensorVoting2::GetNetworkInClass(LDMK_Curves * InputNetwork, LDMK_Points *  Input_NonNullGLs){
  this->Network=InputNetwork;
  this->NonNullGLs=Input_NonNullGLs;
  this->DeltaNetwork.Create_LDMK_N_VoidCurves(10000);
  
  }


//initialize all fields of the class
void TensorVoting2::InitFields(int NX,int NY,int NZ){
    this->TensField.CreateVoidField(NX,NY,NZ);
    this->Functional.CreateVoidField(NX,NY,NZ);
    this->EnhancementField.CreateVoidField(NX,NY,NZ);
    this->NetworkLabels.CreateVoidField(NX,NY,NZ);
    this->denominatorField.CreateVoidField(NX,NY,NZ);
    this->NominatorField.CreateVoidField(NX,NY,NZ);
}

//allocate memory for the lists of the class
void TensorVoting2::InitLists(){
    int i;
    
    
    //1) general list representing filled gaps
    //-> ListFilledGaps[4*i] = segment of node 1
    //-> ListFilledGaps[4*i+1] = element of the segment at node 1
    //-> ListFilledGaps[4*i+2] = segment of node 2
    //-> ListFilledGaps[4*i+3] = element of the segment at node 2
    this->ListFilledGaps = new int[4*this->Network->GetSegNumber()];  
    for (i=0;i<4*this->Network->GetSegNumber();i++) this->ListFilledGaps[i]=-1;
    
    
    //2) lists the number of nodes connected to a segment end
    // NbConnectedSegEnds[2*i]     -> nb of segments connected to 1st extremity of segment i
    // NbConnectedSegEnds[2*i+1]   -> nb of segments connected to 2nd extremity of segment i
    this->NbConnectedSegEnds = new unsigned char[2*this->Network->GetSegNumber()];

    //3) list of the segment ends
    //ListSegmentEnds[i][0] -> ID of the segment
    //ListSegmentEnds[i][1] -> corresponding element in the segment
    this->ListSegmentEnds = new int*[2*this->Network->GetSegNumber()];
    for(i=0;i<2*(this->Network->GetSegNumber());i++) 
      this->ListSegmentEnds[i] = new int[2];
    
    //4) list of isolated islands
    this->ListIslands = new int[this->Network->GetSegNumber()];

    //5) memory for the tested paths
    //Path[i][0] -> X   |  Path[i][1] -> Y   |  Path[i][2] -> Z
    this->Path = new float*[this->MaxPathLength];
    for (i=0;i<this->MaxPathLength;i++)
      this->Path[i] = new float[3];  
    
    this->TmpPath = new float*[this->MaxPathLength];
    for (i=0;i<this->MaxPathLength;i++)
      this->TmpPath[i] = new float[3];  
    
}


//fill the lists which represent the initial newtork
void TensorVoting2::FillStaticLists(){
    int i,j;
    int LocElInSegI,LocElInSegJ,WhichSegEndI,WhichSegEndJ;
    int Test;
    

    //1) count the nunber of ngbh at each segment end
    for (i=0;i<this->Network->GetSegNumber();i++){
        this->NbConnectedSegEnds[2*i]=0;
        this->NbConnectedSegEnds[2*i+1]=0;
    }
    
    for (i=0;i<this->Network->GetSegNumber()-1;i++) if (this->Network->GetElNumber(i)!=0) for(WhichSegEndI=0;WhichSegEndI<2;WhichSegEndI++){
        LocElInSegI=(this->Network->GetElNumber(i)-1)*WhichSegEndI;
        for (j=i+1;j<this->Network->GetSegNumber();j++) if (this->Network->GetElNumber(j)!=0) for(WhichSegEndJ=0;WhichSegEndJ<2;WhichSegEndJ++){
            LocElInSegJ=(this->Network->GetElNumber(j)-1)*WhichSegEndJ;
            if (fabs(this->Network->GetX(i,LocElInSegI)-this->Network->GetX(j,LocElInSegJ))<0.1)
                if (fabs(this->Network->GetY(i,LocElInSegI)-this->Network->GetY(j,LocElInSegJ))<0.1)
                    if (fabs(this->Network->GetZ(i,LocElInSegI)-this->Network->GetZ(j,LocElInSegJ))<0.1){
                        this->NbConnectedSegEnds[2*i+WhichSegEndI]++;
                        this->NbConnectedSegEnds[2*j+WhichSegEndJ]++;
                    }
        }
    }
    
    
    //2) generate the list of islands (used for the enhanced map only!!!)
    this->SizeListIslands=0;
    for (i=0;i<this->Network->GetSegNumber();i++) if ((this->Network->GetElNumber(i)!=0)&&(this->Network->GetElNumber(i)<=this->SizeSmallSegment)) if ((this->NbConnectedSegEnds[2*i]==0)&&(this->NbConnectedSegEnds[2*i+1]==0)){
        this->ListIslands[this->SizeListIslands]=i;
        this->SizeListIslands++;
    }
    
    
    //3) generate the list of segment ends
    this->SizeListSegmentEnds=0;
    for (i=0;i<this->Network->GetSegNumber();i++) if (this->Network->GetElNumber(i)>=this->SizeSmallExtremity) for(WhichSegEndI=0;WhichSegEndI<2;WhichSegEndI++) if(this->NbConnectedSegEnds[2*i+WhichSegEndI]==0){
        Test=0;
        for(j=0;j<this->SizeListIslands;j++) if (this->ListIslands[j]==i) Test=1;
        
        if (Test==0){
            this->ListSegmentEnds[this->SizeListSegmentEnds][0]=i;
            this->ListSegmentEnds[this->SizeListSegmentEnds][1]=(this->Network->GetElNumber(i)-1)*WhichSegEndI;
            this->SizeListSegmentEnds++;
        }
    }
}



//Insert a stick voting field in the tensor field
void TensorVoting2::InsertStickField(int IdSegEnd, float C, float sigma,int BoxSizes){
    int i,j,k;
    float LocVec[3];
    float V_x1,V_y1,V_z1;
    float TempD,Dist;
    int LocX,LocY,LocZ;
    float ArcLength,Curvature,weight,phi;
    int CX,CY,CZ;
    float DX,DY,DZ;
    int x2,y2,z2,z3,y3,x3;
    float AngleIndicator;
    int HalfBoxSizes;
    float ThetaLim,AngleLim;
    
    // 1 ) UPDATE THE TENSOR FIELD

    //1.1) Parameters of the cone in which the tensor field will be written
    ThetaLim=fabs(asin(pow(sigma,2)/(2*sqrt(C))))*180/3.14159;
    if (ThetaLim<90) AngleLim=0.038;
    if (ThetaLim<60) AngleLim=0.270;
    if (ThetaLim<45) AngleLim=0.558;
    if (ThetaLim<30) AngleLim=0.924;

    
    //1.2) define the origin and direction of the stick field
    if (this->Network->GetElNumber(this->ListSegmentEnds[IdSegEnd][0])>=4){ //large segment
        if (this->ListSegmentEnds[IdSegEnd][1]==0){
            CX=static_cast<int>(this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],2)+0.5);
            CY=static_cast<int>(this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],2)+0.5);
            CZ=static_cast<int>(this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],2)+0.5);
            DX=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],3);
            DY=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],3);
            DZ=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],3);
        }
        else{
            CX=static_cast<int>(this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-3)+0.5);
            CY=static_cast<int>(this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-3)+0.5);
            CZ=static_cast<int>(this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-3)+0.5);
            DX=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-4);
            DY=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-4);
            DZ=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-4);
        }
    }
    else{ //small segment
        if (this->ListSegmentEnds[IdSegEnd][1]==0){
            CX=static_cast<int>(this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],0)+0.5);
            CY=static_cast<int>(this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],0)+0.5);
            CZ=static_cast<int>(this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],0)+0.5);
            DX=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],1);
            DY=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],1);
            DZ=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],1);
        }
        else{
            CX=static_cast<int>(this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)+0.5);
            CY=static_cast<int>(this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)+0.5);
            CZ=static_cast<int>(this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)+0.5);
            DX=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-2);
            DY=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-2);
            DZ=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-2);
        }
    }
    
    TempD=sqrt(DX*DX+DY*DY+DZ*DZ);
    DX=DX/TempD;
    DY=DY/TempD;
    DZ=DZ/TempD;
    
    
    //1.3) fill the tensor field
    HalfBoxSizes=static_cast<int>((BoxSizes-1)/2);
    
    for (i=-HalfBoxSizes;i<HalfBoxSizes;i++) for (j=-HalfBoxSizes;j<HalfBoxSizes;j++) for (k=-HalfBoxSizes;k<HalfBoxSizes;k++){
        //1.3.1) Init
        LocX=CX+k;
        LocY=CY+j;
        LocZ=CZ+i;
        
        Dist=sqrt(static_cast<float>((k*k)+(j*j)+(i*i)));
        V_x1=(static_cast<float>(k))/Dist;
        V_y1=(static_cast<float>(j))/Dist;
        V_z1=(static_cast<float>(i))/Dist;
        
        AngleIndicator=DX*V_x1+DY*V_y1+DZ*V_z1;
        
        if ((LocX>0)&&(LocX<this->TensField.NX)&&(LocY>0)&&(LocY<this->TensField.NY)&&(LocZ>0)&&(LocZ<this->TensField.NZ)&&((AngleIndicator>AngleLim)||(Dist<3))){
            // 1.3.2 ) define the vector to add in the tensor field
            TempD=V_x1*DX+V_y1*DY+V_z1*DZ;
            
            LocVec[0]=2*V_x1*TempD-DX;
            LocVec[1]=2*V_y1*TempD-DY;
            LocVec[2]=2*V_z1*TempD-DZ;
            
            // 1.3.3 ) define the local weight
            
            TempD=V_x1*DX+V_y1*DY+V_z1*DZ;
            
            if (TempD<0){weight=0;}
            else{
                if (TempD>1)  TempD=1;
                phi=acos(TempD);
                Curvature=2*sin(phi)/Dist;
                if (phi<0.01) ArcLength=Dist;
                else ArcLength=Dist*phi/sin(phi);
                weight=exp( -(pow(ArcLength,2)+C*pow(Curvature,2))   /  pow(sigma,2)  );
            }
            
            // 1.3.4 ) inject the vector and the weight in the tensor field
            this->TensField.AddTensorisedVector(LocVec,LocX,LocY,LocZ,0,weight);
        }
    }
    
    // 2 ) UPDATE THE FIELD REPRESENTING THE LABELS 
    for (z2=-3;z2<4;z2++) for (y2=-3;y2<4;y2++) for (x2=-3;x2<4;x2++){
        z3=static_cast<int>(this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1])+0.5)+z2;
        y3=static_cast<int>(this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1])+0.5)+y2;
        x3=static_cast<int>(this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1])+0.5)+x2;
        if ((z3>=0)&&(z3<this->TensField.NZ)&&(y3>=0)&&(y3<this->TensField.NY)&&(x3>=0)&&(x3<this->TensField.NX))
            this->NetworkLabels.P(100000000+IdSegEnd,x3,y3,z3);
    }
}



#ifdef COMPILE_WITH_OPENMP

void TensorVoting2::InsertAllStickFields(float C, float sigma,int BoxSizes){
  int i;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i) 
  {
    #pragma omp for
    for(i=0;i<this->SizeListSegmentEnds;i++){
      this->InsertStickField(i,C,sigma,BoxSizes);
    }
  //END FORK FOR THREADS
  }  
}



#else


void TensorVoting2::InsertAllStickFields(float C, float sigma,int BoxSizes){
  int i;
  
  for(i=0;i<this->SizeListSegmentEnds;i++)
    this->InsertStickField(i,C,sigma,BoxSizes);
}

#endif





//Insert the information out of the islands in this->EnhancementField   (related to the 2015 TV algorithm)
//Remark that all islands of 'ListIslands' (from Network) and all points of 'NonNullGLs' are considered for the enhancement field only
void TensorVoting2::GenerateEnhancementField(float sigma,int BoxSizes){
    int i,j,k;
    float Dist,LocWeight;
    int LocX,LocY,LocZ;
    int CX,CY,CZ;
    int IslandId;
    float MaxGL_NonNullGLs;
    float RatioNominatorDenominatorSigmaForNC;
    int HalfBoxSizes;
    float weightForIslands;
    LightFFTconvolver3D LocSmoother;
    
    //init
    RatioNominatorDenominatorSigmaForNC=0.5;   // must be in ]0,1[
    weightForIslands=0.3;  // must be in [0,1]
    HalfBoxSizes=(BoxSizes-1)/2;  //divided by 2 to fit within the for loop
    
    
    this->denominatorField.PutToAllVoxels(0);
    this->NominatorField.PutToAllVoxels(0);
    
    
    //1) generate the field which will be used as a denominator for normalized convolution - part 1: using the islands of ListIslands
    for(IslandId=0;IslandId<this->SizeListIslands;IslandId++){
        CX=static_cast<int>(this->Network->GetX(this->ListIslands[IslandId],this->Network->GetElNumber(this->ListIslands[IslandId])/2)+0.5);
        CY=static_cast<int>(this->Network->GetY(this->ListIslands[IslandId],this->Network->GetElNumber(this->ListIslands[IslandId])/2)+0.5);
        CZ=static_cast<int>(this->Network->GetZ(this->ListIslands[IslandId],this->Network->GetElNumber(this->ListIslands[IslandId])/2)+0.5);
        
        if ((CX>0)&&(CX<this->TensField.NX)&&(CY>0)&&(CY<this->TensField.NY)&&(CZ>0)&&(CZ<this->TensField.NZ))
          this->denominatorField.Add(static_cast<float>(1),CX,CY,CZ);
    }
    
    
    //2) generate the field which will be used as a denominator for normalized convolution - part 2: using the islands of NonNullGLs
    
    //2.1) small preprocessing to make sure that values of are refWeight in [0,1]
    MaxGL_NonNullGLs=fabs(this->NonNullGLs->GetVal(0));
    
    for(IslandId=0;IslandId<this->NonNullGLs->Get_LDMK_PointsNumber();IslandId++)
        if (MaxGL_NonNullGLs<fabs(this->NonNullGLs->GetVal(IslandId)))
            MaxGL_NonNullGLs=fabs(this->NonNullGLs->GetVal(IslandId));
    
    //2.2) fill the field
    for(IslandId=0;IslandId<this->NonNullGLs->Get_LDMK_PointsNumber();IslandId++){
        CX=static_cast<int>(this->NonNullGLs->GetX(IslandId)+0.5);
        CY=static_cast<int>(this->NonNullGLs->GetY(IslandId)+0.5);
        CZ=static_cast<int>(this->NonNullGLs->GetZ(IslandId)+0.5);
        //LocWeight=fabs(this->NonNullGLs->GetVal(IslandId))/MaxGL_NonNullGLs;  // not proper normalized convolution
        LocWeight=1; 
        
        
        if ((CX>0)&&(CX<this->TensField.NX)&&(CY>0)&&(CY<this->TensField.NY)&&(CZ>0)&&(CZ<this->TensField.NZ))
          this->denominatorField.Add(LocWeight,CX,CY,CZ);
    }
    
    //3) convolve denominatorField
    LocSmoother.InitiateConvolver(this->TensField.NX,this->TensField.NY,this->TensField.NZ,1,sigma,sigma,sigma);
    LocSmoother.Convolution(&this->denominatorField);
    
    
    //4) generate the field which will be used as a NominatorField for normalized convolution - part 1: using the islands of ListIslands

    
    for(IslandId=0;IslandId<this->SizeListIslands;IslandId++){
        CX=static_cast<int>(this->Network->GetX(this->ListIslands[IslandId],this->Network->GetElNumber(this->ListIslands[IslandId])/2)+0.5);
        CY=static_cast<int>(this->Network->GetY(this->ListIslands[IslandId],this->Network->GetElNumber(this->ListIslands[IslandId])/2)+0.5);
        CZ=static_cast<int>(this->Network->GetZ(this->ListIslands[IslandId],this->Network->GetElNumber(this->ListIslands[IslandId])/2)+0.5);
        LocWeight=weightForIslands;
        
        if ((CX>0)&&(CX<this->TensField.NX)&&(CY>0)&&(CY<this->TensField.NY)&&(CZ>0)&&(CZ<this->TensField.NZ))
          this->NominatorField.Add(LocWeight,CX,CY,CZ);
    }

    
    
    //5) generate the field which will be used as a NominatorField for normalized convolution - part 2: using the islands of NonNullGLs
    for(IslandId=0;IslandId<this->NonNullGLs->Get_LDMK_PointsNumber();IslandId++){
        CX=static_cast<int>(this->NonNullGLs->GetX(IslandId)+0.5);
        CY=static_cast<int>(this->NonNullGLs->GetY(IslandId)+0.5);
        CZ=static_cast<int>(this->NonNullGLs->GetZ(IslandId)+0.5);
        LocWeight=weightForIslands*fabs(this->NonNullGLs->GetVal(IslandId))/MaxGL_NonNullGLs;   //to make sure values are in [0,weightForIslands]
        
        if ((CX>0)&&(CX<this->TensField.NX)&&(CY>0)&&(CY<this->TensField.NY)&&(CZ>0)&&(CZ<this->TensField.NZ))
          this->NominatorField.Add(LocWeight,CX,CY,CZ);
    }
    
    //6) convolve NominatorField
    LocSmoother.ChangeKernel(sigma*RatioNominatorDenominatorSigmaForNC,sigma*RatioNominatorDenominatorSigmaForNC,sigma*RatioNominatorDenominatorSigmaForNC);
    LocSmoother.Convolution(&this->NominatorField);
    
    
    //7) compute and normalize the EnhancementField
    for (i=0;i<this->TensField.NZ;i++) for (j=0;j<this->TensField.NY;j++) for (k=0;k<this->TensField.NX;k++) if (this->denominatorField.G(k,j,i)>0.001)
        this->EnhancementField.Add(this->NominatorField.G(k,j,i)/this->denominatorField.G(k,j,i),k,j,i);
        
    
    float minGL,maxGL,a,b,tmp;
    minGL=EnhancementField.G(0,0,0);
    maxGL=EnhancementField.G(0,0,0);
    for (i=0;i<EnhancementField.NZ;i++) for (j=0;j<EnhancementField.NY;j++) for (k=0;k<EnhancementField.NX;k++) {
      if (EnhancementField.G(k,j,i)<minGL) minGL=EnhancementField.G(k,j,i);
      if (EnhancementField.G(k,j,i)>maxGL) maxGL=EnhancementField.G(k,j,i);
    }
  
    a=(0-1)/(minGL-maxGL);
    b=1-a*maxGL;
  
    if (minGL==maxGL){
      a=0;
      b=0;
    }
  
    for (i=0;i<EnhancementField.NZ;i++) for (j=0;j<EnhancementField.NY;j++) for (k=0;k<EnhancementField.NX;k++){
      tmp=a*EnhancementField.G(k,j,i)+b;
      EnhancementField.P(tmp,k,j,i);
    }
    
}




#ifdef COMPILE_WITH_OPENMP
 
void TensorVoting2::InsertPlateChainAnchors(float sigma){
  int i,j,k,z2,y2,x2,z3,y3,x3;
  int TempInt4NL,TempInt4EF;
  float TempFl,TempFl2;
  int   radiutInt4NL,radiutInt4EF;
  float radiusFl4NL ,radiusFl4EF;
  float Sx,Sy,Sz,Sxn,Syn,Szn,x2n,y2n,z2n;
  float LocWeight;
  int MaxDist;
  float weightForChainAnchors;
  
  
  //1) initialize general parameters
  MaxDist=7;   //must be higher than 3
  weightForChainAnchors=0.5;   //must be in [0,1]
  
  if (this->Network->GetSegNumber()>=100000000){
      printf("The network contains too much segments\n");
      return;
  }
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j,k,z2,y2,x2,z3,y3,x3,TempInt4NL,TempInt4EF,TempFl,TempFl2,radiutInt4NL,radiutInt4EF,radiusFl4NL,radiusFl4EF,Sx,Sy,Sz,Sxn,Syn,Szn,x2n,y2n,z2n,LocWeight) 
  {
    
    // 2 ) fill the field representing the labels
    #pragma omp for
    for (i=0;i<this->Network->GetSegNumber();i++) if (this->Network->GetElNumber(i)>2*MaxDist+4){
        for (j=MaxDist+2;j<this->Network->GetElNumber(i)-MaxDist-2;j++){
            //2.1) tangent to the skeleton at the current point
            Sx=this->Network->GetX(i,j+3)-this->Network->GetX(i,j-3);
            Sy=this->Network->GetY(i,j+3)-this->Network->GetY(i,j-3);
            Sz=this->Network->GetZ(i,j+3)-this->Network->GetZ(i,j-3);
            Sxn=Sx/(sqrt(pow(Sx,2)+pow(Sy,2)+pow(Sz,2)));
            Syn=Sy/(sqrt(pow(Sx,2)+pow(Sy,2)+pow(Sz,2)));
            Szn=Sz/(sqrt(pow(Sx,2)+pow(Sy,2)+pow(Sz,2)));
            

            //2.2) create two spheres around the current point -- radii are different for the enhancement field and the network labels
            //2.2.1) Spheres radii
            radiutInt4NL=static_cast<int>(this->Network->GetD(i,j)/2+1);  //radius for network labels
            if (radiutInt4NL<1) radiutInt4NL=1;
            radiusFl4NL=static_cast<float>(radiutInt4NL);
            
            radiutInt4EF=static_cast<int>((sigma/2)+(this->Network->GetD(i,j)/2)+1); //radius for enhancement field (larger or equal to the one for the network labels)
            if (radiutInt4EF<1) radiutInt4EF=1;
            radiusFl4EF=static_cast<float>(radiutInt4EF);
            
            for (z2=-radiutInt4EF;z2<radiutInt4EF;z2++)
                for (y2=-radiutInt4EF;y2<radiutInt4EF;y2++)
                    for (x2=-radiutInt4EF;x2<radiutInt4EF;x2++){
                        //2.2.2) define the tested point and the parameters which will check whether they should be written or not in the fields
                        z3=static_cast<int>(this->Network->GetZ(i,j))+z2;
                        y3=static_cast<int>(this->Network->GetY(i,j))+y2;
                        x3=static_cast<int>(this->Network->GetX(i,j))+x2;
                        
                        TempInt4NL=0;
                        TempInt4EF=0;
                        
                        //2.2.3) is the point in the image?
                        if (!((x3>=0)&&(x3<this->TensField.NX)&&(y3>=0)&&(y3<this->TensField.NY)&&(z3>=0)&&(z3<this->TensField.NZ))) {
                          TempInt4NL=1; 
                          TempInt4EF=1;
                          }
                        
                        //is the point in the sphere?
                        if ((TempInt4NL==0)||(TempInt4EF==0)){
                          TempFl=sqrt(pow(static_cast<float>(z2),2)+pow(static_cast<float>(y2),2)+pow(static_cast<float>(x2),2));
                          if (TempFl>radiusFl4NL) TempInt4NL=1;
                          if (TempFl>radiusFl4EF) TempInt4EF=1;
                        }
                        
                        //is the point eligible relative to the tangent to the skeleton
                        if (((TempInt4NL==0)||(TempInt4EF==0))&&(TempFl>3)){
                            TempFl2=(Sxn*static_cast<float>(x2)+Syn*static_cast<float>(y2)+Szn*static_cast<float>(z2))/TempFl;
                            if (TempFl2*TempFl2>0.04){  //angle of about 10 degrees with the skeleton tangent
                              TempInt4NL=1;
                              TempInt4EF=1;
                              }
                        }
                        
                        //point is OK -> can be reached to fill a gap
                        if (TempInt4NL==0){
                          this->NetworkLabels.P(200000000+i,x3,y3,z3);
                        }
                        
                        if (TempInt4EF==0){
                            TempFl-=(this->Network->GetD(i,j)/2);
                            if (TempFl<1) TempFl=1;
                            LocWeight=exp(-(TempFl*TempFl)/(sigma*sigma))*weightForChainAnchors;
                            
                            if (this->EnhancementField.G(x3,y3,z3)<LocWeight)
                                this->EnhancementField.P(LocWeight,x3,y3,z3);
                        }
                    }
        }
    }
  
  //END FORK FOR THREADS
  }

}


#else


void TensorVoting2::InsertPlateChainAnchors(float sigma){
    int i,j,k,z2,y2,x2,z3,y3,x3;
    int TempInt4NL,TempInt4EF;
    float TempFl,TempFl2;
    int MaxDist;
    int   radiutInt4NL,radiutInt4EF;
    float radiusFl4NL ,radiusFl4EF;
    float Sx,Sy,Sz,Sxn,Syn,Szn,x2n,y2n,z2n;
    float LocWeight;
    float weightForChainAnchors;
    
    
    //1) initialize general parameters
    MaxDist=7;   //must be higher than 3
    weightForChainAnchors=0.5;   //must be in [0,1]
    
    if (this->Network->GetSegNumber()>=100000000){
        printf("The network contains too much segments\n");
        return;
    }
    
    // 2 ) fill the field representing the labels
    for (i=0;i<this->Network->GetSegNumber();i++) if (this->Network->GetElNumber(i)>2*MaxDist+4){
        for (j=MaxDist+2;j<this->Network->GetElNumber(i)-MaxDist-2;j++){
            //2.1) tangent to the skeleton at the current point
            Sx=this->Network->GetX(i,j+3)-this->Network->GetX(i,j-3);
            Sy=this->Network->GetY(i,j+3)-this->Network->GetY(i,j-3);
            Sz=this->Network->GetZ(i,j+3)-this->Network->GetZ(i,j-3);
            Sxn=Sx/(sqrt(pow(Sx,2)+pow(Sy,2)+pow(Sz,2)));
            Syn=Sy/(sqrt(pow(Sx,2)+pow(Sy,2)+pow(Sz,2)));
            Szn=Sz/(sqrt(pow(Sx,2)+pow(Sy,2)+pow(Sz,2)));
            

            //2.2) create two spheres around the current point -- radii are different for the enhancement field and the network labels
            //2.2.1) Spheres radii
            radiutInt4NL=static_cast<int>(this->Network->GetD(i,j)/2+1);  //radius for network labels
            if (radiutInt4NL<1) radiutInt4NL=1;
            radiusFl4NL=static_cast<float>(radiutInt4NL);
            
            radiutInt4EF=static_cast<int>((sigma/2)+(this->Network->GetD(i,j)/2)+1); //radius for enhancement field (larger or equal to the one for the network labels)
            if (radiutInt4EF<1) radiutInt4EF=1;
            radiusFl4EF=static_cast<float>(radiutInt4EF);
            
            for (z2=-radiutInt4EF;z2<radiutInt4EF;z2++)
                for (y2=-radiutInt4EF;y2<radiutInt4EF;y2++)
                    for (x2=-radiutInt4EF;x2<radiutInt4EF;x2++){
                        //2.2.2) define the tested point and the parameters which will check whether they should be written or not in the fields
                        z3=static_cast<int>(this->Network->GetZ(i,j))+z2;
                        y3=static_cast<int>(this->Network->GetY(i,j))+y2;
                        x3=static_cast<int>(this->Network->GetX(i,j))+x2;
                        
                        TempInt4NL=0;
                        TempInt4EF=0;
                        
                        //2.2.3) is the point in the image?
                        if (!((x3>=0)&&(x3<this->TensField.NX)&&(y3>=0)&&(y3<this->TensField.NY)&&(z3>=0)&&(z3<this->TensField.NZ))) {
                          TempInt4NL=1; 
                          TempInt4EF=1;
                          }
                        
                        //is the point in the sphere?
                        if ((TempInt4NL==0)||(TempInt4EF==0)){
                          TempFl=sqrt(pow(static_cast<float>(z2),2)+pow(static_cast<float>(y2),2)+pow(static_cast<float>(x2),2));
                          if (TempFl>radiusFl4NL) TempInt4NL=1;
                          if (TempFl>radiusFl4EF) TempInt4EF=1;
                        }
                        
                        //is the point eligible relative to the tangent to the skeleton
                        if (((TempInt4NL==0)||(TempInt4EF==0))&&(TempFl>3)){
                            TempFl2=(Sxn*static_cast<float>(x2)+Syn*static_cast<float>(y2)+Szn*static_cast<float>(z2))/TempFl;
                            if (TempFl2*TempFl2>0.04){  //angle of about 10 degrees with the skeleton tangent
                              TempInt4NL=1;
                              TempInt4EF=1;
                              }
                        }
                        
                        //point is OK -> can be reached to fill a gap
                        if (TempInt4NL==0){
                          this->NetworkLabels.P(200000000+i,x3,y3,z3);
                        }
                        
                        if (TempInt4EF==0){
                            TempFl-=(this->Network->GetD(i,j)/2);
                            if (TempFl<1) TempFl=1;
                            LocWeight=exp(-(TempFl*TempFl)/(sigma*sigma))*weightForChainAnchors;
                            
                            if (this->EnhancementField.G(x3,y3,z3)<LocWeight)
                                this->EnhancementField.P(LocWeight,x3,y3,z3);
                        }
                    }
        }
    }
}



#endif


#ifdef COMPILE_WITH_OPENMP
 
//Compute the saliency map to a curve
void TensorVoting2::CptSaliencyMap(){
  int i,j,k;
  float vec1[3];
  float vec2[3];
  float vec3[3];
  float lambda1,lambda2,lambda3;
  
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(i,j,k,vec1,vec2,vec3,lambda1,lambda2,lambda3) 
  {
    
    // 2 ) fill the field representing the labels
    #pragma omp for
    for (j=0;j<this->TensField.NY;j++) for (i=0;i<this->TensField.NZ;i++) for (k=0;k<this->TensField.NX;k++){
        if ((this->TensField.G(0,0,k,j,i)>0.000001)||(this->TensField.G(1,1,k,j,i)>0.000001)||(this->TensField.G(2,2,k,j,i)>0.000001)){
            
            //perform the PCA
            this->TensField.PCA(vec1,vec2,vec3,&lambda1,&lambda2,&lambda3,k,j,i);
            
            //fill the functional with the saliency map
            this->Functional.P(lambda1-lambda2,k,j,i);
            
            //first eigenvector is also stored
            this->TensField.P(vec1[0],0,0,k,j,i);
            this->TensField.P(vec1[1],1,1,k,j,i);
            this->TensField.P(vec1[2],2,2,k,j,i);
        }
        else{
          this->Functional.P(0,k,j,i);
        }
    }
  
  //END FORK FOR THREADS
  }

} 

#else

 
//Compute the saliency map to a curve
void TensorVoting2::CptSaliencyMap(){
    int i,j,k;
    float vec1[3];
    float vec2[3];
    float vec3[3];
    float lambda1,lambda2,lambda3;
    
    
    for (i=0;i<this->TensField.NZ;i++) for (j=0;j<this->TensField.NY;j++) for (k=0;k<this->TensField.NX;k++){
        if ((this->TensField.G(0,0,k,j,i)>0.000001)||(this->TensField.G(1,1,k,j,i)>0.000001)||(this->TensField.G(2,2,k,j,i)>0.000001)){
            
            //perform the PCA
            this->TensField.PCA(vec1,vec2,vec3,&lambda1,&lambda2,&lambda3,k,j,i);
            
            //fill the functional with the saliency map
            this->Functional.P(lambda1-lambda2,k,j,i);
            
            //first eigenvector is also stored
            this->TensField.P(vec1[0],0,0,k,j,i);
            this->TensField.P(vec1[1],1,1,k,j,i);
            this->TensField.P(vec1[2],2,2,k,j,i);
        }
        else{
          this->Functional.P(0,k,j,i);
        }
    }
} 


#endif


//From a segment-end 'IdSegEnd' different optimal paths of local maxima 
//  -> 'TargetType' is -1 if nothing is found / 1 if a segment end is found / 2 if a segment is found
//  -> LengthPathSteps  not be higher than 1 -- In our tests, 1 seems to give optimal results (probably due to the fact that the path should be sufficiently far away from the origin after a few steps
//  -> ConsideredGLThresh... when defining a path, this search is stopped if reaching this value in the functional
void TensorVoting2::TestPath(int IdSegEnd,int *TargetType, int *ID_target, int *PathLength,float LengthPathSteps, float ConsideredGLThresh,float NoiseLevel){
    float locX,locY,locZ;
    int TXi,TYi,TZi;
    float TXf,TYf,TZf;
    int i,j,k;
    float ValFoncTop;
    int StopSearch,FirstIterations;
    float DirecX,DirecY,DirecZ;
    float CurrentDirection[3];
    float TestedDirection[3];
    float TopTestedDirection[3];
    float tvec1[3];
    float tvec2[3];
    float IniX,IniY,IniZ;
    float theta,phi;
    int TempI;
    float tmpFl,LocAngl;
    int EligiblePoint;
    int IterationsNumber;
    int LXi,LYi,LZi;
    int IniLabel;
    float frandmax;
    float rand1,rand2,rand3;
    float u,v,r,c;
    
    frandmax=static_cast<float>(RAND_MAX);
    
    // 1 : initialisations...
    if (this->ListSegmentEnds[IdSegEnd][1]==0){  // first element of the segment
        IniX=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],0);
        IniY=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],0);
        IniZ=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],0);
        CurrentDirection[0]=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],1);
        CurrentDirection[1]=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],1);
        CurrentDirection[2]=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],0)-this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],1);
    }
    else{   // last element of the segment
        IniX=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1);
        IniY=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1);
        IniZ=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1);
        CurrentDirection[0]=this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetX(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-2);
        CurrentDirection[1]=this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetY(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-2);
        CurrentDirection[2]=this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1)-this->Network->GetZ(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-2);
    }
    
    VecNormalize(CurrentDirection,1);

    if ((IniX<0)||(IniX>this->TensField.NX-1)||(IniY<0)||(IniY>this->TensField.NY-1)||(IniZ<0)||(IniZ>this->TensField.NZ-1)) {*PathLength=0; *TargetType=-1; return;}
    
    *ID_target=-1;
    *TargetType=-1;
    
    //2) follow the path...
    //2.1) two first points
    locZ=IniZ+2*CurrentDirection[2]; locY=IniY+2*CurrentDirection[1]; locX=IniX+2*CurrentDirection[0];
    
    Path[0][0]=IniX; Path[0][1]=IniY; Path[0][2]=IniZ;
    Path[1][0]=locX; Path[1][1]=locY; Path[1][2]=locZ;
    IniLabel=this->NetworkLabels.G(static_cast<int>(IniX+0.5),static_cast<int>(IniY+0.5),static_cast<int>(IniZ+0.5));
    
    //2.2) iterative process for other points
    StopSearch=0;
    FirstIterations=4;
    IterationsNumber=2;
    
    while (((this->Functional.G(locX,locY,locZ)>ConsideredGLThresh)||(this->EnhancementField.G(locX,locY,locZ)>ConsideredGLThresh)||(FirstIterations>0))&&(StopSearch==0)){
        
        //2.2.1) Init parameters for the current point
        
        if ((locX<0)||(locX>this->TensField.NX-1)||(locY<0)||(locY>this->TensField.NY-1)||(locZ<0)||(locZ>this->TensField.NZ-1)) {*PathLength=0; *TargetType=-1; return;}
        
        LXi=static_cast<int>(locX+0.5);
        LYi=static_cast<int>(locY+0.5);
        LZi=static_cast<int>(locZ+0.5);
        
        CptVecsTangentPlane(CurrentDirection,tvec1,tvec2);
        VecNormalize(CurrentDirection,1);
        VecNormalize(tvec1,1);
        VecNormalize(tvec2,1);
        
        ValFoncTop=-1;
        
        //2.2.2) Test potential next points
        for (theta=-24*3.1416/180;theta<24.1*3.1416/180;theta+=4*3.1416/180)
            for (phi=-3.1416;phi<3.1416;phi+=10*3.1416/180){
            
                //2.2.2.1) define the tested direction and point
                TestedDirection[0]=CurrentDirection[0]*cos(theta)+tvec1[0]*sin(theta)*sin(phi)+tvec2[0]*sin(theta)*cos(phi);
                TestedDirection[1]=CurrentDirection[1]*cos(theta)+tvec1[1]*sin(theta)*sin(phi)+tvec2[1]*sin(theta)*cos(phi);
                TestedDirection[2]=CurrentDirection[2]*cos(theta)+tvec1[2]*sin(theta)*sin(phi)+tvec2[2]*sin(theta)*cos(phi);
                
                VecNormalize(TestedDirection,LengthPathSteps);
                TXi=static_cast<int>(locX+TestedDirection[0]+0.5);
                TYi=static_cast<int>(locY+TestedDirection[1]+0.5);
                TZi=static_cast<int>(locZ+TestedDirection[2]+0.5);
                TXf=locX+TestedDirection[0];
                TYf=locY+TestedDirection[1];
                TZf=locZ+TestedDirection[2];
                
                if ((TXi>2)&&(TXi<this->TensField.NX-3)&&(TYi>2)&&(TYi<this->TensField.NY-3)&&(TZi>2)&&(TZi<this->TensField.NZ-3)){
                    //2.2.2.2) does we reach something?  -- if yes the iterative search stops (see the return lines)
                    TempI=this->NetworkLabels.G(TXi,TYi,TZi);
                    for (i=-2;i<3;i++) for (j=-2;j<3;j++) for (k=-2;k<3;k++)
                        if (TempI<this->NetworkLabels.G(TXi+k,TYi+j,TZi+i))
                            TempI=this->NetworkLabels.G(TXi+k,TYi+j,TZi+i);
                    
                    if ((FirstIterations<0)&&(TempI!=IniLabel)){
                        if ((TempI>=100000000)&&(TempI<200000000)){ //a segment end is reached
                            cout << "End of segment " << TempI << " is reached" << " at point (" << TXi << "," << TYi  << "," << TZi << ")" << endl;
                            *TargetType=1;
                            *ID_target=TempI;
                            *PathLength=IterationsNumber;
                            return;
                        }
                        if ((TempI>=200000000)&&(TempI<300000000)){ //a segment is reached
                            cout << "Segment " << TempI << " is reached" << " at point (" << TXi << "," << TYi  << "," << TZi << ")" << endl;
                            *TargetType=2;
                            *ID_target=TempI;
                            *PathLength=IterationsNumber;
                            return;
                        }
                    }
                    
                    //2.2.2.3) if nothing is reached, check whether the tested point is eligible for the path?
                    
                    EligiblePoint=1;
                    VecNormalize(TestedDirection,1);   //normalize with one for computations
                    
                    //...comparison of the tested direction with the average direction of eigenvector 1
                    tmpFl= (this->TensField.G(0,0,TXi,TYi,TZi)+this->TensField.G(0,0,LXi,LYi,LZi))*TestedDirection[0]/2;
                    tmpFl+=(this->TensField.G(1,1,TXi,TYi,TZi)+this->TensField.G(1,1,LXi,LYi,LZi))*TestedDirection[1]/2;
                    tmpFl+=(this->TensField.G(2,2,TXi,TYi,TZi)+this->TensField.G(2,2,LXi,LYi,LZi))*TestedDirection[2]/2;
                    
                    tmpFl=fabs(tmpFl);
                    
                    tmpFl/=sqrt((this->TensField.G(0,0,TXi,TYi,TZi)*this->TensField.G(0,0,TXi,TYi,TZi))+(this->TensField.G(1,1,TXi,TYi,TZi)*this->TensField.G(1,1,TXi,TYi,TZi))+(this->TensField.G(2,2,TXi,TYi,TZi)*this->TensField.G(2,2,TXi,TYi,TZi)));
                    
                    VecNormalize(TestedDirection,LengthPathSteps);  //renormalize with the path step length
                    
                    if (tmpFl>=1) LocAngl=0;
                    else if (tmpFl<=0) LocAngl=90;
                    else LocAngl=acos(tmpFl)*180/3.1416;
                    
                    if (LocAngl>20) EligiblePoint=0;
                    
                    if ((this->Functional.G(TXf,TYf,TZf)<ConsideredGLThresh)&&(this->EnhancementField.G(TXf,TYf,TZf)<ConsideredGLThresh)) EligiblePoint=0;
                    
                    //2.2.2.4) measure the score of the tested point if eligible
                    if (EligiblePoint==1){
                        //rand1=NoiseLevel*(sampleUniform()-0.5);
                        //rand2=NoiseLevel*(sampleUniform()-0.5);
                        //rand3=NoiseLevel*(sampleUniform()-0.5);
                        if (IterationsNumber>=7){
                          rand1=NoiseLevel*sampleNormal();
                          rand2=NoiseLevel*sampleNormal();
                          rand3=NoiseLevel*sampleNormal();
                        }
                        else{
                          rand1=0;
                          rand2=0;
                          rand3=0;
                        }
                        
                        
                        if (ValFoncTop<0){
                            ValFoncTop=this->Functional.G(TXf,TYf,TZf);
                            if (ValFoncTop<this->EnhancementField.G(TXf,TYf,TZf))
                                ValFoncTop=this->EnhancementField.G(TXf,TYf,TZf);
                                
                            TopTestedDirection[0]=TestedDirection[0]+rand1;
                            TopTestedDirection[1]=TestedDirection[1]+rand2;  //expression of the noise (null by default) here
                            TopTestedDirection[2]=TestedDirection[2]+rand3;
                        }
                        if (ValFoncTop<this->EnhancementField.G(TXf,TYf,TZf)){
                            ValFoncTop=this->EnhancementField.G(TXf,TYf,TZf);
                            TopTestedDirection[0]=TestedDirection[0]+rand1;
                            TopTestedDirection[1]=TestedDirection[1]+rand2;  //expression of the noise (null by default) here
                            TopTestedDirection[2]=TestedDirection[2]+rand3;
                        }
                        if (ValFoncTop<this->Functional.G(TXf,TYf,TZf)){
                            ValFoncTop=this->Functional.G(TXf,TYf,TZf);
                            TopTestedDirection[0]=TestedDirection[0]+rand1;
                            TopTestedDirection[1]=TestedDirection[1]+rand2;  //expression of the noise (null by default) here
                            TopTestedDirection[2]=TestedDirection[2]+rand3;
                        }
                    }
                }
            }  //end of for theta phi loop
        
        //update stuffs on the path search
        if (ValFoncTop<0){ *PathLength=IterationsNumber; *ID_target=-1; *TargetType=-1; return;}
        else{
            locZ+=TopTestedDirection[2]; locY+=TopTestedDirection[1]; locX+=TopTestedDirection[0];
            Path[IterationsNumber][0]=locX; Path[IterationsNumber][1]=locY; Path[IterationsNumber][2]=locZ;
            
            CurrentDirection[0]=TopTestedDirection[0];
            CurrentDirection[1]=TopTestedDirection[1];
            CurrentDirection[2]=TopTestedDirection[2];
            //cout <<  "Segment end " << this->ListSegmentEnds[IdSegEnd][0] << " (" << this->ListSegmentEnds[IdSegEnd][1] << "): " << locX  << " " << locY  << " " << locZ  << endl;
            FirstIterations--;
            IterationsNumber++;
            if (IterationsNumber>=this->MaxPathLength) { *PathLength=IterationsNumber; *ID_target=-1; *TargetType=-1; return;}
        }
    }  //end of while loop for path definition
    
    //3) if we're out of the path loop, that means that nothing was reached
    *PathLength=IterationsNumber;
    *ID_target=-1;
    *TargetType=-1;
    return;
}


//Robust version of 'TestPath' where several noisy paths are tested and not a single deterministic one.
//From a segment-end 'IdSegEnd' different optimal paths of local maxima
//Contrary to 'TestPath' the tested paths are noisy
//-> 'TargetType' is -1 if nothing is found / 1 if a segment end is found / 2 if a segment is found
//-> The std dev of the Gaussian noise on the path is 'PropNoise' times the length of the path steps
void TensorVoting2::TestNoisyPaths(int IdSegEnd,int *TargetType, int *ID_target, int *PathLength,float PropNoise){
  float ConsideredGLThresh;
  float LengthPathSteps;
  float NoiseLevel;
  int NbPathsToTest,NbPathsTested;
  int Best_TargetType;
  int Best_ID_target;
  int Best_PathLength;
  int Tst_TargetType;
  int Tst_ID_target;
  int Tst_PathLength;
  int i;
  int EligiblePath;
  float D1[3];
  float D2[3];
  float diam1,diam2;
  int IdReachedSegEnd;
  
  //1) general parameters
  ConsideredGLThresh=0.05;   //when defining a path, this search is stopped if reaching this value in the functional
  LengthPathSteps=1;         // should not be higher than 1 -- In our tests, 1 seems to give optimal results (probably due to the fact that the path should be sufficiently far away from the origin after a few steps)
  NoiseLevel=LengthPathSteps*PropNoise;
  NbPathsToTest=50;

  //2) test different paths
  Best_PathLength=-1;
  srand(time(0));
  
  for (NbPathsTested=0; NbPathsTested<NbPathsToTest;NbPathsTested++){
    Tst_TargetType=-1;
    this->TestPath(IdSegEnd,&Tst_TargetType,&Tst_ID_target,&Tst_PathLength,LengthPathSteps,ConsideredGLThresh,NoiseLevel);
    
    if ((Tst_TargetType==1)||(Tst_TargetType==2)){ //something is reached
      if ((Best_PathLength==-1)||(Best_PathLength>Tst_PathLength)){ //no optimal path yet or shorter than the optimal path
        //Additional tests here to check whether the path is eligible...
        EligiblePath=1;
        
        IdReachedSegEnd=Tst_ID_target-100000000;
        
        //... angle between the path and the reached segment end is not too large
        if (Tst_TargetType==1){
          D1[0]=this->Path[Tst_PathLength-1][0]-this->Path[Tst_PathLength-4][0];
          D1[1]=this->Path[Tst_PathLength-1][1]-this->Path[Tst_PathLength-4][1];
          D1[2]=this->Path[Tst_PathLength-1][2]-this->Path[Tst_PathLength-4][2];
          
          
          
          if (this->ListSegmentEnds[IdReachedSegEnd][1]==0){
            D2[0]=this->Network->GetX(this->ListSegmentEnds[IdReachedSegEnd][0],0)-this->Network->GetX(this->ListSegmentEnds[IdReachedSegEnd][0],3);
            D2[1]=this->Network->GetY(this->ListSegmentEnds[IdReachedSegEnd][0],0)-this->Network->GetY(this->ListSegmentEnds[IdReachedSegEnd][0],3);
            D2[2]=this->Network->GetZ(this->ListSegmentEnds[IdReachedSegEnd][0],0)-this->Network->GetZ(this->ListSegmentEnds[IdReachedSegEnd][0],3);
          }
          else{
            D2[0]=this->Network->GetX(this->ListSegmentEnds[IdReachedSegEnd][0],this->ListSegmentEnds[IdReachedSegEnd][1]-1)-this->Network->GetX(this->ListSegmentEnds[IdReachedSegEnd][0],this->ListSegmentEnds[IdReachedSegEnd][1]-4);
            D2[1]=this->Network->GetY(this->ListSegmentEnds[IdReachedSegEnd][0],this->ListSegmentEnds[IdReachedSegEnd][1]-1)-this->Network->GetY(this->ListSegmentEnds[IdReachedSegEnd][0],this->ListSegmentEnds[IdReachedSegEnd][1]-4);
            D2[2]=this->Network->GetZ(this->ListSegmentEnds[IdReachedSegEnd][0],this->ListSegmentEnds[IdReachedSegEnd][1]-1)-this->Network->GetZ(this->ListSegmentEnds[IdReachedSegEnd][0],this->ListSegmentEnds[IdReachedSegEnd][1]-4);
          }
          
          VecNormalize(D1,1);
          VecNormalize(D2,1);
          if (scalarProd_3vec(D1,D2)>-0.5) EligiblePath=0;   //-0.5 -> angle = Pi/6    --   -0.707 -> angle = Pi/4
        }
        
        //... diameters should not be too different
        if (Tst_TargetType==1){
          if (this->ListSegmentEnds[IdSegEnd][1]==0) diam1=this->Network->GetD(this->ListSegmentEnds[IdSegEnd][0],0);
          else                                       diam1=this->Network->GetD(this->ListSegmentEnds[IdSegEnd][0],this->ListSegmentEnds[IdSegEnd][1]-1);
          
          if (this->ListSegmentEnds[IdReachedSegEnd][1]==0) diam2=this->Network->GetD(this->ListSegmentEnds[IdReachedSegEnd][0],0);
          else                                              diam2=this->Network->GetD(this->ListSegmentEnds[IdReachedSegEnd][0],this->ListSegmentEnds[IdReachedSegEnd][1]-1);
          
          if (diam2<0.33*diam1) EligiblePath=0;
          if (diam2>3*diam1)    EligiblePath=0;
          }
        
        if (EligiblePath==1){
          Best_TargetType=Tst_TargetType;
          Best_ID_target=Tst_ID_target;
          Best_PathLength=Tst_PathLength;
          for (i=0;i<this->MaxPathLength;i++) this->TmpPath[i][0]=this->Path[i][0];
          for (i=0;i<this->MaxPathLength;i++) this->TmpPath[i][1]=this->Path[i][1];
          for (i=0;i<this->MaxPathLength;i++) this->TmpPath[i][2]=this->Path[i][2];
          }
        }
      }
    }
      

  //3) return the best result
  if ((Best_TargetType==1)||(Best_TargetType==2)){
    *TargetType=Best_TargetType;
    *ID_target=Best_ID_target;
    *PathLength=Best_PathLength;
    for (i=0;i<this->MaxPathLength;i++) this->Path[i][0]=this->TmpPath[i][0];
    for (i=0;i<this->MaxPathLength;i++) this->Path[i][1]=this->TmpPath[i][1];
    for (i=0;i<this->MaxPathLength;i++) this->Path[i][2]=this->TmpPath[i][2];
    }
  else{
    *TargetType=-1;
    *ID_target=-1;
    *PathLength=-1;
  }
  
  return;
}




//Identify which segment is reached and at which element. Return -1 in IdSegmentReached if a loop has been computed 
void TensorVoting2::IdentifyAndValidateSegmentAndElement(int TargetId,int PathLength,int * IdSegmentReached,int * IdElementReached,int IdInitSegment){  
    int TopSegment;
    int i,j,TopElem;
    int InitSegment;
    int LoopTst;
    int LocElInSegI,LocElInSegJ;
    float EndX,EndY,EndZ;
    float TmpDist,TopDist;
    
    //init segment 
    InitSegment=this->ListSegmentEnds[IdInitSegment][0];
    
    //a segment is reached
    EndX=this->Path[PathLength-1][0];
    EndY=this->Path[PathLength-1][1];
    EndZ=this->Path[PathLength-1][2];
    
    TopSegment=TargetId-200000000;
    
    TopDist=100000000;
    TopElem=-1;
    for (i=2;i<this->Network->GetElNumber(TopSegment)-1;i++){
        TmpDist=pow(EndX-this->Network->GetX(TopSegment,i),2)+pow(EndY-this->Network->GetY(TopSegment,i),2)+pow(EndZ-this->Network->GetZ(TopSegment,i),2);
        if (TmpDist<TopDist){
            TopElem=i;
            TopDist=TmpDist;
        }
    }
    
    //test whether a loop is generated
    LoopTst=0;
    
    if (InitSegment==TopSegment) LoopTst=0;
    
    LocElInSegI=this->Network->GetElNumber(InitSegment)-1;
    LocElInSegJ=this->Network->GetElNumber(TopSegment)-1;
    
    if (fabs(this->Network->GetX(InitSegment,LocElInSegI)-this->Network->GetX(TopSegment,LocElInSegJ))<0.01)
        if (fabs(this->Network->GetY(InitSegment,LocElInSegI)-this->Network->GetY(TopSegment,LocElInSegJ))<0.01)
            if (fabs(this->Network->GetZ(InitSegment,LocElInSegI)-this->Network->GetZ(TopSegment,LocElInSegJ))<0.01) LoopTst=1;
    
    if (fabs(this->Network->GetX(InitSegment,LocElInSegI)-this->Network->GetX(TopSegment,0))<0.01)
        if (fabs(this->Network->GetY(InitSegment,LocElInSegI)-this->Network->GetY(TopSegment,0))<0.01)
            if (fabs(this->Network->GetZ(InitSegment,LocElInSegI)-this->Network->GetZ(TopSegment,0))<0.01) LoopTst=1;
    
    if (fabs(this->Network->GetX(InitSegment,0)-this->Network->GetX(TopSegment,LocElInSegJ))<0.01)
        if (fabs(this->Network->GetY(InitSegment,0)-this->Network->GetY(TopSegment,LocElInSegJ))<0.01)
            if (fabs(this->Network->GetZ(InitSegment,0)-this->Network->GetZ(TopSegment,LocElInSegJ))<0.01) LoopTst=1;
    
    if (fabs(this->Network->GetX(InitSegment,0)-this->Network->GetX(TopSegment,0))<0.01)
        if (fabs(this->Network->GetY(InitSegment,0)-this->Network->GetY(TopSegment,0))<0.01)
            if (fabs(this->Network->GetZ(InitSegment,0)-this->Network->GetZ(TopSegment,0))<0.01) LoopTst=1;
    
    
    //return the segment and element reached 
    
    if (LoopTst==0){
        *IdSegmentReached=TopSegment;
        *IdElementReached=TopElem;
    }
    else{
        *IdSegmentReached=-1;
        *IdElementReached=-1;
    }
}


//Identify which segment is reached. Return -1 in IdSegmentReached the segment end is already reached or if a loop has been computed 
void TensorVoting2::IdentifyAndValidateSegmentEnd(int NbOfGapsFilled,int ID_target,int InitSegEnd,int * IdSegmentReached){
    int j,k, Test;
    
    *IdSegmentReached=-1;
    
    j=ID_target-100000000;
    
    Test=0;
    //detect potential loops
    if (this->ListSegmentEnds[j][0]==this->ListSegmentEnds[InitSegEnd][0]) Test=1;
    
    //test whether the junction already exists
    if (Test==0)for (k=0;k<NbOfGapsFilled;k++){
        if ((this->ListSegmentEnds[j][0]==this->ListFilledGaps[4*k+2])&&(this->ListSegmentEnds[j][1]==this->ListFilledGaps[4*k+3])&&(this->ListSegmentEnds[InitSegEnd][0]==this->ListFilledGaps[4*k])&&(this->ListSegmentEnds[InitSegEnd][1]==this->ListFilledGaps[4*k+1])){
            Test=1;
        }
        if ((this->ListSegmentEnds[j][0]==this->ListFilledGaps[4*k])&&(this->ListSegmentEnds[j][1]==this->ListFilledGaps[4*k+1])&&(this->ListSegmentEnds[InitSegEnd][0]==this->ListFilledGaps[4*k+2])&&(this->ListSegmentEnds[InitSegEnd][1]==this->ListFilledGaps[4*k+3])){
            Test=1;
        }
    }
    
    if (Test==0) *IdSegmentReached=j;
    
}


//Create a path to fill a gap. Segments diameters are linearly interpolated
//EndType=3 -> segment   | EndType=2 -> segment end
int TensorVoting2::CreatePath(int InitSegment,int InitElement,int ReachedSegment,int ReachedElement,int EndType,int PathLength){
    int i,j,k;
    int CurrentSegment,CurrentDeltaSegment;
    float TpX,TpY,TpZ,TpX2,TpY2,TpZ2;
    float InitDiameter,FinalDiameter;
    
    
    //1) allocate memory for the new segment
    CurrentSegment=this->Network->GenerateVoidSegment(PathLength+1);
    
    if (CurrentSegment==-1)  return -1;
    
    //2) fill the path
    this->Network->PutX(this->Network->GetX(InitSegment,InitElement),CurrentSegment,0);
    this->Network->PutY(this->Network->GetY(InitSegment,InitElement),CurrentSegment,0);
    this->Network->PutZ(this->Network->GetZ(InitSegment,InitElement),CurrentSegment,0);
    this->Network->PutD(1,CurrentSegment,0);
    
    for (i=1;i<PathLength;i++){
        this->Network->PutX(this->Path[i][0],CurrentSegment,i);
        this->Network->PutY(this->Path[i][1],CurrentSegment,i);
        this->Network->PutZ(this->Path[i][2],CurrentSegment,i);
        this->Network->PutD(1,CurrentSegment,i);
    }
    
    this->Network->PutX(this->Network->GetX(ReachedSegment,ReachedElement),CurrentSegment,PathLength);
    this->Network->PutY(this->Network->GetY(ReachedSegment,ReachedElement),CurrentSegment,PathLength);
    this->Network->PutZ(this->Network->GetZ(ReachedSegment,ReachedElement),CurrentSegment,PathLength);
    this->Network->PutD(1,CurrentSegment,PathLength);
    
    //3.1) interpolate the diameters if another segment end is reached
    if (EndType==2){
        InitDiameter=this->Network->GetD(InitSegment,InitElement);
        FinalDiameter=this->Network->GetD(ReachedSegment,ReachedElement);
        for (i=0;i<this->Network->GetElNumber(CurrentSegment);i++){
            this->Network->PutD(InitDiameter+(static_cast<float>(i))/(static_cast<float>(PathLength-1))*(FinalDiameter-InitDiameter),CurrentSegment,i);
        }
    }
    
    //3.2) constant diametre if a segment is reached
    if (EndType==3){
        InitDiameter=this->Network->GetD(InitSegment,InitElement);
        for (i=0;i<this->Network->GetElNumber(CurrentSegment);i++){
            this->Network->PutD(InitDiameter,CurrentSegment,i);
        }
    }
    
    //4) slightly smooth the new segment
    this->Network->SmoothOneSegment(CurrentSegment,1);
    
    
    
    //5) Same thing with the DeltaNetwork network which only contains added segments
    CurrentDeltaSegment=this->DeltaNetwork.GenerateVoidSegment(PathLength+1);
    
    if (CurrentDeltaSegment==-1)  return -2;
    
    this->DeltaNetwork.PutX(this->Network->GetX(InitSegment,InitElement),CurrentDeltaSegment,0);
    this->DeltaNetwork.PutY(this->Network->GetY(InitSegment,InitElement),CurrentDeltaSegment,0);
    this->DeltaNetwork.PutZ(this->Network->GetZ(InitSegment,InitElement),CurrentDeltaSegment,0);
    this->DeltaNetwork.PutD(1,CurrentDeltaSegment,0);
    
    for (i=1;i<PathLength;i++){
        this->DeltaNetwork.PutX(this->Path[i][0],CurrentDeltaSegment,i);
        this->DeltaNetwork.PutY(this->Path[i][1],CurrentDeltaSegment,i);
        this->DeltaNetwork.PutZ(this->Path[i][2],CurrentDeltaSegment,i);
        this->DeltaNetwork.PutD(1,CurrentDeltaSegment,i);
    }
    
    this->DeltaNetwork.PutX(this->Network->GetX(ReachedSegment,ReachedElement),CurrentDeltaSegment,PathLength);
    this->DeltaNetwork.PutY(this->Network->GetY(ReachedSegment,ReachedElement),CurrentDeltaSegment,PathLength);
    this->DeltaNetwork.PutZ(this->Network->GetZ(ReachedSegment,ReachedElement),CurrentDeltaSegment,PathLength);
    this->DeltaNetwork.PutD(1,CurrentDeltaSegment,PathLength);
    
    if (EndType==2){
        InitDiameter=this->Network->GetD(InitSegment,InitElement);
        FinalDiameter=this->Network->GetD(ReachedSegment,ReachedElement);
        for (i=0;i<this->DeltaNetwork.GetElNumber(CurrentDeltaSegment);i++){
            this->DeltaNetwork.PutD(InitDiameter+(static_cast<float>(i))/(static_cast<float>(PathLength-1))*(FinalDiameter-InitDiameter),CurrentDeltaSegment,i);
        }
    }
    
    if (EndType==3){
        InitDiameter=this->Network->GetD(InitSegment,InitElement);
        for (i=0;i<this->DeltaNetwork.GetElNumber(CurrentDeltaSegment);i++){
            this->DeltaNetwork.PutD(InitDiameter,CurrentDeltaSegment,i);
        }
    }
    
    this->DeltaNetwork.SmoothOneSegment(CurrentDeltaSegment,1);
    
    return CurrentSegment;
}



//Computes 'c', 'sigma' et 'WindowSize' as a function of 'dist' and 'angl'.
// -> We lose 1/e of the energy at distance 'dist' from the tensor field origin in the direction pointed by the segment end
// -> We lose (1/e)^2 at the same distance but at an angle 'angl' from  the direction pointed by the segment end
// -> WindowSize is the recommended size of the window in which tensor fields will be computed 
void TensorVoting2::CalcParamStickField(float dist,float angl,float * c,float * sigma,int * WindowSize){
    float RefCurvature;
    
    if (angl<2) angl=2;
    angl=angl*3.14159/180;
    angl=fabs(angl);
    dist=fabs(dist);
    
    RefCurvature=2*sin(angl)/dist;
    
    *c=pow(dist,2)/pow(RefCurvature,2);
    *sigma=dist;
    *WindowSize=2*4.5*dist+1;
}

//fill the gaps after having initiated everything
void TensorVoting2::FillTheGaps(){
  int NbOfGapsFilled;
  int i,j;
  int TargetType,TargetId,PathLength,TopSegment,TopElement;
  
  NbOfGapsFilled=0;
  
  //launch a path from each segment end
  for(i=0;i<this->SizeListSegmentEnds;i++){
      //test a path  (the path is stored in this->Path)
      TargetType=-1;
      //this->TestPath(i,&TargetType,&TargetId,&PathLength);           //to comment for the noisy version of the algorithm
      this->TestNoisyPaths(i,&TargetType,&TargetId,&PathLength,0.2);     //to uncomment for the noisy version of the algorithm
      
      //a segment end is reached
      if (TargetType==1){
          this->IdentifyAndValidateSegmentEnd(NbOfGapsFilled,TargetId,i,&j);
          
          if (j!=-1){
              this->CreatePath(this->ListSegmentEnds[i][0],this->ListSegmentEnds[i][1],this->ListSegmentEnds[j][0],this->ListSegmentEnds[j][1],TargetType,PathLength);
              
              this->ListFilledGaps[4*NbOfGapsFilled]=this->ListSegmentEnds[i][0];
              this->ListFilledGaps[4*NbOfGapsFilled+1]=this->ListSegmentEnds[i][1];
              this->ListFilledGaps[4*NbOfGapsFilled+2]=this->ListSegmentEnds[j][0];
              this->ListFilledGaps[4*NbOfGapsFilled+3]=this->ListSegmentEnds[j][1];
              NbOfGapsFilled++;
          }
      }
      
      //a segment is reached
      if (TargetType==2){
          this->IdentifyAndValidateSegmentAndElement(TargetId,PathLength,&TopSegment,&TopElement,i);
           if(TopElement!=-1){
              this->CreatePath(this->ListSegmentEnds[i][0],this->ListSegmentEnds[i][1],TopSegment,TopElement,TargetType,PathLength);
          }
      }
  }
}



//Treat the lineset 'Network' using the 2015 alternative tensor voting method to  (Risser et al, TMI 2008)
// -> dist and angl:
//  --> An energy loss of 1/e is modeled at a distance 'dista' from the origin (point O) in the direction of the segment end (point A)
//  --> An energy loss of 1/e is modeled at a distance 'dista' from point O in a direction having an angle 'angl' with OA
void TensorVoting2::PerformTensorVoting2(LDMK_Curves * InputNetwork, LDMK_Points *  Input_NonNullGLs,float dista,float angl){
    int i,j;
    float C,sigma;
    int BoxSizes;
    float minX, minY, minZ;
    float WidthX,WidthY,WidthZ;
    
    //1 ) INIT
    printf("Tensor voting launched...\n");
    
    //1.1) COPY THE NETWORK INTO THE CLASS
    this->GetNetworkInClass(InputNetwork,Input_NonNullGLs);
    
    
    //1.2) PARAMETRIZE THE TENSOR VOTING
    this->SizeSmallSegment=5;  //size under which a segment with two isolated segment ends is considered as an island
    this->SizeSmallExtremity=2; //size under which a free segment end can be joined (must be >=2)
    this->MaxPathLength=100;        //max length of a path between two tokens
    
    
    if ((dista>0)&&(angl>=0)) this->CalcParamStickField(dista,angl,&C,&sigma,&BoxSizes); // semi-automatically tuned parameters
    else this->CalcParamStickField(10,45,&C,&sigma,&BoxSizes);
    
    cout <<  "dista=" << dista  <<  "    angl=" << angl  <<  "-->    C=" << C  <<  "   sigma=" << sigma << "  BoxSizes=" << BoxSizes  << endl;

    
    //1.3) ESTIMIATE NETWORK SIZE AND SHIFT IT SO THAT IT HAS A MARGIN OF 5 VOXELS
    this->Network->EstimateROI(&minX,&minY,&minZ,&WidthX,&WidthY,&WidthZ);
    
    WidthX+=10; WidthY+=10; WidthZ+=10; 
    minX-=5;    minY-=5;    minZ-=5;
    
    this->Network->Translate(-minX,-minY,-minZ);
    this->NonNullGLs->Translate(-minX,-minY,-minZ);

    
    //1.4) MEMORY ALLOCATION AND INITATE DIFFERENT LISTS
    
    //1.4.1) Allocate memory for the scalar fields and the tensor field
    this->InitFields(static_cast<int>(WidthX+0.5),static_cast<int>(WidthY+0.5),static_cast<int>(WidthZ+0.5));
    
    //1.4.2) Allocate memory for the different lists
    this->InitLists();
    
    //1.4.3) fill static lists which represent the initial network
    this->FillStaticLists();
    
    
    //2) CREATE THE DIFFERENT FIELDS WHICH WILL HELP TO FILL THE GAPS
    
    printf("Compute different fields and the functional\n");
    
    //2.1) CLEAN-UP EVERYTHING
    cout << "-> Clean-up the fields" << endl;
    this->TensField.PutAllValuesToZero();
    this->Functional.PutToAllVoxels(0);
    this->NetworkLabels.PutToAllVoxels(0);
    
    //2.2) GENERATE THE ENHANCEMENT FIELD USING THE ISLANDS
    cout << "-> GenerateEnhancementField" << endl;
    this->GenerateEnhancementField(sigma,BoxSizes);
    
    //2.3) INSERT ANCHORS FOR THE SEGMENTS
    cout << "-> InsertPlateChainAnchors" << endl;
    this->InsertPlateChainAnchors(sigma);

    //2.4) INSERT TENSOR FIELDS RELATED TO SEGMENT ENDS
    cout << "-> InsertAllStickFields" << endl;
    this->InsertAllStickFields(C,sigma,BoxSizes);
    
    
    //2.5) COMPUTE THE SALIENCY MAP
    cout << "-> CptSaliencyMap" << endl;
    this->CptSaliencyMap();
    
    //this->Functional.Write("SaliencyMap.nii");
    //this->EnhancementField.Write("EnhancementField.nii");
    //this->NetworkLabels.Write("NetworkLabels.nii");
    
    
    //3) SEGMENT END JONCTIONS   --   test the paths from each segment end
    printf("Gap filling\n");
    this->FillTheGaps();

    
    //4) FINALIZE THE COMPUTATIONS --  manage the initial shift to (5,5,5)
    this->Network->Translate(minX,minY,minZ);
    this->NonNullGLs->Translate(minX,minY,minZ);
    this->DeltaNetwork.Translate(minX,minY,minZ);
}




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                            Functions level 1
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



/// read parameters and perform the tensor voting strategy from 2015
void TensorVoting2Launcher(int argc, char **argv){     
    char SkeletonFileName[256];
    char IslandsFileName[256];
    char OutputFile[256];
    LDMK_Curves Network;
    LDMK_Points Islands;
    float CharacteristicLength,CharacteristicAngle;
    TensorVoting2 TV_Manager;
    char OriginalImageName[256];
    char Str_GL[256];
    char Str_Orig[256];
    int KnownOriginalImage;
    int KnownNonNullGL;
    ScalarField OriginalImage;
    
    

    //read parameters
    argc--; argv++;
    strcpy(SkeletonFileName,argv[1]);
    argc--; argv++;
    strcpy(OutputFile,argv[1]);
    argc--; argv++;
    CharacteristicLength=atof(argv[1]);
    argc--; argv++;
    CharacteristicAngle=atof(argv[1]);
    argc--; argv++;
    
    KnownOriginalImage=0;
    strcpy(Str_Orig,"-Orig");
    
    KnownNonNullGL=0;
    strcpy(Str_GL,"-GL");
    
    while (argc>1){
      if (argc>1) if (strcmp(argv[1],Str_Orig)==0) {
        KnownOriginalImage=1;
        argc--; argv++;
        strcpy(OriginalImageName,argv[1]);
        argc--; argv++;
      }
      if (argc>1) if (strcmp(argv[1],Str_GL)==0) {
        KnownNonNullGL=1;
        argc--; argv++;
        strcpy(IslandsFileName,argv[1]);
        argc--; argv++;
      }
    }
    
    
    //Read the mv3d files
    Network.Read(SkeletonFileName);
    
    
    if (KnownNonNullGL==1)
      Islands.Read(IslandsFileName);  //Islands is already initiated as a void LDMK_Points so no need for a 'else'
    
    if (KnownOriginalImage==1){
      OriginalImage.Read(OriginalImageName);
      Network.MillimetersToVoxels(&OriginalImage);
      Islands.MillimetersToVoxels(&OriginalImage);
    }
    
    //perform tensor voting
    TV_Manager.PerformTensorVoting2(&Network,&Islands,CharacteristicLength,CharacteristicAngle);
    
    //Save the transformed network
    if (KnownOriginalImage==1){
      Network.VoxelsToMillimeters(&OriginalImage);
      
      //TV_Manager.DeltaNetwork.VoxelsToMillimeters(&OriginalImage);
      //TV_Manager.DeltaNetwork.Write("DeltaNetwork.mv3d");
    }
  
      

  
    Network.Write(OutputFile);
}




///skeletonize an image
void Skeletonize(int argc, char **argv)
{ 
  ScalarField ImageIn;
  char ImFile[256];
  char OutMV3D[256];
  LDMK_Curves LocNetwork;
  
  //read parameters
  argc--; argv++;
  strcpy(ImFile,argv[1]);
  argc--; argv++;
  strcpy(OutMV3D,argv[1]);
  argc--; argv++;
  
  //read input image
  ImageIn.Read(ImFile);
  
  //skeletonize the shape
  LocNetwork.Skeletonize(&ImageIn);
  
  //save the result
  LocNetwork. Write(OutMV3D);
  
}


///generate a volume from a skeleton
void GenerateVolume(int argc, char **argv)
{ 
  ScalarField ImageIn;
  char ImFile[256];
  char ImFileOut[256];
  char InMV3D[256];
  LDMK_Curves LocNetwork;
  
  //read parameters
  argc--; argv++;
  strcpy(InMV3D,argv[1]);
  argc--; argv++;
  strcpy(ImFile,argv[1]);
  argc--; argv++;
  strcpy(ImFileOut,argv[1]);
  argc--; argv++;
  
  //do the job
  ImageIn.Read(ImFile);
  LocNetwork.Read(InMV3D);
  LocNetwork.Generate3DVolume(ImFileOut,&ImageIn);
  
}

///extract all non-null grey levels of an image and store them in a csv file
void ExtractNonNullGL(int argc, char **argv)
{
    ScalarField ImageIn;
    LDMK_Points tmpPts;
    char ImFile[256];
    char OutMV3D[256];
    
    //read parameters
    argc--; argv++;
    strcpy(ImFile,argv[1]);
    argc--; argv++;
    strcpy(OutMV3D,argv[1]);
    argc--; argv++;
    
    //read input image
    ImageIn.Read(ImFile);
    
    //do the job
    tmpPts.ReadInScalarField(&ImageIn);
    
    //save the network
    tmpPts.Write(OutMV3D);
    
}


///Cleanup a skeleton
void CleanupSkeleton(int argc, char **argv)
{ 
  char InMV3D[256];
  char OutMV3D[256];
  LDMK_Curves SkelToImprove;
  char OriginalImageName[256];
  int KnownOriginalImage;
  ScalarField OriginalImage;

  //read parameters
  argc--; argv++;
  strcpy(InMV3D,argv[1]);
  argc--; argv++;
  strcpy(OutMV3D,argv[1]);
  argc--; argv++;
  
  KnownOriginalImage=0;
  if (argc>1) {
    KnownOriginalImage=1;
    strcpy(OriginalImageName,argv[1]);
    argc--; argv++;
  }
  
  //import the mv3d files in 
  SkelToImprove.Read(InMV3D);
  
  if (KnownOriginalImage==1){
    OriginalImage.Read(OriginalImageName);
    SkelToImprove.MillimetersToVoxels(&OriginalImage);
    }
  
  //cleanup curvesToImprove
  SkelToImprove.CleanUp(); 

  //save the result
  if (KnownOriginalImage==1)
    SkelToImprove.VoxelsToMillimeters(&OriginalImage);
    
  SkelToImprove.Write(OutMV3D,0);
}







void RawCurve2Mv3d(int argc, char **argv)
{ 
  char InRawFile[256];
  char OutMV3D[256];
  LDMK_Curves SkelLoc;
  int SubsamplingFactor;
  
  //read parameters
  argc--; argv++;
  strcpy(InRawFile,argv[1]);
  argc--; argv++;
  SubsamplingFactor=atoi(argv[1]);
  argc--; argv++;
  strcpy(OutMV3D,argv[1]);
  argc--; argv++;
  
  //read the mv3d files
  SkelLoc.ReadInRawFile(InRawFile,SubsamplingFactor);
  
  //save the result
  SkelLoc.Write(OutMV3D);
}


void mv3d2vtk(int argc, char **argv)
{ 
  char InMV3D[256];
  char OutVTK[256];
  LDMK_Curves SkelLoc;
  
  
  //read parameters
  argc--; argv++;
  strcpy(InMV3D,argv[1]);
  argc--; argv++;
  strcpy(OutVTK,argv[1]);
  argc--; argv++;
  
  //read the mv3d files
  SkelLoc.Read(InMV3D);
  
  //save the result
  SkelLoc.ExportAsVtk(OutVTK,0);
}



///transform the [Mv3dFileIn] coordinates and diameters from voxels to mm according to the image 2 world properties of [RefScalField]
void AffTransfoPointSet(int argc, char **argv)
{ 
  char InMV3D[256];
  char OutMV3D[256];
  char TransfoMatrix[256];
  LDMK_Curves SkelLoc;
  float TransfoMat[4][4];
  
  
  //read parameters
  argc--; argv++;
  strcpy(InMV3D,argv[1]);
  argc--; argv++;
  strcpy(TransfoMatrix,argv[1]);
  argc--; argv++;
  strcpy(OutMV3D,argv[1]);
  argc--; argv++;
  
  //import the mv3d files in 
  SkelLoc.Read(InMV3D);
  
  //read the 4*4 matrix
  Read_quat4t4mat(TransfoMatrix,TransfoMat);
  
  //do the job
  SkelLoc.AffineTransfo(TransfoMat);
  
  //save the result
  SkelLoc.Write(OutMV3D);
}




///Mask a skeleton
void MaskSkeleton(int argc, char **argv)
{ 
  char InMV3D[256];
  char ImFile[256];
  char OutMV3D[256];
  ScalarField Mask;
  LDMK_Curves SkelToImprove;
  float x_sk,y_sk,z_sk;
  int SegRemove;
  int i,j;
  
  //read parameters
  argc--; argv++;
  strcpy(InMV3D,argv[1]);
  argc--; argv++;
  strcpy(ImFile,argv[1]);
  argc--; argv++;
  strcpy(OutMV3D,argv[1]);
  argc--; argv++;
  
  //import the mv3d files in 
  SkelToImprove.Read(InMV3D);
  
  //read the mask
  Mask.Read(ImFile);

  //mask the skeleton
  for (i=0;i<SkelToImprove.GetSegNumber();i++) if (SkelToImprove.GetElNumber(i)>0){
    SegRemove=0;
    for (j=0;j<SkelToImprove.GetElNumber(i);j++){
      x_sk=SkelToImprove.GetX(i,j)*Mask.World2Image[0][0]+SkelToImprove.GetY(i,j)*Mask.World2Image[0][1]+SkelToImprove.GetZ(i,j)*Mask.World2Image[0][2]+Mask.World2Image[0][3];
      y_sk=SkelToImprove.GetX(i,j)*Mask.World2Image[1][0]+SkelToImprove.GetY(i,j)*Mask.World2Image[1][1]+SkelToImprove.GetZ(i,j)*Mask.World2Image[1][2]+Mask.World2Image[1][3];
      z_sk=SkelToImprove.GetX(i,j)*Mask.World2Image[2][0]+SkelToImprove.GetY(i,j)*Mask.World2Image[2][1]+SkelToImprove.GetZ(i,j)*Mask.World2Image[2][2]+Mask.World2Image[2][3];
    
      if (fabs(Mask.G(x_sk,y_sk,z_sk))<0.1) SegRemove++;
    }
    
    if (SegRemove>1+SkelToImprove.GetElNumber(i)/5) SkelToImprove.DeleteSegment(i);
    if ((SkelToImprove.GetElNumber(i)<4)&&(SegRemove>1)) SkelToImprove.DeleteSegment(i);
  }

  //save the result
  SkelToImprove.Write(OutMV3D,0);
  
}



///Smooth a skeleton
void SmoothSkeleton(int argc, char **argv)
{ 
  char InMV3D[256];
  char OutMV3D[256];
  int NbIt;
  LDMK_Curves SkelToImprove;
  
  //read parameters
  argc--; argv++;
  strcpy(InMV3D,argv[1]);
  argc--; argv++;
  NbIt=atoi(argv[1]);
  argc--; argv++;
  strcpy(OutMV3D,argv[1]);
  argc--; argv++;
  
  //import the mv3d file
  SkelToImprove.Read(InMV3D);
  
  //do the job
  SkelToImprove.Smooth(NbIt);
  
  //save the result
  SkelToImprove.Write(OutMV3D,0);
  
}





/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                        INPUTS MANAGEMENT
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void usage(){
  cerr << "Usage: CurvesAna <options>\n";
  cerr << "\n";
  cerr << "Where <options> are the following:\n";
  cerr << "\n";
  cerr << "\n";
  cerr << "-RawCurve2Mv3d [RawFileIn][SubsamplingFactor][Mv3dFileOut]\n";
  cerr << "    -> export a raw file containing a 3D curve into a mv3d file.\n";
  cerr << "    -> One point every [SubsamplingFactor] points is loaded... all points are loaded if equal to 1.\n";
  cerr << "\n";
  cerr << "-mv3d2vtk [Mv3dFileIn] [VtkFileOut]\n";
  cerr << "    -> export a mv3d file into a vtk file.\n";
  cerr << "\n";
  cerr << "-AffTransfoPointSet [Mv3dFileIn] [TransfoMatrix] [Mv3dFileOut]\n";
  cerr << "    -> transform the [Mv3dFileIn] using the 4*4 transformation matrix [TransfoMatrix]. Result is in [Mv3dFileOut].\n";
  cerr << "\n";
  cerr << "-MaskSkeleton [Mv3dFileIn] [ImageMask] [Mv3dFileOut]\n";
  cerr << "    -> Mask [Mv3dFileIn] (in world coord.) with [ImageMask]. Result is put in [Mv3dFileOut].\n";
  cerr << "\n";
  cerr << "-SmoothSkeleton [Mv3dFileIn] [NbIt] [Mv3dFileOut]\n";
  cerr << "    -> Smooth a skeleton by averaging the neighbors [NbIt] times.\n";
  cerr << "\n";
  cerr << "-Skeletonize [SegImag] [Mv3dFile]\n";
  cerr << "    -> Skeletonize [SegImag] into the mv3d file [Mv3dFile]. Coordinates are given in the world domain\n";
  cerr << "       and not the image domain (expect if the image 2 world matrix is identity)\n";
  cerr << "\n";
  cerr << "-GenerateVolume [Mv3dFile] [RefImag] [GeneratedVolume]\n";
  cerr << "    -> Generate the 3D volume of [Mv3dFile] in the image domain of [RefImag]. Image to world properties of \n";
  cerr << "       [RefImag] are taken into account. Volume is saved in the nifti file [GeneratedVolume].\n";
  cerr << "\n";
  cerr << "-ExtractNonNullGL [RefImag] [CSVFile]\n";
  cerr << "    -> Detect all non null values of [RefImag] and store them in a csv file [CSVFile].\n";
  cerr << "\n";
  cerr << "-CleanupSkeleton [Mv3dFileIn] [Mv3dFileOut] <[OriginalImage]>\n";
  cerr << "    -> Cleanup [Mv3dFileIn] and put the result in [Mv3dFileOut].\n";
  cerr << "    -> If the coordinates of the mv3d are in mm and not in voxels (eg. using -Skeletonize), please\n";
  cerr << "       give the name of the original image so that voxel resolution is known.\n";
  cerr << "\n";
  cerr << "-TensorVoting [InputLineset][OutputLineset] [RefDist][RefAngl] <-GL [NonNullGL]><-Orig [OriginalImage]>\n";
  cerr << "    -> Close the gaps in the mv3d lineset [InputLineset] using tensor voting and save the results in\n";
  cerr << "       the mv3d file [OutputLineset].\n";
  cerr << "    -> [RefDist] and [RefAngl] are the caracteristic distance (in voxels) and angle (in degrees).\n";
  cerr << "    -> If the non-null grey levels are known, please give the csv file [NonNullGL]. Data of this file\n";
  cerr << "       are separated by tabs. Each row of the file gives the x, y, z coordinates of the point in the same\n";
  cerr << "       coordinate system as [InputLineset], plus the corresponding intensity.\n";
  cerr << "    -> If the coordinates of the mv3d are in mm and not in voxels, which is the case when using -Skeletonize,\n";
  cerr << "       please give the name of the original image [OriginalImage] so that voxel resolution is known.\n";
  cerr << "\n";
  cerr << "\n";
  
  exit(1);
}

int main(int argc, char **argv){
  bool done;
  done=false;
  

  // Check command line
  if (argc < 2) {
    usage();
  }
  else {
     if (done == false) if (strcmp(argv[1], "-RawCurve2Mv3d") == 0) {
      RawCurve2Mv3d(argc,argv);
      done = true;
    }
    if (done == false) if (strcmp(argv[1], "-mv3d2vtk") == 0) {
      mv3d2vtk(argc,argv);
      done = true;
    }
     if (done == false) if (strcmp(argv[1], "-AffTransfoPointSet") == 0) {
      AffTransfoPointSet(argc,argv);
      done = true;
    }
     if (done == false) if (strcmp(argv[1], "-MaskSkeleton") == 0) {
      MaskSkeleton(argc,argv);
      done = true;
    }
     if (done == false) if (strcmp(argv[1], "-SmoothSkeleton") == 0) {
      SmoothSkeleton(argc,argv);
      done = true;
    }
    if (done == false) if (strcmp(argv[1], "-Skeletonize") == 0) {
          Skeletonize(argc,argv);
          done = true;
    }
    if (done == false) if (strcmp(argv[1], "-GenerateVolume") == 0) {
          GenerateVolume(argc,argv);
          done = true;
    }
    if (done == false) if (strcmp(argv[1], "-ExtractNonNullGL") == 0) {
          ExtractNonNullGL(argc,argv);
          done = true;
    }
    if (done == false) if (strcmp(argv[1], "-CleanupSkeleton") == 0) {
          CleanupSkeleton(argc,argv);
          done = true;
    }
    if (done == false) if (strcmp(argv[1], "-TensorVoting") == 0) {
          TensorVoting2Launcher(argc,argv);
          done = true;
    }
  }
  
  return 0;
}
