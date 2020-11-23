//mex transpose image stack
#include "Windows.h"
#include "mex.h"
#include <iostream>
#include <time.h> 


using namespace std;

typedef unsigned short uint16;

struct _Sinf
{
    char* SrcFn;
    char* SaveFn;
    char* StrPath;
    double* Dimensions;
    double Freq;

} Sinf;

void gettime(char* Datestr, const uint16 Bufsz){
    time_t rawtime;
    struct tm * timeinfo;
    
    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (Datestr, Bufsz, "%c", timeinfo);

}

void mexFunction(
				 int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
    char *Filename;
    mxArray *tmp;
    FILE* infile;
    FILE* Outfile;
  
    
    if (nrhs != 1) {
        mexErrMsgTxt("Input required!....");
        return;
    } 
           
    tmp = mxGetField(prhs[0], 0, "Source" );
    if( tmp == NULL )
        mexErrMsgTxt("Source not found");
    else {
        Filename = mxArrayToString(tmp);
        if(strlen(Filename) == 0)
            Sinf.SrcFn = NULL;
        else {
           Sinf.SrcFn =  Filename;
           infile = fopen(Filename, "rb");
           if (infile == NULL){
                mexErrMsgTxt("Error opening input file"); 
                return;
           }
        }
    } 
    
        tmp = mxGetField(prhs[0], 0, "Save" );
    if( tmp == NULL )
        mexErrMsgTxt("Save not found");
    else {
        Filename = mxArrayToString(tmp);
        if(strlen(Filename) == 0)
            Sinf.SaveFn = NULL;
        else
           Sinf.SaveFn =  Filename;
           Outfile = fopen(Filename, "wb");
        if (Outfile==NULL){
            mexErrMsgTxt("Error opening output file"); 
            return;
        }
        fclose(Outfile);
    } 
        
    tmp = mxGetField(prhs[0], 0, "StrPath" );
    if( tmp == NULL )
        mexErrMsgTxt("Path string not found");
    else {
        Filename = mxArrayToString(tmp);
        if(strlen(Filename) != 0){
            Sinf.StrPath = Filename;
        }else{
           mexErrMsgTxt("Path string empty");
           return;
        }
    }
        
    tmp = mxGetField(prhs[0],0, "Freq" );
    if( tmp == NULL ){
        mexErrMsgTxt("Freq not found");
        Sinf.Freq = 0.;
    }
    else {   
        Sinf.Freq = mxGetScalar(tmp);
    }
        
    tmp = mxGetField(prhs[0],0, "Dimensions" );
    if( tmp == NULL )
        mexErrMsgTxt("DisplayDim not found");
    else {   
        Sinf.Dimensions = mxGetPr(tmp);
         mexPrintf("LineSz : %d, Lines : %d, Slices : %d, StackSz : %d, Frequency : %e   \n", 
        (UINT)Sinf.Dimensions[0], (UINT)Sinf.Dimensions[1], (UINT)Sinf.Dimensions[2], (UINT)Sinf.Dimensions[3], Sinf.Freq);
         if(Sinf.Dimensions[2] > 1)
             mexPrintf("Stack contains both green and red channel. Red channel will be removed....\n");
    }
    mexEvalString("drawnow");

    
    UINT lines = (UINT)Sinf.Dimensions[1];
    FILE **  ppfile = new FILE * [lines];
    char* strp = Sinf.StrPath;
    
    //stdio has a limit on open files < 512
    if(lines > 500) _setmaxstdio(lines+12);
    
    char** Fn = new char * [lines];    
    for(UINT i = 0; i < lines; i++){
        Fn[i] = new char[200];
        sprintf(Fn[i], "%s\\line_%d.tmp", strp, i);
        ppfile[i] = fopen(Fn[i], "wb");
        char ErrStr[200];
        if (ppfile[i] == NULL){
            sprintf(ErrStr, "Error opening: %s\\line_%d.tmp", strp, i);
            mexErrMsgTxt(ErrStr); 
        }
    }
    
    UINT StackSz = (UINT)Sinf.Dimensions[3];
    UINT LineSz = (UINT)Sinf.Dimensions[0];
    UINT ChanSz = (UINT)Sinf.Dimensions[2];
    size_t Bsz = LineSz * lines * 2 * ChanSz;
    char* BUFO = new char[Bsz];
    
    printf("\nProcessing Frames:       ");
    for(UINT i = 0; i < StackSz; i++){
        UINT cnt = (UINT)fread(BUFO, sizeof(char), Bsz, infile);
        if( cnt == Bsz){
            for(UINT j = 0; j < lines; j++){
                fwrite( (void*)(BUFO + j*LineSz*2), sizeof(char), LineSz*2, ppfile[j]);    
            }
            if(i%100 == 0){
                printf("\b\b\b\b\b%5i", i);
                mexEvalString("drawnow");
            }
             
        } else {
            char ErrStr[200];
            sprintf(ErrStr, "Not enough data for Frame: %d \n", i);
            mexErrMsgTxt(ErrStr);
        }
    }
    for(UINT i = 0; i < lines; i++){
        fclose(ppfile[i]);       
    } 

    delete [] ppfile;
    delete [] BUFO;
    fclose(infile);
    
    //finished writing data to line files
    char Header[500];
    const uint16 Bufsz = 100;
    char strDatetime[Bufsz];
    gettime(strDatetime, Bufsz);
    sprintf(Header, "%d %d %d %e %s", StackSz, LineSz, lines, Sinf.Freq, strDatetime);
    
    Outfile = fopen(Sinf.SaveFn, "wb");
    fwrite(Header, 1, 500, Outfile);
    printf("\nAll Frames processed. \nTransposing lines: ");
    
    Bsz = LineSz * StackSz * 2;
    char* BUFI = new char[Bsz];
    BUFO = new char[StackSz * 2];
    
    for(UINT i = 0; i < lines; i++){
        infile = fopen(Fn[i], "rb+");  
        UINT cnt = (UINT)fread(BUFI, 1, Bsz, infile);
        if(cnt == Bsz){
            for(UINT j = 0; j < LineSz; j++){
                for(UINT k = 0; k < StackSz; k++){ 
                    BUFO[k*2] = BUFI[k*LineSz*2 + j*2];
                    BUFO[k*2+1] = BUFI[k*LineSz*2+1 + j*2];
                }
                fwrite(BUFO, 1, StackSz*2, Outfile);
            }
        } else {
                char ErrStr[200];
                sprintf(ErrStr, "Not enough data for Line: %d \n", i);
                mexErrMsgTxt(ErrStr); 
        }
            
        fclose(infile);       
        remove(Fn[i]);
        delete [] Fn[i];
        if(i%10 == 0){
            printf(".");
            mexEvalString("drawnow");
        }
    }
    printf("\nDone!\n");
    
    delete [] Fn;
    delete [] BUFI;
    delete [] BUFO;
    
    fclose(Outfile);
    
}