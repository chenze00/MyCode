//argv[1] is full path  root file name
//argv[2] is output path
using namespace std;
#include "MyCode/TupleProducer/interface/MyTauClass.h"

int main(int argc, char * argv[]){
    
    TString input_filename=TString(argv[1]);
    TString outputpath=TString(argv[2]);
    TFile *file=new TFile(input_filename,"READ");
    TTree *tree=(TTree*)file->Get("taus");
    
    MyTauClass tau(tree,input_filename);
    
    tau.Loop(outputpath);
    
    return(EXIT_SUCCESS);
    
}
