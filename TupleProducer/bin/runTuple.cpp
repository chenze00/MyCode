using namespace std;
#include "MyCode/TupleProducer/interface/MyTauClass.h"


int main(int argc, char * argv[]){
    
    TString input_filename=TString(argv[1]);
    TFile *file=new TFile(input_filename,"READ");
    TTree *tree=(TTree*)file->Get("taus");
    
    MyTauClass tau(tree,input_filename);
    
    tau.Loop();
    
    return(EXIT_SUCCESS);
    
}
