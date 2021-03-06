#include "Analysis.h"
#include "makeHists.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;
int main(int argc, char** argv)
{
    //usage output
    if (argc < 4)
    {
        std::cout << "usage : " << argv[0] << " ntuple_file_name ntuple_version output_tag EXTRAOPT=(either SYST or SKIM)" << std::endl;
        return 0;
    }

    //
    // read in filelist
    // filelist will contain a list of file names to run over without .root at the end
    //
    // e.g.  If we have,
    //
    // ${WVZ_DATA_PATH}/babies/WVZ2018_v0.0.5/
    //      dy_m1050_madgraph_1.root
    //      dy_m50_madgraph_1.root
    //      ggh_hzz4l_powheg_1.root
    //      sts_4f_leptonic_madgraph_1.root
    //      ttbar_amcatnlo_1.root
    //      ...
    //
    // The input file list should have
    //
    // $ cat  configs/Filelist_bkg.list
    // dy_m1050_madgraph_1
    // dy_m50_madgraph_1
    // ggh_hzz4l_powheg_1
    // sts_4f_leptonic_madgraph_1
    // ttbar_amcatnlo_1
    //
    // N.B. without .root
    //

    TString InputOptionString = argv[1];

    // Input file name list
    std::vector<TString> filenamelist;

    // If the input option string for the first argument
    // contains ".root" it means it's not a text file of
    // list of file names but rather the literal file
    // name to be run over in the input directory
    if (InputOptionString.Contains(".root"))
    {
        TString filename = InputOptionString;
        filename.ReplaceAll(".root","");
        filenamelist.push_back(filename);
    }
    else
    {
        ifstream infile;
        infile.open(argv[1], ios::in);
        // parsing the input file and obtaining file name list
        TString InputRoot;
        while (infile >> InputRoot)
        {

            // Form full path to the input root file
            TString RootName = InputRoot;
            filenamelist.push_back(RootName);

        }
    }

    int count = 0; // To keep track of how many input files were ran over

    cout << "*************************" << endl;
    cout << "*** Analysis Begin ****" << endl;
    cout << "*************************" << endl;
    cout << endl;

    // if argv[1] is of the format "file.root" with a ".root" extension
    // take the argument as the file name to be run over

    const char * userDataPath = std::getenv("WVZ_DATA_PATH");
    if(!userDataPath) {
        throw std::runtime_error("WVZ_DATA_PATH not set!");
    }

    for (auto& InputRoot : filenamelist)
    {

        // Form full path to the input root file
        TString RootName = InputRoot;
        TString RootAdd = (TString)(std::string(userDataPath) + "/babies/" + argv[2] + "/" + InputRoot + ".root");

        // Increase # of input files ran over
        count++;

        // Create the analysis instance
        Analysis Run(argv[1], InputRoot);

        cout << " Initial Begin: " << endl;

        // Set the input root file and the "index" of the input root file
        Run.Initial(RootAdd, count);

        cout << " Loop Begin:" << endl;

        //
        // Loop!
        //
        bool dosyst = false;
        bool doskim = false;
        TString filtercuts = "";
        if (argc>=5)
        {
            std::string extraopt = argv[4];
            if (TString(extraopt).EqualTo("SYST"))
                dosyst = true;
            else if (TString(extraopt).EqualTo("SKIM"))
                doskim = true;
            else
                filtercuts = argv[4];
        }
        if (argc>=6)
        {
            filtercuts = argv[5];
        }

        Run.Loop(argv[2], argv[3], dosyst, doskim, filtercuts);

        // Done
        Run.End(count);

        Run.Output();

        Run.Finish(count);
    }

    cout << endl;
    cout << "**********************" << endl;
    cout << "**** Analysis End ****" << endl;
    cout << "**********************" << endl;


    return 1;
}
