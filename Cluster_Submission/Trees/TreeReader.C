#include <time.h> // Needed for "clock" / tree-reading time

clock_t start = clock();

void TreeReader()
{

    // Here create ROOT file to WRITE to, and all histograms/tree
    TFile* f = new TFile("TestS450.root", "RECREATE");

    // Here define the histograms you're going to want to see. ("Name", "Title", Bins, Start, End)
    TH2F* ID_Z1_AoQ = new TH2F("ID_Z1_AoQ", "Z1 vs A/Q", 500, 1.8, 2.8, 500, 40, 90);
    TH2F* ID_dBrho_dE4 = new TH2F("ID_dBrho_dE4", "dE vs dBRho", 500, -1000, 1000, 500, 0, 4000);
    TH2F* ID_dE_x2 = new TH2F("ID_dE_x2", "dE vs X2", 500, 0, 4000, 500, -200, 100);
    TH2F* ID_dE_x4 = new TH2F("ID_dE_x4", "dE vs X4", 500, 0, 4000, 500, -200, 100);
    TH2F* ID_x2_x4 = new TH2F("ID_x2_x4", "X2 vs X4", 500, -200, 100, 500, -200, 100);

    /* 
    Here define filepath to ROOT tree that we want to READ. 
    I'm looping over multiple trees in the following way -- PROBABLY NOT OPTIMAL / "correct"(?)
    Will send you a better way to do it if I find one. 
    */
    TString path = "~/../../lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Trees/Run_75/";
    vector<TString> file_array;
    TString file;

    /*
    Define `i` based on how many files you want to read, and how they're named.
    In this example I have 3 files called S452_Tree_10_1, S452_Tree_10_2, S452_Tree_10_3
    But yours may be different.
    */
    for (int i = 0; i < 127; i++)
    {   

        TString file = Form("tree75_%d.root",i*2); 
        file_array.push_back(file);
    }

    // Declare any variables you might need later in the reader here.

    // Loop through the dimension of the `file_array`. I had 3 files so I go from j = 0, to j = 2
    for (int j = 0; j < 127; j++)
    {
        TString totpath = path + file_array[j];

        // Open the file (or try) and quit if it fails (generally because it can't find due to a mistake here)
        TFile* MyFile = TFile::Open(totpath);
        if (MyFile == 0)
        {
            cout << "Error: The Tree is Dead; Plant Another... " << endl;
            gROOT->ProcessLine(".q");
        }

        // Just for checking which file is currently being processed
        cout << "Processing file: " << file_array[j] << endl;

        /*
        Here the tree-reader starts. All trees I've seen are called `AnalysisxTree`,
        So this probably doesn't need changing. Can be checked by opening the tree in ROOT though
        */
        TTreeReader MyReader("AnalysisxTree", MyFile);

        /*
        Here you'll point to BRANCHES of the Tree, so you can access them to plot in histograms.
        You can (should) open the tree in ROOT to check how each branch is named, below are common examples
        
        If you need an array use TTreeReaderArray, if you need a value use TTreeReaderValue etc.
        This bit is something I always fuck up just because its annoying to remember what to use and when, 
        You'll get the hang of it playing around.
        */
        TTreeReaderValue<float> FRS_Z(MyReader, "AnlEvent.pFRS_z");
        TTreeReaderValue<float> FRS_AoQ(MyReader, "AnlEvent.pFRS_AoQ");
        TTreeReaderValue<float> FRS_X2(MyReader, "AnlEvent.pFRS_ID_x2");
        TTreeReaderValue<float> FRS_X4(MyReader, "AnlEvent.pFRS_ID_x4");
        TTreeReaderValue<float> FRS_dE(MyReader, "AnlEvent.pFRS_Music_dE");



        /*
        Here I'm defining some cuts. You'll need the X and Y variables to exist before doing this,
        I.e. You'll need to hav pointed to the correct branches and named them things like FRS_Z
        */
        TCutG *Poly1 = new TCutG("Poly1",10);
        Poly1->SetVarX("FRS_AoQ");
        Poly1->SetVarY("FRS_Z");
        Poly1->SetPoint(0,2.66481,79.0649);
        Poly1->SetPoint(1,2.66481,79.0649);
        Poly1->SetPoint(2,2.64949,78.8895);
        Poly1->SetPoint(3,2.64182,78.6736);
        Poly1->SetPoint(4,2.64582,78.4441);
        Poly1->SetPoint(5,2.66714,78.2957);
        Poly1->SetPoint(6,2.6838,78.4037);
        Poly1->SetPoint(7,2.68913,78.7073);
        Poly1->SetPoint(8,2.67414,78.9705);
        Poly1->SetPoint(9,2.66481,79.0649);

        TCutG *Poly2 = new TCutG("Poly2",10);
        Poly2->SetVarX("FRS_AoQ");
        Poly2->SetVarY("FRS_Z");
        Poly2->SetPoint(0,2.65915,78.1472);
        Poly2->SetPoint(1,2.65915,78.1472);
        Poly2->SetPoint(2,2.64515,78.0325);
        Poly2->SetPoint(3,2.63949,77.7896);
        Poly2->SetPoint(4,2.64715,77.5399);
        Poly2->SetPoint(5,2.66048,77.486);
        Poly2->SetPoint(6,2.67847,77.5939);
        Poly2->SetPoint(7,2.6828,77.8166);
        Poly2->SetPoint(8,2.67314,78.0393);
        Poly2->SetPoint(9,2.65915,78.1472);

        // Start of reading process
        while (MyReader.Next())
        {    
            
            /*
            This bits kinda tricky and I don't fully understand. 
            The gist, I think, is that FRS_Z is POINTING to the value of z,
            and so you need to create a `float z` variable to then assign it this value using the asterisk operator.
            All I know for sure is, this is the only way it ever works. You can't do ->Fill(FRS_AoQ, FRS_Z)
            */
            float z = *FRS_Z;
            float aoq = *FRS_AoQ;
            float x2 = *FRS_X2;
            float x4 = *FRS_X4;
            float dE = *FRS_dE;

            // Filling
            ID_Z1_AoQ->Fill(aoq, z);
            ID_dE_x2->Fill(x2,dE);
            ID_dE_x4->Fill(x4,dE);
            ID_x2_x4->Fill(x2,x4);


            float brho1 = (1 - (x2 - 1.829)/12.397) *13.3486;
            float brho2 = (1 - (x4 - 1.829*x2)/12.397) * 9.7927;
            float dbrho = brho2-brho1;

            ID_dBrho_dE4->Fill(dbrho, dE);

        } // reader loop
    } // file loop

    f->Write();
    f->Close();

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    cout << "Time taken to process: " << time_spent << " seconds." << endl;

    gROOT->ProcessLine(".q");

} // end of ImplantDecayMPC
