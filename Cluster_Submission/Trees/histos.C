#include <iostream>
#include <TTree.h>
#include <TTreeReader.h>

using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void histos(){
	cout<< "time to plot some histograms :)" <<endl;

	}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void aoq(){

	//TCHAIN START
	TChain* ch1 = new TChain("AnalysisxTree");
    	TString path = "";
   		TString file;
    	TString totpath;

    for (int i = 0; i < 2; i++){   
        file = Form("tree75_%d.root",i*2);
        totpath = path + file;
        ch1->Add(totpath);
        cout << "file: " << i << " added to chain." << endl;
    }
	//TCHAIN END

	//TCUTG START
    TFile *cut=new TFile("blob.root"); 
    TCutG *mycutg;
    mycutg=(TCutG*)cut->Get("CUTG"); 
    mycutg->SetName("mycutg");
	//TCUTG END

	//CUT PLOT START
	auto *c1=new TCanvas("c1","c1",800,600);
	TH2D* h =new TH2D("h","",500, 2.5, 2.8, 500, 70, 85);
	ch1->Draw("AnlEvent.pFRS_z:AnlEvent.pFRS_AoQ>>h","","colz");
	//ch1->Draw("AnlEvent.pFRS_z:AnlEvent.pFRS_AoQ>>h");
	//h->Draw("colz [mycutg]");
    //mycutg->Draw("same");
	//CUT PLOT END	
	
	//REMOVE STATS BOX
	gStyle->SetOptStat(0);
	gPad->Update();	

	}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void dEdBRho(){

	//TCHAIN START
	TChain* ch1 = new TChain("AnalysisxTree");
    	TString path = "";
   		TString file;
    	TString totpath;

    for (int i = 0; i < 2; i++){   
        file = Form("tree75_%d.root",i*2);
        totpath = path + file;
        ch1->Add(totpath);
        cout << "file: " << i << " added to chain." << endl;
    }
	//TCHAIN END

	//LIST OF DE VARIABLES START
	//AnlEvent.pFRS_Music_dE
	//AnlEvent.pFRS_dEdeg
	//AnlEvent.pFRS_dEdegoQ
	//LIST OF DE VARIABLES END


	//CUTS START 
	TCutG *middle = new TCutG("middle",13);
	middle->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	middle->SetVarY("AnlEvent.pFRS_Music_dE");
	middle->SetPoint(0,-3.67869,2340.82);
	middle->SetPoint(1,-3.67755,2417.59);
	middle->SetPoint(2,-3.65017,2430.82);
	middle->SetPoint(3,-3.43795,2191.26);
	middle->SetPoint(4,-3.42996,2110.53);
	middle->SetPoint(5,-3.42426,1996.71);
	middle->SetPoint(6,-3.42426,1946.41);
	middle->SetPoint(7,-3.44936,1950.38);
	middle->SetPoint(8,-3.45278,2066.85);
	middle->SetPoint(9,-3.45963,2089.35);
	middle->SetPoint(10,-3.67926,2302.44);
	middle->SetPoint(11,-3.67926,2346.12);
	middle->SetPoint(12,-3.67869,2340.82);

	TCutG *top = new TCutG("top",13);
	top->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	top->SetVarY("AnlEvent.pFRS_Music_dE");
    top->SetPoint(0,-3.42825,2152.88);
    top->SetPoint(1,-3.50127,2268.03);
    top->SetPoint(2,-3.54292,2323.62);
    top->SetPoint(3,-3.59369,2358.03);
    top->SetPoint(4,-3.61651,2367.29);
    top->SetPoint(5,-3.63248,2410.97);
    top->SetPoint(6,-3.67641,2412.29);
    top->SetPoint(7,-3.67755,2511.56);
    top->SetPoint(8,-3.54292,2446.71);
    top->SetPoint(9,-3.45563,2348.76);
    top->SetPoint(10,-3.42483,2290.53);
    top->SetPoint(11,-3.42996,2151.56);
    top->SetPoint(12,-3.42825,2152.88);

	TCutG *bottom = new TCutG("bottom",11);
	bottom->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	bottom->SetVarY("AnlEvent.pFRS_Music_dE");
    bottom->SetPoint(0,-3.45849,2008.62);
    bottom->SetPoint(1,-3.5161,2113.18);
    bottom->SetPoint(2,-3.60852,2213.76);
    bottom->SetPoint(3,-3.66215,2273.32);
    bottom->SetPoint(4,-3.68611,2287.88);
    bottom->SetPoint(5,-3.67926,2163.47);
    bottom->SetPoint(6,-3.58856,2072.15);
    bottom->SetPoint(7,-3.5298,2005.97);
    bottom->SetPoint(8,-3.46248,1926.56);
    bottom->SetPoint(9,-3.45849,2015.24);
    bottom->SetPoint(10,-3.45849,2008.62);
   
	TCutG *blur = new TCutG("blur",10);
	blur->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	blur->SetVarY("AnlEvent.pFRS_Music_dE");
    blur->SetPoint(0,-3.9993,2481.12);
    blur->SetPoint(1,-3.99417,2013.91);
    blur->SetPoint(2,-3.66785,1994.06);
    blur->SetPoint(3,-3.67869,2432.15);
    blur->SetPoint(4,-3.86581,2457.29);
    blur->SetPoint(5,-3.91943,2477.15);
    blur->SetPoint(6,-3.95537,2478.47);
    blur->SetPoint(7,-3.9839,2481.12);
    blur->SetPoint(8,-3.99987,2479.79);
    blur->SetPoint(9,-3.9993,2481.12);
	//CUTS END

	//NO CUT PLOT START
	auto *c2=new TCanvas("c2","dEvsdBRho",800,600);
	TH2D* dEvsdBRho=new TH2D("dEvsdBRho","",500, -6, 0, 1000, -1000, 1500);
	//ch1->Draw("AnlEvent.pFRS_Music_dE:((1 - (AnlEvent.pFRS_ID_x4/7981)) * 9.7927)>>h");
	ch1->Draw("AnlEvent.pFRS_Music_dE:((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)>>dEvsdBRho");
	blur->Draw("same");	
	dEvsdBRho->Draw("colz");
	//top->Draw("same");
	//middle->Draw("same");
	//bottom->Draw("same");
	//NO CUT PLOT END


	
	//CUT PLOTS START
	TCanvas *c3=new TCanvas("c3","Middle Cut",800,600);
	ch1->Draw("AnlEvent.pFRS_z:AnlEvent.pFRS_AoQ>>MiddleCut(500, 2.5, 2.8, 500, 70, 85)","middle","colz");

	TCanvas *c4=new TCanvas("c4","Top Cut",800,600);
	ch1->Draw("AnlEvent.pFRS_z:AnlEvent.pFRS_AoQ>>TopCut(500, 2.5, 2.8, 500, 70, 85)","top","colz");

	TCanvas *c5=new TCanvas("c5","Bottom Cut",800,600);
	ch1->Draw("AnlEvent.pFRS_z:AnlEvent.pFRS_AoQ>>BottomCut(500, 2.5, 2.8, 500, 70, 85)","bottom","colz");

	TCanvas *c6=new TCanvas("c6","Blur Cut",800,600);
	ch1->Draw("AnlEvent.pFRS_z:AnlEvent.pFRS_AoQ>>BlurCut(500, 2.5, 2.8, 500, 70, 85)","blur","colz");
	//CUTS PLOT END
	

}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void BRho(){

	//TCHAIN START
	TChain* ch1 = new TChain("AnalysisxTree");
    	TString path = "Run_75/";
   		TString file;
    	TString totpath;

    for (int i = 0; i < 128; i++){   
        file = Form("tree75_%d.root",i*2);
        totpath = path + file;
        ch1->Add(totpath);
        cout << "file: " << i << " added to chain." << endl;
    }
	//TCHAIN END

    //FRS_brho1 = ((1 - (FRS_ID_x2/7981)) *13.3486);
    //FRS_brho2 = ((1 - (FRS_ID_x4 - 1.268*FRS_ID_x2)/7981) * 9.7927);

	//PLOTTING START
	auto *c7=new TCanvas("c7","BRho1vsBRho2",800,600);
	TH2D *BRho1vsBRho2=new TH2D("BRho1vsBRho2","",5000, 8, 12, 5000, 11, 14);
	ch1->Draw("((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486):((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)>>BRho1vsBRho2", "", "colz");
	BRho1vsBRho2->GetYaxis()->SetTitle("BRho1 (TA-S2)");
	BRho1vsBRho2->GetXaxis()->SetTitle("BRho2 (S2-S4)");

	//PLOTTING END

}
	
	
void stats_dEdBRho(){

	//TCHAIN START
	TChain* ch1 = new TChain("AnalysisxTree");
    	TString path = "";
   		TString file;
    	TString totpath;

    for (int i = 0; i < 2; i++){   
        file = Form("tree75_%d.root",i*2);
        totpath = path + file;
        ch1->Add(totpath);
        cout << "file: " << i << " added to chain." << endl;
    }
	//TCHAIN END


	//CUTS START
	//CORRECT DATA (ALL)
	TCutG *correct = new TCutG("correct",6);
	correct->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	correct->SetVarY("AnlEvent.pFRS_Music_dE");
	correct->SetPoint(0,-4.3471,3100.24);
	correct->SetPoint(1,-4.23085,454.923);
	correct->SetPoint(2,-2.77434,425.267);
	correct->SetPoint(3,-3.02735,3052.79);
	correct->SetPoint(4,-4.34026,3106.17);
	correct->SetPoint(5,-4.3471,3100.24);

	//INCORRECT DATA (ALL)
	TCutG *wrong = new TCutG("wrong",6);
	wrong->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	wrong->SetVarY("AnlEvent.pFRS_Music_dE");
	wrong->SetPoint(0,-1.19475,3106.17);
	wrong->SetPoint(1,-1.0785,460.854);
	wrong->SetPoint(2,0.378009,431.198);
	wrong->SetPoint(3,0.125,3058.72);
	wrong->SetPoint(4,-1.18791,3112.1);
	wrong->SetPoint(5,-1.19475,3106.17);

	//CORRECT DATA MINUS FISSION FRAGS
	TCutG *goodblob = new TCutG("goodblob",6);
	goodblob->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	goodblob->SetVarY("AnlEvent.pFRS_Music_dE");
	goodblob->SetPoint(0,-4.21718,2655.4);
	goodblob->SetPoint(1,-4.19666,1741.99);
	goodblob->SetPoint(2,-2.93846,1724.2);
	goodblob->SetPoint(3,-2.98632,2655.4);
	goodblob->SetPoint(4,-3.94365,2696.92);
	goodblob->SetPoint(5,-4.21718,2655.4);

	//INCORRECT DATA MINUS FISSION FRAGS
	TCutG* badblob = new TCutG("badblob",6);
	badblob->SetVarX("((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)");
	badblob->SetVarY("AnlEvent.pFRS_Music_dE");
	badblob->SetPoint(0,-1.05115,2679.12);
	badblob->SetPoint(1,-1.03063,1765.72);
	badblob->SetPoint(2,0.227571,1747.92);
	badblob->SetPoint(3,0.179705,2679.12);
	badblob->SetPoint(4,-0.777626,2720.64);
	badblob->SetPoint(5,-1.05115,2679.12);


	auto *c8=new TCanvas("c7","dEvsdBRho",800,600);
	//TH2D* dEvsdBRho=new TH2D("dEvsdBRho","",500, -6, 0, 1000, -1000, 1500);
	//ch1->Draw("AnlEvent.pFRS_Music_dE:((1 - (AnlEvent.pFRS_ID_x4/7981)) * 9.7927)>>h");
	ch1->Draw("AnlEvent.pFRS_Music_dE:((1 - (AnlEvent.pFRS_ID_x4 - 1.268*AnlEvent.pFRS_ID_x2)/7981) * 9.7927)-((1 - (AnlEvent.pFRS_ID_x2/7981)) *13.3486)>>dEvsdBRho(500, -6, 0, 1000, 0, 3000)", "", "colz");
	//dEvsdBRho->Draw("[goodblob] colz");

}	


void dE_MUSIC_deg(){

		//TCHAIN START
	TChain* ch1 = new TChain("AnalysisxTree");
    	TString path = "Run_75/";
   		TString file;
    	TString totpath;

    for (int i = 0; i < 2; i++){   
        file = Form("tree75_%d.root",i*2);
        totpath = path + file;
        ch1->Add(totpath);
        cout << "file: " << i << " added to chain." << endl;
    }
	//TCHAIN END


	auto *c7=new TCanvas("c7","dE_MUSIC_deg",800,600);
	TH2D *dE_MUSIC_deg=new TH2D("dE_MUSIC_deg","",5000, 0, 100, 5000, 0, 3000);
	ch1->Draw("AnlEvent.pFRS_Music_dE:AnlEvent.pFRS_dEdeg>>dE_MUSIC_deg", "", "colz");
	dE_MUSIC_deg->GetYaxis()->SetTitle("FRS_Music_dE");
	dE_MUSIC_deg->GetXaxis()->SetTitle("FRS_dEdeg");

}
