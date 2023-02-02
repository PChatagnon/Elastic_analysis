#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "bib/Particleclass.h"
//#include "bib/TCSfunc.h"
#include "bib/ElasticPlotClass.h"
//#include "bib/TCSMomentumCorrection.h"
#include "bib/LeptonIDClass.h"
#include "bib/FiducialCuts.h"
#include "bib/ElasticEvent.h"
#include "bib/ElasticMCEvent.h"
#include "bib/RunSelectorClass.h"
#include "bib/InputParser.h"

#include "reader.h"

#include <ctime> // time_t
#include <cstdio>
using namespace std;

#define ADDVAR(x, name, t, tree) tree->Branch(name, x, TString(name) + TString(t))

int analysisElastic()
{

	time_t begin, intermediate, end; // time_t is a datatype to store time values.

	time(&begin); // note time before execution

	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(111);
	gStyle->SetPalette(55);
	gStyle->SetLabelSize(.05, "xyz");
	gStyle->SetTitleSize(.05, "xyz");
	gStyle->SetTitleSize(.07, "t");
	gStyle->SetMarkerStyle(13);
	gStyle->SetOptFit(1);


	Int_t argc = gApplication->Argc();
	char **argv = gApplication->Argv();
	Input input(argc, argv);

	bool IsData = true;
	bool IsHipo = true;

	bool IsSimu = false;
	IsSimu = input.cmdOptionExists("-IsSimu");

	bool RGA_Fall2018 = false; // inbending or outbending in the end
	RGA_Fall2018 = input.cmdOptionExists("-RGA_Fall2018");

	bool inbending = true;
	inbending = !input.cmdOptionExists("-outbending");

	if (input.cmdOptionExists("-usage"))
	{
		cout << "Use as : clas12root -l analysisElastic.C -o ouputname -f files -ef -inbending\n";
	}

	cout << "////////////////////////////////////////////"
		 << "\n";
	cout << "Run with the following options : "
		 << "\n";
	cout << "Inbending : " << inbending << "\n";
	cout << "RGA_Fall2018 : " << RGA_Fall2018 << "\n";
	cout << "IsSimu : " << IsSimu << "\n";
	cout << "////////////////////////////////////////////"
		 << "\n";

	double nbrecEvent = 0;
	int nbf = 0;
	double nEventTCS = 0;
	double denom = 0;
	double nCD = 0;
	double nFD = 0;

	int nbevent_after_EC = 0;
	int nbevent_after_posi = 0;

	int AfterCuts = 0;

	TString nameFiles = "";

	TString type = "REC";

	////////////////////////////////////////////
	// Init run selector
	////////////////////////////////////////////

	RunSelector Run_Selector;

	//////////////////////////////////////////////
	// Momentum Correction
	//////////////////////////////////////////////
	//MomentumCorrection MomCorr;

	///////////////////////////////////////////
	// Setup the TTree output
	TString output_file = (TString)(input.getCmdOption("-o")); // argv[4]);
	TFile *outFile = new TFile(Form("outputElastic_%s.root", output_file.Data()), "recreate");
	TTree *outT = new TTree("tree", "tree");
	//TTree *outT_Gen = new TTree("tree_Gen", "tree_Gen");

	TLorentzVector tree_Electron, tree_Proton;
	outT->Branch("Electron", "TLorentzVector", &tree_Electron);
	outT->Branch("Proton", "TLorentzVector", &tree_Proton);

	TLorentzVector tree_Electron_Gen, tree_Proton_Gen;
	outT->Branch("Electron_Gen", "TLorentzVector", &tree_Electron_Gen);
	outT->Branch("Proton_Gen", "TLorentzVector", &tree_Proton_Gen);

	int trigger_bit;
	outT->Branch("trigger_bit", &trigger_bit, "trigger_bit/I");

	TString fvars[] = {
		"evt_num", "Q2", "W", "weight", "run",
		"electron_Nphe", "electron_SF", "electron_score",
		"status_elec", "status_prot",
		"chi2_proton",
		"vx_elec", "vy_elec", "vz_elec",
		"vx_prot", "vy_prot", "vz_prot",
		"PCAL_x_elec", "PCAL_y_elec",
		"PCAL_sector_elec",
		"PCAL_energy_elec", "ECIN_energy_elec",
		"electron_HTCC_ECAL_match", 
		"PCAL_U_elec", "PCAL_V_elec", "PCAL_W_elec", 
		
		"evt_num",
		"Q2_Gen",
		"W_Gen",
		"vz_elec_Gen",
		"vz_prot_Gen",
		
		};

	std::map<TString, Float_t> outVars;
	for (size_t i = 0; i < sizeof(fvars) / sizeof(TString); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}

	TString fvars_Gen[] = {
		"evt_num",
		"Q2_Gen",
		"W_Gen",
		"vz_elec_Gen",
		"vz_prot_Gen",
	};

	std::map<TString, Float_t> outVars_Gen;
	if (IsSimu)
	{
		for (size_t i = 0; i < sizeof(fvars_Gen) / sizeof(TString); i++)
		{
			outVars_Gen[fvars_Gen[i]] = 0.;
			ADDVAR(&(outVars_Gen[fvars_Gen[i]]), fvars_Gen[i], "/F", outT_Gen);
		}

		outT_Gen->Branch("Electron_Gen", "TLorentzVector", &tree_Electron_Gen);
		outT_Gen->Branch("Proton_Gen", "TLorentzVector", &tree_Proton_Gen);
	}
	///////////////////////////////////////////

	int nbtc = 0;
	int nbJPSI = 0;

	///////////////////////////////////////////
	// TMVA PID for Positron
	///////////////////////////////////////////
	// PositronIdentification PositronPID("MLP method", "TMVAClassification_MLP6D.weights.xml", InputParameters.MLPscoreCut, InputParameters.MLPMomCut);
	// PositronPID.InitializePositronIdentification();

	TString electron_bdt_weights;

	if (inbending)
	{
		electron_bdt_weights = "ML_weights/TMVAClassification_BDT_neg_inbending.weights.xml";
	}
	else
	{
		electron_bdt_weights = "ML_weights/TMVAClassification_BDT_neg_outbending.weights.xml";
	}

	
	LeptonIdentification ElectronPID("BDT", electron_bdt_weights, 0.0, 4.0);
	ElectronPID.InitializeBDT_new_PositronIdentification();

	///////////////////////////////////////////
	// Plots
	///////////////////////////////////////////
	ElasticPlots Plots;
	Plots.Initialize_1D();
	Plots.Initialize_2D();
	Plots.SetOutputFolder("Plots");

	

	int corrrad = 0;
	int nbEvent = 0;

	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	for (Int_t i = input.getCmdIndex("-f") + 2; i < input.getCmdIndex("-ef") + 1; i++)
	{
		if (TString(argv[i]).Contains("MC") || IsSimu)
		{
			IsData = false;
		}

		if (TString(argv[i]).Contains(".hipo"))
		{
			nbf++;
			nameFiles = TString(argv[i]);
		}

		if (TString(argv[5]).Contains(".root"))
		{
			IsHipo = false;
			nameFiles = TString(argv[i]);
		}

		////////////////////////////////////////////
		cout << "////////////////////////////////////////////"
			 << "\n";
		if (IsData)
			cout << "Running on Data"
				 << "\n";
		else if (IsSimu)
			cout << "Running on Simulation"
				 << "\n";

		cout << TString(argv[i]) << "\n";
		cout << "Is hipo ? " << IsHipo << "\n";
		cout << "////////////////////////////////////////////"
			 << "\n";
		////////////////////////////////////////////

		////////////////////////////////////////////
		// hipo reader
		hipo::reader reader;
		hipo::dictionary factory;
		hipo::event hipo_event;
		////////////////////////////////////////////

		reader.open(nameFiles);
		reader.readDictionary(factory);
		factory.show(); // might be useful to remove this

		hipo::bank EVENT(factory.getSchema("REC::Event"));
		hipo::bank PART(factory.getSchema("REC::Particle"));
		hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
		hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
		hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
		hipo::bank RUN(factory.getSchema("RUN::config"));
		hipo::bank MCPART(factory.getSchema("MC::Particle"));
		hipo::bank MCEVENT(factory.getSchema("MC::Event"));
		hipo::bank TRACK(factory.getSchema("REC::Track"));
		hipo::bank TRAJ(factory.getSchema("REC::Traj"));

		outFile->cd();

		while (((reader.next() && IsHipo)) /*&& nbEvent < 100000*/)
		{

			nbEvent++;
			if (nbEvent % 30000 == 0)
			{
				time(&intermediate);
				double intermediate_time = difftime(intermediate, begin);

				cout << nbEvent << " events processed in " << intermediate_time << "s"
					 << "\n";
			}

			double w = 1;
			
			Event ev;
			MCEvent MC_ev;

			int run = 0;

			// Get banks
			reader.read(hipo_event);
			hipo_event.getStructure(MCPART);
			hipo_event.getStructure(MCEVENT);
			hipo_event.getStructure(RUN);
			hipo_event.getStructure(PART);
			hipo_event.getStructure(SCIN);
			hipo_event.getStructure(CHE);
			hipo_event.getStructure(CALO);
			hipo_event.getStructure(EVENT);
			hipo_event.getStructure(TRAJ);
			hipo_event.getStructure(TRACK);

			// Number of total event
			Plots.Fill_1D("evt_count", 0, 1);

			if (MCPART.getSize() < 1 && (!IsData))
				continue;

			run = RUN.getInt("run", 0);
			trigger_bit = RUN.getLong("trigger", 0);
			int np_input = PART.getRows();
			ev.Set_trigger_bit(trigger_bit);

			if (IsSimu)
			{

				MC_ev.Set_MC_Particles(MCEVENT, MCPART);
				MC_ev.Get_Kinematics();

				ev.Set_Weight(w);

				outVars_Gen["Q2_Gen"] = MC_ev.Q2_Gen;
				outVars_Gen["W_Gen"] = MC_ev.W_Gen;
				outVars_Gen["vz_elec_Gen"] = MC_ev.vz_elec_Gen;
				outVars_Gen["vz_prot_Gen"] = MC_ev.vz_prot_Gen;

				tree_Electron_Gen = MC_ev.Electron.Vector;
				tree_Proton_Gen = MC_ev.Proton.Vector;

				outT_Gen->Fill();
			}

			///////////////////////////////////////////
			// Filter good runs for data only
			///////////////////////////////////////////
			if (!Run_Selector.Is_Good_Run(run) && IsData && RGA_Fall2018)
				continue;

			///////////////////////////////////////////
			// Get Particles and cut on event topology
			///////////////////////////////////////////
			ev.Set_Particles(PART);

			if (!ev.pass_topology_cut())
			{
				continue;
			}
			///////////////////////////////////////////

			///////////////////////////////////////////
			// Associate detector responses and do EC cuts
			///////////////////////////////////////////
			ev.Apply_EC_Cuts(CALO); //Just assign CALO response, no cuts
			ev.Associate_detector_resp(CHE, SCIN);
			ev.Set_Nphe_HTCC();
			///////////////////////////////////////////

			///////////////////////////////////////////
			// TMVA
			///////////////////////////////////////////
			ElectronPID.Evaluate_BDT_new(ev.Electron);
			ev.Set_Elec_score(ElectronPID.score);
			///////////////////////////////////////////

			///////////////////////////////////////////
			// Radiative correction
			///////////////////////////////////////////
			//ev.Apply_Radiative_Correction(InputParameters.RadCorr);
			ev.Compute_SF();
			///////////////////////////////////////////

			///////////////////////////////////////////
			// Momentum MC correction
			///////////////////////////////////////////
			//ev.Apply_MC_Correction(MomCorr);
			///////////////////////////////////////////

			ev.Set_Run_Number(run);

			if ((ev.pass_EC_cut()))
			{

				///////////////////////////////////////////
				// Compute kinematics
				///////////////////////////////////////////
				ev.Get_Kinematics();

				outVars["evt_num"] = nbEvent;
				outVars["Q2"] = ev.Q2;
				outVars["W"] = ev.W;
				outVars["weight"] = ev.weight;
				outVars["run"] = ev.run;
				outVars["electron_Nphe"] = ev.electron_Nphe;
				outVars["electron_SF"] = ev.electron_SF;
				outVars["electron_score"] = ev.electron_score;
				outVars["status_elec"] = ev.Electron.status;
				outVars["status_prot"] = ev.Proton.status;
				outVars["chi2_proton"] = ev.Proton.chi2;
				outVars["vz_elec"] = ev.Electron.vertex.z;
				outVars["vz_prot"] = ev.Proton.vertex.z;
				outVars["electron_HTCC_ECAL_match"] = (ev.Electron.SectorCalo(ECAL, PCAL) == ev.Electron.SectorChe(HTCC)) ? 1. : 0.0;
				outVars["PCAL_sector_elec"] = ev.Electron.SECTOR_CALO(PCAL);
				outVars["PCAL_x_elec"] = ev.Electron.X_CALO(PCAL);
				outVars["PCAL_y_elec"] = ev.Electron.Y_CALO(PCAL);
				outVars["PCAL_U_elec"] = ev.Electron.U_CALO(PCAL);
				outVars["PCAL_V_elec"] = ev.Electron.V_CALO(PCAL);
				outVars["PCAL_W_elec"] = ev.Electron.W_CALO(PCAL);
				outVars["PCAL_energy_elec"] = ev.Electron.Energy(ECAL, PCAL);
				outVars["ECIN_energy_elec"] = ev.Electron.Energy(ECAL, ECIN);

				outVars["Q2_Gen"] = MC_ev.Q2_Gen;
				outVars["W_Gen"] = MC_ev.W_Gen;
				outVars["vz_elec_Gen"] = MC_ev.vz_elec_Gen;
				outVars["vz_prot_Gen"] = MC_ev.vz_prot_Gen;


				tree_Electron = ev.Electron.Vector;
				tree_Proton = ev.Proton.Vector;

				tree_Electron_Gen = MC_ev.Electron.Vector;
				tree_Proton_Gen = MC_ev.Proton.Vector;

				outT->Fill();
			}
		}
	}


	outT->Write();
	outFile->Write();
	outFile->Close();

	cout << "nb of file " << nbf << "\n";
	cout << "nb of events " << nbEvent << "\n";
	
	time(&end); 
	double difference = difftime(end, begin);
	printf("All this work done in only %.2lf seconds. Congratulations !\n", difference);

	gApplication->Terminate();

	return 0;
}
