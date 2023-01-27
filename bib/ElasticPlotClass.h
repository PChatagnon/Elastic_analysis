#ifndef ElasticPlotClass
#define ElasticPlotClass

const std::vector<vector<TString>> ListPlot1D{

    // Format : { "name", "title;tile_x;title_y" , bin_x , min_x , max_x }
    {"evt_count", "evt_count", "6", "0", "6"},
    {"efficiency", "efficiency", "7", "0", "7"},
   };

const std::vector<vector<TString>> ListPlot2D{

   };

class ElasticPlots
{
public:
        TString outputFolder = "Plots";
        std::unordered_map<string, TH1F *> Plot1D;
        std::unordered_map<string, TH2F *> Plot2D;

        void SetOutputFolder(TString input_outputFolder)
        {
                outputFolder = input_outputFolder;
        }

        void SavePlots(TString format)
        {
                for (auto plot : Plot1D)
                {
                        TCanvas *can = new TCanvas("", "can", 1500, 1000);
                        can->cd();
                        plot.second->Draw();
                        can->SaveAs(outputFolder + (plot.second->GetName()) + "." + format);
                }

                for (auto plot : Plot2D)
                {
                        TCanvas *can = new TCanvas("", "can", 1500, 1000);
                        can->cd();
                        plot.second->Draw("col");
                        can->SaveAs(outputFolder + (plot.second->GetName()) + "." + format);
                }
        }

        void Initialize_1D()
        {
                for (int i = 0; i < ListPlot1D.size(); i++)
                {
                        TString name = ListPlot1D[i][0];
                        TString title = ListPlot1D[i][1];
                        TString bin_x = ListPlot1D[i][2];
                        TString min_x = ListPlot1D[i][3];
                        TString max_x = ListPlot1D[i][4];

                        Plot1D[(string)name.Data()] = new TH1F(name, title, stoi((string)bin_x.Data()), stof((string)min_x.Data()), stof((string)max_x.Data()));
                }
        }

        void Initialize_2D()
        {
                for (int i = 0; i < ListPlot2D.size(); i++)
                {
                        TString name = ListPlot2D[i][0];
                        TString title = ListPlot2D[i][1];
                        TString bin_x = ListPlot2D[i][2];
                        TString min_x = ListPlot2D[i][3];
                        TString max_x = ListPlot2D[i][4];
                        TString bin_y = ListPlot2D[i][5];
                        TString min_y = ListPlot2D[i][6];
                        TString max_y = ListPlot2D[i][7];

                        Plot2D[(string)name.Data()] = new TH2F(name, title, stoi((string)bin_x.Data()), stof((string)min_x.Data()), stof((string)max_x.Data()), stoi((string)bin_y.Data()), stof((string)min_y.Data()), stof((string)max_y.Data()));
                }
        }

        void Fill_1D(TString name, double value, double weight)
        {
                Plot1D[(string)name.Data()]->Fill(value, weight);
        }

        void Fill_2D(TString name, double value_x, double value_y, double weight)
        {
                Plot2D[(string)name.Data()]->Fill(value_x, value_y, weight);
        }

        void Add_Hist_1D(TString name, TString title, int bin, float min_x, float max_x)
        {
                Plot1D[(string)name.Data()] = new TH1F(name, title, bin, min_x, max_x);
        }

        void Add_Hist_2D(TString name, TString title, int bin_x, float min_x, float max_x, int bin_y, float min_y, float max_y)
        {
                Plot2D[(string)name.Data()] = new TH2F(name, title, bin_x, min_x, max_x, bin_y, min_y, max_y);
        }

        /*   void Draw_Hist_1D(TString name)
           {
                   TCanvas *can = new TCanvas("canChi2", "canChi2", 4000, 2000);
           Chi2ElectronAvant->Draw();
           can->SaveAs("canChi2.pdf");
           can->SaveAs("canChi2.root");
           Chi2ProtonAvant->SaveAs("Chi2ProtonAvant.root");
           }

           void Draw_Hist_2D(TString name)
           {
                    TCanvas *can = new TCanvas("canChi2", "canChi2", 4000, 2000);
           Chi2ElectronAvant->Draw();
           can->SaveAs("canChi2.pdf");
           can->SaveAs("canChi2.root");
           Chi2ProtonAvant->SaveAs("Chi2ProtonAvant.root");
           }
   */
        void Draw_All_1D()
        {
        }

        void Draw_All_2D()
        {
        }
};

#endif
