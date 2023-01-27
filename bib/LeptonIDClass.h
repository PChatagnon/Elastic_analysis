#ifndef LeptonIDClass
#define LeptonIDClass

class LeptonIdentification
{
public:
        TString methodName;
        TString weightFile;
        double scoreCut;
        double momentumThreshold;
        Float_t P, Theta, Phi;
        Float_t SFPCAL, SFECIN, SFECOUT;
        Float_t m2PCAL, m2ECIN, m2ECOUT;
        TMVA::Reader *PositronIdentificationReader;
        double score;

        LeptonIdentification(TString input_methodName, TString input_weightFile, double input_scoreCut, double input_momentumThreshold)
        {
                methodName = input_methodName;
                weightFile = input_weightFile;
                PositronIdentificationReader = new TMVA::Reader("!Color:!Silent");
                scoreCut = input_scoreCut;
                momentumThreshold = input_momentumThreshold;

                score = -1.;
        }

        void AddVariable(Float_t input_variable, TString input_variable_string)
        {
                PositronIdentificationReader->AddVariable(input_variable_string, &input_variable);
        }

        void InitializePositronIdentification()
        {
                PositronIdentificationReader->AddVariable("SFPCAL", &SFPCAL);
                PositronIdentificationReader->AddVariable("SFECIN", &SFECIN);
                PositronIdentificationReader->AddVariable("SFECOUT", &SFECOUT);
                PositronIdentificationReader->AddVariable("m2PCAL", &m2PCAL);
                PositronIdentificationReader->AddVariable("m2ECIN", &m2ECIN);
                PositronIdentificationReader->AddVariable("m2ECOUT", &m2ECOUT);
                PositronIdentificationReader->BookMVA(methodName, weightFile);
        }

        void InitializeBDT_new_PositronIdentification()
        {
                PositronIdentificationReader->AddVariable("P", &P);
                PositronIdentificationReader->AddVariable("Theta", &Theta);
                PositronIdentificationReader->AddVariable("Phi", &Phi);
                PositronIdentificationReader->AddVariable("SFPCAL", &SFPCAL);
                PositronIdentificationReader->AddVariable("SFECIN", &SFECIN);
                PositronIdentificationReader->AddVariable("SFECOUT", &SFECOUT);
                PositronIdentificationReader->AddVariable("m2PCAL", &m2PCAL);
                PositronIdentificationReader->AddVariable("m2ECIN", &m2ECIN);
                PositronIdentificationReader->AddVariable("m2ECOUT", &m2ECOUT);
                PositronIdentificationReader->BookMVA(methodName, weightFile);
        }

        void Evaluate(Particle vLepton)
        {

                SFPCAL = (vLepton.Energy(ECAL, PCAL)) / vLepton.Vector.P();
                SFECIN = (vLepton.Energy(ECAL, ECIN)) / vLepton.Vector.P();
                SFECOUT = (vLepton.Energy(ECAL, ECOUT)) / vLepton.Vector.P();

                double M2PCAL = -1.;
                double M2ECIN = -1.;
                double M2ECOUT = -1.;

                for (int i = 0; i < vLepton.Calorimeter.size(); i++)
                {
                        if (PCAL == vLepton.Calorimeter[i].layer)
                        {
                                M2PCAL = (vLepton.Calorimeter[i].m2u + vLepton.Calorimeter[i].m2v + vLepton.Calorimeter[i].m2w) / 3.;
                        }
                        if (ECIN == vLepton.Calorimeter[i].layer)
                        {
                                M2ECIN = (vLepton.Calorimeter[i].m2u + vLepton.Calorimeter[i].m2v + vLepton.Calorimeter[i].m2w) / 3.;
                        }
                        if (ECOUT == vLepton.Calorimeter[i].layer)
                        {
                                M2ECOUT = (vLepton.Calorimeter[i].m2u + vLepton.Calorimeter[i].m2v + vLepton.Calorimeter[i].m2w) / 3.;
                        }
                }

                m2PCAL = M2PCAL;
                m2ECIN = M2ECIN;
                m2ECOUT = M2ECOUT;
                
                score = PositronIdentificationReader->EvaluateMVA(methodName);
        }

        void Evaluate_BDT_new(Particle vLepton)
        {

                SFPCAL = (vLepton.Energy(ECAL, PCAL)) / vLepton.Vector.P();
                SFECIN = (vLepton.Energy(ECAL, ECIN)) / vLepton.Vector.P();
                SFECOUT = (vLepton.Energy(ECAL, ECOUT)) / vLepton.Vector.P();

                double M2PCAL = -1.;
                double M2ECIN = -1.;
                double M2ECOUT = -1.;

                for (int i = 0; i < vLepton.Calorimeter.size(); i++)
                {
                        if (PCAL == vLepton.Calorimeter[i].layer)
                        {
                                M2PCAL = (vLepton.Calorimeter[i].m2u + vLepton.Calorimeter[i].m2v + vLepton.Calorimeter[i].m2w) / 3.;
                        }
                        if (ECIN == vLepton.Calorimeter[i].layer)
                        {
                                M2ECIN = (vLepton.Calorimeter[i].m2u + vLepton.Calorimeter[i].m2v + vLepton.Calorimeter[i].m2w) / 3.;
                        }
                        if (ECOUT == vLepton.Calorimeter[i].layer)
                        {
                                M2ECOUT = (vLepton.Calorimeter[i].m2u + vLepton.Calorimeter[i].m2v + vLepton.Calorimeter[i].m2w) / 3.;
                        }
                }

                m2PCAL = M2PCAL;
                m2ECIN = M2ECIN;
                m2ECOUT = M2ECOUT;

                P = vLepton.Vector.P();
                Theta = vLepton.Vector.Theta();
                Phi = vLepton.Vector.Phi();
                
                score = PositronIdentificationReader->EvaluateMVA(methodName);
        }

        bool Accept(Particle vPositron)
        {

                //cout<<vPositron.Vector.P()<<" "<<momentumThreshold<<" "<<score<<" "<<scoreCut<<endl;
                if ((float)vPositron.Vector.P() > (float)momentumThreshold)
                        return ((float)score > (float)scoreCut);
                else return true;
        }
};

#endif
