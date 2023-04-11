#ifndef ElasticEvent
#define ElasticEvent

class Event
{

public:
        // Particle and Lorentz vectors
        Particle Electron;
        Particle Proton;
        TLorentzVector vBeam;
        TLorentzVector vRestProton;

        // Event weight
        float weight = 1;

        // Number of reconstructed particle of each type
        int np = 2; // total number of particles
        int recem = 0;
        int recp = 0;

        // Kinematic variables
        float W;
        float Q2;

        // Electron ID
        float electron_SF;
        float electron_score;
        float electron_Nphe;

        // run number and trigger bit
        float run;
        int trigger_bit;

        Event()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_trigger_bit(long input_trigger_bit)
        {
                trigger_bit = (long)input_trigger_bit;
        }

        void Set_Weight(float w_in)
        {
                weight = w_in;
        }

        void Set_Particles(hipo::bank PART)
        {

                for (int i = 0; i < PART.getRows(); i++)
                {

                        int pid = PART.getInt("pid", i);
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float beta = PART.getFloat("beta", i);
                        int status = abs(PART.getInt("status", i));
                        int charge = PART.getInt("charge", i);
                        float chi2 = PART.getFloat("chi2pid", i);
                        float vx = PART.getFloat("vx", i);
                        float vy = PART.getFloat("vy", i);
                        float vz = PART.getFloat("vz", i);
                        float vt = PART.getFloat("vt", i);

                        if (pid == 11)
                        {
                                if (status > 2000 && recem==0)
                                {
                                        Electron.Vector.SetXYZM(px, py, pz, me);
                                        Electron.index = i;
                                        Electron.pid = 11;
                                        Electron.beta = beta;
                                        Electron.status = status;
                                        Electron.chi2 = chi2;
                                        Electron.vertex.x = vx;
                                        Electron.vertex.y = vy;
                                        Electron.vertex.z = vz;
                                        Electron.vt = vt;
                                        recem++;
                                }
                        }

                        if (pid == 2212)
                        {
                                Proton.Vector.SetXYZM(px, py, pz, mp);
                                Proton.index = i;
                                Proton.pid = 2212;
                                Proton.beta = beta;
                                Proton.status = status;
                                Proton.chi2 = chi2;
                                Proton.vertex.x = vx;
                                Proton.vertex.y = vy;
                                Proton.vertex.z = vz;
                                Proton.vt = vt;
                                recp++;
                        }
                }
        }

        bool pass_topology_cut()
        {
                return (recem > 0);
        }

        void Apply_EC_Cuts(hipo::bank CALO)
        {
                Electron = ApplyECcuts(Electron, CALO);
                Proton = ApplyECcuts(Proton, CALO);
        }

        void Associate_DC_traj(hipo::bank TRAJ)
        {
                Electron.Associate_DC_traj_to_Particle(TRAJ);
                //Proton.Associate_DC_traj_to_Particle(TRAJ);  
        }

        bool pass_EC_cut()
        {
                return (Electron.passEC && Proton.passEC);
        }

        void Compute_SF()
        {
                electron_SF = ((Electron.Energy(ECAL, PCAL) + Electron.Energy(ECAL, ECIN) + Electron.Energy(ECAL, ECOUT))) / Electron.Vector.P();
        }

        void Set_Elec_score(float input_elec_score)
        {
                electron_score = input_elec_score;
        }

        void Set_Nphe_HTCC()
        {
                electron_Nphe = Electron.nphe(15);
        }

        void Associate_detector_resp(hipo::bank CHE, hipo::bank SCIN)
        {
                vector<Particle> Particles = {Electron, Proton};

                CalorimeterResp Calo;
                CheResp Che;
                ScinResp Scin;
                for (int i = 0; i < 2; i++)
                {

                        for (int c = 0; c < CHE.getRows(); c++)
                        {
                                int Chepindex = CHE.getInt("pindex", c);
                                int Chedetector = CHE.getInt("detector", c);
                                int Chesector = CHE.getInt("sector", c);
                                float Chenphe = CHE.getFloat("nphe", c);
                                float Chetime = CHE.getFloat("time", c);
                                float Chechi2 = CHE.getFloat("chi2", c);
                                float Chex = CHE.getFloat("x", c);
                                float Chey = CHE.getFloat("y", c);
                                float Chez = CHE.getFloat("z", c);

                                if (Chepindex == (Particles[i].index))
                                {
                                        Che.detector = Chedetector;
                                        Che.pindex = Chepindex;
                                        Che.sector = Chesector;
                                        Che.nphe = Chenphe;
                                        Che.time = Chetime;
                                        Che.chi2 = Chechi2;
                                        Che.x = Chex;
                                        Che.y = Chey;
                                        Che.z = Chez;
                                        Particles[i].Cherenkov.push_back(Che);
                                }
                        }
                        for (int c = 0; c < SCIN.getRows(); c++)
                        {
                                int Scindetector = SCIN.getInt("detector", c);
                                int Scinpindex = SCIN.getInt("pindex", c);
                                float Scintime = SCIN.getFloat("time", c);
                                float Scinpath = SCIN.getFloat("path", c);
                                float Scinenergy = SCIN.getFloat("energy", c);
                                int Scinsector = SCIN.getInt("sector", c);

                                if (Scinpindex == (Particles[i].index))
                                {
                                        Scin.detector = Scindetector;
                                        Scin.pindex = Scinpindex;
                                        Scin.t = Scintime;
                                        Scin.path = Scinpath;
                                        Scin.energy = Scinenergy;
                                        Scin.sector = Scinsector;

                                        if (Particles[i].Scintillator.energy < Scinenergy)
                                        {
                                                Particles[i].Scintillator = Scin;
                                        };
                                }
                        }
                }

                Electron = Particles[0];
                Proton = Particles[1];

                //Apply momentum corrections to electron
                double ex = Electron.Vector.Px();
                double ey = Electron.Vector.Py();
                double ey = Electron.Vector.Pz();
                int esec = Electron.SECTOR_CALO(PCAL);
                double fe = dppC(ex, ey, ez, esec, 0) + 1;
                TLorentzVector eleC = ROOT::Math::PxPyPzMVector(ex*fe, ey*fe, ez*fe, 0);
                Electron.Vector = eleC;
        }

        /*void Apply_MC_Correction(MomentumCorrection MomCorr)
        {
                Proton = MomCorr.Apply_MC_Correction(Proton);
        }*/

        void Get_Kinematics()
        {
                W = (vBeam + vRestProton - Electron.Vector).M();
                Q2 = (vBeam - Electron.Vector).M2();
        }

        void Set_SF(float input_electron_SF)
        {
                electron_SF = input_electron_SF;
        }

        void Set_Run_Number(int input_run)
        {
                run = (float)(input_run);
        }
};

#endif
