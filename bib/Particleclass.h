#ifndef Particleclass
#define Particleclass

const float me = 0.000511;
const float mp = 0.938;
float ebeam = 10.604;

const float Pi = 3.14159269;

const bool PIDflag = true;
const bool SFdist = true;

const int LTCC = 16;
const int HTCC = 15;
const int ECAL = 7;
const int DC = 6;

const int PCAL = 1;
const int ECIN = 4;
const int ECOUT = 7;

class Track
{

public:
        int detector = 0;
        int pindex = 0;
        int sector = 0;
        int layer = 0;
        float x = 0.;
        float y = 0.;
        float z = 0.;
        float chi2 = 0.0;
};

class Traj
{

public:
        int detector = 0;
        int pindex = 0;
        int layer = 0;
        float x = -10000.;
        float y = -10000.;
        float z = -10000.;
        float edge = -10000.;
};

class CalorimeterResp
{

public:
        int detector, pindex, sector, layer;
        float energy, x, y, z, u, v, w;
        float du, dv, dw;
        float m2u, m2v, m2w;
        float m3u, m3v, m3w;
};

class CheResp
{

public:
        int detector, pindex, sector;
        float nphe, x, y, z, time, chi2;
};

class ScinResp
{

public:
        int detector, pindex, sector;
        float t, path;
        float energy = 0.0;
};

class vertex
{
public:
        float x;
        float y;
        float z;
};

class Particle
{

public:
        int index = -1;

        int pid = -5000;
        float beta = -5000.;
        float vt = -5000.;
        int status = 0;
        float chi2 = -5000.;
        bool passEC = true;
        vertex vertex;
        TLorentzVector Vector;
        ScinResp Scintillator;

        vector<CalorimeterResp> Calorimeter;
        vector<CheResp> Cherenkov;
        vector<Traj> Trajs;

        float nphe(int det)
        {
                float np = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                np += Cherenkov[i].nphe;
                }
                return np;
        }

        float TimeChe(int det)
        {
                float time = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                time = Cherenkov[i].time;
                                
                }
                return time;
        }

        float SectorChe(int det)
        {
                int sector = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                        {
                                sector = Cherenkov[i].sector;
                                break;
                        }
                }
                return sector;
        }

        float SectorCalo(int det, int layer)
        {
                int sector = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer)
                        {
                                if (det == Calorimeter[i].detector)
                                {
                                        sector = Calorimeter[i].sector;
                                        break;
                                }
                        }
                }
                return sector;
        }

        float Chi2Che(int det)
        {
                float chi2 = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                chi2 = Cherenkov[i].chi2;
                }
                return chi2;
        }

        float XChe(int det)
        {
                float X = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                X = Cherenkov[i].x;
                }
                return X;
        }

        float YChe(int det)
        {
                float Y = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                Y = Cherenkov[i].y;
                }
                return Y;
        }

        float X_CALO(int layer)
        {
                float X = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                X = Calorimeter[i].x;
                                break;
                }
                return X;
        }

        float Y_CALO(int layer)
        {
                float Y = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].y;
                                break;
                }
                return Y;
        }

        float U_CALO(int layer)
        {
                float Y = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].u;
                                break;
                }
                return Y;
        }

        float V_CALO(int layer)
        {
                float Y = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].v;
                                break;
                }
                return Y;
        }

        float W_CALO(int layer)
        {
                float Y = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].w;
                                break;
                }
                return Y;
        }

        float SECTOR_CALO(int layer)
        {
                float Y = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].sector;
                                break;
                }
                return Y;
        }

        float Energy(int det, int layer)
        {
                float en = 0.;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer)
                        {
                                if (det == Calorimeter[i].detector)
                                        en = en + Calorimeter[i].energy;
                        }
                }
                return en;
        }

        void Associate_DC_traj_to_Particle(hipo::bank TRAJ)
        {
                
                for (int t = 0; t < TRAJ.getRows(); t++)
                {
                        Traj new_Traj;
                        int TRAJ_detector = TRAJ.getInt("detector", t);
                        int TRAJ_pindex = TRAJ.getInt("pindex", t);
                        int TRAJ_layer = TRAJ.getInt("layer", t);
                        float TRAJ_X = TRAJ.getFloat("x", t);
                        float TRAJ_Y = TRAJ.getFloat("y", t);
                        float TRAJ_Z = TRAJ.getFloat("z", t);
                        float TRAJ_edge = TRAJ.getFloat("edge", t);

                        if (TRAJ_pindex == index && TRAJ_detector == DC)
                        { 
                                new_Traj.detector = TRAJ_detector;
                                new_Traj.pindex = TRAJ_pindex;
                                new_Traj.layer = TRAJ_layer;
                                new_Traj.x = TRAJ_X;
                                new_Traj.y = TRAJ_Y;
                                new_Traj.z = TRAJ_Z; 
                                new_Traj.edge = TRAJ_edge; 
                                Trajs.push_back(new_Traj);
                        }
                }     
        }
};

class ThetaPhi
{
public:
        float theta;
        float phi;
};


#endif
