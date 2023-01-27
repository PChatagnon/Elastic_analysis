#ifndef ElasticMCEvent
#define ElasticMCEvent

class MCEvent
{

public:
        // Particle and Lorentz vectors
        TLorentzVector Electron;
        TLorentzVector Proton;
        TLorentzVector vBeam;
        TLorentzVector vRestProton;

        // Kinematic variables
        float Q2_Gen;
        float W_Gen;

        // MC vertex
        float vz_elec_Gen;
        float vz_prot_Gen;

        // MC weight
        float w;

        MCEvent()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_MC_Particles(hipo::bank MCEVENT, hipo::bank MCPART)
        {

                Electron.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), me);
                Proton.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), mp);

                vz_elec_Gen = MCPART.getFloat("vz", 0);
                vz_prot_Gen = MCPART.getFloat("vz", 1);
        }

        void Get_Kinematics()
        {
                W_Gen = (vBeam + vRestProton - Electron).M();
                Q2_Gen = (vBeam - Electron).M2();
        }
};

#endif
