#ifndef TCSFiducialCuts
#define TCSFiducialCuts

// EC cuts
const float PCALlimits[6][2][3] = {

	// sector1
	{
		// down
		{7.81688, 9.78924, 9.47359},
		// up
		{403.388, 402.06, 393.895}},

	// sector2
	{
		// down
		{6.55955, 8.62768, 8.57818},
		// up
		{403.687, 402.389, 402.064}},

	// sector3
	{
		// down
		{6.88656, 9.23112, 8.23956},
		// up
		{403.566, 403.875, 403.622}},

	// sector4
	{
		// down
		{7.14218, 19.2814, 8.26354},
		// up
		{404.295, 403.021, 392.355}},

	// sector5
	{
		// down
		{6.6288, 8.73336, 9.28017},
		// up
		{402.824, 402.915, 403.634}},

	// sector6
	{
		// down
		{7.68173, 9.12088, 8.13996},
		// up
		{402.915, 403.681, 403.886}}};

const float ECINlimits[6][2][3] = {

	// sector1
	{
		// down
		{13.1302, 10.1434, 13.4929},
		// up
		{397.865, 398.946, 363.28}},

	// sector2
	{
		// down
		{24.8122, 23.3301, 15.0611},
		// up
		{385.095, 398.682, 353.012}},

	// sector3
	{
		// down
		{35.1522, 18.2694, 22.4406},
		// up
		{396.237, 397.406, 353.071}},

	// sector4
	{
		// down
		{10.7544, 15.3729, 30.4228},
		// up
		{394.161, 398.805, 363.549}},

	// sector5
	{
		// down
		{10.912, 35.8515, 24.0282},
		// up
		{396.647, 398.616, 363.371}},

	// sector6
	{
		// down
		{22.6877, 25.3906, 29.905},
		// up
		{396.942, 397.33, 363.586}}};

const float ECOUTlimits[6][2][3] = {

	// sector1
	{
		// down
		{10.4308, 56.8701, 35.1741},
		// up
		{399.925, 400.357, 378.102}},

	// sector2
	{
		// down
		{10.0306, 18.2873, 24.7615},
		// up
		{410.416, 412.219, 366.346}},

	// sector3
	{
		// down
		{11.1158, 15.9105, 23.3541},
		// up
		{410.511, 412.472, 377.301}},

	// sector4
	{
		// down
		{9.35355, 17.7095, 45.6757},
		// up
		{411.27, 412.996, 377.179}},

	// sector5
	{
		// down
		{14.9061, 34.8015, 16.3359},
		// up
		{408.767, 412.057, 377.174}},

	// sector6
	{
		// down
		{9.43871, 28.0928, 43.4728},
		// up
		{409.363, 410.544, 378.023}}};


bool pass_Dead_Paddle_PCAL(int sector, float U, float V, float W)
{

	if (sector == 1)
	{
		return !( (W>72. && W<93.) ||  (W>210. && W<231.) );
	}
	else if (sector == 2)
	{
		return !( (V>100. && V<115.) );
	}
	else if (sector == 4)
	{
		return !( (V>228. && V<242.) ||  (V<15.) );
	}
	else if (sector == 6)
	{
		return !( (W>170. && V<194.) );
	}
	else
		return true;

}



Particle ApplyECcuts(Particle particle, hipo::bank CALO)
{

	bool pass = true;
	CalorimeterResp Calo;
	for (int c = 0; c < CALO.getRows(); c++)
	{

		int Calopindex = CALO.getInt("pindex", c);
		int Calosector = CALO.getInt("sector", c);
		int Calolayer = CALO.getInt("layer", c);
		int Calodetector = CALO.getInt("detector", c);
		float Caloenergy = CALO.getFloat("energy", c);
		float Calox = CALO.getFloat("x", c);
		float Caloy = CALO.getFloat("y", c);
		float Caloz = CALO.getFloat("z", c);
		float Calou = CALO.getFloat("lu", c);
		float Calov = CALO.getFloat("lv", c);
		float Calow = CALO.getFloat("lw", c);
		float Calodu = CALO.getFloat("du", c);
		float Calodv = CALO.getFloat("dv", c);
		float Calodw = CALO.getFloat("dw", c);
		float Calom2u = CALO.getFloat("m2u", c);
		float Calom2v = CALO.getFloat("m2v", c);
		float Calom2w = CALO.getFloat("m2w", c);
		float Calom3u = CALO.getFloat("m3u", c);
		float Calom3v = CALO.getFloat("m3v", c);
		float Calom3w = CALO.getFloat("m3w", c);

		if (Calopindex == (particle.index))
		{
			Calo.detector = Calodetector;
			Calo.pindex = Calopindex;
			Calo.sector = Calosector;
			Calo.layer = Calolayer;
			Calo.energy = Caloenergy;
			Calo.x = Calox;
			Calo.y = Caloy;
			Calo.z = Caloz;
			Calo.u = Calou;
			Calo.v = Calov;
			Calo.w = Calow;
			Calo.du = Calodu;
			Calo.dv = Calodv;
			Calo.dw = Calodw;
			Calo.m2u = Calom2u;
			Calo.m2v = Calom2v;
			Calo.m2w = Calom2w;
			Calo.m3u = Calom3u;
			Calo.m3v = Calom3v;
			Calo.m3w = Calom3w;

			if (Calolayer == 1 && pass_Dead_Paddle_PCAL(Calosector, Calou, Calov, Calow) /*&& Calov > PCALlimits[Calosector - 1][0][1] && Calov < PCALlimits[Calosector - 1][1][1] && Calow > PCALlimits[Calosector - 1][0][2] && Calow < PCALlimits[Calosector - 1][1][2]*/)
			{
				particle.Calorimeter.push_back(Calo);
			}

			else if (Calolayer == 4 /*&& Calou>ECINlimits[Calosector-1][0][0] && Calou<ECINlimits[Calosector-1][1][0] && Calov>ECINlimits[Calosector-1][0][1] && Calov<ECINlimits[Calosector-1][1][1] && Calow>ECINlimits[Calosector-1][0][2] && Calow<ECINlimits[Calosector-1][1][2]*/)
			{
				particle.Calorimeter.push_back(Calo);
			}

			else if (Calolayer == 7 /* && Calou>ECOUTlimits[Calosector-1][0][0] && Calou<ECOUTlimits[Calosector-1][1][0] && Calov>ECOUTlimits[Calosector-1][0][1] && Calov<ECOUTlimits[Calosector-1][1][1] && Calow>ECOUTlimits[Calosector-1][0][2] && Calow<ECOUTlimits[Calosector-1][1][2]*/)
			{
				particle.Calorimeter.push_back(Calo);
			}
			else
				pass = false;
		}
	}
	particle.passEC = pass;
	return particle;
}

#endif
