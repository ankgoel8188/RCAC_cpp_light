// #include <cfloat>
#include <px4_platform_common/defines.h>

struct RCACTuneParams {
    float p0            = 0;
    float N_nf          = 1;
    float Ru            = 1;
    float alpha_PID     = 1;
    float RCAC_EN       = false;
};

struct RCACInitParams {
    float lambda        = 1;
    float Rz            = 1;
    int   errorNormMode = 0;
    float lim_int       = FLT_MAX;
    bool  RBlock_EN     = false;
};

struct RCACParams {
    RCACTuneParams tuneParams;
    RCACInitParams initParams;
};


// Contains setters for tunable parameters; Getters for all Parameters
class RCACParams_IO {
    public:

    RCACParams_IO() : RCACParams_IO(nullptr) {}
    RCACParams_IO(RCACParams * paramsPtr_in) : paramsPtr(paramsPtr_in) {}

    // Struct Block Access to the RCAC Params
    const RCACParams & get_RCACParams() { return *paramsPtr; }

    const RCACTuneParams & get_tuneParams() { return paramsPtr->tuneParams; }

    const RCACInitParams & get_initParams() { return paramsPtr->initParams; }

    void set_tuneParams(const RCACTuneParams & tuneParams_in ) { paramsPtr->tuneParams = tuneParams_in; }

    void set_RCACParams(const RCACParams & RCACParams_in ) { *paramsPtr = RCACParams_in; }

    // Individual Variable Access to the RCAC Params

	float get_Ru() { return paramsPtr->tuneParams.Ru; }
	float get_Rz() { return paramsPtr->initParams.Rz; }
	float get_P0() { return paramsPtr->tuneParams.p0; }
	float get_alpha() { return paramsPtr->tuneParams.alpha_PID; }
	bool  get_switch() { return paramsPtr->tuneParams.RCAC_EN; }

    	void set_Ru(float Ru_in)
	{
		// if (pitch_Ru_in != rcac_pitch_Ru)
		PX4_INFO("[RCAC] Ru: %6.4f", (double)Ru_in);
		paramsPtr->tuneParams.Ru = Ru_in;
	}

	void set_P0(float P0_in)
	{
		// if (rcac_pitch_P0 != pitch_P0_in)
		PX4_INFO("[RCAC] P0: %6.4f", (double)P0_in);
		paramsPtr->tuneParams.p0 = P0_in;
	}

	void set_alpha(float alpha_in)
	{
		PX4_INFO("[RCAC] Alpha: %6.4f", (double)alpha_in);
		paramsPtr->tuneParams.alpha_PID = alpha_in;
	}

	void set_RCAC_EN(bool RCAC_EN_in)
	{
		PX4_INFO("[RCAC] RCAC Switch: %s", pitch_SW_in ? "true" : "false");
		paramsPtr->tuneParams.RCAC_EN = RCAC_EN_in;
	}

	void set_N_nf(bool N_nf_in)
	{
		PX4_INFO("[RCAC] N_nf: %6.4f", (double)N_nf_in);
		paramsPtr->tuneParams.N_nf = N_nf_in;
	}

    private:
    RCACParams * paramsPtr;
};
