#include <cfloat>

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