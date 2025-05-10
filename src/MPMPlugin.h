

#ifndef __MPM_PLUGIN_h__
#define __MPM_PLUGIN_h__

//#include <GEO/GEO_Point.h>
//
#include <SOP/SOP_Node.h>
#include <UT/UT_Vector3.h>
#include <Eigen/Dense>
#include "MPMSolver.h"

namespace HDK_Sample {
    class SOP_SnowSim : public SOP_Node
    {
    public:
        static OP_Node* myConstructor(OP_Network*, const char*,
            OP_Operator*);

        void readSDFFromVDB(const GU_Detail* sdfGdp);

        /// Stores the description of the interface of the SOP in Houdini.
        /// Each parm template refers to a parameter.
        static PRM_Template		 myTemplateList[];

        /// This optional data stores the list of local variables.
        static CH_LocalVariable	 myVariables[];

        MPMSolver solver;

        const GU_PrimVDB* vdbPrimSDF;

        


    protected:

        SOP_SnowSim(OP_Network* net, const char* name, OP_Operator* op);
        virtual ~SOP_SnowSim();

        bool isTimeDependent() const;



        /// Disable parameters according to other parameters.
        virtual unsigned		 disableParms();


        /// cookMySop does the actual work of the SOP computing
        virtual OP_ERROR		 cookMySop(OP_Context& context);

        /// This function is used to lookup local variables that you have
        /// defined specific to your SOP.
        virtual bool		 evalVariableValue(
            fpreal& val,
            int index,
            int thread);
        // Add virtual overload that delegates to the super class to avoid
        // shadow warnings.
        virtual bool		 evalVariableValue(
            UT_String& v,
            int i,
            int thread)
        {
            return evalVariableValue(v, i, thread);
        }

    private:
         ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        int		myCurrPoint;
        int		myTotalPoints;
        int     prevFrame;


        // MATERIAL PROPERTIES TAB
        fpreal PARTICLE_SEP(fpreal t) { return evalFloat("particle_sep", 0, t); }
        fpreal CRIT_COMPRESSION(fpreal t) { return evalFloat("crit_compression", 0, t); }
        fpreal CRIT_STRETCH(fpreal t) { return evalFloat("crit_stretch", 0, t); }
        fpreal HARDENING(fpreal t) { return evalFloat("hardening", 0, t); }
        fpreal INIT_DENSITY(fpreal t) { return evalFloat("init_density", 0, t); }
        fpreal YOUNG_MODULUS(fpreal t) { return evalFloat("young_modulus", 0, t); }
        fpreal POISSON(fpreal t) { return evalFloat("poisson", 0, t); }

		// SIMULATION TAB
		fpreal TIME_STEP(fpreal t) { return evalFloat("dt", 0, t); }

        UT_Vector3 GRAVITY_FORCE(fpreal t) {
            return UT_Vector3(evalFloat("grav", 0, t), evalFloat("grav", 1, t), evalFloat("grav", 2, t));
        }

        fpreal GROUND_PLANE(fpreal t) { return evalFloat("ground_plane", 0, t); }

        int FRAME_START(fpreal t) {
            return evalInt("sim_range", 0, t);
        }

        int FRAME_END(fpreal t) {
            return evalInt("sim_range", 1, t);
        }

        UT_Vector3 BOUNDS_SIZE(fpreal t) {
            return UT_Vector3(evalFloat("sim_bounds", 0, t), evalFloat("sim_bounds", 1, t), evalFloat("sim_bounds", 2, t));
        }

        UT_Vector3 BOUNDS_POS(fpreal t) {
            return UT_Vector3(evalFloat("sim_pos", 0, t), evalFloat("sim_pos", 1, t), evalFloat("sim_pos", 2, t));
        }



    };
}

#endif
