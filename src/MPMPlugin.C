


#include <UT/UT_DSOVersion.h>
//#include <RE/RE_EGLServer.h>


#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>


#include <limits.h>
#include "MPMPlugin.h"
#include "MPMSolver.h"
using namespace HDK_Sample;

//
// Help is stored in a "wiki" style text file. 
//
// See the sample_install.sh file for an example.
//
// NOTE : Follow this tutorial if you have any problems setting up your visual studio 2008 for Houdini 
//  http://www.apileofgrains.nl/setting-up-the-hdk-for-houdini-12-with-visual-studio-2008/


///
/// newSopOperator is the hook that Houdini grabs from this dll
/// and invokes to register the SOP.  In this case we add ourselves
/// to the specified operator table.
///
void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(
	    new OP_Operator("CusSnowSim",			// Internal name
			    "SnowSim",			// UI name
				 SOP_SnowSim::myConstructor,	// How to build the SOP
				 SOP_SnowSim::myTemplateList,	// My parameters
			     1,				// Min # of sources
			     2,				// Max # of sources
				 SOP_SnowSim::myVariables,	// Local variables
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//PUT YOUR CODE HERE
//You need to declare your parameters here
//Example to declare a variable for angle you can do like this :
//static PRM_Name		angleName("angle", "Angle");

//// Material Properties Tab
//static PRM_Name separationName("particle_sep", "Particle Separation");
//static PRM_Name compressionName("crit_compression", "Critical Compression");
//static PRM_Name stretchName("crit_stretch", "Critical Stretch");
//static PRM_Name hardeningName("hardening", "Hardening Coefficient");
//static PRM_Name densityName("init_density", "Initial Density");
//static PRM_Name youngName("young_modulus", "Initial Young's Modulus");
//static PRM_Name poissonName("poisson", "Poisson's Ratio");
//
//// Simulation Tab
//static PRM_Name gravityName("gravity", "Gravity Force");
//static PRM_Name groundName("ground_plane", "Ground Plane");
//static PRM_Name resetName("reset_cache", "Reset Cache");




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version


// PUT YOUR CODE HERE
// You need to setup the initial/default values for your parameters here
// For example : If you are declaring the inital value for the angle parameter
// static PRM_Default angleDefault(30.0);	


//// DEFAULT VALUES
//static PRM_Default separationDefault(0.1);
//static PRM_Default compressionDefault(2.5);
//static PRM_Default stretchDefault(1.5);
//static PRM_Default hardeningDefault(10.0);
//static PRM_Default densityDefault(1000.0);
//static PRM_Default youngDefault(1e5);
//static PRM_Default poissonDefault(0.2);
//
//// Default for vector parameters
//static PRM_Default gravityDefault[] = {
//	PRM_Default(0.0), PRM_Default(-9.8), PRM_Default(0.0)
//};
//static PRM_Default groundDefault[] = {
//	PRM_Default(0.0), PRM_Default(0.0), PRM_Default(0.0)
//};
//
//
//// RANGES
//static PRM_Range positiveRange(PRM_RANGE_PRM, 0.0, PRM_RANGE_UI, 10000.0);
//static PRM_Range poissonRange(PRM_RANGE_PRM, 0.0, PRM_RANGE_UI, 0.5);









////////////////////////////////////////////////////////////////////////////////////////



PRM_Template
SOP_SnowSim::myTemplateList[] = {
	// PUT YOUR CODE HERE
	// You now need to fill this template with your parameter name and their default value
	// EXAMPLE : For the angle parameter this is how you should add into the template
	// PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0),
	// Similarly add all the other parameters in the template format here



	//// Material Properties Tab
	//PRM_Template(PRM_SEPARATOR, 1, 0), // tab separator
	////PRM_Template(PRM_LABEL, 1, 0, &PRM_Name("material_tab", "Material Properties")),

	//PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &separationName, &separationDefault, 0, &positiveRange),
	//PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &compressionName, &compressionDefault, 0, &positiveRange),
	//PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &stretchName, &stretchDefault, 0, &positiveRange),
	//PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &hardeningName, &hardeningDefault, 0, &positiveRange),
	//PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &densityName, &densityDefault, 0, &positiveRange),
	//PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &youngName, &youngDefault, 0, &positiveRange),
	//PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &poissonName, &poissonDefault, 0, &poissonRange),

	//// Simulation Tab
	//PRM_Template(PRM_SEPARATOR, 1, 0),
	////PRM_Template(PRM_LABEL, 1, 0, &PRM_Name("sim_tab", "Simulation")),

	//PRM_Template(PRM_FLT_J, PRM_Template::PRM_EXPORT_MIN, 3, &gravityName, gravityDefault),
	//PRM_Template(PRM_FLT_J, PRM_Template::PRM_EXPORT_MIN, 3, &groundName, groundDefault),

	//PRM_Template(PRM_TOGGLE, 1, &resetName, 0), // reset cache button




	/////////////////////////////////////////////////////////////////////////////////////////////

	PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
SOP_SnowSim::myVariables[] = {
    { "PT",	VAR_PT, 0 },		// The table provides a mapping
    { "NPT",	VAR_NPT, 0 },		// from text string to integer token
    { 0, 0, 0 },
};

bool
SOP_SnowSim::evalVariableValue(fpreal &val, int index, int thread)
{
    // myCurrPoint will be negative when we're not cooking so only try to
    // handle the local variables when we have a valid myCurrPoint index.
    if (myCurrPoint >= 0)
    {
	// Note that "gdp" may be null here, so we do the safe thing
	// and cache values we are interested in.
	switch (index)
	{
	    case VAR_PT:
		val = (fpreal) myCurrPoint;
		return true;
	    case VAR_NPT:
		val = (fpreal) myTotalPoints;
		return true;
	    default:
		/* do nothing */;
	}
    }
    // Not one of our variables, must delegate to the base class.
    return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node *
SOP_SnowSim::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_SnowSim(net, name, op);
}

SOP_SnowSim::SOP_SnowSim(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op)
{
    myCurrPoint = -1;	// To prevent garbage values from being returned
}

SOP_SnowSim::~SOP_SnowSim() {}

unsigned
SOP_SnowSim::disableParms()
{
    return 0;
}

OP_ERROR
SOP_SnowSim::cookMySop(OP_Context &context)
{
	fpreal		 now = context.getTime();

	// PUT YOUR CODE HERE
	// Decare the necessary variables and get always keep getting the current value in the node
	// For example to always get the current angle thats set in the node ,you need to :
	//    float angle;
	//    angle = ANGLE(now)       
    //    NOTE : ANGLE is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
	//LSystem myplant;
	
	MPMSolver MPMSolver(Eigen::Vector3f(2.0f, 2.0f, 2.0f), 0.05f, Eigen::Vector3f(0.0f, 0.0f, 0.0f), 0.00001f,
		0.025f, 0.0075f, 10.f, 400.f, 140000.f, 0.2f);

	/*float separation = PARTICLE_SEP(now);
	float critCompression = CRIT_COMPRESSION(now);
	float critStretch = CRIT_STRETCH(now);
	float hardening = HARDENING(now);
	float initDensity = INIT_DENSITY(now);
	float youngModulus = YOUNG_MODULUS(now);
	float poisson = POISSON(now);

	Eigen::Vector3f gravity = GRAVITY(now);
	Eigen::Vector3f ground_plane = GROUND_PLANE(now);

	int resetCache = RESET_CACHE(now);*/







	//if (const GU_Detail* inputGeometry = inputGeo(0, context)) {
	//	GA_Offset ptoff;
	//	GA_FOR_ALL_PTOFF(inputGeometry, ptoff) {
	//		const UT_Vector3& pos = inputGeometry->getPos3(ptoff);

	//		Eigen::Vector3f position(pos.x(), pos.y(), pos.z());
	//		Eigen::Vector3f velocity(0.0f, -5.0f, 0.0f);  

	//		MPMSolver.addParticle(particle(position, velocity, 1.0f));
	//	}
	//}


	// THIS IS JUST A TEST TO SEE IF I CAN DO SOMETHING TO INCOMMING POINTS ****


    UT_Interrupt	*boss;


	if (error() < UT_ERROR_ABORT)
	{
		UT_Interrupt* boss = UTgetInterrupt();
		
		gdp->clearAndDestroy();

		// Start interrupt monitor
		if (boss->opStart("Building SnowSim"))
		{
			gdp->clearAndDestroy();

			if (!lockInput(0, context)) {
				addWarning(SOP_MESSAGE, "Failed to lock input 0.");
				return error();
			}
			const GU_Detail* input = inputGeo(0, context);
			if (!input) {
				addWarning(SOP_MESSAGE, "Input geometry is null after locking.");
				unlockInput(0);
				return error();
			}


			
			GA_Offset ptoff;
			GA_FOR_ALL_PTOFF(input, ptoff) {
				const UT_Vector3& pos = input->getPos3(ptoff);

				Eigen::Vector3f p(pos.x(), pos.y(), pos.z());
				p.y() += 0.5f;

				GA_Offset new_pt = gdp->appendPoint();
				gdp->setPos3(new_pt, UT_Vector3(p.x(), p.y(), p.z()));
			}
			

			unlockInput(0);
	
		}

		boss->opEnd();  // Always call this after opStart
	}

	return error();
}

