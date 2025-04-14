
#if 0

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
#include "LSYSTEMPlugin.h"
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
	    new OP_Operator("CusLsystem",			// Internal name
			    "MyLsystem",			// UI name
			     SOP_Lsystem::myConstructor,	// How to build the SOP
			     SOP_Lsystem::myTemplateList,	// My parameters
			     0,				// Min # of sources
			     0,				// Max # of sources
			     SOP_Lsystem::myVariables,	// Local variables
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//PUT YOUR CODE HERE
//You need to declare your parameters here
//Example to declare a variable for angle you can do like this :
//static PRM_Name		angleName("angle", "Angle");

// Material Properties Tab
static PRM_Name separationName("particle_sep", "Particle Separation");
static PRM_Name compressionName("crit_compression", "Critical Compression");
static PRM_Name stretchName("crit_stretch", "Critical Stretch");
static PRM_Name hardeningName("hardening", "Hardening Coefficient");
static PRM_Name densityName("init_density", "Initial Density");
static PRM_Name youngName("young_modulus", "Initial Young's Modulus");
static PRM_Name poissonName("poisson", "Poisson's Ratio");

// Simulation Tab
static PRM_Name gravityName("gravity", "Gravity Force");
static PRM_Name groundName("ground_plane", "Ground Plane");
static PRM_Name resetName("reset_cache", "Reset Cache");








//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version


// PUT YOUR CODE HERE
// You need to setup the initial/default values for your parameters here
// For example : If you are declaring the inital value for the angle parameter
// static PRM_Default angleDefault(30.0);	



// DEFAULT VALUES
static PRM_Default separationDefault(0.1);
static PRM_Default compressionDefault(2.5);
static PRM_Default stretchDefault(1.5);
static PRM_Default hardeningDefault(10.0);
static PRM_Default densityDefault(1000.0);
static PRM_Default youngDefault(1e5);
static PRM_Default poissonDefault(0.2);

// Default for vector parameters
static PRM_Default gravityDefault[] = {
    PRM_Default(0.0), PRM_Default(-9.8), PRM_Default(0.0)
};
static PRM_Default groundDefault[] = {
    PRM_Default(0.0), PRM_Default(0.0), PRM_Default(0.0)
};


// RANGES
static PRM_Range positiveRange(PRM_RANGE_PRM, 0.0, PRM_RANGE_UI, 10000.0);
static PRM_Range poissonRange(PRM_RANGE_PRM, 0.0, PRM_RANGE_UI, 0.5);






////////////////////////////////////////////////////////////////////////////////////////

PRM_Template
SOP_Lsystem::myTemplateList[] = {
// PUT YOUR CODE HERE
// You now need to fill this template with your parameter name and their default value
// EXAMPLE : For the angle parameter this is how you should add into the template
// PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0),
// Similarly add all the other parameters in the template format here

// Material Properties Tab
PRM_Template(PRM_SEPARATOR, 1, 0), // tab separator
PRM_Template(PRM_LABEL, 1, 0, &PRM_Name("material_tab", "Material Properties")),

PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &separationName, &separationDefault, 0, &positiveRange),
PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &compressionName, &compressionDefault, 0, &positiveRange),
PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &stretchName, &stretchDefault, 0, &positiveRange),
PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &hardeningName, &hardeningDefault, 0, &positiveRange),
PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &densityName, &densityDefault, 0, &positiveRange),
PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &youngName, &youngDefault, 0, &positiveRange),
PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &poissonName, &poissonDefault, 0, &poissonRange),

// Simulation Tab
PRM_Template(PRM_SEPARATOR, 1, 0),
PRM_Template(PRM_LABEL, 1, 0, &PRM_Name("sim_tab", "Simulation")),

PRM_Template(PRM_FLT_J, PRM_Template::PRM_EXPORT_MIN, 3, &gravityName, gravityDefault),
PRM_Template(PRM_FLT_J, PRM_Template::PRM_EXPORT_MIN, 3, &groundName, groundDefault),

PRM_Template(PRM_TOGGLE, 1, &resetName, 0), // reset cache button




/////////////////////////////////////////////////////////////////////////////////////////////

    PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
SOP_Lsystem::myVariables[] = {
    { "PT",	VAR_PT, 0 },		// The table provides a mapping
    { "NPT",	VAR_NPT, 0 },		// from text string to integer token
    { 0, 0, 0 },
};

bool
SOP_Lsystem::evalVariableValue(fpreal &val, int index, int thread)
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
SOP_Lsystem::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_Lsystem(net, name, op);
}

SOP_Lsystem::SOP_Lsystem(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op)
{
    myCurrPoint = -1;	// To prevent garbage values from being returned
}

SOP_Lsystem::~SOP_Lsystem() {}

unsigned
SOP_Lsystem::disableParms()
{
    return 0;
}

OP_ERROR
SOP_Lsystem::cookMySop(OP_Context &context)
{
	fpreal		 now = context.getTime();

	// PUT YOUR CODE HERE
	// Decare the necessary variables and get always keep getting the current value in the node
	// For example to always get the current angle thats set in the node ,you need to :
	//    float angle;
	//    angle = ANGLE(now)       
    //    NOTE : ANGLE is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
	LSystem myplant;

	float separation = PARTICLE_SEP(now);
	float critCompression = CRIT_COMPRESSION(now);
	float critStretch = CRIT_STRETCH(now);
	float hardening = HARDENING(now);
	float initDensity = INIT_DENSITY(now);
	float youngModulus = YOUNG_MODULUS(now);
	float poisson = POISSON(now);
	
	UT_Vector3 gravity;
	GRAVITY(gravity, now);
	
	UT_Vector3 groundPlane;
	GROUND_PLANE(groundPlane, now);
	
	int resetCache = RESET_CACHE(now);










	///////////////////////////////////////////////////////////////////////////

	//PUT YOUR CODE HERE
	// Next you need to call your Lystem cpp functions 
	// Below is an example , you need to call the same functions based on the variables you declare
    // myplant.loadProgramFromString("F\nF->F[+F]F[-F]";  
    // myplant.setDefaultAngle(30.0f);
    // myplant.setDefaultStep(1.0f);





	///////////////////////////////////////////////////////////////////////////////

	// PUT YOUR CODE HERE
	// You the need call the below function for all the genrations ,so that the end points points will be
	// stored in the branches vector , you need to declare them first

	//for (int i = 0; i < generations ; i++)
	//{
	//	  myplant.process(i, branches);
	//}





	///////////////////////////////////////////////////////////////////////////////////


	// Now that you have all the branches ,which is the start and end point of each point ,its time to render 
	// these branches into Houdini 
    

	// PUT YOUR CODE HERE
	// Declare all the necessary variables for drawing cylinders for each branch 
    float		 rad, tx, ty, tz;
    int			 divisions, plane;
    int			 xcoord =0, ycoord = 1, zcoord =2;
    float		 tmp;
    UT_Vector4		 pos;
    GU_PrimPoly		*poly;
    int			 i;
    UT_Interrupt	*boss;

    // Since we don't have inputs, we don't need to lock them.

    divisions  = 5;	// We need twice our divisions of points
    myTotalPoints = divisions;		// Set the NPT local variable value
    myCurrPoint   = 0;			// Initialize the PT local variable



    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT)
    {
	boss = UTgetInterrupt();
	if (divisions < 4)
	{
	    // With the range restriction we have on the divisions, this
	    //	is actually impossible, but it shows how to add an error
	    //	message or warning to the SOP.
	    addWarning(SOP_MESSAGE, "Invalid divisions");
	    divisions = 4;
	}
	gdp->clearAndDestroy();

	// Start the interrupt server
	if (boss->opStart("Building LSYSTEM"))
	{
        // PUT YOUR CODE HERE
	    // Build a polygon
	    // You need to build your cylinders inside Houdini from here
		// TIPS:
		// Use GU_PrimPoly poly = GU_PrimPoly::build(see what values it can take)
		// Also use GA_Offset ptoff = poly->getPointOffset()
		// and gdp->setPos3(ptoff,YOUR_POSITION_VECTOR) to build geometry.






























		////////////////////////////////////////////////////////////////////////////////////////////

	    // Highlight the star which we have just generated.  This routine
	    // call clears any currently highlighted geometry, and then it
	    // highlights every primitive for this SOP. 
	    select(GU_SPrimitive);
	}

	// Tell the interrupt server that we've completed. Must do this
	// regardless of what opStart() returns.
	boss->opEnd();
    }

    myCurrPoint = -1;
    return error();
}

#endif