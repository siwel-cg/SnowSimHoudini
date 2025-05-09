#include <UT/UT_DSOVersion.h>
//#include <RE/RE_EGLServer.h>


#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_VDBUtils.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <GU/GU_PrimVDB.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/GridTransformer.h>

#include <iostream>
#include <limits.h>
#include "MPMPlugin.h"
#include "MPMSolver.h"
#include "grid.h"
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
				 OP_FLAG_GENERATOR)
	    );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Material Properties Tab
static PRM_Name compressionName("crit_compression", "CritCompression");
static PRM_Name stretchName("crit_stretch", "CritStretch");
static PRM_Name hardeningName("hardening", "Hardening");
static PRM_Name densityName("init_density", "Density");
static PRM_Name youngName("young_modulus", "YoungsMod");
static PRM_Name poissonName("poisson", "Poisson");

// Simulation Tab
//static PRM_Name gravityName("gravity", "Gravity Force");
//static PRM_Name timeStepName("dt", "Time Step");
//static PRM_Name groundName("ground_plane", "Ground Plane");
//static PRM_Name boundsSizeName("sim_bounds", "Simulation Bounds");
//static PRM_Name boundsPosName("sim_pos", "Simulation Position");
//static PRM_Name resetName("reset_cache", "Reset Cache");


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//// DEFAULT VALUES
static PRM_Default compressionDefault(0.05);
static PRM_Default stretchDefault(0.005);
static PRM_Default hardeningDefault(10.0);
static PRM_Default densityDefault(600.0);
static PRM_Default youngDefault(180000);
static PRM_Default poissonDefault(0.4);



////////////////////////////////////////////

//static PRM_Default timestepDefault(0.001);
//static PRM_Default groundDefault(-2.5);
////Default for vector parameters
//static PRM_Default boundSizeDefault[] = {
//	PRM_Default(2.5), PRM_Default(2.5), PRM_Default(2.5)
//};
//static PRM_Default boundPosDefault[] = {
//	PRM_Default(0.0), PRM_Default(0.0), PRM_Default(0.0)
//};
//static PRM_Default gravityDefaults[] = {
//    PRM_Default(0.0f), // X
//    PRM_Default(-9.8f), // Y
//    PRM_Default(0.0f)  // Z
//};


// RANGES
static PRM_Range compressionRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 1.0);
static PRM_Range stretchRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 1.0);
static PRM_Range hardeningRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 100.0);
static PRM_Range densityRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 1000.0);
static PRM_Range youngRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 1e6);
static PRM_Range poissonRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 0.499);

//static PRM_Range timeRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 0.01);
//static PRM_Range groundRange(PRM_RANGE_UI, -100.0, PRM_RANGE_RESTRICTED, 100.0);
//static PRM_Range boundsSizeRange(PRM_RANGE_UI, 0.0, PRM_RANGE_RESTRICTED, 100.0);
//static PRM_Range boundsPosRange(PRM_RANGE_UI, -100.0, PRM_RANGE_RESTRICTED, 100.0);
//static PRM_Range gravityRange(PRM_RANGE_UI, -100.0, PRM_RANGE_RESTRICTED, 100.0);



////////////////////////////////////////////////////////////////////////////////////////


PRM_Template
SOP_SnowSim::myTemplateList[] = {


	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &compressionName, &compressionDefault, 0, &compressionRange),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &stretchName, &stretchDefault, 0, &stretchRange),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &hardeningName, &hardeningDefault, 0, &hardeningRange),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &densityName, &densityDefault, 0, &densityRange),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &youngName, &youngDefault, 0, &youngRange),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &poissonName, &poissonDefault, 0, &poissonRange),

	//PRM_Template(PRM_FLT_J,  PRM_Template::PRM_EXPORT_MIN, 3, &gravityName,     gravityDefaults,     0, &gravityRange),
	//PRM_Template(PRM_FLT,    PRM_Template::PRM_EXPORT_MIN, 1, &timeStepName,    &timestepDefault,    0, &timeRange),
	//PRM_Template(PRM_FLT,    PRM_Template::PRM_EXPORT_MIN, 1, &groundName,      &groundDefault,      0, &groundRange),
	//PRM_Template(PRM_FLT_J,  PRM_Template::PRM_EXPORT_MIN, 3, &boundsSizeName,  boundSizeDefault,    0, &boundsSizeRange),
	//PRM_Template(PRM_FLT_J,  PRM_Template::PRM_EXPORT_MIN, 3, &boundsPosName,   boundPosDefault,     0, &boundsPosRange),





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
	prevFrame = -1; 


	// INITILIZE MPM SOLVER
	solver = MPMSolver();
	solver.computeInitialDensity();
}

SOP_SnowSim::~SOP_SnowSim() {}

unsigned
SOP_SnowSim::disableParms()
{
    return 0;
}

bool SOP_SnowSim::isTimeDependent() const
{
	return true;
}


void SOP_SnowSim::readSDFFromVDB(const GU_Detail* sdfGdp)
{
	if (!sdfGdp) {
		UTprintf("No SDF input provided\n");
		return;
	}

	vdbPrimSDF = nullptr;
	for (GA_Iterator it(sdfGdp->getPrimitiveRange()); !it.atEnd(); ++it)
	{
		const GA_Primitive* prim = sdfGdp->getPrimitive(*it);
		if (prim->getTypeId() == GA_PRIMVDB)
		{
			const GU_PrimVDB* vdbPrim = static_cast<const GU_PrimVDB*>(prim);
			const openvdb::GridBase::ConstPtr gridBase = vdbPrim->getConstGridPtr();
			if (gridBase->isType<openvdb::FloatGrid>())
			{
				vdbPrimSDF = vdbPrim;
				break;
			}
		}
	}

	if (!vdbPrimSDF) {
		UTprintf("No valid float VDB found for SDF\n");
	}
}

#if 1
OP_ERROR
SOP_SnowSim::cookMySop(OP_Context& context)
{
	flags().setTimeDep(true);
	fpreal now = context.getTime();
	int frame = context.getFrame();


	// TEST INPUT SOURCE
	lockInputs(context);
	const GU_Detail* inGdp = inputGeo(0, context);
	if (!inGdp) {
		UTprintf("SOP_SnowSim: INVALID INPUT\n");
		unlockInputs();
		return error();
	}

	const GU_Detail* sdfGdp = inputGeo(1, context);
	// Read the SDF VDB if available
	if (sdfGdp) {
		readSDFFromVDB(sdfGdp);
	}

	//float separation = PARTICLE_SEP(now);
	float critCompression = CRIT_COMPRESSION(now);
	float critStretch = CRIT_STRETCH(now);
	float hardening = HARDENING(now);
	float initDensity = INIT_DENSITY(now);
	float youngModulus = YOUNG_MODULUS(now);
	float poisson = POISSON(now);

	//float timeStep = TIME_STEP(now);
	//float groundPlane = GROUND_PLANE(now);
	//Eigen::Vector3f boundsSize = Eigen::Vector3f(BOUNDS_SIZE(now)[0], BOUNDS_SIZE(now)[1], BOUNDS_SIZE(now)[2]);
	//Eigen::Vector3f boundsPos = Eigen::Vector3f(BOUNDS_POS(now)[0], BOUNDS_POS(now)[1], BOUNDS_POS(now)[2]);

	if (frame <= 1) {

		// RESET SOLVER 

		/*MPMSolver(Eigen::Vector3f gridDim, float spacing, Eigen::Vector3f gridOrigin, float dt,
			float critCompression, float critStretch,
			float hardeningCoeff, float initialDensity, float youngsMod,
			float poissonRatio);*/

		solver = MPMSolver(Eigen::Vector3f(2.5, 2.5, 2.5), 0.1, Eigen::Vector3f(0.0f, 0.0f, 0.0f), -2.5, 0.001,
					0.05f, 0.005f, 10.f, 600.f, 180000.f, 0.35, nullptr);

		//solver = MPMSolver(boundsSize, 0.1, boundsPos, groundPlane, timeStep,
		//	critCompression, critStretch, hardening, initDensity, youngModulus, poisson, vdbPrimSDF);

		const GU_Detail* inGdp = inputGeo(0, context);
		if (inGdp) {
			GA_Offset ptoff;
			GA_FOR_ALL_PTOFF(inGdp, ptoff) {
				const UT_Vector3 pos = inGdp->getPos3(ptoff);
				Eigen::Vector3f position(pos.x(), pos.y(), pos.z());
				Eigen::Vector3f velocity(0.0f, 0.0f, 0.0f);
				float mass = 1.0f;

				solver.addParticle(MPMParticle(position, velocity, mass));
			}
		}


		solver.computeInitialDensity();
		prevFrame = 1;
	}


	for (int f = prevFrame + 1; f <= frame; ++f) {
		solver.step();
	}
	prevFrame = frame;

	unlockInputs();

	gdp->clearAndDestroy();

	// INSTATIATE THE POINTS IN SPACE
	for each(MPMParticle p in solver.getParticles())
	{
		Eigen::Vector3f pos = p.position;
		GA_Offset pt = gdp->appendPoint();
		gdp->setPos3(pt, UT_Vector3(pos.x(), pos.y(), pos.z()));
	}


	// GROUND PLANE

	//if (frame == 1) {
	//	float size = 5.0f;

	//	UT_Vector3 p0(-size, groundPlane, -size);
	//	UT_Vector3 p1(size, groundPlane, -size);
	//	UT_Vector3 p2(size, groundPlane, size);
	//	UT_Vector3 p3(-size, groundPlane, size);

	//	GA_Offset pt0 = gdp->appendPointOffset();
	//	GA_Offset pt1 = gdp->appendPointOffset();
	//	GA_Offset pt2 = gdp->appendPointOffset();
	//	GA_Offset pt3 = gdp->appendPointOffset();

	//	gdp->setPos3(pt0, p0);
	//	gdp->setPos3(pt1, p1);
	//	gdp->setPos3(pt2, p2);
	//	gdp->setPos3(pt3, p3);

	//	GEO_PrimPoly* quad = (GEO_PrimPoly*)gdp->appendPrimitive(GA_PRIMPOLY);
	//	quad->appendVertex(pt0);
	//	quad->appendVertex(pt1);
	//	quad->appendVertex(pt2);
	//	quad->appendVertex(pt3);
	//	quad->close();

	//	// Optional: color the plane gray
	//	//GA_RWHandleV3 cd(gdp->addAttrib("Cd", GA_ATTRIB_POINT, UT_Vector3(0.4f, 0.4f, 0.4f)));
	//	//cd.set(pt0, UT_Vector3(0.4f));
	//	//cd.set(pt1, UT_Vector3(0.4f));
	//	//cd.set(pt2, UT_Vector3(0.4f));
	//	//cd.set(pt3, UT_Vector3(0.4f));
	//}
	


	// BOUNDING BOX ( ONLY SHOW AT INITIALIZATION )
	//UT_Vector3 center(boundsPos.x(), boundsPos.y(), boundsPos.z());
	//UT_Vector3 half(boundsSize.x() * 0.5f, boundsSize.y() * 0.5f, boundsSize.z() * 0.5f);

	//if (frame == 1) {
	//	GA_Offset bb[8];
	//	for (int k = 0; k < 2; ++k) {
	//		for (int j = 0; j < 2; ++j) {
	//			for (int i = 0; i < 2; ++i) {
	//				int idx = i + 2 * (j + 2 * k);
	//				bb[idx] = gdp->appendPoint();
	//				UT_Vector3 p(
	//					center.x() + (i ? half.x() : -half.x()),
	//					center.y() + (j ? half.y() : -half.y()),
	//					center.z() + (k ? half.z() : -half.z())
	//				);
	//				gdp->setPos3(bb[idx], p);
	//			}
	//		}
	//	}

	//	auto drawEdge = [&](int a, int b) {
	//		GU_PrimPoly* line = GU_PrimPoly::build(gdp, false);
	//		line->appendVertex(bb[a]);
	//		line->appendVertex(bb[b]);
	//		};

	//	drawEdge(0, 1); drawEdge(1, 3); drawEdge(3, 2); drawEdge(2, 0); // bottom
	//	drawEdge(4, 5); drawEdge(5, 7); drawEdge(7, 6); drawEdge(6, 4); // top
	//	drawEdge(0, 4); drawEdge(1, 5); drawEdge(2, 6); drawEdge(3, 7); // sides
	//}





	myCurrPoint = -1;
	return error();
}
#endif
