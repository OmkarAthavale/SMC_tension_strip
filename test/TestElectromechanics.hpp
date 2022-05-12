#include <cxxtest/TestSuite.h>
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "FileComparison.hpp"
#include "FileFinder.hpp"

class TestTensionGenerationStrip : public CxxTest::TestSuite
{
public:
    void TestTensionStrip()
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000);

        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01/*stepsize*/, 0.1/*length*/, 0.1/*width*/, 0.1/*depth*/);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1, 0.1 /*as above with a different stepsize*/);

        HeartConfig::Instance()->SetSimulationDuration(40.0);

        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.0); // all the X=0.0 nodes

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,0.01/*contraction model ODE timestep*/);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestCardiacElectroMechanicsExample2");

        problem.Solve();
        CompressibleMooneyRivlinMaterialLaw<2> law(2.0,1.0); // random (non-cardiac) material law
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, true /*deformation affects cell models*/);
        problem.SetNoElectricsOutput();

        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[5](0), 0.090464, 1e-4);

        FileFinder finder1("TestCardiacElectroMechanicsExample/deformation/solution_40.nodes", RelativeTo::ChasteTestOutput);
        FileFinder finder2("TestCardiacElectroMechanicsExample2/deformation/solution_40.nodes", RelativeTo::ChasteTestOutput);
        FileComparison comparer(finder1,finder2);
        TS_ASSERT(comparer.CompareFiles());
    }
};