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

class TestCardiacElectroMechanicsTutorial : public CxxTest::TestSuite
{
public:
    void TestCardiacElectroMechanicsExample()
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000);

        HeartConfig::Instance()->SetSimulationDuration(40.0);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     0.1,  // width of square (cm)
                                                     5,    // Number mechanics elements in each direction
                                                     10,   // Number electrics elements in each direction
                                                     &cell_factory,
                                                     KERCHOFFS2003,  // The contraction model (see below)
                                                     1.0,  // mechanics solve timestep
                                                     0.01, // contraction model ode timestep
                                                     "TestCardiacElectroMechanicsExample" /* output directory */);
        problem.Solve();

    }

    void TestCardiacElectroMechanicsExampleAgain()
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

    void dontTestTwistingCube()
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-1000*1000);

        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01/*stepsize*/, 0.1/*length*/, 0.1/*width*/, 0.1/*depth*/);

        QuadraticMesh<3> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1, 0.1 /*as above with a different stepsize*/);

        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 2, 0.0);

        HeartConfig::Instance()->SetSimulationDuration(50.0);

        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        OutputFileHandler handler("TutorialFibreFiles");
        out_stream p_file = handler.OpenOutputFile("5by5by5_fibres.ortho");

        *p_file << mechanics_mesh.GetNumElements() << "\n"; // first line is number of entries
        for (unsigned i=0; i<mechanics_mesh.GetNumElements(); i++)
        {
            double X = mechanics_mesh.GetElement(i)->CalculateCentroid()(0);
            double theta = M_PI/3 - 10*X*2*M_PI/3; // 60 degrees when X=0, -60 when X=0.1;
            *p_file <<  "0 " << cos(theta)  << " " << sin(theta)  // first three entries are fibre direction
                    << " 0 " << -sin(theta) << " " << cos(theta)  // next three are sheet direction
                    << " 1 0 0\n";                                // then normal to sheet direction
        }
        p_file->close();
        out_stream p_file2 = handler.OpenOutputFile("5by5by5_fibres.orthoquad");

        // Mechanics deformation solvers use 3rd order quadrature rules
        GaussianQuadratureRule<3> quad_rule(3);
        QuadraturePointsGroup<3> quad_points(mechanics_mesh, quad_rule);

        *p_file2 << quad_points.Size() << "\n";
        for (unsigned i=0; i<quad_points.Size(); i++)
        {
            double X = quad_points.rGet(i)(0);
            double theta = M_PI/3 - 10*X*2*M_PI/3;
            *p_file2 <<  "0 " << cos(theta)  << " " << sin(theta)
                     << " 0 " << -sin(theta) << " " << cos(theta)
                     << " 1 0 0\n";
        }
        p_file2->close();

        FileFinder finder = handler.FindFile("5by5by5_fibres.orthoquad");
        problem_defn.SetVariableFibreSheetDirectionsFile(finder, true);

        CardiacElectroMechanicsProblem<3,1> problem(COMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestCardiacElectroMech3dTwistingCube");

        problem.Solve();
    }
};