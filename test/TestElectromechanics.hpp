#include <cxxtest/TestSuite.h>
#include "../src/ICCFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "NonlinearElasticityTools.hpp"
#include "FileComparison.hpp"
#include "FileFinder.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "Debug.hpp"

#define PROBLEM_SPACE_DIM 2
class TestTensionGenerationStrip : public CxxTest::TestSuite
{
public:
    void TestTensionStrip() throw(Exception)
    {        
        TetrahedralMesh<PROBLEM_SPACE_DIM,PROBLEM_SPACE_DIM> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.025/*stepsize*/, 0.2/*length*/, 0.9/*width*/, 0.1/*depth*/);
        QuadraticMesh<PROBLEM_SPACE_DIM> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.05/*stepsize*/, 0.2/*length*/, 0.9/*width*/, 0.1/*depth*/);

        std::set<unsigned> iccNodes;
        for (unsigned i=0; i < electrics_mesh.GetNumAllNodes() ; ++i) iccNodes.insert(i);
        ICCFactory<PROBLEM_SPACE_DIM> cell_factory(iccNodes);

        // Print mesh summary
        TRACE("Electrics nodes: " << electrics_mesh.GetNumAllNodes());
        TRACE("Mechanics nodes: " << mechanics_mesh.GetNumAllNodes());

        std::vector<unsigned> fixed_nodes_LHS
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.0); // all the X=0.0 nodes
        std::vector<unsigned> fixed_nodes_RHS
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.9); // all the X=0.0 nodes

        // set fixed nodes

        int n = sizeof(fixed_nodes_LHS) / sizeof(fixed_nodes_LHS[0]);
        std::vector<unsigned> fixed_nodes;
        std::vector<unsigned>::iterator it, st;

        std::sort(fixed_nodes_LHS.begin(), fixed_nodes_LHS.end());
        std::sort(fixed_nodes_RHS.begin(), fixed_nodes_RHS.end());
        it = std::set_union(fixed_nodes_LHS.begin(), fixed_nodes_LHS.end(),fixed_nodes_RHS.begin(), fixed_nodes_RHS.end(), fixed_nodes.begin());

        for (st = fixed_nodes.begin(); st != it; ++st)
            std::cout << ' ' << *st;
        std::cout << '\n';
        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,0.01/*contraction model ODE timestep*/);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        HeartConfig::Instance()->SetSimulationDuration(10000.0);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12)); // these are quite smaller than cm values
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.2, 0.2)); // these are quite smaller than cm values
    
        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestSMCTensionStrip");

        problem.Solve();
        // NashHunterPoleZeroLaw<2> law(); // random (non-cardiac) material law
        // problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        // problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, false /*deformation affects cell models*/);
        // problem.SetNoElectricsOutput();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};