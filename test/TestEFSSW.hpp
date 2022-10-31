#ifndef TESTMINIMAL_HPP_
#define TESTMINIMAL_HPP_

#define PROBLEM_SPACE_DIM 2
#define PROBLEM_ELEMENT_DIM 2
/**
 * @file
 * This test runs a minimal simulation with Simplified Imtiaz cells in a 2D mesh
 */

#include <cxxtest/TestSuite.h>
/* Most Chaste code uses PETSc to solve linear algebra problems.  This involves starting PETSc at the beginning of a test-suite
 * and closing it at the end.  (If you never run code in parallel then it is safe to replace PetscSetupAndFinalize.hpp with FakePetscSetup.hpp)
 */

#include "Debug.hpp"

#include "ChasteEllipsoid.hpp"
#include "ChastePoint.hpp"

#include "../src/DummyDerivedCa.hpp"
#include "../src/ICCSMC.hpp"

#include "AbstractCardiacCellFactory.hpp"
#include "../src/BidomainProblemNeural.hpp"
#include "../src/ICCFactory.hpp"

#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include "../src/CardiacSimulationArchiverNeural.hpp"


#include "PetscSetupAndFinalize.hpp"

class TestEFS : public CxxTest::TestSuite
{
  private:
  double Beta_Baker2018(double f_EFS)
  {
      double a = 0.00506829502520040;
      double b = 0.0208263263233946;
      double c = -0.00489642388727810;
      double d = -0.289000341617805;

      double toHz = 1.0;
      double endings_per_ICC = 1.0;

      double eval_val = a*exp(b * f_EFS * toHz * endings_per_ICC) + c*exp(d * f_EFS * toHz * endings_per_ICC);

      if (eval_val > 0.007)
      {
          return 0.007;
      } else if (eval_val < 0.0001)
      {
          return 0.0001;
      } else
      {
          return eval_val;
      }
  }

  double GBKmax_Kim2003(double f_EFS)
  {
      double a = 0.147525397773565;
      double b = -0.175087725001323;
      double c = 1.00016775704983;
      double d = 0.0285035244071441;

      double toHz = 1.0;
      double endings_per_ICC = 1.0;

      double eval_val = a*exp(b * f_EFS * toHz * endings_per_ICC) + c*exp(d * f_EFS * toHz * endings_per_ICC);

      if (eval_val > 2.5)
      {
          return 2.5;
      } else if (eval_val < 1.15)
      {
          return 1.15;
      } else
      {
          return eval_val;
      }
  }

  public:
  void TestBaseline() throw(Exception)
  {

    // -------------- OPTIONS ----------------- //
    std::string mesh_ident = "EFS_problem_0-5_0-025";
    std::string output_dir = mesh_ident + "-BaselineCheckpoint";
    unsigned bath_attr = 0;
    unsigned icc_attr = 1;
    double duration = 10.0;      // ms
    double print_step = 5.0;        // ms
    // ---------------------------------------- //

    // Mesh location
    std::string mesh_dir = "projects/mesh/EFS_problem/" + mesh_ident;
    TrianglesMeshReader<PROBLEM_ELEMENT_DIM,PROBLEM_SPACE_DIM> mesh_reader(mesh_dir.c_str());

    // Initialise mesh variables
    std::set<unsigned> iccNodes;
    unsigned nElements = 0;
    DistributedTetrahedralMesh<PROBLEM_ELEMENT_DIM,PROBLEM_SPACE_DIM> mesh;

    // Cell labels
    std::set<unsigned> ICC_id;
    ICC_id.insert(icc_attr);
    std::set<unsigned> bath_id;
    bath_id.insert(bath_attr);

    // Construct ICC mesh network from mesh file
    mesh.ConstructFromMeshReader(mesh_reader);
    nElements = mesh.GetNumLocalElements();

    // Define boundary nodes as bath
    double eleIdentify = 0;
    for (DistributedTetrahedralMesh<PROBLEM_ELEMENT_DIM,PROBLEM_SPACE_DIM>::ElementIterator iter = mesh.GetElementIteratorBegin(); iter != mesh.GetElementIteratorEnd(); ++iter)
    {
      eleIdentify = iter->GetAttribute();
      if (eleIdentify == icc_attr) // ICC=1 and Bath=0
      {
        if(!iter->GetNode(0)->IsBoundaryNode())
        {
          iccNodes.insert(iter->GetNodeGlobalIndex(0));
        }
        if(!iter->GetNode(1)->IsBoundaryNode())
        {
          iccNodes.insert(iter->GetNodeGlobalIndex(1));
        }
        // 2D has two nodes per element for line elements?
        if(!iter->GetNode(2)->IsBoundaryNode())
        {
          iccNodes.insert(iter->GetNodeGlobalIndex(2));
        }
      }
    }

    // Print mesh summary
    TRACE("Number of elements: " << nElements);
    TRACE("Number of ICC nodes: " << iccNodes.size());
    TRACE("Total number of nodes: " << mesh.GetNumAllNodes());

    // Initialise problem with cells
    ICCFactory<PROBLEM_SPACE_DIM> network_cells(iccNodes);
    BidomainProblemNeural<PROBLEM_SPACE_DIM> bidomain_problem(&network_cells, true);
    bidomain_problem.SetMesh( &mesh );

    // Modify simulation config
    HeartConfig::Instance()->Reset();
    HeartConfig::Instance()->SetSimulationDuration(duration);
    HeartConfig::Instance()->SetOutputDirectory(output_dir.c_str());
    HeartConfig::Instance()->SetOutputFilenamePrefix("results");
    HeartConfig::Instance()->SetTissueAndBathIdentifiers(ICC_id, bath_id);
    HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12)); // these are quite smaller than cm values
    HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.2, 0.2)); // these are quite smaller than cm values
    HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
    HeartConfig::Instance()->SetCapacitance(2.5);
    HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
    HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, print_step);

    // Update problem from config
    bidomain_problem.SetWriteInfo();
    bidomain_problem.Initialise();    // resets initial conditions and time to 0.0 ms

    TRACE("Starting Solve");
    // Solve problem
    bidomain_problem.Solve();

    CardiacSimulationArchiverNeural< BidomainProblemNeural<PROBLEM_SPACE_DIM> >::Save(bidomain_problem, output_dir + "/checkpoint_problem");

    // Print summary to terminal
    HeartEventHandler::Headings();
    HeartEventHandler::Report();
  };

};

#endif /*TESTMINIMAL_HPP_*/
