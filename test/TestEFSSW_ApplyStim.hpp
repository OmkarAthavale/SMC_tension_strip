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
  void TestRestartingEFS() throw(Exception)
  {

    // -------------- OPTIONS ----------------- //
    std::string mesh_ident = "EFS_problem_0-5_0-025";
    std::string chkpt_dir = mesh_ident + "-BaselineCheckpoint";
    double added_duration = 10000.0;      // ms
    std::string freq_file = "/home/chaste/input1.txt";                    // Hz
    std::string output_dir = chkpt_dir + "_EFS";
    // ---------------------------------------- //

    // Read input file
    std::ifstream in(freq_file.c_str());
    std::string line;
    std::getline(in,line);  
    double freq = std::stod(line);

    BidomainProblemNeural<PROBLEM_SPACE_DIM>* p_bidomain_problem = CardiacSimulationArchiverNeural< BidomainProblemNeural<PROBLEM_SPACE_DIM> >::Load(chkpt_dir + "/checkpoint_problem");

    std::set<unsigned> ICC_id;
    ICC_id.insert(1);
    std::set<unsigned> iccNodes;

    AbstractTetrahedralMesh<PROBLEM_ELEMENT_DIM,PROBLEM_SPACE_DIM>& mesh = p_bidomain_problem->rGetMesh();
    AbstractCardiacTissue<PROBLEM_ELEMENT_DIM,PROBLEM_SPACE_DIM>* tissue = p_bidomain_problem->GetTissue();
    double eleIdentify = 0;
    for (AbstractTetrahedralMesh<PROBLEM_ELEMENT_DIM,PROBLEM_SPACE_DIM>::ElementIterator iter = mesh.GetElementIteratorBegin(); iter != mesh.GetElementIteratorEnd(); ++iter)
    {
      eleIdentify = iter->GetAttribute();
      if (eleIdentify == 1) // ICC=1 and Bath=0
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

    double ex_val = freq;
    double in_val = freq;

    TRACE("Excitatory (Hz): " << ex_val);
    TRACE("Inhibitory (Hz): " << in_val);
	int iterate=0;
	
    for (std::set<unsigned>::iterator it = iccNodes.begin(); it != iccNodes.end(); ++it){
      tissue->GetCardiacCell(*it)->SetParameter("excitatory_neural", ex_val);
      tissue->GetCardiacCell(*it)->SetParameter("inhibitory_neural", in_val);
    }
	
    HeartConfig::Instance()->SetSimulationDuration(p_bidomain_problem->GetCurrentTime() + added_duration); //ms
    HeartConfig::Instance()->SetOutputDirectory(output_dir);

    p_bidomain_problem->Solve();

    // Print summary to terminal
    HeartEventHandler::Headings();
    HeartEventHandler::Report();

    delete p_bidomain_problem;

  };

};

#endif /*TESTMINIMAL_HPP_*/
