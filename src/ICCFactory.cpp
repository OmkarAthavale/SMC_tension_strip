#include "ICCFactory.hpp"

template<unsigned DIM>
AbstractCardiacCell* ICCFactory<DIM>::CreateCardiacCellForTissueNode(Node<DIM>* pNode)
{
  unsigned index = pNode->GetIndex();
  double y = pNode->GetPoint()[1];
  if(setICCNode.find(index) != setICCNode.end())
  {
    CellICCSMCFromCellML* cell = new CellICCSMCFromCellML(this->mpSolver, this->mpZeroStimulus);
    
    cell->SetParameter("E_K", -70.0-4.0*y);

    return cell;
  }

  return new DummyDerivedCa(this->mpSolver, this->mpZeroStimulus);

}

// Explicit instantiation
template class ICCFactory<1>;
template class ICCFactory<2>;
template class ICCFactory<3>;