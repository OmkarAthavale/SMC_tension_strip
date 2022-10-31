#ifndef ICCFACTORY_HPP_
#define ICCFACTORY_HPP_

#include <set>

#include "AbstractCardiacCell.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "../src/DummyDerivedCa.hpp"
#include "../src/ICCSMC.hpp"

template<unsigned DIM>
class ICCFactory : public AbstractCardiacCellFactory<DIM>
{
  private:
  std::set<unsigned> setICCNode;

  public:
  ICCFactory(std::set<unsigned> iccNodes) : 
  AbstractCardiacCellFactory<DIM>(), 
  setICCNode(iccNodes)
  {};

  // Destructor
  virtual ~ICCFactory(){};

  AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode);
};

#endif