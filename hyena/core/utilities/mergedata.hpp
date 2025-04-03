// Copyright (C) 2009-2010 Matthias Messner, Michael Messner, Franz
// Rammerstorfer, Peter Urthaler
//
// This file is part of HyENA - a C++ boundary element methods library.
//
// HyENA is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// HyENA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for more
// details.
//
// You should have received a copy of the GNU Lesser Public License along with
// HyENA. If not, see <http://www.gnu.org/licenses/>.

/**
 * @file   mergedata.hpp
 * @author TT
 * @geistigervater Rf
 * @date   created:     14.07.11
 *         last change: 14.07.11
 */

#ifndef MERGEDATA_HPP
#define MERGEDATA_HPP


namespace hyena
{

  class MergeData
  {
  public:

    template<typename LDOF_ARRAY,
             typename IN_DATA,
             typename OUT_DATA>
    void operator()(const LDOF_ARRAY& ldofs,
                    const IN_DATA& in_data,
                    OUT_DATA& out_data,
                    int idx_offset=0, int node_base_idx=0) const
    {
		for(unsigned int cnt=0;cnt<ldofs.size();++cnt) {
			if(idx_offset==0){ // default case, use original version
				out_data[ldofs[cnt]->getID()] = in_data[cnt];
			}else{ // using idx_offset for cu_ldofs
				//out_data[ldofs[cnt]->getID()] = in_data[cnt+idx_offset]; // this is not good as we cannot prevent cu_ldofs being out of order wrt. idx so simple counting won't work!
				out_data[ldofs[cnt]->getID()] = in_data[ldofs[cnt]->getIDX()+idx_offset];

				// for debugging ...
				unsigned int node = ldofs[cnt]->getGDof()->getSuperElement(0)->getElement()->getNode( ldofs[cnt]->getGDof()->getReferenceElementIdx(0) )->getInputId();
				unsigned int target = (node-node_base_idx)*3+ldofs[cnt]->getLIDX();
				if(target != ldofs[cnt]->getID()) printf("\n target ID mismatch nd=%d 3(nd-b)+d=%d, id=%d",node,target,ldofs[cnt]->getID()); //indicates that a node is not numbered according to the elements list
				//out_data[target] = in_data[ldofs[cnt]->getIDX()+idx_offset];
				//printf("\n%d %d %d %% %d==%d is %d",ldofs[cnt]->getID(),node, ldofs[cnt]->getIDX()+idx_offset, (node-1)*3+ldofs[cnt]->getLIDX(), ldofs[cnt]->getID(), ((node-1)*3+ldofs[cnt]->getLIDX()) == ldofs[cnt]->getID() );
			}
		}
    }

  };

} // end namespace hyena
//se=gdof->getSuperElement(support_cnt); se->getElement()->getNode((*ldof_it)->getGDof()->getReferenceElementIdx(support_cnt))->getInputId();
#endif
