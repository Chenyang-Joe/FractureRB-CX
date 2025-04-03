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
 * @file    mesh.hpp
 * @ingroup mesh
 * @author  Peter, Rf
 * @date    created:     14.12.09
 *          last change: 14.12.09
 */
#ifndef mesh_tpl
#define mesh_tpl


template<size_t N>
class EdgeKey
{

public:
  unsigned int ids[N];

  EdgeKey(const unsigned int (ids_)[N])
  {
    for(unsigned int cnt =0;cnt<N;++cnt)
          ids[cnt] = ids_[cnt];
    // TODO faster sort
    std::sort(ids, ids+N);

  }

  EdgeKey(const EdgeKey & edge_key_) {
    for(unsigned int cnt =0;cnt<N;++cnt)
      this->ids[cnt] = edge_key_.ids[cnt];
  }

  bool operator<(const EdgeKey& key) const {
    for(unsigned int cnt=0;cnt<N;++cnt) {
      if (ids[cnt] > key.ids[cnt])
        break;
      if (ids[cnt] < key.ids[cnt])
        return true;
    }
    return false;
  }

};


//cctor
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
Mesh<SHAPE, ORDER>::
Mesh(const Mesh& m)
{
  //std::cout << "Copy Constructor" << std::endl;
  nodes_.reserve(m.nodes_.size());
  const_nodeIter nodeit;
  for (nodeit = m.beginNode(); nodeit != m.endNode(); ++nodeit){
    node_type* copied_node = new node_type(**nodeit);
    nodes_.push_back(copied_node);
  }

  elements_.reserve(m.elements_.size());
  const_elementIter elementit;
  for (elementit = m.beginElement(); elementit != m.endElement(); ++elementit){
    const unsigned int id =  (*elementit)->getId();
    const unsigned int input_id = (*elementit)->getInputId();
    boost::array<node_type*, num_nodes> nodes_ptr;
		for(unsigned int n=0; n<num_nodes; ++n) {
      node_type* node = nodes_[(*elementit)->getNode(n)->getId()] ;
      nodes_ptr[n] = node;
    }
    elements_.push_back( new element_type(id, input_id, nodes_ptr) );
  }
  node_id_map_ = m.node_id_map_;
  element_id_map_ = m.element_id_map_;

  initEdges();
  initSupport();
}


template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
Mesh<SHAPE, ORDER>& Mesh<SHAPE, ORDER>::
operator= (const Mesh& m)
{
  this->clear();
  nodes_.reserve(m.nodes_.size());
  const_nodeIter nodeit;
  for (nodeit = m.beginNode(); nodeit != m.endNode(); ++nodeit){
    node_type* copied_node = new node_type(**nodeit);
    nodes_.push_back(copied_node);
  }

  elements_.reserve(m.elements_.size());
  const_elementIter elementit;
  for (elementit = m.beginElement(); elementit != m.endElement(); ++elementit){
    const unsigned int id =  (*elementit)->getId();
    const unsigned int input_id = (*elementit)->getInputId();
    boost::array<node_type*, num_nodes> nodes_ptr;
    for(unsigned int n=0; n<num_nodes; ++n) {
      node_type* node = nodes_[(*elementit)->getNode(n)->getId()] ;
      nodes_ptr[n] = node;
    }
    elements_.push_back( new element_type(id, input_id, nodes_ptr) );
  }
  node_id_map_ = m.node_id_map_;
  element_id_map_ = m.element_id_map_;

  initEdges();
  initSupport();
  return *this;
}


// init nodes
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
void Mesh<SHAPE, ORDER>::
initNodes(const std::map<unsigned int, std::vector<double> >& inp_nodes)
{
  nodes_.reserve(inp_nodes.size());
  // iterator
  std::map<unsigned int, std::vector<double> >::const_iterator iter
    = inp_nodes.begin();
  const std::map<unsigned int, std::vector<double> >::const_iterator
    iter_end = inp_nodes.end();

  // loop inp_nodes
  unsigned int id;
  for( id=0; iter!=iter_end; ++iter, ++id) {
    HYENA_ASSERT( iter->second.size() == dim );
    // create temporary point
    point_type p;
    // loop over all coords per node
		for(unsigned int c=0; c<dim; ++c)
      p[c] = iter->second[c];
    // insert node in array (heap)
    nodes_.push_back( new node_type(id,iter->first,p));
    node_id_map_[iter->first]=id;
  }
}



// init element
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
void Mesh<SHAPE, ORDER>::
initElements(const std::map<unsigned int, std::vector<unsigned int> >&
             inp_elements)
{
  // iterator
  std::map<unsigned int, std::vector<unsigned int> >::const_iterator iter
    = inp_elements.begin();
  const std::map<unsigned int, std::vector<unsigned int> >::const_iterator
    iter_end = inp_elements.end();
  const_nodeIter nodeit;

  // loop inp_elements
  unsigned int id=0;
  for( ; iter!=iter_end; ++iter, ++id ) {
    HYENA_ASSERT( iter->second.size() == num_nodes );

    // loop over all nodes per element
    boost::array<node_type*, num_nodes > nodes_ptr;
		for(unsigned int n=0; n<num_nodes; ++n) {
		  unsigned int node_input_id = iter->second[n];
		  id_map_iter_type node_it = node_id_map_.find(node_input_id);
		  if (node_it == node_id_map_.end())
		    HYENA_ERROR_MSG("Mesh:initElements node id not found");
      unsigned int node_id = node_it->second;
      nodes_ptr[n] = nodes_.at(node_id);
    }

    // insert element in array (heap)
    elements_.push_back(new element_type(id, iter->first, nodes_ptr) );
    element_id_map_[iter->first]=id;
  }
}


// init edges
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
void Mesh<SHAPE, ORDER>::initEdges()
{

  // stores the edge key
  typedef EdgeKey<shape_dim> edge_key_type;
  // stores element id and local edge id
  typedef std::pair<unsigned int, unsigned int> edge_el_type;
  typedef typename std::map<edge_key_type,
    std::vector<edge_el_type > > edge_map_type;

  edge_map_type edge_map;


  // iterate over all elements
  const_elementIter eIt = elements_.begin();
  for(;eIt!=elements_.end();++eIt) {
    // iterate over all edges of element
    for (unsigned int edge=0; edge<num_edges;++edge) {
      // got local node ids on edge


      // get node ids and store in array
      unsigned int node_ids[shape_dim];
      for (unsigned int cnt=0;cnt<shape_dim;++cnt) {
        node_ids[cnt]=(*eIt)->getNode((edge+cnt)%num_topo_nodes)->getId();
      }
      // create edge key with node ids
      edge_key_type key(node_ids);
      unsigned int el_id = (*eIt)->getId();
      // append element and local edge id to vector stored at key
      edge_map[key].push_back(edge_el_type(el_id,edge));
    }
  }

  edges_.reserve(edge_map.size());

  // iterate over all elements of edge map
  typename edge_map_type::const_iterator it = edge_map.begin();
  // init edge id counter
  unsigned int edge_id=0;
  for(;it!=edge_map.end();++it,++edge_id) {
    // temporary arrays for element ptr and local edge id
    typename edge_type::element_array_type curr_elements;
    typename edge_type::local_id_array_type curr_local_edge_id;

    HYENA_ASSERT(it->second.size() == edge_type::sides);
    // iterate over both sides
    for (unsigned int cnt=0;cnt < it->second.size();++cnt) {
      edge_el_type edge_element = it->second[cnt];
      // get and store element ptr
      curr_elements.push_back(elements_[edge_element.first]);
      // store local edge id
      curr_local_edge_id.push_back(edge_element.second);
    }


    // create edge and add to vector
    edges_.push_back(new edge_type( edge_id,
                                  curr_elements,
                                  curr_local_edge_id));

    // set edge at element
    for (unsigned int cnt=0;cnt < it->second.size();++cnt) {
      //add edge id to element
      unsigned int element_id = it->second[cnt].first;
      unsigned int local_edge_id = it->second[cnt].second;
      elements_.at(element_id)->setEdge(edges_[edge_id],local_edge_id);
    }

  }

}

// find elements with have the same member point.
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
void Mesh<SHAPE, ORDER>::
initSupport()
{

  // determine unordered supports
  typedef std::vector<element_type*> element_vector_type;
  std::map<unsigned int, element_vector_type> supports;
  const_elementIter eIter = elements_.begin();
  for (; eIter!=elements_.end(); ++eIter)
    for (unsigned int node=0; node<num_nodes; ++node)
      supports[(*eIter)->getNode(node)->getId()].push_back(*eIter);

  /**
   * Reorder all neighboring elements of each corner node in positive
   * order (right hand rule), this is needed for the computation of the C_ij
   * matrix on non smooth parts of the boundary
   *
   * \code
   *
   *                   <---.
   *              _______   \
   *             /\     /\   |
   *            /  \ 0 /  \  |
   *           / 10 \ /  3 \
   *          |----- * -----|
   *           \ 12 / \  7 /
   *            \  / 8 \  /
   *             \/_____\/
   *
   * \endcode
   *
   * The picure above shows one node (*) with all its neighboring elements and
   * their respective ids. It does not matter where to start with the first
   * element, but only to list all its neighbors in positive order, f.e here:
   * 3, 0, 10, 12, 8, 7.
   **/
  // loop over all nodes (for which the support is to be determined in
  // reordered position)
  typename std::map<unsigned int, element_vector_type >::
    const_iterator sIter = supports.begin();

  for( ; sIter!=supports.end(); ++sIter ) {
    // number of supports for this node
    unsigned int num_support = sIter->second.size();

    // dummy vector to collect ordered support
    element_vector_type dummy;
    dummy.resize(num_support);

    // flag if neighbor element is found
    bool found;

    // only the orientation, not the first neighbor matters, so leave the
    // first neighbor of each node unchanged and start from the second
    dummy[0] = sIter->second[0];

    for(unsigned int m=1; m<num_support; ++m) {
      // the support of the actual element is not found yet
      found = false;
      // the actual element is the previous one in the reordered support list
      // (dummy)
      element_type* actual_elem = dummy[m-1];

      // counting forward from the actual node, take the last node on the
      // current element - this node is only shared with the neighbor element
      // in poitive orientation
      unsigned int last_node(-1);

      for(unsigned int count=0; count<num_topo_nodes; ++count)
        if( actual_elem->getNode(count)->getId() == sIter->first )
          last_node = actual_elem->
            getNode((count+num_topo_nodes-1)%num_topo_nodes)->getId();

      // check which of the other supports shares this node, that one is the
      // next. this can only be checked by going trough the UNORDERED
      // support list, doing that, the actual element has to be excluded,
      // because it is already listed in the actual list (supports).
      for(unsigned int rest=1; rest<num_support; ++rest) {
        for(unsigned int counter=0; counter<num_topo_nodes; ++counter) {
          unsigned int curr_node =
              sIter->second[rest]->getNode(counter)->getId();

          if( (curr_node==last_node) && (sIter->second[rest]!=actual_elem) ) {
            if (m>1) {
              // check if the element found is not already in list
              // if so wrong orientation of element
              if (dummy[m-2] == sIter->second[rest]) {
                std::string id1 = boost::lexical_cast<std::string>(dummy[m-2]->getId());
                std::string id2 = boost::lexical_cast<std::string>(dummy[m-1]->getId());
                HYENA_ERROR_MSG("Wrong orientation of elements "+id1+" "+id2);
              }
            }
            dummy[m] = sIter->second[rest];
            found = true;
          }
          // if the element which shares the same geometry node is found the
          // loop over all gdofs per element does need to be finished
          if(found) break;
        }
        // if the neighbouring element is found the loop over the remaining
        // elements does not need to be finished either
        if(found) break;
      }
      // if a neighboring element cannot be found throw out all remaining
      // elements (rest) there is no need to continue with the reordering,
      // i.e. the boundary is probably open
      if(!found) break;
    }

    if (found == false)
      dummy = sIter->second;

    nodes_.at(sIter->first)->initElements(dummy);
  }
}






// transform the geometry
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
const double Mesh<SHAPE, ORDER>::
transformGeometry()
{
	// the maximal radius is to determine
	double max_r = getMaxDistance();
 
  // node iterator
  nodeIter nIter = beginNode();
  const_nodeIter nIend = endNode();
 
  //transform all nodes with max_r
  for( ; nIter!=nIend; ++nIter ) {
    point_type p = (*nIter)->getPoint();
    p /= max_r;
    (*nIter)->setPoint(p);
  }
  return max_r;
}

template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
void Mesh<SHAPE, ORDER>::getMaxEdgeRatio(double& max_ratio,
                                         unsigned int& element_id) const
{
  // initialize 
  double first_edge, second_edge, tmp_ratio;
  max_ratio = 0.0;
  element_id = -1;
  
  // element iterator
  const_elementIter eIter = beginElement();
  const_elementIter eIend = endElement();

  const unsigned int num_nodes = element_type::num_nodes;

  // iterate over all elements
  for(; eIter!=eIend; ++eIter) {
    for(unsigned int m=0; m<num_nodes; ++m) {
      first_edge = 
        ((*eIter)->getNode((m+1)%num_nodes)->getPoint() - 
         (*eIter)->getNode((m  )%num_nodes)->getPoint()).norm(); 
      for(unsigned int n=m+1; n<num_nodes; ++n) {
        second_edge = 
          ((*eIter)->getNode((n+1)%num_nodes)->getPoint() - 
           (*eIter)->getNode((n  )%num_nodes)->getPoint()).norm(); 
        
        // temporary edge ratio
        tmp_ratio = first_edge/second_edge;

        // invert if ratio smaller than one
        if(tmp_ratio<1.0) tmp_ratio = 1.0/tmp_ratio;

        // increase max_ratio and set element_id if temporary ratio is bigger
        if(tmp_ratio>max_ratio) {
          max_ratio = tmp_ratio;
          element_id = (*eIter)->getInputId();
        }
      } //m
    } // n
  } // eIter
}






// get maximal distance in all space directions
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
const double Mesh<SHAPE, ORDER>::getMaxDistance() const
{
	// the maximal radius is to determine
	double max_r = 0.;
  
  // node iterator
  const_nodeIter nIter = beginNode();
  const_nodeIter nIend = endNode();

	switch(dim)
		{
		case 2:
			HYENA_ERROR_MSG("geometry transformation for 2D not implemented");
			break;

		case 3:
			{
        // initialize all min/max coordinates with the first nodal values
        nIter = beginNode();
				Point3 max_point( (*nIter)->getPoint() );
				Point3 min_point( (*nIter)->getPoint() );

				// coordinates of the actual node
				Point3 actual_point;

				// compare with all other nodes to find min/max coordinates
        for( ; nIter!=nIend; ++nIter ) {
					actual_point = (*nIter)->getPoint();

          // find out max coordinates
					if( actual_point[0] > max_point[0] ) max_point[0] = actual_point[0];
					if( actual_point[1] > max_point[1] ) max_point[1] = actual_point[1];
					if( actual_point[2] > max_point[2] ) max_point[2] = actual_point[2];

					// find out min coordinates
					if( actual_point[0] < min_point[0] ) min_point[0] = actual_point[0];
					if( actual_point[1] < min_point[1] ) min_point[1] = actual_point[1];
					if( actual_point[2] < min_point[2] ) min_point[2] = actual_point[2];
				}

				// compute the maximal radius of the cuboid (this can only be the
				// spatial diagonal)
				max_r = max_point.distance(min_point);
				break;
			}

		default: HYENA_ERROR_MSG("unknown space DIMENSION");
		}
  return max_r;
}


//transform the geometry back
template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
void Mesh<SHAPE, ORDER>::
transformGeometryBack(double max_r)
{
  nodeIter nIter = beginNode();
  const_nodeIter nIend = endNode();

  nIter = beginNode();
  for( ; nIter!=nIend; ++nIter ) {
    point_type p = (*nIter)->getPoint();
    p *= max_r;
    (*nIter)->setPoint(p);
  }
}

//////////////////////////////////////////////////////////////////////
//   OUTPUT   ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
std::ostream& operator << (std::ostream& out,
                           const Mesh<SHAPE, ORDER>& m)
{
  typename Mesh<SHAPE, ORDER>::const_nodeIter nodeit;
  for (nodeit = m.beginNode(); nodeit != m.endNode(); nodeit++)
    out << nodeit->first << " " << *nodeit->second << std::endl;
  out << std::endl;

  typename Mesh<SHAPE, ORDER>::const_elementIter elementit;
  for (elementit = m.beginElement(); elementit != m.endElement(); elementit++)
    out << elementit->first << " " <<  *elementit->second << std::endl;

  return out;
}



#endif // include guard
