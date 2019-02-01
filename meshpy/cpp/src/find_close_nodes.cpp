

#include "find_close_nodes.H"

namespace py = pybind11;


std::tuple<py::array_t<int>, int> find_close_nodes(
    py::array_t<double> coordinates_numpy, double eps)
{
  // Get the data from the numpy array.
  auto coordinates = coordinates_numpy.unchecked<2>();

  // Number of nodes.
  unsigned int n_nodes = coordinates.shape(0);
  unsigned int n_dim = coordinates.shape(1);

  // Return array with index of partner for each node. -1 means the node does not have a partner.
  auto partner_list = py::array_t<int>(n_nodes);
  py::buffer_info partner_list_buffer = partner_list.request();
  int *partner_list_ptr = (int *)partner_list_buffer.ptr;
  for (unsigned int i = 0; i < n_nodes; i++) partner_list_ptr[i] = -1;

  // Loop over nodes and check if it matches with each other node.
  unsigned int partner = 0;
  bool this_is_partner;
  double distance;
  for (unsigned int i = 0; i < n_nodes; i++)
  {
    // Only check this node if it does not already have partners.
    if (partner_list_ptr[i] == -1)
    {
      this_is_partner = false;

      // Check with all nodes after this one.
      for (unsigned int j = i + 1; j < n_nodes; j++)
      {
        // Calculate the distance between the two nodes.
        distance = 0.;
        for (unsigned int k = 0; k < n_dim; k++)
          distance += pow(coordinates(i, k) - coordinates(j, k), 2);
        distance = sqrt(distance);

        // Check if nodes are within a tolerance.
        if (distance < eps)
        {
          // If the nodes are partner, set the flag for node i.
          this_is_partner = true;

          // The case where node j node already has a partner is not yet implemented.
          if (partner_list_ptr[j] != -1)
            throw std::runtime_error(
                "The case where a node connects two other nodes that are more than eps apart is "
                "not yet implemented.");

          // Set the current partner index for node j.
          partner_list_ptr[j] = partner;
        }
      }


      // Set the partner values for node i too.
      if (this_is_partner)
      {
        partner_list_ptr[i] = partner;
        partner++;
      }
    }
  }

  // Return the partner list, as well as the number of partners.
  return std::make_tuple(partner_list, partner);
}
