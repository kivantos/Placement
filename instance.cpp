#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <limits>

#include "instance.h"

using namespace std;

void 
Instance::read_file(const char* filename)
{
   fstream file(filename);             // open file
   if (! file)
   {
         cerr << "Could not open file";
         return;
   }

   string line;
   getline(file, line);                // get first line of file
   stringstream ss(line);              // convert line to a stringstream
   ss >> _num_cells;                         // for which we can use << as with cin
   _cells.clear();

   while (getline(file, line))
   {
         stringstream ss(line);
         Cell c;
         ss >> c.width;
         if (! ss)
         {
               cerr << "Invalid file format: Could not read cell. Width=" << c.width << endl;
               return;
         }
         ss >> c.height;
         _cells.push_back(c);
   }

   if (_cells.size() != _num_cells)
   {
      cerr << "Invalid file format: _num_cells=" << _num_cells << ", vector has "
      << _cells.size() << " elements." << endl;
      return;
   }
   
   _best_placement.resize(_num_cells);
}


void Instance::minimum_perimeter()
{
   vector<size_t> pi;
   vector<size_t> sigma;
   vector<size_t> free_indices;
   vector<pair<x_coord, y_coord> > placement;
   pi.resize(_num_cells);
   sigma.resize(_num_cells);
   placement.resize(_num_cells);
   
   size_t idx = _num_cells;
   
   for (size_t i = 0; i<_num_cells; i++)
   {
      free_indices.push_back(i);
   }
   
   find_perimeter_for_all_permutations(pi, sigma, free_indices, placement, idx, true);
   
}


void Instance::find_perimeter_for_all_permutations(std::vector<size_t> & pi, 
                                            std::vector<size_t> & sigma, 
                                            std::vector<size_t> & free_indices,
                                            vector<pair<x_coord, y_coord> > & placement,
                                            size_t idx,
                                            bool pi_or_sigma)
{
   if (pi_or_sigma)
   {
      for (size_t i = 0; i < idx; i++)
      {
         size_t tmp = free_indices[i];
         pi[tmp] = idx;
         free_indices[i] = free_indices[idx - 1];
         find_perimeter_for_all_permutations(pi, sigma, free_indices, placement, idx-1, true);
         free_indices[i] = tmp;
      }
      if (idx == 0)
      {
         vector<size_t> free_indices_sigma;
         
         for (size_t i = 0; i<_num_cells; i++)
         {
            free_indices_sigma.push_back(i);
         }
         find_perimeter_for_all_permutations(pi, sigma, free_indices_sigma, placement, _num_cells, false);
      }
   }
   else
   {
      for (size_t i = 0; i < idx; i++)
      {
         size_t tmp = free_indices[i];
         sigma[tmp] = idx;
         free_indices[i] = free_indices[idx - 1];
         find_perimeter_for_all_permutations(pi, sigma, free_indices, placement, idx-1, false);
         free_indices[i] = tmp;
      }
      if (idx == 0)
      {
         if (placement_for_pi_sigma(pi, sigma, placement))
         {
            for (size_t i = 0; i<_num_cells; i++)
            {
               _best_placement[i] = placement[i];
            }
         }
//          cout << "Tiefe erreicht: pi             sigma" << endl;
//          for (size_t j = 0; j < pi.size(); j++)
//          {
//             cout << pi[j] << " ";
//          }
//          cout << "               ";
//          for (size_t j = 0; j < pi.size(); j++)
//          {
//             cout << sigma[j] << " ";
//          }
//          cout << endl;
      }
   }
   
}


bool Instance::placement_for_pi_sigma(vector<size_t> const & pi, 
                                      vector<size_t> const & sigma, 
                                      vector<pair<x_coord, y_coord> > & placement)
{
   map<pair<size_t, long int> > Q;
   Q[0] = 0;
   Q[_num_cells + 1] = 100000; //TODO: INF
   
   for (size_t i = 0; i < _num_cells; i++)
   {
      size_t c = pi[i];
      size_t p = sigma[c];
      long int l_p;
      
      map<pair<size_t, long int> >::iterator it;
      it = Q.upper_bound(p+1);
      it--;
      placement[c].first = (*it).second;
      l_p = placement[c].first + _cells[c].width;
      Q[p+1] = l_p;
      
      if (l_p >= _best_perimeter)
      {
         return false;
      }
      
      it = Q.upper_bound(p+1);
      while ((*it).second <= l_p)
      {
         Q.erase(it);
         it = Q.upper_bound(p+1);
      }
   }
   
   for (size_t i = 0; i < _num_cells; i++)
   {
      size_t c = sigma[i];
      size_t p = pi[c];
      long int l_p;
      
      map<pair<size_t, long int> >::iterator it;
      it = Q.upper_bound(p+1);
      it--;
      placement[c].first = (*it).second;
      l_p = placement[c].first + _cells[c].width;
      Q[p+1] = l_p;
      
      if (l_p >= _best_perimeter)
      {
         return false;
      }
      
      it = Q.upper_bound(p+1);
      while ((*it).second <= l_p)
      {
         Q.erase(it);
         it = Q.upper_bound(p+1);
      }
   }
   
}


void
Instance::print()
{
   cout << "INSTANCE: (" << _num_cells << " cells)" << endl;
   for (size_t i=0; i<_cells.size(); i++)
   {
      cout << _cells[i].width << " " << _cells[i].height << endl;
   }
   cout << endl;
}
