#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <climits>

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
   _total_cell_size = 0;

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
         _total_cell_size += c.width * c.height;
   }

   if (_cells.size() != _num_cells)
   {
      cerr << "Invalid file format: _num_cells=" << _num_cells << ", vector has "
      << _cells.size() << " elements." << endl;
      return;
   }

   best_placement.resize(_num_cells);
   best_perimeter = INT_MAX;
}


void Instance::minimum_perimeter()
{
   vector<size_t>                  pi;                 //Permutation pi
   vector<size_t>                  pi_inverse;         //Inverse of permutation pi
   vector<size_t>                  free_indices_pi;    //Uninitialized indices of pi
   vector<size_t>                  sigma;              //Permutation sigma
   vector<size_t>                  sigma_inverse;      //Inverse of permutation sigma
   vector<size_t>                  free_indices_sigma; //Uninitialized indices of sigma
   vector<pair<x_coord, y_coord> > placement;          //Coordinates of the cells in a placement
   vector<pair<x_coord, y_coord> > best_placement;     //Coordinates of current best placement
   map<size_t, long int>           Q;                  //Priority queue needed for the placement alg.
   long int best_perimeter = INT_MAX;                  //Current best perimeter of placement


   //Reserve memory for all vectors.
   //We only allocate memory once for every vector.
   pi.resize(_num_cells);
   pi_inverse.resize(_num_cells);
   sigma.resize(_num_cells);
   sigma_inverse.resize(_num_cells);
   placement.resize(_num_cells);
   best_placement.resize(_num_cells);


   //At the beginning, all indices of pi and sigma are uninitialized.
   for (size_t i = 0; i<_num_cells; i++)
   {
      free_indices_pi.push_back(i);
      free_indices_sigma.push_back(i);
   }

   find_perimeter_for_all_permutations(pi,
                                       pi_inverse,
                                       free_indices_pi,
                                       sigma,
                                       sigma_inverse,
                                       free_indices_sigma,
                                       placement,
                                       best_placement,
                                       best_perimeter,
                                       Q,
                                       _num_cells,  //Index to be fixed first.
                                       true);

   //Output
   cout << endl << "BBX=" << best_perimeter << endl;
   cout << "Placement:" << endl;
   for (size_t i = 0; i < _num_cells; i++)
   {
      cout << "Circuit " << i << ": (" << best_placement[i].first
      << ", " << best_placement[i].second << ")" << endl;
   }
}


void Instance::find_perimeter_for_all_permutations(vector<size_t> & pi,
                                                   vector<size_t> & pi_inverse,
                                                   vector<size_t> & free_indices_pi,
                                                   vector<size_t> & sigma,
                                                   vector<size_t> & sigma_inverse,
                                                   vector<size_t> & free_indices_sigma,
                                                   vector<pair<x_coord, y_coord> > & placement,
                                                   vector<pair<x_coord, y_coord> > & best_placement,
                                                   long int & best_perimeter,
                                                   map<size_t, long int> & Q,
                                                   size_t idx,
                                                   bool pi_or_sigma)
{
   if (pi_or_sigma)
   //The permutation pi is not completely determined yet.
   //Fix the next value in pi, and then recursively find the next values for the permutation.
   {
      for (size_t i = 0; i < idx; i++)
      {
         size_t tmp = free_indices_pi[i];
         pi[tmp] = idx-1;
         pi_inverse[idx-1] = tmp;
         free_indices_pi[i] = free_indices_pi[idx - 1];

         find_perimeter_for_all_permutations(pi,
                                             pi_inverse,
                                             free_indices_pi,
                                             sigma,
                                             sigma_inverse,
                                             free_indices_sigma,
                                             placement,
                                             best_placement,
                                             best_perimeter,
                                             Q,
                                             idx-1,   //Index to be fixed next.
                                             true);

         free_indices_pi[i] = tmp;
      }
      if (idx == 0)
      //Permutation pi is now completely determined.
      //Now recursively determine sigma.
      //Remark: At this point the vector free_indices_sigma contains all entries
      //from 0 to _num_cells-1. However, they are not ordered.
      {
         find_perimeter_for_all_permutations(pi,
                                             pi_inverse,
                                             free_indices_pi,
                                             sigma,
                                             sigma_inverse,
                                             free_indices_sigma,
                                             placement,
                                             best_placement,
                                             best_perimeter,
                                             Q,
                                             _num_cells,
                                             false);
      }
   }
   else
   //The permutation pi is fully determined. The permutation sigma is not yet
   //fully determined. Fix the next value and recurse on the rest.
   {
      for (size_t i = 0; i < idx; i++)
      {
         size_t tmp = free_indices_sigma[i];
         sigma[tmp] = idx-1;
         sigma_inverse[idx-1] = tmp;
         free_indices_sigma[i] = free_indices_sigma[idx - 1];

         find_perimeter_for_all_permutations(pi,
                                             pi_inverse,
                                             free_indices_pi,
                                             sigma,
                                             sigma_inverse,
                                             free_indices_sigma,
                                             placement,
                                             best_placement,
                                             best_perimeter,
                                             Q,
                                             idx-1,
                                             false);

         free_indices_sigma[i] = tmp;
      }
      if (idx == 0)
      //At this point both, pi and sigma are fully determined.
      //Now the placement corresponding to pi and sigma can be calculated.
      {
         long int bbx_area = placement_for_pi_sigma(pi, pi_inverse, sigma, sigma_inverse, placement);
         if (bbx_area > 0)
         {
            cout << "IMPROVED BBX_length to " << best_perimeter
            << ", bbx_area=" << bbx_area << " / " << _total_cell_size << endl;
            for (size_t i = 0; i<_num_cells; i++)
            {
               best_placement[i] = placement[i];
            }
//             if (bbx_area <= _total_cell_size + 1)
//             {
//                cout << "Already found optimum!" << endl;
//                return;
//             }
         }
         else
         {
//             cout << "        not improved......" << endl;
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


long int Instance::placement_for_pi_sigma(vector<size_t> const & pi,
                                      vector<size_t> const & pi_inverse,
                                      vector<size_t> const & sigma,
                                      vector<size_t> const & sigma_inverse,
                                      vector<pair<x_coord, y_coord> > & placement)
{
   long int total_width(0);
   long int total_height(0);

   Q.clear();
   Q[0] = 0;
   Q[_num_cells + 1] = INT_MAX;

   for (size_t i = 0; i < _num_cells; i++)
   {
      size_t c = pi[i];
      size_t p = sigma_inverse[c];
      long int l_p;

      map<size_t, long int>::iterator it;
      it = Q.upper_bound(p+1);
      it--;
      placement[c].first = (*it).second;
      l_p = placement[c].first + _cells[c].width;

      if (l_p >= best_perimeter)
      {
//          cout << "   Break calculation, in x dimension. Looked at "
//          << i+1 << " from " << _num_cells << " cells." << endl;
         return 0;
      }

      Q[p+1] = l_p;
      total_width = max(total_width, l_p);

      it = Q.upper_bound(p+1);
      while ((*it).second <= l_p)
      {
         Q.erase(it);
         it = Q.upper_bound(p+1);
      }
   }


   Q.clear();
   Q[0] = 0;
   Q[_num_cells + 1] = 100000; //TODO: INF
   for (size_t i = _num_cells; i > 0;)
   {
      i--;
      size_t c = sigma[i];
      size_t p = pi_inverse[c];
      long int l_p;

      map<size_t, long int>::iterator it;
      it = Q.upper_bound(p+1);
      it--;
      placement[c].second = (*it).second;
      l_p = placement[c].second + _cells[c].height;

      if (l_p + total_width >= best_perimeter)
      {
//          cout << "        Break calculation, in y dimension. Looked at "
//          << _num_cells-i << " from " << _num_cells << " cells." << endl;
         return 0;
      }

      Q[p+1] = l_p;
      total_height = max(total_height, l_p);

      it = Q.upper_bound(p+1);
      while ((*it).second <= l_p)
      {
         Q.erase(it);
         it = Q.upper_bound(p+1);
      }
   }
   best_perimeter = total_width + total_height;
   return total_width*total_height;
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
