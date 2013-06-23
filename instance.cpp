#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <climits>

#include "instance.h"

using namespace std;

void Instance::read_file(const char* filename)
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
   ss >> _num_cells;                   // for which we can use << as with cin
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
}


void Instance::minimum_perimeter()
{
   vector<size_t>                  pi;                 //Permutation pi
   vector<size_t>                  pi_inverse;         //Inverse of permutation pi
   vector<size_t>                  free_indices_pi;    //Uninitialized indices of pi
   vector<size_t>                  best_pi;            //Permutation with best result
   vector<size_t>                  old_pi;
   vector<size_t>                  sigma;              //Permutation sigma
   vector<size_t>                  sigma_inverse;      //Inverse of permutation sigma
   vector<size_t>                  free_indices_sigma; //Uninitialized indices of sigma
   vector<size_t>                  best_sigma;         //Permutation with best result
   vector<size_t>                  old_sigma;
   vector<pair<x_coord, y_coord> > placement;          //Coordinates of the cells in a placement
   vector<pair<x_coord, y_coord> > best_placement;     //Coordinates of current best placement
   map<size_t, long int>           Q;                  //Priority queue needed for the placement alg.
   long int best_perimeter = INT_MAX;                  //Current best perimeter of placement;
   long int lower_bound_perimeter = 0;


   //Reserve memory for all vectors.
   //We only allocate memory once for every vector.
   pi.resize(_num_cells);
   pi_inverse.resize(_num_cells);
   free_indices_pi.resize(_num_cells);
   best_pi.resize(_num_cells);
   old_pi.resize(_num_cells);

   sigma.resize(_num_cells);
   sigma_inverse.resize(_num_cells);
   free_indices_sigma.resize(_num_cells);
   best_sigma.resize(_num_cells);
   old_sigma.resize(_num_cells);

   placement.resize(_num_cells);
   best_placement.resize(_num_cells);

   for (size_t k = 1; k <= _num_cells; k++)
   {
      //At the beginning, all indices of pi and sigma are uninitialized.
      for (size_t i = 0; i<k; i++)
      {
         free_indices_pi[i] = i;
         free_indices_sigma[i] = i;
         if (i!=k-1) old_pi[i] = best_pi[i];
         if (i!=k-1) old_sigma[i] = best_sigma[i];
      }
      best_perimeter = INT_MAX;

      cout << "Calculate for k=" << k << "..." << endl;

      bool success;
      success = find_placement_one_free_cell(old_pi,
                                             best_pi,
                                             old_sigma,
                                             best_sigma,
                                             placement,
                                             best_placement,
                                             best_perimeter,
                                             lower_bound_perimeter,
                                             k,
                                             Q);

      if (!success)
      {
         find_perimeter_for_all_permutations(pi,
                                             pi_inverse,
                                             free_indices_pi,
                                             best_pi,
                                             sigma,
                                             sigma_inverse,
                                             free_indices_sigma,
                                             best_sigma,
                                             placement,
                                             best_placement,
                                             best_perimeter,
                                             lower_bound_perimeter,
                                             k,  //Number of cells we look at.
                                             Q,
                                             k,  //Index to be fixed first.
                                             true);

         lower_bound_perimeter = best_perimeter;
      }
      else
      {
         cout << "Could use old placement to insert one new cell!" << endl;
      }

      cout << "...done. Lower bound improved to " << lower_bound_perimeter << endl;
   }

   //Output
   //TODO: Correct format.
   cout << endl << "BBX=" << best_perimeter << endl;
   cout << "Placement:" << endl;
   for (size_t i = 0; i < _num_cells; i++)
   {
      cout << "Circuit " << i << ": (" << best_placement[i].first
      << ", " << best_placement[i].second << ")" << endl;
   }
}


bool Instance::find_perimeter_for_all_permutations(vector<size_t> & pi,
                                                   vector<size_t> & pi_inverse,
                                                   vector<size_t> & free_indices_pi,
                                                   vector<size_t> & best_pi,
                                                   vector<size_t> & sigma,
                                                   vector<size_t> & sigma_inverse,
                                                   vector<size_t> & free_indices_sigma,
                                                   vector<size_t> & best_sigma,
                                                   vector<pair<x_coord, y_coord> > & placement,
                                                   vector<pair<x_coord, y_coord> > & best_placement,
                                                   long int & best_perimeter,
                                                   long int & lower_bound_perimeter,
                                                   size_t k,
                                                   map<size_t, long int> & Q,
                                                   size_t idx,
                                                   bool pi_or_sigma)
{
   bool found_optimum = false;

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


         found_optimum = find_perimeter_for_all_permutations(pi,
                                                            pi_inverse,
                                                            free_indices_pi,
                                                            best_pi,
                                                            sigma,
                                                            sigma_inverse,
                                                            free_indices_sigma,
                                                            best_sigma,
                                                            placement,
                                                            best_placement,
                                                            best_perimeter,
                                                            lower_bound_perimeter,
                                                            k,
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
         found_optimum = find_perimeter_for_all_permutations(pi,
                                                            pi_inverse,
                                                            free_indices_pi,
                                                            best_pi,
                                                            sigma,
                                                            sigma_inverse,
                                                            free_indices_sigma,
                                                            best_sigma,
                                                            placement,
                                                            best_placement,
                                                            best_perimeter,
                                                            lower_bound_perimeter,
                                                            k,
                                                            Q,
                                                            k,
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

         found_optimum = find_perimeter_for_all_permutations(pi,
                                                            pi_inverse,
                                                            free_indices_pi,
                                                            best_pi,
                                                            sigma,
                                                            sigma_inverse,
                                                            free_indices_sigma,
                                                            best_sigma,
                                                            placement,
                                                            best_placement,
                                                            best_perimeter,
                                                            lower_bound_perimeter,
                                                            k,
                                                            Q,
                                                            idx-1,
                                                            false);

         free_indices_sigma[i] = tmp;
      }
      if (idx == 0)
      //At this point both, pi and sigma are fully determined.
      //Now the placement corresponding to pi and sigma can be calculated.
      {
         found_optimum = placement_for_pi_sigma(pi,
                                                pi_inverse,
                                                best_pi,
                                                sigma,
                                                sigma_inverse,
                                                best_sigma,
                                                placement,
                                                best_placement,
                                                best_perimeter,
                                                lower_bound_perimeter,
                                                k,
                                                Q);
      }
   }

   if (found_optimum) return true;   //Optimum placement found.
   else return false;                //Optimum placement not yet found.
}


bool Instance::placement_for_pi_sigma(vector<size_t> const & pi,
                                      vector<size_t> const & pi_inverse,
                                      vector<size_t> & best_pi,
                                      vector<size_t> const & sigma,
                                      vector<size_t> const & sigma_inverse,
                                      vector<size_t> & best_sigma,
                                      vector<pair<x_coord, y_coord> > & placement,
                                      vector<pair<x_coord, y_coord> > & best_placement,
                                      long int & best_perimeter,
                                      long int & lower_bound_perimeter,
                                      size_t k,
                                      map<size_t, long int> & Q)
{
   //Find x-coordinates
   long int total_width(0);
   Q.clear();
   Q[0] = 0;
   Q[k + 1] = INT_MAX;

   for (size_t i = 0; i < k; i++)
   {
      size_t c = pi[i];
      size_t p = sigma_inverse[c];
      long int l_p;

      map<size_t, long int>::iterator it;
      it = Q.upper_bound(p+1);
      it--;
      placement[c].first = (*it).second;
      l_p = placement[c].first + _cells[c].width;

      total_width = max(total_width, l_p);
      if (total_width >= best_perimeter) return false;   //Total width is already too long.

      Q[p+1] = l_p;

      it = Q.upper_bound(p+1);
      while ((*it).second <= l_p)
      {
         Q.erase(it);
         it = Q.upper_bound(p+1);
      }
   }


   //Find y-coordinates
   long int total_height(0);
   Q.clear();
   Q[0] = 0;
   Q[k + 1] = INT_MAX;
   for (size_t i = k; i > 0;)
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

      total_height = max(total_height, l_p);
      if (total_height + total_width >= best_perimeter) return false;

      Q[p+1] = l_p;

      it = Q.upper_bound(p+1);
      while ((*it).second <= l_p)
      {
         Q.erase(it);
         it = Q.upper_bound(p+1);
      }
   }


   //The currently calculated placement is the best so far.
   //Store the values.
   best_perimeter = total_width + total_height;
   for (size_t i = 0; i < k; i++)
   {
      best_placement[i].first  = placement[i].first;
      best_placement[i].second = placement[i].second;
      best_pi[i] = pi[i];
      best_sigma[i] = sigma[i];
   }

   if (best_perimeter == lower_bound_perimeter) return true;  //Optimum placement found.
   else return false;
}


bool Instance::find_placement_one_free_cell(vector< size_t > const & pi_old,
                                            vector< size_t > & best_pi,
                                            vector< size_t > const & sigma_old,
                                            vector< size_t > & best_sigma,
                                            vector< pair< x_coord, y_coord > >& placement,
                                            vector< pair< x_coord, y_coord > >& best_placement,
                                            long int& best_perimeter,
                                            long int& lower_bound_perimeter,
                                            size_t k,
                                            std::map< size_t,
                                            long int >& Q)
{
   if (k <= 1) return false;

   bool found_optimum;
   vector<size_t> pi;
   vector<size_t> pi_inverse;
   vector<size_t> sigma;
   vector<size_t> sigma_inverse;

   pi.resize(k);
   pi_inverse.resize(k);
   sigma.resize(k);
   sigma_inverse.resize(k);

   for (size_t i_pi = 0; i_pi < k; i_pi++)
   {
      for (size_t t = 0; t < k; t++)
      {
         if (t<i_pi)       pi[t] = pi_old[t];
         else if (t==i_pi) pi[t] = k-1;
         else if (t>i_pi)  pi[t] = pi_old[t-1];
      }
      for (size_t t = 0; t < k; t++)
      {
         pi_inverse[pi[t]] = t;
      }

      for (size_t i_sigma = 0; i_sigma < k; i_sigma++)
      {
         for (size_t t = 0; t < k; t++)
         {
            if (t<i_sigma)       sigma[t] = sigma_old[t];
            else if (t==i_sigma) sigma[t] = k-1;
            else if (t>i_sigma)  sigma[t] = sigma_old[t-1];
         }
         for (size_t t = 0; t < k; t++)
         {
            sigma_inverse[sigma[t]] = t;
         }

//          cout << "pi old: ";
//          for (size_t t = 0; t < k-1; t++) cout << pi_old[t] << " ";
//          cout << endl << "pi new: ";
//          for (size_t t = 0; t < k; t++) cout << pi[t] << " ";
//          cout << endl << "sigma old: ";
//          for (size_t t = 0; t < k-1; t++) cout << sigma_old[t] << " ";
//          cout << endl << "sigma new: ";
//          for (size_t t = 0; t < k; t++) cout << sigma[t] << " ";
//          cout << endl << endl;

         found_optimum = placement_for_pi_sigma(pi,
                                                pi_inverse,
                                                best_pi,
                                                sigma,
                                                sigma_inverse,
                                                best_sigma,
                                                placement,
                                                best_placement,
                                                best_perimeter,
                                                lower_bound_perimeter,
                                                k,
                                                Q);

         if (found_optimum) return true;  //Optimum placement found!
      }
   }
   return false;
}


void Instance::print()
{
   cout << "INSTANCE: (" << _num_cells << " cells)" << endl;
   for (size_t i=0; i<_cells.size(); i++)
   {
      cout << "Cell " << i << ": " << _cells[i].width << " " << _cells[i].height << endl;
   }
   cout << endl;
}
