#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <climits>
#include <sstream>

#include "instance.h"
#include "plotter.h"

using namespace std;



/**This class is used as a priority queue.
 * All elements of the queue have a key and a value.
 * They are sorted by their key.
 *
 * There will always be at least two elements using the keys
 * 0 and n+1, where n is the number of cells to place in the algorithm.
 * The values are 0 respectively INT_MAX.
 * All other elements will have keys inbetween 0 and n+1.
 * Note that the cell with index i will always have the key i+1.
 * Hereby the indices of the cells are 0,1,...,n-1.
 *
 * Remark: We are aware of the fact that the current implementation
 * may not satisfy the theoretical runtime. However, in practice
 * this is much faster than using maps, since we do not need to
 * allocate memery every time we insert elements.
 * This will not be a good solution for large instances.
 * But the implementation does not run in any realistic time on large
 * instances anyway.
 **/
class Queue
{
public:
   Queue(unsigned int size);

   /**Set the value of the element which has key _key.
    **/
   void set_value(size_t _key, long int value);

   /**Delete all elements but the first and the last.
    **/
   void reset();

   /**Return the value of the element with largest key smaller than _key.
    **/
   long int val_of_pred(size_t _key);

   /**Delete all elements which have a key greater than _key and
    * value smaller than the value with key _key.
    **/
   void delete_smaller_succ(size_t _key);

   /**Print the queue.
    **/
   void print();

private:
   vector<pair<size_t, long int> > _q;        //The queue
   size_t                          _length;   //Length of queue for internal use only.
};

Queue::Queue(unsigned int size)
{
   _q.resize(size+2);
   _q[0].first = 0;
   _q[0].second = 0;
   _q[1].first = size+1;
   _q[1].second = INT_MAX;
   _length = 2;
}

void Queue::set_value(size_t key, long int value)
{
   size_t cur_key = 1;
   while (_q[cur_key].first < key)
   {
      cur_key++;
   }
   for (size_t i = _length; i > cur_key; i--)
   {
      _q[i].first = _q[i-1].first;
      _q[i].second = _q[i-1].second;
   }
   _q[cur_key].first = key;
   _q[cur_key].second = value;
   _length += 1;
}

void Queue::reset()
{
   _q[0].first = 0;
   _q[0].second = 0;
   _q[1].first = _q.size()-1;
   _q[1].second = INT_MAX;
   _length = 2;
}


long int Queue::val_of_pred(size_t key)
{
   size_t cur_key = 1;
   while (_q[cur_key].first < key)
   {
      cur_key++;
   }
   return _q[cur_key-1].second;
}

void Queue::delete_smaller_succ(size_t key)
{
   size_t cur_key = 1;
   while (_q[cur_key].first < key)
   {
      cur_key++;
   }
   if (_q[cur_key].first != key)
   {
      cout << "FEHLER::::::" << endl;
   }
   if (cur_key >= _length-1)
   {
      cout << "Fehler mit cur_key" << endl;
      print();
      cout << "key=" << key << endl << endl;
   }
   if (_q[cur_key+1].second >= _q[cur_key].second)
   {
      return;
   }

   long int compare_value = _q[cur_key].second;
   size_t write = cur_key+1;
   size_t read = cur_key+2;

   while (read < _length)
   {
      if (_q[read].second > compare_value)
      {
         _q[write].first = _q[read].first;
         _q[write].second = _q[read].second;
         write++;
         read++;
      }
      else
      {
         read++;
      }
   }
   _length = write;
}


void Queue::print()
{
   cout << "Q " << _length << "/" << _q.size() << ": ";
   for (size_t i = 0; i<_q.size(); i++)
   {
      cout << "(" << _q[i].first << "," << _q[i].second << ") ";
   }
   cout << endl;
}




void Instance::read_file(const char* filename)
{
   vector<Cell> cells_puffer;
   _cells.clear();

   _logfile.open("logfile.txt", ios::out|ios::trunc);

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
         cells_puffer.push_back(c);
         cells_puffer[cells_puffer.size()-1].idx = cells_puffer.size()-1;
         _total_cell_size += c.width * c.height;
   }

   if (cells_puffer.size() != _num_cells)
   {
      cerr << "Invalid file format: _num_cells=" << _num_cells << ", vector has "
      << cells_puffer.size() << " elements." << endl;
      return;
   }

   bool find_longest = true;
   size_t best_idx;
   long int highest_value;

   while (!cells_puffer.empty())
   {
      highest_value = 0;
      for (size_t i = 0; i < cells_puffer.size(); i++)
      {
         if (find_longest)
         {
            if (cells_puffer[i].width > highest_value)
            {
               highest_value = cells_puffer[i].width;
               best_idx = i;
            }
         }
         else
         {
            if (cells_puffer[i].height > highest_value)
            {
               highest_value = cells_puffer[i].width;
               best_idx = i;
            }
         }
      }
      _cells.push_back(cells_puffer[best_idx]);
      cells_puffer[best_idx] = cells_puffer[cells_puffer.size()-1];
      cells_puffer.pop_back();

      find_longest = find_longest ? false : true;
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
   Queue Q(_num_cells);                                  //Priority queue needed for the placement alg.
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

      _logfile << "Calculate for k=" << k << "..." << endl;

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
         _logfile << "Could use old placement to insert one new cell!" << endl;
      }

      _logfile << "...done. Lower bound improved to " << lower_bound_perimeter << endl;

      string plot_filename;
      if (k==0) plot_filename = "placement_0.eps";
      else if (k==1) plot_filename = "placement_1.eps";
      else if (k==2) plot_filename = "placement_2.eps";
      else if (k==3) plot_filename = "placement_3.eps";
      else if (k==4) plot_filename = "placement_4.eps";
      else if (k==5) plot_filename = "placement_5.eps";
      else if (k==6) plot_filename = "placement_6.eps";
      else if (k==7) plot_filename = "placement_7.eps";
      else if (k==8) plot_filename = "placement_8.eps";
      else if (k==9) plot_filename = "placement_9.eps";
      else if (k==10) plot_filename = "placement_10.eps";
      else if (k==11) plot_filename = "placement_11.eps";
      else plot_filename = "placement_last.eps";

      plot_placement(best_placement, plot_filename, k);

   }

   output(best_placement);

   plot_placement(best_placement, "placement_final.eps", best_placement.size());
}

void
Instance::plot_placement(std::vector< std::pair< x_coord, y_coord > > const & placement,
                         string filename_string,
                         size_t size)
{
   const char *filename = filename_string.c_str();
   long int x_min = 0;
   long int x_max = 0;
   long int y_min = 0;
   long int y_max = 0;

   for (size_t i = 0; i < size; i++)
   {
      x_max = max(x_max, placement[i].first + _cells[i].width);
      y_max = max(y_max, placement[i].second + _cells[i].height);
   }

   PLOTTER plot;
   plot.initialize(filename,
                   (double) x_min,
                   (double) y_min,
                   (double) x_max,
                   (double) y_max);

   plot.rectangle_empty((double) x_min,
                        (double) y_min,
                        (double) x_max,
                        (double) y_max,
                        "black");

   string color;

   for (size_t i = 0; i < size; i++)
   {
      if (i==0) color = "blue";
      else if (i==1) color = "red";
      else if (i==2) color = "yellow";
      else if (i==3) color = "forrestgreen";
      else if (i==4) color = "orange";
      else if (i==5) color = "violet";
      else if (i==6) color = "magenta";
      else if (i==7) color = "pink";
      else if (i==8) color = "brown";
      else if (i==9) color = "grey";
      else if (i==10) color = "black";

      plot.rectangle((double) placement[i].first,
                     (double) placement[i].second,
                     (double) (placement[i].first + _cells[i].width),
                     (double) placement[i].second + _cells[i].height,
                     color);
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
                                                   Queue & Q,
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
                                      Queue & Q)
{
   //Find x-coordinates
   long int total_width(0);
   Q.reset();

   for (size_t i = 0; i < k; i++)
   {
      size_t c = pi[i];
      size_t p = sigma_inverse[c];
      placement[c].first = Q.val_of_pred(p+1);
      long int l_p = placement[c].first + _cells[c].width;

      total_width = max(total_width, l_p);
      if (total_width >= best_perimeter) return false;   //Total width is already too long.

      Q.set_value(p+1, l_p);
      Q.delete_smaller_succ(p+1);
   }


   //Find y-coordinates
   long int total_height(0);
   Q.reset();

   for (size_t i = k; i > 0;)
   {
      i--;
      size_t c = sigma[i];
      size_t p = pi_inverse[c];

      placement[c].second = Q.val_of_pred(p+1);
      long int l_p = placement[c].second + _cells[c].height;

      total_height = max(total_height, l_p);
      if (total_height + total_width >= best_perimeter) return false;

      Q.set_value(p+1, l_p);
      Q.delete_smaller_succ(p+1);
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
                                            Queue & Q)
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


void Instance::output(const std::vector< std::pair< x_coord, y_coord > >& placement)
{
   long int width = 0;
   long int height = 0;
   for (size_t i = 0; i < _num_cells; i++)
   {
      width = max(width, placement[i].first + _cells[i].width);
      height = max(height, placement[i].second + _cells[i].height);
   }

   _logfile << endl << "BBX=" << width << " * " << height << endl
   << "Perimeter = " << width+height << endl;
   _logfile << "Placement:" << endl;
   cout << width << " " << height << endl;

   for (size_t i = 0; i < _num_cells; i++)
   {
      for (size_t j = 0; j < _num_cells; j++)
      {
         if (_cells[j].idx == i)
         {
            _logfile << "Circuit " << i << ": (" << placement[j].first
            << ", " << placement[j].second << ")" << endl;
            cout << placement[j].first << " " << placement[j].second << endl;
            j = _num_cells;
         }
      }
   }
}


void Instance::print()
{
   _logfile << "INSTANCE: (" << _num_cells << " cells)" << endl;
   for (size_t i=0; i<_cells.size(); i++)
   {
      _logfile << "Cell " << i << ": " << _cells[i].width << " " << _cells[i].height << endl;
   }
   _logfile << endl;
}
