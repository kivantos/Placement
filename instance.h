#ifndef INSTANCE_H
#define INSTANCE_H
#include <map>


typedef long int x_coord;
typedef long int y_coord;

class Cell
{
public:
   long int height;
   long int width;
};


class Instance
{
public:

   void read_file(char const * filename);

   /**Finds a placement minimizing the bounding box.
    * Output the placement to the standard output.
    **/
   void minimum_perimeter();


   //DEBUG FUNCTIONS
   /**Prints Instance information.
    **/
   void print();

private:

   /**Recursively creates all permutations for pi and sigma.
    * The values are stored in pi respectively sigma.
    * At the lowest level of the recursion, the placement corresponding
    * to the permutations pi and sigma is determined.
    * The current best perimeter during the calculation is stored
    * at best_perimeter. The current best placement is stored in the vector best_placement.
    **/
   void find_perimeter_for_all_permutations(std::vector<size_t> & pi,
                                            std::vector<size_t> & pi_inverse,
                                            std::vector<size_t> & free_indices_pi,
                                            std::vector<size_t> & sigma,
                                            std::vector<size_t> & sigma_inverse,
                                            std::vector<size_t> & free_indices_sigma,
                                            std::vector<std::pair<x_coord, y_coord> > & placement,
                                            std::vector<std::pair<x_coord, y_coord> > & best_placement,
                                            long int & best_perimeter,
                                            std::map<size_t, long int> & Q,
                                            size_t idx,
                                            bool pi_or_sigma);


   long int placement_for_pi_sigma(std::vector<size_t> const & pi,
                               std::vector<size_t> const & pi_inverse,
                               std::vector<size_t> const & sigma,
                               std::vector<size_t> const & sigma_inverse,
                               std::vector<std::pair<x_coord, y_coord> > & placement);

   unsigned int      _num_cells;        //Dimension
   long int          _total_cell_size;  //Sum of all cell sizes.
   std::vector<Cell> _cells;
};

#endif