#ifndef INSTANCE_H
#define INSTANCE_H


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
   void find_perimeter_for_all_permutations(std::vector<size_t> & pi, 
                                            std::vector<size_t> & pi_inverse, 
                                            std::vector<size_t> & sigma, 
                                            std::vector<size_t> & sigma_inverse, 
                                            std::vector<size_t> & free_indices, 
                                            std::vector<std::pair<x_coord, y_coord> > & placement,
                                            size_t idx,
                                            bool pi_or_sigma);
   
   long int placement_for_pi_sigma(std::vector<size_t> const & pi, 
                               std::vector<size_t> const & pi_inverse, 
                               std::vector<size_t> const & sigma, 
                               std::vector<size_t> const & sigma_inverse, 
                               std::vector<std::pair<x_coord, y_coord> > & placement);
   
   unsigned int      _num_cells; //Dimension
   std::vector<Cell> _cells;
   
   std::vector<std::pair<x_coord, y_coord> > _best_placement;
   long int          _best_perimeter;
   long int          _total_cell_size;
   std::map<size_t, long int> Q;
};

#endif