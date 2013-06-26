/**
 * Authors: Christoph Hunkenschröder, Rasmus Schroeder
 * Date: June 2013
 *
 * Programme description:
 * The programme calculates for a given list of cells with widths and heights
 * a placement minimizing the bounding box length.
 * The placement is written to the standard output.
 *
 * Blala von Christoph.
 **/

#include <iostream>
#include <vector>

#include "instance.h"

using namespace std;

int main(int argc, char **argv)
{
   if (! argc == 2)
   {
      cout << "Programme call: ./placement <INSTANCE>" << endl;
      return 1;
   }

   Instance Inst;
   Inst.read_file(argv[1]);
   Inst.print();
   Inst.minimum_perimeter();

   return 0;
}
