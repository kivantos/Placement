#include <iostream>
#include <vector>
#include <map>

#include "instance.h"

using namespace std;

int main(int argc, char **argv) 
{
   if (! argc == 2)
   {
      cout << "Programme call: ./placement <INSTANCE>" << endl;
      return 1;
   }
   
   vector<int> test;

   Instance Inst;
   Inst.read_file(argv[1]);
   Inst.minimum_perimeter();
   
   Inst.print();
   return 0;
   
}
