#ifndef DEFINE_PLOTTER_H
#define DEFINE_PLOTTER_H

#include <iostream>
#include <fstream>

using namespace std;


class PLOTTER
{
public:

   ////////////////////////////
   //   DRAWING SECTION      //
   ////////////////////////////

   void
   line(double x1,
        double y1,
        double x2,
        double y2,
        string color);

   void
   rectangle(double x1,
            double y1,
            double x2,
            double y2,
             string color);

   void
   rectangle_empty(double x1,
                  double y1,
                  double x2,
                  double y2,
                  string color);

   void
   point(double x,
         double y,
         string color);

   //TODO: Waagerechte Linie im Kreis entfernen.
   void
   circle(double x,
          double y,
          double radius,
          string color);




   ////////////////////////////
   //      META SECTION      //
   ////////////////////////////

   inline double
   px(double x)
   {
      return ((x - _x_min)*_scale);
   }

   inline double
   py(double y)
   {
      return ((y - _y_min)*_scale);
   }

   bool
   initialize(const char *filename,
              double x_min,
              double y_min,
              double x_max,
              double y_max);

   void
   close();


   void
   print_data();




private:

   bool
   write_file_head();

   unsigned int
   color(string color);

   string
   color_code(unsigned int color);


   fstream      _file;

   const char        *_filename;
   char        *_author;

   double       _x_min;
   double       _y_min;
   double       _x_max;
   double       _y_max;


   // Decodes the last color used.
   // 0 - black
   // 1 - blue
   // 2 - red
   // 3 - green
   // 4 - yellow
   // 5 - violet
   unsigned int _current_color;
   double       _scale;


};


#endif