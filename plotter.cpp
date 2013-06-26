#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "plotter.h"

using namespace std;

bool
PLOTTER::initialize ( const char *filename,
                      double x_min,
                      double y_min,
                      double x_max,
                      double y_max )
{
   _filename = filename;
   //_author   = "Rasmus Schroeder";

   //strcat(_author, "Rasmus Schroeder");

   _x_min         = x_min;
   _y_min         = y_min;
   _x_max         = x_max;
   _y_max         = y_max;
   _current_color = 0;
   _scale         = min(550.0/(_x_max - _x_min), 730.0/(_y_max - _y_min));

   _file.open(_filename, ios::out | ios::trunc);

   if (!_file.is_open())
   {
      cout << "ERROR: Failed to open plotfile." << endl;
      return false;
   }

   if (!write_file_head())
   {
      cout << "ERROR: Failed to write head of plotfile." << endl;
      return false;
   }

   if (_file.bad())
   {
      cout << "ERROR: Plotfile is corrupt." << endl;
      return false;
   }

   return true;
}



bool
PLOTTER::write_file_head()
{
   if (!_file.is_open())
   {
      _file.open(_filename, ios::out | ios::trunc);
   }

   if (_file.bad()) return false;

   //cout << "Writing file head..." << endl;

   _file << "% Created by " << endl// << _author << endl
   << "% Date: " << endl << endl

   << "/VLine{dup 0 moveto Height lineto} def" << endl
   << "/HLine{dup 0 exch moveto Width exch lineto} def" << endl
   << "/c{gsave currentpoint dup 4 3 roll dup true charpath pathbbox 3 -1 roll "
   << "add 2 div 5 4 roll exch sub 3 1 roll add 2 div 5 -1 roll exch sub exch "
   << "grestore rmoveto show pop} def" << endl
   << "/ur{4 copy 0 0 0 setrgbcolor rectstroke" << endl
   << "/RectHeight exch def /RectWidth exch def /RectY exch def /RectX exch def "
   << "/Text exch def Text () ne {" << endl
   << "gsave /Helvetica findfont FontScale scalefont setfont newpath 0 0 moveto "
   << "Text false charpath flattenpath pathbbox grestore 2 index sub /TextHeight "
   << "exch store 2 index sub /TextWidth exch store pop pop" << endl
   << "/Helvetica findfont RectWidth TextWidth div RectHeight TextHeight div lt "
   << "{RectWidth TextWidth div}{RectHeight TextHeight div} ifelse 2 div "
   << "FontScale mul scalefont setfont" << endl
   << "Text RectX RectWidth 2 div add RectY RectHeight 2 div add moveto c} "
   << "if} def" << endl
   << "/r{1 sub NumRect 1 sub div 3 div 2 mul 1 1 sethsbcolor 4 copy rectfill 4 "
   << "copy 0 0 0 setrgbcolor rectstroke" << endl
   << "/RectHeight exch def /RectWidth exch def /RectY exch def /RectX exch "
   << "def /Text exch def Text () ne {" << endl
   << "gsave /Helvetica findfont FontScale scalefont setfont newpath 0 0 "
   << "moveto Text false charpath flattenpath pathbbox grestore 2 index sub "
   << "/TextHeight exch store 2 index sub /TextWidth exch store pop pop" << endl
   << "/Helvetica findfont RectWidth TextWidth div RectHeight TextHeight div "
   << "lt {RectWidth TextWidth div}{RectHeight TextHeight div} ifelse 2 div  "
   << "FontScale mul scalefont setfont" << endl
   << "Text RectX RectWidth 2 div add RectY RectHeight 2 div add moveto c} "
   << "if} def" << endl
   << "/rect{1 sub NumRect 1 sub div 3 div 2 mul 1 1 sethsbcolor 4 copy "
   << "rectfill  0 0 0 setrgbcolor rectstroke } def" << endl << endl

   << "0 setlinewidth" << endl
   << "20 20 translate" << endl
   //<< _scale << " " << _scale << " scale" << endl
   //<< -_x_min << " " << -_y_min << " translate" << endl
   << "0 0 0 setrgbcolor" << endl << endl;

   _current_color = 0;

   if (_file.bad()) return false;

   return true;
}



void
PLOTTER::line(double x1,
              double y1,
              double x2,
              double y2,
              string color_name)
{
   unsigned int col = color(color_name);

   if (col != _current_color)
   {
      _current_color = col;

      _file << color_code(col) << " setrgbcolor" << endl;
   }

   _file << px(x1) << " " << py(y1) << " moveto " << px(x2) << " " << py(y2)
   << " lineto stroke" << endl;
}

void
PLOTTER::rectangle(double x1,
                   double y1,
                   double x2,
                   double y2,
                   string color_name)
{
   unsigned int col = color(color_name);

   if (col != _current_color)
   {
      _current_color = col;

      _file << color_code(col) << " setrgbcolor" << endl;
   }

   _file << px(x1) << " " << py(y1) << " moveto "
   << px(x2) << " " << py(y1) << " lineto "
   << px(x2) << " " << py(y2) << " lineto "
   << px(x1) << " " << py(y2) << " lineto closepath fill" << endl;
}

void
PLOTTER::rectangle_empty(double x1,
                        double y1,
                        double x2,
                        double y2,
                        string color_name)
{
   unsigned int col = color(color_name);

   if (col != _current_color)
   {
      _current_color = col;

      _file << color_code(col) << " setrgbcolor" << endl;
   }

   _file << px(x1) << " " << py(y1) << " moveto "
   << px(x2) << " " << py(y1) << " lineto "
   << px(x2) << " " << py(y2) << " lineto "
   << px(x1) << " " << py(y2) << " lineto closepath stroke" << endl;
}



void
PLOTTER::point(double x,
               double y,
               string color_name)
{
   unsigned int col = color(color_name);

   if (col != _current_color)
   {
      _current_color = col;

      _file << color_code(col) << " setrgbcolor" << endl;
   }

   _file << px(x) << " " << py(y) << " moveto currentpoint "
   << _scale << " 0 360 arc fill" << endl;
}


void
PLOTTER::circle(double x,
                double y,
                double radius,
                string color_name)
{
   unsigned int col = color(color_name);

   if (col != _current_color)
   {
      _current_color = col;
      _file << color_code(col) << " setrgbcolor" << endl;
   }

   //TODO: Waagerechte Linie im Kreis noch entfernen.
   _file << px(x) << " " << py(y) << " moveto currentpoint "
   << radius*_scale << " 0 360 arc stroke" << endl;
}


unsigned int
PLOTTER::color(string color_name)
{
//    string str = trim(color_name);
   string str = color_name;


   if (str.compare("black") == 0)        return 0;
   else if (str.compare("Black") == 0)   return 0;
   else if (str.compare("BLACK") == 0)   return 0;
   else if (str.compare("schwarz") == 0) return 0;
   else if (str.compare("Schwarz") == 0) return 0;
   else if (str.compare("SCHWARZ") == 0) return 0;

   else if (str.compare("blue") == 0)    return 1;
   else if (str.compare("Blue") == 0)    return 1;
   else if (str.compare("BLUE") == 0)    return 1;
   else if (str.compare("blau") == 0)    return 1;
   else if (str.compare("Blau") == 0)    return 1;
   else if (str.compare("BLAU") == 0)    return 1;

   else if (str.compare("red") == 0)     return 2;
   else if (str.compare("Red") == 0)     return 2;
   else if (str.compare("RED") == 0)     return 2;
   else if (str.compare("rot") == 0)     return 2;
   else if (str.compare("Rot") == 0)     return 2;
   else if (str.compare("ROT") == 0)     return 2;

   else if (str.compare("green") == 0)   return 3;
   else if (str.compare("Green") == 0)   return 3;
   else if (str.compare("GREEN") == 0)   return 3;
   else if (str.compare("gruen") == 0)   return 3;
   else if (str.compare("Gruen") == 0)   return 3;
   else if (str.compare("GRUEN") == 0)   return 3;

   else if (str.compare("yellow") == 0)   return 4;
   else if (str.compare("grey") == 0)   return 5;
   else if (str.compare("violet") == 0)   return 6;
   else if (str.compare("magenta") == 0)   return 7;
   else if (str.compare("pink") == 0)   return 8;
   else if (str.compare("brown") == 0)   return 9;
   else if (str.compare("orange") == 0)   return 10;
   else if (str.compare("forrestgreen") == 0)   return 11;

   else return 12;
}


string
PLOTTER::color_code(unsigned int color)
{
   if (color == 0) return "0 0 0";          //black
   else if (color == 1) return "0 0 1";     //blue
   else if (color == 2) return "1 0 0";     //red
   else if (color == 3) return "0 1 0";     //green
   else if (color == 4) return "1 1 0";     //yellow
   else if (color == 5) return "0.6 0.6 0.6";     //grey
   else if (color == 6) return "0.5 0.2 0.5";     //violet
   else if (color == 7) return "1 0 1";     //magenta
   else if (color == 8) return "1 0.7 0.8";     //pink
   else if (color == 9) return "0.5 0.3 0.07";     //brown
   else if (color == 10) return "1 0.55 0";     //orange
   else if (color == 11) return "0.133 0.54 0.133";     //forrestgreen
   else return "0 0 0";
}


void
PLOTTER::close()
{
   _file.close();
}


void
PLOTTER::print_data()
{
   cout << "Plotter info:" << endl
   << "Filename: " << _filename << endl
   //<< "Author: " << _author << endl
   << "x_min=" << _x_min << endl
   << "y_min=" << _y_min << endl
   << "x_max=" << _x_max << endl
   << "y_max=" << _y_max << endl << endl;
}
