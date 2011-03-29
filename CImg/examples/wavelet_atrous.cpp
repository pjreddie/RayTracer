/*
 #
 #  File        : wavelet_atrous.cpp
 #                ( C++ source file )
 #
 #  Description : Performs a 2D or 3D 'a trous' wavelet transform
 #                (using a cubic spline) on an image or a video sequence.
 #                This file is a part of the CImg Library project.
 #                ( http://cimg.sourceforge.net )
 #
 #  Author      : Renaud Peteri
 #                ( Renaud.Peteri(at)mines-paris.org )
 #
 #  Institution : CWI, Amsterdam
 #
 #  Date        : February 2005
 #
 #  References  : Starck, J.-L., Murtagh, F. and Bijaoui, A.,
 #                Image Processing and Data Analysis: The Multiscale Approach,
 #                Cambridge University Press, 1998.
 #                (Hardback and softback, ISBN 0-521-59084-1 and 0-521-59914-8.)
 #
 #  License     : CeCILL v2.0
 #                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed by the CeCILL  license under French law and
 #  abiding by the rules of distribution of free software.  You can  use,
 #  modify and/ or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
 #
*/

#include "CImg.h"
using namespace cimg_library;
#ifndef cimg_imagepath
#define cimg_imagepath "img/"
#endif

// Define convolution mask for the X-Axis.
CImg<float> mask_x(const unsigned char scale) {
  unsigned char d1 = (unsigned char)std::pow(2.0,(double)(scale-1));
  unsigned char d2 = (unsigned char)std::pow(2.0,(double)(scale));
  unsigned char cx = (unsigned char)std::pow(2.0,(double)(scale));
  unsigned char res = (unsigned char)std::pow(2.0,(double)scale);
  CImg<float> m(2*res +1,1,1);m.fill(0);
  m(cx) = 6.0;
  m(cx-d1) =  m(cx+d1) =4.0;
  m(cx-d2) =  m(cx+d2) =1.0;
  m /= 16.0;
  return m;
}

// Define convolution mask for the Y-Axis.
CImg<float> mask_y(const unsigned char scale) {
  unsigned char d1 = (unsigned char)std::pow(2.0,(double)(scale-1));
  unsigned char d2 = (unsigned char)std::pow(2.0,(double)(scale));
  unsigned char cy = (unsigned char)std::pow(2.0,(double)(scale));
  unsigned char res = (unsigned char)std::pow(2.0,(double)scale);
  CImg<float> m(1,2*res +1);m.fill(0);
  m(0,cy) = 6.0;
  m(0,cy-d1) =  m(0,cy+d1) =4.0;
  m(0,cy-d2) =  m(0,cy+d2) =1.0;
  m /= 16.0;
  return m;
}

// Define convolution mask for the T-Axis.
CImg<float> mask_t(const unsigned char scale) {
  unsigned char d1 = (unsigned char)std::pow(2.0,(double)(scale-1));
  unsigned char d2 = (unsigned char)std::pow(2.0,(double)(scale));
  unsigned char ct = (unsigned char)std::pow(2.0,(double)(scale));
  unsigned char res = (unsigned char)std::pow(2.0,(double)scale);
  CImg<float> m(1,1,2*res +1);m.fill(0);
  m(0,0,ct) = 6.0;
  m(0,0,ct-d1) =  m(0,0,ct+d1) =4.0;
  m(0,0,ct-d2) =  m(0,0,ct+d2) =1.0;
  m /= 16.0;
  return m;
}

/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {

  cimg_usage("Perform an 'a trous' wavelet transform (using a cubic spline) on an image or on a video sequence.\n"
             "This wavelet transform is undecimated and produces 2 images/videos at each scale. For an example of\n"
             "decomposition on a video, try -i img/trees.inr (sequence from the MIT).\n"
             "\t(Type -h for help)");

  // Read command line parameters
  const char
    *name_i  = cimg_option("-i",cimg_imagepath "lena.pgm","Input image or video"),
    *name_o  = cimg_option("-o","","Name of the multiscale analysis output"),
    *axe_dec = cimg_option("-axe",(char*)NULL,"Perform the multiscale decomposition in just one direction ('x', 'y' or 't')");
  const unsigned int
    s = cimg_option("-s",3,"Scale of decomposition");

  const bool help = cimg_option("-h",false,"Display Help");
  if(help) exit(0);

  // Initialize Image Data
  std::fprintf(stderr," - Load image sequence '%s'...\n",cimg::basename(name_i));
  const CImg<float> texture_in(name_i);
  CImg<float> mask_conv_x, mask_conv_y, mask_conv_t;
  CImgList<float> res(s,texture_in.width(),texture_in.height(),texture_in.depth());
  CImgList<float> wav(s,texture_in.width(),texture_in.height(),texture_in.depth());
  cimglist_for(res,l) { res(l).fill(0.0); wav(l).fill(0.0);}
  unsigned int i;

  if (!axe_dec){
    // Perform the multiscale decomposition in all directions
    for(i=0;i<s;i++){
      std::fprintf(stderr," - Performing scale %u ...\n",i+1);
      if(i==0){ res(i) =  texture_in;} else {  res(i) = res(i-1);}
      mask_conv_x = mask_x(i+1);
      res(i) = res(i).get_convolve(mask_conv_x);
      mask_conv_y = mask_y(i+1);
      res(i) = res(i).get_convolve(mask_conv_y);
      mask_conv_t = mask_t(i+1);
      res(i) = res(i).get_convolve(mask_conv_t);
      if(i==0){wav(i) = texture_in - res(i);}  // res(0) and wav(0) are the 1st scale of decompostion
      else {wav(i) = res(i-1) - res(i);}
    } }

  if (axe_dec) {
    // Perform the multiscale decomposition in just one direction
    char c;
    c = cimg::uncase(axe_dec[0]);
    fprintf(stderr," - Decompose the image along axe '%c'\n",c); fflush(stdout);

    switch(c) {
    case 'x': {
      for(i=0;i<s;i++) {
        std::fprintf(stderr," - Performing scale %u ...\n",i+1);
        if(i==0){ res(i) =  texture_in;} else {  res(i) = res(i-1);}
        mask_conv_x = mask_x(i+1);
        res(i) = res(i).get_convolve(mask_conv_x);
        if(i==0){wav(i) = texture_in - res(i);}
        else {wav(i) = res(i-1) - res(i);}}}
      break;

    case 'y': {
      for(i=0;i<s;i++) {
        std::fprintf(stderr," - Performing scale %u ...\n",i+1);
        if(i==0){ res(i) =  texture_in;} else {  res(i) = res(i-1);}
        mask_conv_y = mask_y(i+1);
        res(i) = res(i).get_convolve(mask_conv_y);
        if(i==0){wav(i) = texture_in - res(i);}
        else {wav(i) = res(i-1) - res(i);}}}
      break;

    case 't': {
      for(i=0;i<s;i++) {
        std::fprintf(stderr," - Performing scale %u ...\n",i+1);
        if(i==0){ res(i) =  texture_in;} else {  res(i) = res(i-1);}
        mask_conv_t = mask_t(i+1);
        res(i) = res(i).get_convolve(mask_conv_t);
        if(i==0){wav(i) = texture_in - res(i);}
        else {wav(i) = res(i-1) - res(i);}}}
      break;

    default: throw CImgException("Error, unknow decompostion axe '%c', try 'x', 'y' or 't'",c);
    }
    fputc('\n',stderr);
  }

  if (*name_o){
    // Save the Multi-Scale Analysis
    std::fprintf(stderr," - Saving of all output sequences : %s in the msa/ directory... \n",cimg::basename(name_o));
    int count = 1; // res0 = original image
    char filename[256] = "", filename_wav[256] = "";
    char STmp[3] = "";
    const int err = std::system("mkdir msa");
    if (!err) for(i=0;i<s;i++) {
      std::strcpy( filename, "msa/res" );
      std::strcpy( filename_wav, "msa/wav" );
      if( count < 10 )
        { std::strcat( filename, "0" );std::strcat( filename_wav, "0" );}
      std::sprintf( STmp, "%d_", count );
      std::strcat( filename, STmp ); std::strcat( filename_wav, STmp );
      std::strcat( filename,name_o);std::strcat( filename_wav,name_o);
      res(i).save(filename);
      wav(i).save(filename_wav);
      count++;
    }
  }

  // Result visualization
  const float value = 255;
  for(i=0;i<s;i++) {
    res[i].normalize(0,255).draw_text(2,2,"Scale %d",&value,0,1,13,i);
    wav[i].normalize(0,255).draw_text(2,2,"Scale %d",&value,0,1,13,i);
  }

  CImgDisplay disp(res,"Approximations levels by increasing scale",0);
  CImgDisplay disp2(wav,"Wavelet coefficients by increasing scale",0);
  while (!disp.is_closed() && !disp.is_keyQ() && !disp.is_keyESC() &&
         !disp2.is_closed() && !disp2.is_keyQ() && !disp2.is_keyESC()) {
    if (disp.is_resized()) disp.resize().display(res);
    if (disp2.is_resized()) disp2.resize().display(wav);
    CImgDisplay::wait(disp,disp2);
  }

  return 0;
}
