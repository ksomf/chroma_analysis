//  ==============
//  QCDSF/ANALYSIS
//  ==============
//! Fourier transform phase factor support taken from Chroma
//! \file

#include "sftmom.h"

namespace qa
{


  multi1d<int> SftMom::crtesn(int ipos, const multi1d<int>& latt_size)
  {
    multi1d<int> coord(latt_size.size());

    for(int i=0; i < latt_size.size(); ++i)
    {
      coord[i] = ipos % latt_size[i];
      ipos = ipos / latt_size[i];
    }

    return coord;
  }


  void SftMom::init(int mom2_max, multi1d<int> mom_off, bool avg_mom, int j_decay)
  {
    decay_dir     = j_decay;  
    mom_offset    = mom_off;  
    avg_equiv_mom = avg_mom;  

    multi1d<int> mom_size ;
    if ((j_decay<0)||(j_decay>=Nd)) {
      mom_size.resize(Nd) ;
    } else {
      mom_size.resize(Nd-1) ;
    }

    int L;
    int mom_vol = 1;

    for (L=1; L*L <= mom2_max; ++L) ;

    for(int mu=0; mu < mom_size.size(); ++mu) {
      if (avg_equiv_mom) {  
        mom_vol      *= L;
        mom_size[mu]  = L;
      } else {
        mom_vol      *= (2*L) + 1;
        mom_size[mu]  = (2*L) + 1;
      }
    }

    num_mom = 0;


    for(int n=0; n < mom_vol; ++n) {
      multi1d<int> mom = crtesn(n, mom_size);

      int mom2 = 0 ;

      for(int mu=0; mu < mom_size.size(); ++mu) {
        if (!avg_equiv_mom) mom[mu] -= L ;
        mom2 += mom[mu]*mom[mu];
      }

      if (mom2 > mom2_max) {
        continue;
      } else if (avg_equiv_mom) {
        bool skip = false ;
        for(int mu=0; mu < mom_size.size()-1; ++mu)
  	for(int nu=mu+1; nu < mom_size.size(); ++nu)
  	  if (mom[nu] > mom[mu]) skip=true;

        if (!skip) ++num_mom ;
      } else {
        ++num_mom ;
      }
    }


    mom_list.resize(num_mom, mom_size.size());

    int mom_num = 0;

    for(int n=0; n < mom_vol; ++n)
      {
        multi1d<int> mom = crtesn(n, mom_size);

        int mom2 = 0 ;

        for(int mu=0; mu < mom_size.size(); ++mu) {
  	if (!avg_equiv_mom) mom[mu] -= L ;
  	mom2 += mom[mu]*mom[mu];
        }

        if (mom2 > mom2_max) {
  	continue;
        } else if (avg_equiv_mom) {
  	bool skip = false ;
  	for(int mu=0; mu < mom_size.size()-1; ++mu)
  	  for(int nu=mu+1; nu < mom_size.size(); ++nu)
  	    if (mom[nu] > mom[mu]) skip = true ;

  	if (!skip) mom_list[mom_num++] = mom ;
        } else {
  	for (int mu=0; mu < mom_size.size(); ++mu) {
  	  mom_list[mom_num][mu] = mom_offset[mu] + mom[mu]  ;
  	}
  	++mom_num ;
        }
      }


  }


  int SftMom::momToNum(const multi1d<int>& mom_in) const
  {
    multi1d<int> mom;

    mom = mom_in;

    for(int mom_num=0; mom_num < num_mom; ++mom_num)
      {
        bool match = true ;
        for (int mu=0; mu < mom.size(); ++mu)
  	{
  	  if (mom_list[mom_num][mu] != mom[mu])
  	    {
  	      match = false ;
  	      break;
  	    }
  	}
        if (match) return mom_num ;
      }
    return -1;
  }

}
