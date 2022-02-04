#ifndef RESULT_SFTMOM_H
#define RESULT_SFTMOM_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stack>
#include "multi.h"

#define Nd 4

using namespace std;

namespace qa
{

  class SftMom
  {
    public:
      SftMom(int mom2_max, multi1d<int> mom_offset_,
           bool avg_equiv_mom_ = false, int j_decay = -1)
      { init(mom2_max, mom_offset_, avg_equiv_mom_, j_decay); }


      int numMom() const { return num_mom; }
      multi1d<int> numToMom(int mom_num) const { return mom_list[mom_num]; }
      int momToNum(const multi1d<int>& mom_in) const;

    private:
      SftMom() {}

      multi1d<int> crtesn(int ipos, const multi1d<int>& latt_size);

      void init(int mom2_max, multi1d<int> mom_offset,
              bool avg_mom_ = false, int j_decay = -1);

      multi2d<int> mom_list;
      bool         avg_equiv_mom;
      int          decay_dir;
      int          num_mom;
      multi1d<int> origin_offset;
      multi1d<int> mom_offset;
      multi1d<int> mom_degen;
  };

}

#endif
