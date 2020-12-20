//////////////////////////////////////////////////////////////////////////////////////
/// \file profoundoutput.h
/// \brief Output module for COST PROFOUND
///
/// \author Joerg Steinkamp
/// $Date: Mon Jan 30 11:06:14 2017 $
///
///////////////////////////////////////////////////////////////////////////////////////
#ifndef LPJ_GUESS_PROFOUND_OUTPUT_H
#define LPJ_GUESS_PROFOUND_OUTPUT_H
#include "outputmodule.h"
#include "outputchannel.h"
#include "gutil.h"
namespace GuessOutput {
  class ProfoundOutput : public OutputModule {
  public:
    ProfoundOutput();
    ~ProfoundOutput();
    void init();
    void outannual(Gridcell& gridcell);
    void outdaily(Gridcell& gridcell);
  private:
    void define_output_tables();
    xtring file_diamclass_dbh;
    Table out_diamclass_dbh;
    xtring file_dbh_domhei;
    Table out_dbh_domhei;
    xtring file_dom_height;
    Table out_dom_height;
    xtring file_diamclass_mai;
    Table out_diamclass_mai;
    xtring file_diamclass_ba;
    Table out_diamclass_ba;
    xtring file_diamclass_harv;
    Table out_diamclass_harv;
    xtring file_diamclass_stemno;
    Table out_diamclass_stemno;
    xtring file_diamclass_age;
    Table out_diamclass_age;
    xtring file_diamclass_mortstemno;
    Table out_diamclass_mortstemno;
    xtring file_diamclass_harvstemno;
    Table out_diamclass_harvstemno;
    xtring file_diamclass_dens;
    Table out_diamclass_dens;
    xtring file_diamclass_vol;
    Table out_diamclass_vol;
    xtring file_diamclass_cveg;
    Table out_diamclass_cveg;

    xtring file_daily_gpp;
    Table out_daily_gpp;
    xtring file_daily_npp;
    Table out_daily_npp;
    xtring file_daily_ra;
    Table out_daily_ra;
    xtring file_daily_fpar;
    Table out_daily_fpar;
    xtring file_daily_cflux;
    Table out_daily_cflux;
  };
}
#endif
