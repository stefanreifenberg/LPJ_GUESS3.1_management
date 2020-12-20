//////////////////////////////////////////////////////////////////////////////////////
/// \file profoundoutput.cpp
/// \brief Output module for COST PROFOUND
///
/// \author Joerg Steinkamp
/// $Date: Mon Jan 30 11:06:14 2017 $
///
///////////////////////////////////////////////////////////////////////////////////////
#include "config.h"
#include "profoundoutput.h"
#include "parameters.h"
#include "guess.h"

#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace GuessOutput {
  REGISTER_OUTPUT_MODULE("profound", ProfoundOutput)
  ProfoundOutput::ProfoundOutput() {
    declare_parameter("file_dbh_domhei", &file_dbh_domhei, 300, "Mean dbh of 100/ha talles trees");
    declare_parameter("file_dom_height", &file_dom_height, 300, "Mean height of 100/ha tallest trees");
    declare_parameter("file_diamclass_dbh", &file_diamclass_dbh, 300, "Mean diameter at breastheight");
    declare_parameter("file_diamclass_mai", &file_diamclass_mai, 300, "Mean annual stand volume increment");
    declare_parameter("file_diamclass_ba", &file_diamclass_ba, 300, "Basal Area per diameter class");
    declare_parameter("file_diamclass_harv", &file_diamclass_harv, 300, "Harvested volume per diameter class");
    declare_parameter("file_diamclass_stemno", &file_diamclass_stemno, 300, "Remaining stem no. per diameter class");
    declare_parameter("file_diamclass_age", &file_diamclass_age, 300, "Avg. age per diameter class");
    declare_parameter("file_diamclass_mortstemno", &file_diamclass_mortstemno, 300, "Stem no. removed by natural mortality per diameter class");
    declare_parameter("file_diamclass_harvstemno", &file_diamclass_harvstemno, 300, "Harvested stem no. per diameter class");
    declare_parameter("file_diamclass_dens", &file_diamclass_dens, 300, "Tree density per diameter class");
    declare_parameter("file_diamclass_vol", &file_diamclass_vol, 300, "Stand volume per diameter class");
    declare_parameter("file_diamclass_cveg", &file_diamclass_cveg, 300, "Carbon mass per diameter class");

    // daily
    declare_parameter("file_daily_gpp", &file_daily_gpp, 300, "Daily gross primary production");
    declare_parameter("file_daily_npp", &file_daily_npp, 300, "Daily net primary production");
    declare_parameter("file_daily_ra", &file_daily_ra, 300, "Daily autotrophic respiration");
    declare_parameter("file_daily_cflux", &file_daily_cflux, 300, "Daily heterotrophic respiration");
    declare_parameter("file_daily_fpar", &file_daily_fpar, 300, "Daily fraction of photosynthetic active radiation");
  }

  ProfoundOutput::~ProfoundOutput() {
  }
  void ProfoundOutput::init() {
    define_output_tables();
  }
  void ProfoundOutput::define_output_tables() {
    // create a vector with the pft names
    std::vector<std::string> pfts;
    pftlist.firstobj();
    while (pftlist.isobj) {
      Pft& pft=pftlist.getobj();

      pfts.push_back((char*)pft.name);
      pftlist.nextobj();
    }
    // pfts
    ColumnDescriptors pft_columns_diam;
    pft_columns_diam += ColumnDescriptor("Diam", 10, 2);
    pft_columns_diam += ColumnDescriptors(pfts, 10, 3);
    pft_columns_diam += ColumnDescriptor("Total", 10, 3);

    ColumnDescriptors pft_columns;
    pft_columns += ColumnDescriptors(pfts, 10, 3);
    pft_columns += ColumnDescriptor("Total", 10, 3);

    ColumnDescriptors pft_columns_wide;
    pft_columns_wide += ColumnDescriptors(pfts, 13, 6);
    pft_columns_wide += ColumnDescriptor("Total", 13, 6);

    ColumnDescriptors cflux_columns_wide;
    cflux_columns_wide += ColumnDescriptor("GPP", 13, 6);
    cflux_columns_wide += ColumnDescriptor("NPP", 13, 6);
    cflux_columns_wide += ColumnDescriptor("RA", 13, 6);
    cflux_columns_wide += ColumnDescriptor("RH", 13, 6);
    cflux_columns_wide += ColumnDescriptor("NEE", 13, 6);

    create_output_table(out_dbh_domhei, file_dbh_domhei, pft_columns);
    create_output_table(out_dom_height, file_dom_height, pft_columns);

    create_output_table(out_diamclass_dbh, file_diamclass_dbh, pft_columns_diam);
    create_output_table(out_diamclass_mai, file_diamclass_mai, pft_columns_diam);
    create_output_table(out_diamclass_ba, file_diamclass_ba, pft_columns_diam);
    create_output_table(out_diamclass_harv, file_diamclass_harv, pft_columns_diam);
    create_output_table(out_diamclass_stemno, file_diamclass_stemno, pft_columns_diam);
    create_output_table(out_diamclass_age, file_diamclass_age, pft_columns_diam);
    create_output_table(out_diamclass_mortstemno, file_diamclass_mortstemno, pft_columns_diam);
    create_output_table(out_diamclass_harvstemno, file_diamclass_harvstemno, pft_columns_diam);
    create_output_table(out_diamclass_dens, file_diamclass_dens, pft_columns_diam);
    create_output_table(out_diamclass_vol, file_diamclass_vol, pft_columns_diam);
    create_output_table(out_diamclass_cveg, file_diamclass_cveg, pft_columns_diam);

    // daily
    create_output_table(out_daily_gpp, file_daily_gpp, pft_columns_wide);
    create_output_table(out_daily_npp, file_daily_npp, pft_columns_wide);
    create_output_table(out_daily_ra, file_daily_ra, pft_columns_wide);
    create_output_table(out_daily_cflux, file_daily_cflux, cflux_columns_wide);
    create_output_table(out_daily_fpar, file_daily_fpar, pft_columns);
  }
  void ProfoundOutput::outdaily(Gridcell& gridcell) {
    if (vegmode == POPULATION)
      return;

    //if (date.year >= nyear_spinup) {
    if (date.get_calendar_year() >= plant_year) {
      
      double lon = gridcell.get_lon();
      double lat = gridcell.get_lat();

      double gc_gpp = 0.0;
      double gc_npp = 0.0;
      double gc_ra = 0.0;
      double gc_fpar = 0.0;

      OutputRows out(output_channel, lon, lat, date.get_calendar_year(), date.day + 1);

      pftlist.firstobj();
      while (pftlist.isobj) {
        Pft& pft = pftlist.getobj();

        double pft_gpp = 0.0;
        double pft_npp = 0.0;
        double pft_ra = 0.0;
        double pft_fpar = 0.0;

        Gridcellpft& gridcellpft = gridcell.pft[pft.id];

        Gridcell::iterator gc_itr = gridcell.begin();
        // Loop through Stands
        while (gc_itr != gridcell.end()) {
          Stand& stand = *gc_itr;

          // Loop through Patches
          stand.firstobj();
          while (stand.isobj) {
            Patch& patch = stand.getobj();

            Vegetation& vegetation = patch.vegetation;
            vegetation.firstobj();
            while (vegetation.isobj) {
              Individual& indiv = vegetation.getobj();

              if (indiv.id != -1 && indiv.alive) {
                if (indiv.pft.id==pft.id) {
                  pft_gpp  += indiv.gpp * 1000.0 / (double)stand.npatch(); // g
                  pft_npp  += indiv.npp * 1000.0 / (double)stand.npatch(); // g
                  pft_ra   += indiv.ra * 1000.0 / (double)stand.npatch();   // g
                  pft_fpar += indiv.fpar / (double)stand.npatch(); // g
                }
              }
              vegetation.nextobj();
            } // individuals
            stand.nextobj();
          } // patches
          ++gc_itr;
        } // stands

        out.add_value(out_daily_gpp, pft_gpp);
        out.add_value(out_daily_npp, pft_npp);
        out.add_value(out_daily_ra, pft_ra);
        out.add_value(out_daily_fpar, pft_fpar);

        gc_gpp += pft_gpp;
        gc_npp += pft_npp;
        gc_ra += pft_ra;
        gc_fpar += pft_fpar;

        pftlist.nextobj();
      } // pftlist
      //Totals
      out.add_value(out_daily_gpp, gc_gpp);
      out.add_value(out_daily_npp, gc_npp);
      out.add_value(out_daily_ra, gc_ra);
      out.add_value(out_daily_fpar, gc_fpar);

      out.add_value(out_daily_cflux, gc_gpp);
      out.add_value(out_daily_cflux, gc_npp);
      out.add_value(out_daily_cflux, gc_ra);

      double gc_rh = 0.0;

      Gridcell::iterator gc_itr = gridcell.begin();
      // Loop through Stands
      while (gc_itr != gridcell.end()) {
        Stand& stand = *gc_itr;

        // Loop through Patches
        stand.firstobj();
        while (stand.isobj) {
          Patch& patch = stand.getobj();

          gc_rh += patch.rh * 1000.0 / (double)stand.npatch(); // g

          stand.nextobj();
        }
        gc_itr++;
      }
      out.add_value(out_daily_cflux, gc_rh);
      out.add_value(out_daily_cflux, gc_npp - gc_rh);
    }
    return;
  }
  void ProfoundOutput::outannual(Gridcell& gridcell) {
    if (vegmode == POPULATION)
      return;

    //if (date.year >= nyear_spinup) {
    if (date.get_calendar_year() >= plant_year) {
    
      double lon = gridcell.get_lon();
      double lat = gridcell.get_lat();

      // list of growing 1D vectors pftlist.nobj nindiv
      std::map< std::string, std::vector<double> > dbh_domhei;
      std::map< std::string, std::vector<double> > dom_height; 

      // list of fixed vectors pftlist.nobj ndiamclass
      std::map< std::string, std::vector<double> > diamclass_dbh; 
      std::map< std::string, std::vector<double> > diamclass_mai;
      std::map< std::string, std::vector<double> > diamclass_ba;
      std::map< std::string, std::vector<double> > diamclass_harv;
      std::map< std::string, std::vector<double> > diamclass_stemno;
      std::map< std::string, std::vector<double> > diamclass_age;
      std::map< std::string, std::vector<double> > diamclass_mortstemno;
      std::map< std::string, std::vector<double> > diamclass_harvstemno;
      std::map< std::string, std::vector<double> > diamclass_dens;
      std::map< std::string, std::vector<double> > diamclass_vol;
      std::map< std::string, std::vector<double> > diamclass_cveg;

      OutputRows out(output_channel, lon, lat, date.get_calendar_year());
      std::vector<double> dbh_domhei_total(0);
      std::vector<double> dom_height_total(0);

      pftlist.firstobj();
      while (pftlist.isobj) {
        Pft& pft = pftlist.getobj();
        char *pftname = (char*) pft.name;
        dbh_domhei[pftname] = std::vector<double>(0);
        dom_height[pftname] = std::vector<double>(0);

        diamclass_dbh[pftname].assign(ndiamclass, 0.0);
        diamclass_mai[pftname].assign(ndiamclass, 0.0);
        diamclass_ba[pftname].assign(ndiamclass, 0.0);
        diamclass_harv[pftname].assign(ndiamclass, 0.0);
        diamclass_stemno[pftname].assign(ndiamclass, 0.0);
        diamclass_age[pftname].assign(ndiamclass, 0.0);
        diamclass_mortstemno[pftname].assign(ndiamclass, 0.0);
        diamclass_harvstemno[pftname].assign(ndiamclass, 0.0);
        diamclass_dens[pftname].assign(ndiamclass, 0.0);
        diamclass_vol[pftname].assign(ndiamclass, 0.0);
        diamclass_cveg[pftname].assign(ndiamclass, 0.0);

        Gridcellpft& gridcellpft = gridcell.pft[pft.id];

        Gridcell::iterator gc_itr = gridcell.begin();
        // Loop through Stands
        while (gc_itr != gridcell.end()) {
          Stand& stand = *gc_itr;

          // Loop through Patches
          stand.firstobj();
          while (stand.isobj) {
            Patch& patch = stand.getobj();
            Patchpft& patchpft = patch.pft[pft.id];

            std::transform(patchpft.diamclass_harvstemno.begin(), patchpft.diamclass_harvstemno.end(), patchpft.diamclass_harvstemno.begin(),
                           std::bind1st(std::multiplies<double>(), 10000.0 / patcharea / (double)stand.npatch())); // per ha
            std::transform(diamclass_harvstemno[pftname].begin(), diamclass_harvstemno[pftname].end(),
                           patchpft.diamclass_harvstemno.begin(), diamclass_harvstemno[pftname].begin(), std::plus<double>());

            std::transform(patchpft.diamclass_mortstemno.begin(), patchpft.diamclass_mortstemno.end(), patchpft.diamclass_mortstemno.begin(),
                           std::bind1st(std::multiplies<double>(), 10000.0 / patcharea / (double)stand.npatch())); // per ha
            std::transform(diamclass_mortstemno[pftname].begin(), diamclass_mortstemno[pftname].end(),
                           patchpft.diamclass_mortstemno.begin(), diamclass_mortstemno[pftname].begin(), std::plus<double>());

            // Loop through individuals
            Vegetation& vegetation = patch.vegetation;
            vegetation.firstobj();
            while (vegetation.isobj) {
              Individual& indiv = vegetation.getobj();
              if (indiv.id != -1 && indiv.alive) {
                if (indiv.pft.id == pft.id && indiv.age > 0) {
                  int nindiv    = (int)(ceil(indiv.densindiv * patcharea));
                  double diam   = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);
                  int diamclass = (int)(floor(diam / diamclasswidth));
                  double ba     = pow(diam / 2.0, 2.0) * PI * indiv.densindiv;
                  if (diamclass >= ndiamclass)
                    diamclass = ndiamclass - 1;
                  std::vector<double> indiv_values;
                  indiv_values.assign(nindiv, diam * 100.0); // cm
                  dbh_domhei[pftname].insert(dbh_domhei[pftname].end(), indiv_values.begin(), indiv_values.end());
                  indiv_values.assign(nindiv, indiv.height);
                  dom_height[pftname].insert(dom_height[pftname].end(), indiv_values.begin(), indiv_values.end());

                  diamclass_dbh[pftname][diamclass]    += 100.0 * diam * (double)nindiv * 10000.0 / patcharea / (double)stand.npatch(); // cm
                  diamclass_stemno[pftname][diamclass] += (double)nindiv / (double)stand.npatch() * 10.0;  // per ha
                  diamclass_ba[pftname][diamclass]     += ba / (double)stand.npatch() * 10000.0; // per ha
                  diamclass_dens[pftname][diamclass]   += indiv.densindiv / (double)stand.npatch() * 10000.0;  // per ha
                  diamclass_cveg[pftname][diamclass]   += indiv.cmass_wood() / (double)stand.npatch();
                  diamclass_age[pftname][diamclass]    += indiv.age * (double)nindiv / (double)stand.npatch() * 10.0;
                }
              } // alive check
              vegetation.nextobj();
            } // individuals
            stand.nextobj();
          } // patches
          ++gc_itr;
        } // stands
        std::sort(dbh_domhei[pftname].begin(), dbh_domhei[pftname].end());
        std::sort(dom_height[pftname].begin(), dom_height[pftname].end());
        std::reverse(dbh_domhei[pftname].begin(), dbh_domhei[pftname].end());
        std::reverse(dom_height[pftname].begin(), dom_height[pftname].end());
        dbh_domhei_total.insert(dbh_domhei_total.end(), dbh_domhei[pftname].begin(), dbh_domhei[pftname].end());
        dom_height_total.insert(dom_height_total.end(), dom_height[pftname].begin(), dom_height[pftname].end());

        if (dbh_domhei[pftname].size() > 0) {
          if (dbh_domhei[pftname].size() > 10)
            dbh_domhei[pftname].resize(10);
          out.add_value(out_dbh_domhei, std::accumulate(dbh_domhei[pftname].begin(), dbh_domhei[pftname].end(), 0.0) / dbh_domhei[pftname].size());
        } else {
          out.add_value(out_dbh_domhei, -999.999);
        }
        if (dom_height[pftname].size() > 0) { 
          if (dom_height[pftname].size() > 10)
            dom_height[pftname].resize(10);
          out.add_value(out_dom_height, std::accumulate(dom_height[pftname].begin(), dom_height[pftname].end(), 0.0) / dom_height[pftname].size());
        } else {
          out.add_value(out_dom_height, -999.999);
        }

        pftlist.nextobj();
      } // pftlist
      std::sort(dbh_domhei_total.begin(), dbh_domhei_total.end());
      std::sort(dom_height_total.begin(), dom_height_total.end());
      std::reverse(dbh_domhei_total.begin(), dbh_domhei_total.end());
      std::reverse(dom_height_total.begin(), dom_height_total.end());
      if (dbh_domhei_total.size() > 10)
        dbh_domhei_total.resize(10);
      if (dom_height_total.size() > 10)
        dom_height_total.resize(10);

      out.add_value(out_dbh_domhei, std::accumulate(dbh_domhei_total.begin(), dbh_domhei_total.end(), 0.0) / dbh_domhei_total.size());
      out.add_value(out_dom_height, std::accumulate(dom_height_total.begin(), dom_height_total.end(), 0.0) / dom_height_total.size());

      for (int i = 0; i < ndiamclass; i++) {
        OutputRows out_diam(output_channel, lon, lat, date.get_calendar_year());
        out_diam.add_value(out_diamclass_dbh,        100.0 * (i + 1) * diamclasswidth); // cm        
        out_diam.add_value(out_diamclass_ba,         100.0 * (i + 1) * diamclasswidth); // cm        
        out_diam.add_value(out_diamclass_stemno,     100.0 * (i + 1) * diamclasswidth); // cm
        out_diam.add_value(out_diamclass_age,        100.0 * (i + 1) * diamclasswidth); // cm
        out_diam.add_value(out_diamclass_mortstemno, 100.0 * (i + 1) * diamclasswidth); // cm
        out_diam.add_value(out_diamclass_harvstemno, 100.0 * (i + 1) * diamclasswidth); // cm
        out_diam.add_value(out_diamclass_dens,       100.0 * (i + 1) * diamclasswidth); // cm       
        out_diam.add_value(out_diamclass_cveg,       100.0 * (i + 1) * diamclasswidth); // cm

        double dbh;
        double dbh_total        = 0.0;
        double ba_total         = 0.0;
        double stemno_total     = 0.0;
        double harvstemno_total = 0.0;
        double mortstemno_total = 0.0;
        double dens_total       = 0.0;
        double cveg_total       = 0.0;

        pftlist.firstobj();
        while (pftlist.isobj) {
          Pft& pft = pftlist.getobj();
          char* pftname = (char*) pft.name;
  
          if (diamclass_stemno[pftname][i] > 0) {
            
            dbh = diamclass_dbh[pftname][i] / diamclass_stemno[pftname][i];
            dbh_total += diamclass_dbh[pftname][i];
          } else {
            dbh = -999.999;
          }
          out_diam.add_value(out_diamclass_dbh, dbh);
        //out.add_value(out_diamclass_mai, -999.9);
          out_diam.add_value(out_diamclass_ba, diamclass_ba[pftname][i]);
        //out.add_value(out_diamclass_harv, diamclass_harv[i]);
          out_diam.add_value(out_diamclass_stemno, diamclass_stemno[pftname][i]);
          out_diam.add_value(out_diamclass_age, diamclass_age[pftname][i] / diamclass_stemno[pftname][i]);
          out_diam.add_value(out_diamclass_harvstemno, diamclass_harvstemno[pftname][i]);
          out_diam.add_value(out_diamclass_mortstemno, diamclass_mortstemno[pftname][i]);
          out_diam.add_value(out_diamclass_dens, diamclass_dens[pftname][i]);
        //out.add_value(out_diamclass_vol, -999.9);
          out_diam.add_value(out_diamclass_cveg, diamclass_cveg[pftname][i]);

          ba_total         += diamclass_ba[pftname][i];
          stemno_total     += diamclass_stemno[pftname][i];
          harvstemno_total += diamclass_harvstemno[pftname][i];
          mortstemno_total += diamclass_mortstemno[pftname][i];
          dens_total       += diamclass_dens[pftname][i];
          cveg_total       += diamclass_cveg[pftname][i];

        pftlist.nextobj();
      } // pftlist
      out_diam.add_value(out_diamclass_dbh, dbh_total / stemno_total);
      out_diam.add_value(out_diamclass_ba, ba_total);
      out_diam.add_value(out_diamclass_stemno, stemno_total);
      out_diam.add_value(out_diamclass_harvstemno, harvstemno_total);
      out_diam.add_value(out_diamclass_mortstemno, mortstemno_total);
      out_diam.add_value(out_diamclass_age, -999.999);
      out_diam.add_value(out_diamclass_dens, dens_total);
      out_diam.add_value(out_diamclass_cveg, cveg_total);
    }

    }
    return;
  }
}
