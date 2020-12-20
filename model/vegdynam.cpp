///////////////////////////////////////////////////////////////////////////////////////
/// \file vegdynam.cpp
/// \brief Vegetation dynamics and disturbance
///
/// \author Ben Smith
/// $Date: 2014-06-23 15:50:25 +0200 (Mo, 23 Jun 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
// DISTURBANCE
// Generic patch-destroying disturbance with a prescribed probability
// Internal function - do not call from framework

void disturbance(Patch& patch, double disturb_prob) {

	// DESCRIPTION
	// Destroys all biomass in a patch with a certain stochastic probability. 
	// Biomass enters the litter, which is not affected by the disturbance.
	// NB: cohort and individual mode only

	// INPUT PARAMETER
	// disturb_prob = the probability of a disturbance this year

	if (randfrac(patch.stand.seed)<disturb_prob) {

		Vegetation& vegetation = patch.vegetation;

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv = vegetation.getobj();

      // f_js_20170207+
      if (indiv.pft.lifeform == TREE) {
        Patchpft& patchpft = patch.pft[indiv.pft.id];
        double diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);
        int diamclass = (int)(floor(diam / diamclasswidth));
        if (diamclass >= ndiamclass)
          diamclass = ndiamclass - 1;
        patchpft.diamclass_mortstemno[diamclass] += indiv.densindiv * patcharea ;
        patchpft.cmort += indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_sap + indiv.cmass_heart - indiv.cmass_debt; // f_js_20170705
      }
      // f_js_20170207-
      
			indiv.kill();

			vegetation.killobj();
		}

		patch.disturbed = true;
		patch.age = 0;
	}

	else patch.disturbed = false;
}

// Thinning based on Diameter

void thinning(Stand& stand, Patch& patch) {
	if (plant_year == 0)
		return;
	if (date.get_calendar_year() <= plant_year)
		return;

	Vegetation& vegetation = patch.vegetation;

	const int npft = pftlist.nobj;
	std::vector<double> ba_total;
	std::vector<std::vector<double> > height;
	ba_total.resize(npft);
	height.resize(npft);
	std::fill(ba_total.begin(), ba_total.end(), 0);

	// calculate the total basal area per pft first and save all indiv.heights
	int ipft = 0;
	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();
		if ((date.get_calendar_year() - plant_year) % pft.thinning_interv == 0) {
			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();
				if (indiv.id != -1 && indiv.alive) {
					if (indiv.pft.id == pft.id) {
						double diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);
						ba_total[ipft] += indiv.densindiv * pow(diam / 2.0, 2.0) * PI;
						height[ipft].push_back(indiv.height);
					}
				}
				vegetation.nextobj();
			}

		}
		ipft++;
		pftlist.nextobj();
	}

	// now perform the thinning
	ipft = 0;
	double charvest_flux = 0.0; // SR_ needed for harvest flux
	double frac_left; // SR_
	bool killed;
	
	pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();
		if (pft.lifeform == TREE && (date.get_calendar_year() - plant_year) % pft.thinning_interv == 0) {
			vegetation.firstobj();
			while (vegetation.isobj) {
				Individual& indiv = vegetation.getobj();
				if (indiv.id != -1 && indiv.alive) {
					if (indiv.pft.id == pft.id) {
						// distinguish between thinning regimes
						
						double diam = pow(indiv.height / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3);

            // f_js_20170206+
            int diamclass = (int)(floor(diam / diamclasswidth));
            if (diamclass >= ndiamclass)
              diamclass = ndiamclass - 1;
            Patchpft& patchpft = indiv.patchpft();
            // f_js_20170206-	

						
						if (pft.thinning_regime == ABOVE && diam >= pft.min_diam) {
							
							// Target diameter type of continous cover forestry
							if (diam > pft.max_diam) {  // if diameter is above the max diameter = kill the object

								patch.pft[indiv.pft.id].litter_leaf += indiv.cmass_leaf*pft.harvest_frac;
								patch.pft[indiv.pft.id].litter_root += indiv.cmass_root*pft.harvest_frac;
								charvest_flux = (indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_sap + indiv.cmass_heart -
									indiv.cmass_debt)*pft.harvest_frac;
								

                // f_js_20170206
                patchpft.diamclass_harvstemno[diamclass] += indiv.densindiv * patcharea ;

								frac_left = 1 - pft.harvest_frac;

								indiv.densindiv *= frac_left;
								indiv.cmass_leaf *= frac_left;
								indiv.cmass_root *= frac_left;
								indiv.cmass_sap *= frac_left;
								indiv.cmass_debt *= frac_left;
								indiv.cmass_heart *= frac_left;

								allometry(indiv);
					
							}
							else {

								

								patch.pft[indiv.pft.id].litter_leaf += indiv.cmass_leaf*pft.thinning_frac;
								patch.pft[indiv.pft.id].litter_root += indiv.cmass_root*pft.thinning_frac;
								charvest_flux = (indiv.cmass_leaf + indiv.cmass_root + indiv.cmass_sap + indiv.cmass_heart -
									indiv.cmass_debt)*pft.thinning_frac;

                // f_js_20170206
                patchpft.diamclass_harvstemno[diamclass] += indiv.densindiv * patcharea * pft.thinning_frac;

								frac_left = 1- pft.thinning_frac;								

								indiv.densindiv *= frac_left;
								indiv.cmass_leaf *= frac_left;
								indiv.cmass_root *= frac_left;
								indiv.cmass_sap *= frac_left;
								indiv.cmass_debt *= frac_left;
								indiv.cmass_heart *= frac_left;
																
								allometry(indiv);
														

							}
							
						}
						else if (pft.thinning_regime == BELOW && diam <= pft.max_diam) {
							
							if (diam < pft.min_diam) {
								patch.pft[indiv.pft.id].litter_leaf += indiv.cmass_leaf*pft.thinning_frac;
								patch.pft[indiv.pft.id].litter_root += indiv.cmass_root*pft.thinning_frac;
								charvest_flux += (indiv.cmass_sap + indiv.cmass_heart -
									indiv.cmass_debt)*pft.thinning_frac;
								patch.pft[indiv.pft.id].litter_wood += (indiv.cmass_sap + indiv.cmass_heart -
									indiv.cmass_debt)*pft.thinning_frac; // Harvest

                // f_js_20170206
                patchpft.diamclass_harvstemno[diamclass] += indiv.densindiv * patcharea * pft.thinning_frac;

								frac_left = 1 - pft.thinning_frac;

								indiv.densindiv *= frac_left;
								indiv.cmass_leaf *= frac_left;
								indiv.cmass_root *= frac_left;
								indiv.cmass_sap *= frac_left;
								indiv.cmass_debt *= frac_left;
								indiv.cmass_heart *= frac_left;

								
								allometry(indiv);

								
								
							}
							
						}

						else if (pft.thinning_regime == NONE) { 					
							
							
						}
					}
				}
				if (negligible(indiv.densindiv)) {
					vegetation.killobj();
					killed = true;
				}
				vegetation.nextobj();
				
			}
			
		}
		ipft++;
		pftlist.nextobj();
	}
	
	patch.fluxes.report_flux(Fluxes::HARVESTC, charvest_flux);

	return;
}



void finalharvest(Stand& stand, Patch& patch) {
  int rotation_interv = 9999;
	pftlist.firstobj();
	while (pftlist.isobj) {
    Pft& pft = pftlist.getobj();
    rotation_interv = min(rotation_interv, pft.rotation_interv);
    pftlist.nextobj();
	}

  pftlist.firstobj();
	while (pftlist.isobj) {
		Pft& pft = pftlist.getobj();
		if (pft.lifeform == TREE && (date.get_calendar_year() - plant_year) % rotation_interv == 0) {
			disturbance(patch, 1.0);
		}
    pftlist.nextobj();
	}
  return;
}


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2001
// Fulton, MR 1991 Adult recruitment rate as a function of juvenile growth in size-
//   structured plant populations. Oikos 61: 102-105.
// Pacala SW, Canham, CD & Silander JA Jr 1993 Forest models defined by field
//   measurements: I. The design of a northeastern forest simulator. Canadian Journal
//   of Forest Research 23: 1980-1988.
// Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65:
//   51-70.
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Smith, B, Prentice, IC & Sykes, M (2001) Representation of vegetation dynamics in
//   the modelling of terrestrial ecosystems: comparing two contrasting approaches
//   within European climate space. Global Ecology and Biogeography 10: 621-637
// Thonicke, K, Venevsky, S, Sitch, S & Cramer, W (2001) The role of fire disturbance
//   for global vegetation dynamics: coupling fire into a Dynamic Global Vegetation
//   Model. Global Ecology and Biogeography 10: 661-677.
// Zwillinger, D 1996 CRC Standard Mathematical Tables and Formulae, 30th ed. CRC
//   Press, Boca Raton, Florida.