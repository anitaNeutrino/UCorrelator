#ifndef UCORRELATOR_MAP_STAT_TOOLS_H
#define UCORRELATOR_MAP_STAT_TOOLS_H

/*  Tools written by John Russell <jwruss@hawaii.edu> to account for the implicit spherical geometry
 *  of ANITA interferometric maps.
 */


class TH2;


namespace UCorrelator {

	/*  Implementation of interferometric map mean which accounts for spherical geometry by default.
	 *  The parameter "sphWeight" is assuming spherical geometry by default, otherwise flat.
	 *  The parameter "antOffset" is to account for antenna offset in degrees.
	 */
	double getMapMean(const TH2 * interfMap, bool sphWeight = true, double antOffset = 10);

	/*  Implementation of interferometric map RMS which accounts for spherical geometry by default.
	 */
	double getMapRMS(const TH2 * interfMap, bool sphWeight = true, double antOffset = 10);

	/*  Implementation of interferometric map SNR which accounts for spherical geometry by default.
	 */
	double getMapSNR(const TH2 * interfMap, bool sphWeight = true, double antOffset = 10);

	/*  Implementation of interferometric map peak Z-score which accounts for spherical geometry by default.
	 */
	double getMapPeakZScore(const TH2 * interfMap, bool sphWeight = true, double antOffset = 10);
}


#endif
