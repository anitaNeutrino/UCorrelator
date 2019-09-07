#include "UCMapStatTools.h"
#include "TH2.h"
#include <cmath>

#ifndef DEG2RAD
#define DEG2RAD M_PI / 180
#endif


double UCorrelator::getMapMean(const TH2 * interfMap, bool sphWeight, double antOffset) {

	int numBinsX = interfMap -> GetNbinsX();
	int numBinsY = interfMap -> GetNbinsY();

	double binYWtSum = 0, mapMean = 0;
	for (int binY = 1; binY <= numBinsY; ++binY) {  //  With histograms, bin = 0 is underflow, bin = nbins + 1 is overflow.

		double binXSum = 0;
		for (int binX = 1; binX <= numBinsX; ++binX) binXSum += interfMap -> GetBinContent(binX, binY);
		double binYWt = sphWeight ? cos(DEG2RAD * (interfMap -> GetYaxis() -> GetBinCenter(binY) + antOffset)) : 1;
		binYWtSum += binYWt;
		mapMean += binXSum * binYWt;
	}
	mapMean /= numBinsX * binYWtSum;

	return mapMean;
}


double UCorrelator::getMapRMS(const TH2 * interfMap, bool sphWeight, double antOffset) {

	int numBinsX = interfMap -> GetNbinsX();
	int numBinsY = interfMap -> GetNbinsY();

	double binYWtSum = 0, mapMean = 0, mapSqMean = 0;
	for (int binY = 1; binY <= numBinsY; ++binY) {

		double binXSum = 0, binXSqSum = 0;
		for (int binX = 1; binX <= numBinsX; ++binX) {

			double binValue = interfMap -> GetBinContent(binX, binY);
			binXSum += binValue;
			binXSqSum += binValue * binValue;
		}
		double binYWt = sphWeight ? cos(DEG2RAD * (interfMap -> GetYaxis() -> GetBinCenter(binY) + antOffset)) : 1;
		binYWtSum += binYWt;
		mapMean += binXSum * binYWt;
		mapSqMean += binXSqSum * binYWt;
	}
	double binArea = numBinsX * binYWtSum;
	mapMean /= binArea;
	mapSqMean /= binArea;

	return sqrt(mapSqMean - mapMean * mapMean);
}


double UCorrelator::getMapSNR(const TH2 * interfMap, bool sphWeight, double antOffset) {

	double mapPeak = interfMap -> GetBinContent(interfMap -> GetMaximumBin());
	double mapRMS = UCorrelator::getMapRMS(interfMap, sphWeight, antOffset);

	return mapPeak / mapRMS;
}


double UCorrelator::getMapPeakZScore(const TH2 * interfMap, bool sphWeight, double antOffset) {

	double mapPeak = interfMap -> GetBinContent(interfMap -> GetMaximumBin());

	int numBinsX = interfMap -> GetNbinsX();
	int numBinsY = interfMap -> GetNbinsY();

	double binYWtSum = 0, mapMean = 0, mapSqMean = 0;
	for (int binY = 1; binY <= numBinsY; ++binY) {

		double binXSum = 0, binXSqSum = 0;
		for (int binX = 1; binX <= numBinsX; ++binX) {

			double binValue = interfMap -> GetBinContent(binX, binY);
			binXSum += binValue;
			binXSqSum += binValue * binValue;
		}
		double binYWt = sphWeight ? cos(DEG2RAD * (interfMap -> GetYaxis() -> GetBinCenter(binY) + antOffset)) : 1;
		binYWtSum += binYWt;
		mapMean += binXSum * binYWt;
		mapSqMean += binXSqSum * binYWt;
	}
	double binArea = numBinsX * binYWtSum;
	mapMean /= binArea;
	mapSqMean /= binArea;

	return (mapPeak - mapMean) / sqrt(mapSqMean - mapMean * mapMean);
}
