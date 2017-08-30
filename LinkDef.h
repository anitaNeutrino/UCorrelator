#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all namespaces;

#pragma link C++ namespace UCorrelator+;
#pragma link C++ namespace UCorrelator::gui+;
#pragma link C++ namespace UCorrelator::flags;
#pragma link C++ namespace UCorrelator::peakfinder;
#pragma link C++ namespace UCorrelator::shape;
#pragma link C++ namespace UCorrelator::spectrum;
#pragma link C++ namespace UCorrelator::image;

#pragma link C++ class UCorrelator::peakfinder::FineMaximum;
#pragma link C++ class UCorrelator::peakfinder::RoughMaximum;
#pragma link C++ class UCorrelator::AntennaPositions;
#pragma link C++ class UCorrelator::Baseline;
#pragma link C++ class UCorrelator::TimeDependentAverage;
#pragma link C++ class UCorrelator::Correlator;
#pragma link C++ class UCorrelator::ComplicatedNotchFilter;
#pragma link C++ class UCorrelator::SineSubtractFilter;
#pragma link C++ class UCorrelator::AdaptiveFilterAbby;
#pragma link C++ class UCorrelator::CombinedSineSubtractFilter;
#pragma link C++ class UCorrelator::AdaptiveMinimumPhaseFilter;
#pragma link C++ class UCorrelator::AdaptiveBrickWallFilter; 
#pragma link C++ class UCorrelator::AdaptiveButterworthFilter;
#pragma link C++ class UCorrelator::WaveformCombiner;
#pragma link C++ class UCorrelator::Analyzer;
#pragma link C++ class UCorrelator::AnalysisConfig;
#pragma link C++ class UCorrelator::PointingResolution;
#pragma link C++ class UCorrelator::ProbabilityMap;
#pragma link C++ class UCorrelator::PointingResolutionModel;
#pragma link C++ class UCorrelator::ConstantPointingResolutionModel;

#pragma link C++ class UCorrelator::gui::Map+; 
#pragma link C++ class UCorrelator::gui::SummaryText; 

#pragma link C++ class AnitaNoiseMachine+;



#endif

