void drawGroupings(const UCorrelator::ProbabilityMap * m) 
{
  TCanvas * c = 0;
  for (int level = 0; level < m->NLevels(); level++)
  {
    c = new TCanvas(TString::Format("c%d\n", level),TString::Format("Level = %g", m->getLevel(level)),800,800); 
    std::vector<double> counts(m->segmentationScheme()->NSegments());

    m->groupAdjacent(m->getWgtAboveLevel(level), 0, &counts[0]);


    for (int i = 0; i < counts.size(); i++) 
    {
      if (counts[i] <= 1.1) 
      {
        counts[i] = 0; 
      }
    }
    m->segmentationScheme()->Draw("colz",&counts[0]); 

    TH2 * h = (TH2*)  c->FindObject("tmp_copy"); 
    h->SetTitle(TString::Format("Multiplicity MahalanobisDistance = %g\n", m->getLevel(level))); 
    h->SetMaximum(100); 
    c->SetLogz(); 
    c->Print("groupings.gif+100"); 
  }

   c->Print("groupings.gif++"); 

}

