
double sum(int n,const double * W)
{
  double answer = 0;
  for (int i = 0; i < n; i++) answer += W[i];

  return answer;
}


void doCutTable(int N, TCut base, TChain & c, TCut * cuts, const char ** names, FILE * fout = stdout) 
{

  c.SetEstimate(c.GetEntries() * 10); 
  //need a dummy variable 
  int n = c.Draw("theta",base,"goff"); 

  double total = sum(n, c.GetW()); 

  fprintf(fout,"\\begin{tabular}{l|l|l|l}\n"); 
  fprintf(fout,"Cut & N After (Cumulative) & N If Only Cut & N If All But \\\\\n"); 
  fprintf(fout,"None & %0.2f & %0.2f &  \\\\\n", total,total); 
  

  TCut cumulative = base;

  for (int i = 0; i < N; i++) 
  {
    cumulative *= cuts[i]; 

    //only
    int nonly = c.Draw("theta", base * cuts[i],"goff"); 
    double only = sum(nonly, c.GetW()); 

    //cumulative
    int ncum = c.Draw("theta",cumulative,"goff"); 
    double cum = sum(ncum, c.GetW()); 

    //all but 
    
    TCut cut_all_but = base; 

    for (int j = 0; j < N; j++)
    {
      if (i==j) continue; 
      cut_all_but *= cuts[j]; 
    }
    
    int nbut = c.Draw("theta",cut_all_but,"goff"); 
    double but = sum(nbut, c.GetW()); 
    fprintf(fout,"%s & %0.2f & %0.2f & %0.2f \\\\\n",names[i], cum,only,but); 

  }

  fprintf(fout,"\\end{tabular}\n"); 

}

void makeMCCutTable(const char * addstr = "thermalTrees/simulated*.root", FILE * fout = stdout)
{

  TChain c("simulation"); 
  c.Add(addstr); 


  const char * names[] = { "MostImpulsivePointsMC", "TrigMostImpulsive","MaxPeak$<$1000","Blast","Thermal" }; 

  TCut base ("(isMostImpulsive * weight)"); 
  TCut points ( "pointsToMC"); 
  TCut trig (" ((iteration < 5 && HPolTrigger) || (iteration > 4 && VPolTrigger))" ); 
  TCut maxpeak (" (MaxPeak < 1000)" ); 
  TCut blast (" (!payloadBlast)" ); 
  TCut thermal (" (F > 3)" ); 

  TCut  cuts[] = {points,trig,maxpeak,blast, thermal}; 

  doCutTable(5, base, c, cuts, names,fout); 





}

void makeDataCutTable(FILE * fout = stdout)
{

  TChain c("anita3"); 
  c.Add("thermalTrees/a3all*.root"); 


  const char * names[] = {"TrigMostImpulsive","MaxPeak$<$1000","$\\chi^2 < 0.01$","Blast","Below Horizontal" }; 

  TCut base ("(isMostImpulsive)"); 
  TCut trig (" ((iteration < 5 && HPolTrigger) || (iteration > 4 && VPolTrigger))" ); 
  TCut maxpeak (" (MaxPeak < 1000)" ); 
  TCut blast (" (!payloadBlast)" ); 
  TCut down (" theta > 0" ); 
  TCut thermal (" (F > 3)" ); 

  TCut quality = trig * maxpeak * blast; 
  TCut  cuts[] = {quality, down, thermal}; 

  doCutTable(3, base, c, cuts, names, fout); 

}

void makeCutTable(bool mc = true) 
{
  if (mc) makeMCCutTable() ; 
  else makeDataCutTable(); 

}



