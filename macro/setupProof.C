TProof * setupProof(int nworkers = 8) 
{

  TProof * proof = TProof::Open("",TString::Format("workers=%d",nworkers)); 

  // proof->Load("proofloader.C",true); 
  // proof->Load("cuts.C"); 
  
  return proof; 
}
