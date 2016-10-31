TProof * setupProof(int nworkers = 8) 
{

  TProof * proof = TProof::Open("",TString::Format("workers=%d",nworkers)); 

  proof->Load("macro/proofloader.C",true); 
  proof->Load("macro/cuts.C"); 
  
  return proof; 
}
