
{
   ROOT::Math::GSLRandomEngine rng; 
   rng.Initialize(); 


   TH2D hist("hist","hist",50,-5,5,50,-5,5); 
   TFile f("gausquad.root","RECREATE"); 

   TTree * tree = new TTree("gausquad","gausquad"); 
   
   double sigma_x, sigma_y, x0, y0, rho, dx,dy; 
   dx = hist->GetXaxis()->GetBinWidth(1); 
   dy = hist->GetYaxis()->GetBinWidth(1); 

   double reco_sigma_x, reco_sigma_y, reco_x0, reco_y0, reco_rho, reco_val, val;; 


   tree->Branch("sigma_x",&sigma_x); 
   tree->Branch("sigma_y",&sigma_y); 
   tree->Branch("x0",&x0); 
   tree->Branch("y0",&y0); 
   tree->Branch("rho",&rho); 
   tree->Branch("dx",&dx); 
   tree->Branch("dy",&dy); 
   tree->Branch("val",&val); 
   tree->Branch("reco_sigma_x",&reco_sigma_x); 
   tree->Branch("reco_sigma_y",&reco_sigma_y); 
   tree->Branch("reco_x0",&reco_x0); 
   tree->Branch("reco_y0",&reco_y0); 
   tree->Branch("reco_rho",&reco_rho); 
   tree->Branch("reco_val",&reco_val); 
   
   for (int j = 0; j < 100; j++)
   {

     sigma_x = gRandom->Uniform(1.5,5); 
     sigma_y = gRandom->Uniform(1.5,5); 
     x0 = gRandom->Uniform(-1,1); 
     y0 = gRandom->Uniform(-1,1); 
     rho = gRandom->Uniform(-0.3,0.3); 
     val = gRandom->Uniform(100,500); 
     hist.Reset(); 

     for (int xi = 1;xi <= hist.GetNbinsX(); xi++) 
     {
       for (int yi = 1; yi <= hist.GetNbinsY(); yi++)
       {

         double x = hist.GetXaxis()->GetBinCenter(xi); 
         double y = hist.GetYaxis()->GetBinCenter(yi); 
         hist.SetBinContent(xi,yi,    val * exp(- (x-x0)*(x-x0)/(2*sigma_x*sigma_x * (1-rho*rho)) - (y-y0)*(y-y0) / (2*sigma_y * sigma_y * (1-rho*rho))  + rho * (x-x0)* (y-y0) / (sigma_x * sigma_y * (1-rho*rho)))); 
       }
     }


     UCorrelator::peakfinder::FineMaximum max; 
     UCorrelator::peakfinder::doPeakFindingQuadratic25(&hist,&max); 

     reco_sigma_x = max.sigma_x; 
     reco_sigma_y = max.sigma_y; 
     reco_x0 = max.x; 
     reco_y0 = max.y; 
     reco_rho = max.covar / (max.sigma_x * max.sigma_y); 
     reco_val = max.val; 
     /*
     printf("x: %f\n", max.x); 
     printf("y: %f\n", max.y); 
     printf("sigma_x: %f\n", max.sigma_x); 
     printf("sigma_y: %f\n", max.sigma_y); 
     printf("rho: %f\n", max.covar / (max.sigma_x * max.sigma_y)); 
     printf("val: : %f\n", max.val); 
     */

     tree->Fill(); 
   }

   tree->Write(); 

   f.Close(); 


}
