//Peng Cao
// Some temp inactive commands that are useful in ROOT interactive mode.
TCanvas *c1 = new TCanvas("c1","multipads",900,700);
//c1->Divide(1,1,0,0);
c1->Divide(2,2,0.02,0.02);

TH2F *h1 = new TH2F("h1","theta:phi",3000,0,360,3000,-90,90);
TH2F *h2 = new TH2F("h2","imagePeak:hilbertPeak",3000,0,500,3000,0,0.25);
TH1F *h3 = new TH1F("h3","peakPower histogram",108,-5,25);
TH1F *h4 = new TH1F("h4","peakFrequency histogram",108,0,1.3);
//Payload Blast:
c1->cd(1);
anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>h1","summary->flags.isPayloadBlast == 1","colz")

c1->cd(2);
anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h2","summary->flags.isPayloadBlast == 1","colz")


c1->cd(3);
anita4->Draw("summary->coherent[0][0].peakPower[0]>>h3","summary->flags.isPayloadBlast == 1","")
anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]","summary->flags.isPayloadBlast == 1","same*H")

c1->cd(4);
anita4->Draw("summary->coherent_filtered[0][0].peakFrequency[0]","summary->flags.isPayloadBlast == 1","*H")
anita4->Draw("summary->coherent[0][0].peakFrequency[0]>>h4","summary->flags.isPayloadBlast == 1","same")


//Satellite cut:

TCanvas *c1 = new TCanvas("c1","multipads",900,700);
//c1->Divide(1,1,0,0);
c1->Divide(2,2,0.02,0.02);

TH2F *h1 = new TH2F("h1","theta:phi",3000,0,360,3000,-90,90);
TH2F *h2 = new TH2F("h2","imagePeak:hilbertPeak",3000,0,500,3000,0,0.25);
TH1F *h3 = new TH1F("h3","peakPower histogram",108,-5,25);
TH1F *h4 = new TH1F("h4","peakFrequency histogram",108,0,1.3);
c1->cd(1);
anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>h1","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","colz");

c1->cd(2);
anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h2","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","colz");

c1->cd(3);
anita4->Draw("summary->coherent[0][0].peakPower[0]>>h3","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","");

anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","same*H");

c1->cd(4);

anita4->Draw("summary->coherent_filtered[0][0].peakFrequency[0]>>h4","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","*H");

anita4->Draw("summary->coherent[0][0].peakFrequency[0]","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","same");





sun cut:
c1->cd(1);
anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>h1","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","colz")

c1->cd(2);
anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h2","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","colz")

c1->cd(3);
anita4->Draw("summary->coherent[0][0].peakPower[0]>>h3","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","")

anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","same*H")

c1->cd(4);
anita4->Draw("summary->coherent_filtered[0][0].peakFrequency[0]>>h4(100,0,1.3)","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","*H")

anita4->Draw("summary->coherent[0][0].peakFrequency[0]","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","same")




wais cuts:

c1->cd(1);
anita4->Draw("-1*summary->peak[0][0].theta + summary->wais.theta:FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>>h1","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")
c1->cd(2);
anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h2","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")
c1->cd(3);

anita4->Draw("summary->coherent[0][0].peakPower[0]>>h3","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","")

anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","same*H")

c1->cd(4);

anita4->Draw("summary->coherent_filtered[0][0].peakFrequency[0]>>h4","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","*H")


anita4->Draw("summary->coherent[0][0].peakFrequency[0]","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","same")



Other:
not sun && not wais && not satellite:

c1->cd(1);
anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>h1","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1)","colz")



c1->cd(2);
anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h2","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1)","colz")
c1->cd(3);
anita4->Draw("summary->coherent[0][0].peakPower[0]>>h3","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1)","")

anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1)","same*H")

c1->cd(4);

anita4->Draw("summary->coherent_filtered[0][0].peakFrequency[0]>>h4","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1)","*H")

anita4->Draw("summary->coherent[0][0].peakFrequency[0]","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1)","same")



anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h(100,,,100,,)","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1) && (summary->coherent_filtered[0][0].peakPower>12)","colz")

anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>h(300,,,300,,)","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90)) &&(summary->flags.isPayloadBlast != 1)&& (summary->coherent_filtered[0][0].peakPower>12)","colz")


all events:

anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h(300,,,300,,)","(summary->coherent_filtered[0][0].peakPower>12)&&(summary->flags.isPayloadBlast != 1)","colz")

anita4->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].peakHilbert>>h(3000,,,3000,,)&&(summary->flags.isPayloadBlast != 1)","","colz")

anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>h(300,,,300,,)","(summary->coherent_filtered[0][0].peakPower>12)&&(summary->flags.isPayloadBlast != 1)","colz")



p2/p1:
satelite:
anita4->Draw("summary->peak[0][1].value/summary->peak[0][0].value>>h1(1000,,)","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","same")

sun:
anita4->Draw("summary->peak[0][1].value/summary->peak[0][0].value>>h2(1000,,)","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","same")

wais:
anita4->Draw("summary->peak[0][1].value/summary->peak[0][0].value>>h3(1000,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","same")

other:
anita4->Draw("summary->peak[0][1].value/summary->peak[0][0].value>>h4(1000,,)","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90))&&(summary->flags.isPayloadBlast != 1)","same")


all:
anita4->Draw("summary->peak[0][1].value/summary->peak[0][0].value>>h5(1000,,)","&&(summary->flags.isPayloadBlast != 1)","same")



PeakPower:
satelite:
anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]>>h1(100,,)","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","same")

sun:
anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]>>h2(100,,)","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","same")

wais:
anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]>>h3(100,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","same")

other:
anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]>>h4(100,,)","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90))&&(summary->flags.isPayloadBlast != 1)","same")


all:
anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]>>h5(100,,)","(summary->flags.isPayloadBlast != 1)","")




anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,,,108,,)","(summary->flags.isPayloadBlast != 1)","colz")
anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,,,108,,)","(summary->flags.isPayloadBlast != 1)","colz")
anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakPower[0]>>h(108,,,108,,)","(summary->flags.isPayloadBlast != 1)","colz")
anita4->Draw("summary->coherent_filtered[0][0].peakPower[0]:summary->coherent_filtered[0][0].peakFrequency[0]:summary->coherent_filtered[0][0].bandwidth[0]>>h(108,,,108,,)","(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,,,108,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")



bandwidth vs peakFrequency:

satelite:
anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,0,1.3,50,0,0.5)","(-1*summary->peak[0][0].theta>-5) &&(-1*summary->peak[0][0].theta<3) &&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>331)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<345)&&(summary->flags.isPayloadBlast != 1)","colz")

sun:
anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,0,1.3,50,0,0.5)","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","colz")

wais:
anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,0,1.3,50,0,0.5)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")

other:
anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,0,1.3,50,0,0.5)","!((FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-20)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<20))&&!((-1*summary->peak[0][0].theta + summary->wais.theta>-0.8)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.8)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<2)) &&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>324)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<352))&&!((FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>40)&&(FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<90))&&(summary->flags.isPayloadBlast != 1)","colz")


all:
anita4->Draw("summary->coherent_filtered[0][0].bandwidth[0]:summary->coherent_filtered[0][0].peakFrequency[0]>>h(108,0,1.3,50,0,0.5)","(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->realTime>>h(300,,,300,,)","(summary->flags.pulser == 1 ) && (FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-5)&& (FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<5) && (FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0)>-2) && (FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0)<2)","colz")

anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->realTime - 1481362282>>h(300,,,300,,)","(summary->flags[1].pulser == 1 )","colz")



anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")

anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0):summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0): FFTtools::wrap(pat->heading - summary->wais.phi,360,180)>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("FFTtools::wrap(pat->heading - summary->wais.phi,360,180):summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("pat->pitch:summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("pat->roll:summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("header->phiTrigMaskHOffline","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0):header->l3TrigPatternH>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("summary->peak[0][1].value:summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")


anita4->Draw("summary->peak[0][0].value:summary->realTime>>h(300,,,300,,)","summary->flags.pulser == 0","colz")
anita4->Draw("summary->eventNumber:summary->realTime>>h(300,,,300,,)","summary->flags.pulser == 0","colz")
anita4->DrawDerivative("summary->eventNumber:summary->realTime>>h(300,,,300,,)","summary->flags.pulser == 0","colz")
anita4->Draw("summary->realTime>>h(300,,)","summary->flags.pulser == 0","colz")
anita4->Draw("pat->heading:summary->realTime>>h(300,,,300,,)","summary->flags.pulser == 0","colz")
anita4->Draw("summary->wais.phi:summary->realTime>>h(300,,,300,,)","summary->flags.pulser == 0","colz")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","")


rough:
anita4->Draw("FFTtools::wrap(summary->peak[0][0].dphi_rough,360,0):summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")

anita4->Draw("FFTtools::wrap(summary->peak[0][0].dtheta_rough,360,0):summary->realTime>>h(300,,,300,,)","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")

anita4->Draw("summary->peak[0][0].dtheta_rough:summary->wais.phi","(-1*summary->peak[0][0].theta + summary->wais.theta>-0.4)&&(-1*summary->peak[0][0].theta + summary->wais.theta<0.4)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>-1)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)<1)&&(summary->flags.isPayloadBlast != 1)","colz")

anita4->Draw("summary->wais.phi:summary->realTime>>h(300,,,300,,)","","")

anita4->Draw("FFTtools::wrap(-1*summary->peak[0][0].theta + summary->wais.theta,360,0):FFTtools::wrap(summary->peak[0][0].phi-summary.wais.phi,360,0)>>h(300,,,300,,)","","")

anita4->Draw("FFTtools::wrap(-1*summary->peak[0][0].theta + summary->wais.theta,360,0):summary->realTime>>h(300,,,300,,)","","colz")
anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary.wais.phi,360,0):summary->realTime>>h(300,,,300,,)","","colz")

anita4->Draw("summary->channels[0][11].phi:summary->realTime>>h(300,,,300,,)","","")


anita4->Draw("summary->channels[0][20].rms:FFTtools::wrap(summary->wais.phi,360,0)>>h(300,,,300,,)","","colz")
anita4->Draw("summary->channels[0][20].rms:summary->realTime>>h(300,,,300,,)","","")

anita4->Draw("summary->channels[0][1].rms:FFTtools::wrap(summary->wais.phi,360,0)>>h(300,,,300,,)","","colz")
anita4->Draw("summary->channels[0][1].rms:summary->realTime>>h(300,,,300,,)","","")

anita4->Draw("summary->channels[0][1].avgPower:summary->realTime>>h(300,,,300,,)","","")
anita4->Draw("summary->channels[0][1].peakHilbert:summary->realTime>>h(300,,,300,,)","","")

anita4->Draw("summary->channels[0][].rms:FFTtools::wrap((2 -summary->channels[0][].phi)*22.5 + summary->wais.phi,360,0)>>h(300,,,300,,)","summary->channels[0][].channelNumber<16","colz")

anita4->Draw("summary->channels[0][].rms:FFTtools::wrap((2 -summary->channels[0][].phi)*22.5 + summary->wais.phi,360,0)>>h(300,,,300,,)","summary->channels[0][].channelNumber<16","",1,4)
anita4->Draw("summary->channels[0][].rms:FFTtools::wrap((2 -summary->channels[0][].phi)*22.5 + summary->wais.phi,360,0)","summary->channels[0][].channelNumber>=16 && summary->channels[0][].channelNumber<32","same",1,4)
anita4->Draw("summary->channels[0][].rms:FFTtools::wrap((2 -summary->channels[0][].phi)*22.5 + summary->wais.phi,360,0)","summary->channels[0][].channelNumber>=32 && summary->channels[0][].channelNumber<48","same",1,4)

anita4->Draw("summary->channels[0][].rms:summary->channels[0][].channelNumber>>h(300,,,300,,)","summary->channels[0][].channelNumber<16","", 1,2)

anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)","","");

root doTMVA_A4.C'(140,140,0,"max_300000_sinsub_10_3_ad_2",1)'


TestTree->Draw("Fisher>>h(300,-90,67)","classID==1","")
TestTree->Draw("Fisher","classID==0","same")

TrainTree->Draw("Fisher>>h(300,-90,67)","classID==1","")
TrainTree->Draw("Fisher","classID==0","same")


TestTree->Draw("kNN>>h(14,-0.2,1.2)","classID==1","")
TestTree->Draw("kNN","classID==0","same")

TrainTree->Draw("kNN>>h(14,-0.2,1.2)","classID==1","")
TrainTree->Draw("kNN","classID==0","same")



TestTree->Draw("MLP>>h(14,-0.2,1.2)","classID==1","")
TestTree->Draw("MLP","classID==0","same")

TrainTree->Draw("MLP>>h(14,-0.2,1.2)","classID==1","")
TrainTree->Draw("MLP","classID==0","same")

anita4->Draw("summary->realTime:summary->run:summary->eventNumber:pat->realTime:pat->longitude","","para")
anita4->Draw("summary->channels[][].ant:summary->channels[][].getPhi(summary->channels[][].ant)","","")

anita4->Draw("header->triggerTimeNs:summary->realTime","","")
anita4->Draw("header->triggerTimeNs - summary->wais.distance/3.0*10:summary->realTime","","")
anita4->Draw("header->triggerTimeNs - summary->wais.distance/3.0*10:summary->realTime","abs(header->triggerTimeNs - summary->wais.distance/3.0*10)<20000","")
anita4->Draw("header->triggerTimeNs - sqrt(pow(summary->wais.distance,2)+pow(pat->altitude,2))/3.0*10:sqrt(pow(summary->wais.distance,2)+pow(pat->altitude,2))","","")

anita4->Draw("header->triggerTimeNs - summary->wais.distance/3.0*10:summary->wais.distance/3.0*10","","")



anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->realTime - 1481362282>>h1(300,,,300,,)","(summary->flags.pulser == 1 )","colz")
h1->GetXaxis()->SetTitle("time/[s]")
h1->GetYaxis()->SetTitle("deltaPhi/[degree]")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0):summary->realTime - 1481362282>>h2(300,,,300,,)","(summary->flags.pulser == 1 )","colz")
h2->GetXaxis()->SetTitle("time/[s]")
h2->GetYaxis()->SetTitle("deltaTheta/[degree]")

anita4->Draw("FFTtools::wrap(summary->peak[1][0].theta-summary->wais.theta,360,0):summary->realTime - 1481362282>>h2(300,,,300,,)","(summary->flags.pulser == 5 )","colz")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0):summary->wais.phi","(summary->flags.pulser == 1 ) && (summary->realTime - 1481362282)<5000","colz")



anita4->Draw("FFTtools::wrap(summary->peak[1][0].phi-summary->wais.phi,360,0):summary->realTime - 1481362282>>h1","(summary->flags.pulser == 5 )","colz")
h1->GetXaxis()->SetTitle("time/[s]")
h1->GetYaxis()->SetTitle("deltaPhi/[degree]")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta + summary->peak[1][0].theta-2*summary->wais.theta,360,0):summary->realTime - 1481362282>>h2(300,,,300,,)","(summary->flags.pulser == 5 )","colz")


anita4->Draw("FFTtools::wrap(summary->peak[1][0].theta-summary->wais.theta,360,0):summary->realTime - 1481362282>>h2(300,,,300,,)","(summary->flags.pulser == 1 )","colz")
h2->GetXaxis()->SetTitle("time/[s]")
h2->GetYaxis()->SetTitle("deltaTheta/[degree]")


anita4->Draw("FFTtools::wrap(summary->wais.phi,360,0):summary->realTime - 1481362282>>h(300,,,300,,)","(summary->flags.pulser == 1 )","colz")
h->GetXaxis()->SetTitle("time/[s]")
h->GetYaxis()->SetTitle("deltaPhi/[degree]")


anita4->Draw("pat->pitch:summary->realTime - 1481362282>>pitch(300,,,300,,)","(summary->flags.pulser == 1 )","colz")
pitch->GetXaxis()->SetTitle("time/[s]")
pitch->GetYaxis()->SetTitle("pitch/[degree]")


anita4->Draw("pat->roll:summary->realTime - 1481362282>>roll(300,,,300,,)","(summary->flags.pulser == 1 )","colz")
roll->GetXaxis()->SetTitle("time/[s]")
roll->GetYaxis()->SetTitle("roll/[degree]")
anita4->Draw("pat->roll:pat->payloadTime - 1481362282 + payloadTimeUs/1000000>>roll(300,,,300,,)","(summary->flags.pulser == 1 )","colz")


anita4->Draw("summary->peak[0][0].value:summary->realTime - 1481362282","(summary->flags.pulser == 1 )","colz")

anita4->Draw("summary->realTime - 1481362282>>eventRate(7260,,)","(summary->flags.pulser == 1 )","")
anita4->Scan("summary->realTime - 1481362266>>eventRate(726,,)","","")

anita4->Draw("summary->wais.distance:summary->realTime - 1481362282","(summary->flags.pulser == 1 )","colz")
anita4->Draw("FFTtools::wrap(summary->wais.phi - pat->heading,360,0):summary->realTime - 1481362282","(summary->flags.pulser == 1 )","colz")
anita4->Draw("summary->wais.theta:summary->realTime - 1481362282","(summary->flags.pulser == 1 )","colz")
anita4->Draw("pat->altitude:summary->realTime - 1481362282","(summary->flags.pulser == 1 )","colz")


deconvImpulsivity @ deconvolved.impulsivityMeasure[][] 
deconvLinearPolFraction @ sqrt(TMath::Power(deconvolved[][].Q,2) + TMath::Power(deconvolved[][].U,2))/deconvolved[][].I  
deconvLinearPolAngle @ TMath::ATan2(deconvolved[][].U,deconvolved[][].Q)/90/TMath::Pi() 


anita4->Draw("summary->deconvolved.impulsivityMeasure[0][0]:summary->realTime - 1481362282","","colz")
anita4->Draw("sqrt(TMath::Power(summary->deconvolved[][].Q,2) + TMath::Power(summary->deconvolved[][].U,2))/summary->deconvolved[][].I:summary->realTime - 1481362282","","colz")


anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>skymap(1000,,,1000,,)","","colz");
skymap->GetXaxis()->SetTitle("phi - heading/[degree]")
skymap->GetYaxis()->SetTitle("theta/[degree]")


anita4->Draw("summary->coherent[0][0].peakPower[0]:summary->coherent[0][0].peakFrequency[0]","summary->flags.isPayloadBlast != 1","colz")
anita4->Draw("summary->coherent[0][0].peakPower[0]:summary->coherent[0][0].peakFrequency[0]","abs(-1*summary->peak[0][0].theta + 14)<5 && abs(FFTtools::wrap(summary->peak[0][0].phi-pat.heading) - 184)<9","colz");

anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>skymap(1000,,,1000,,)","abs(-1*summary->peak[0][0].theta + 14)<5 && abs(FFTtools::wrap(summary->peak[0][0].phi-pat.heading) - 184)<9","colz");


from south pole:
anita4->Draw("summary->peak[0][0].longitude:summary->peak[0][0].latitude>>skymap(1000,,,1000,,)","abs(-1*summary->peak[0][0].theta + 14)<5 && abs(FFTtools::wrap(summary->peak[0][0].phi-pat.heading) - 184)<9 && FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<182","colz");
skymap->GetXaxis()->SetTitle("latitude/[degree]")
skymap->GetYaxis()->SetTitle("longitude/[degree]")
anita4->Draw("summary->coherent[0][0].peakPower[0]:summary->coherent[0][0].peakFrequency[0]>>hist(300,,,300,,)","abs(-1*summary->peak[0][0].theta + 14)<5 && abs(FFTtools::wrap(summary->peak[0][0].phi-pat.heading) - 184)<9 && FFTtools::wrap(summary->peak[0][0].phi-pat.heading)<182","colz")
hist->GetXaxis()->SetTitle("peakFrequency/[GHz]")
hist->GetYaxis()->SetTitle("peakPower/[dB]")

near south pole
anita4->Draw("summary->peak[0][0].longitude:summary->peak[0][0].latitude>>skymap(1000,,,1000,,)","abs(-1*summary->peak[0][0].theta + 14)<5 && abs(FFTtools::wrap(summary->peak[0][0].phi-pat.heading) - 184)<9 && FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>182","colz");
anita4->Draw("summary->coherent[0][0].peakPower[0]:summary->coherent[0][0].peakFrequency[0]","abs(-1*summary->peak[0][0].theta + 14)<5 && abs(FFTtools::wrap(summary->peak[0][0].phi-pat.heading) - 184)<9 && FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>182","colz")

anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->realTime - 1481362282>>h1(300,,,300,,)","(summary->flags.pulser == 1 ) && abs(FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0))<4","colz")

anita4->Draw("-1*summary->peak[0][0].theta + summary->wais.theta:FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0)>>h1","","colz")



anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->wais.phi>>h1(300,,,300,,)","","colz")
anita4->Draw("summary->wais.phi:summary->realTime","","")
anita4->Draw("pat->heading:summary->realTime + pat->payloadTimeUs/1000000 - 1481362282","","")
anita4->Draw("summary->peak[0][0].value:summary->realTime + summary->payloadTimeUs/1000000 - 1481362282","","")

anita4->Draw("summary->realTime + pat->payloadTimeUs/1000000 -  1481362266>>eventRate(7260,,)","","")
anita4->Draw("summary->eventNumber:summary->realTime + pat->payloadTimeUs/1000000 -  1481362266","","")


anita4->Draw("summary->coherent[1][0].Q:summary->realTime","abs(summary->coherent[1][0].Q)<200000","")

anita4->Draw("summary->coherent[0][0].Q:summary->realTime","abs(summary->coherent[0][0].Q)<200000","")


adu5PatTree->Draw("pat->pitch:pat->payloadTime + pat->payloadTimeUs/1000000 - 1481362282","","")

adu5PatTree->Draw("pat->heading:pat->payloadTime + pat->payloadTimeUs/1000000 - 1481362282","","")
adu5PatTree->Draw("pat->roll:pat->payloadTime + pat->payloadTimeUs/1000000 - 1481362282","","")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->realTime - 1481362282 >>h1(300,,,1000,-0.6,0.6)","summary->flags.pulser == 1","colz")
anita4->Draw("FFTtools::wrap(summary->peak[1][0].phi-summary->wais.phi,360,0):summary->realTime - 1481362282 >>h1(300,,,1000,-0.6,0.6)","summary->flags.pulser == 5","colz")


anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0):summary->realTime - 1481362282 >>h1(300,,,1000,-0.6,0.6)","summary->flags.pulser == 1","colz")

anita4->Draw("FFTtools::wrap(summary->peak[1][0].theta-summary->wais.theta,360,0):summary->realTime - 1481362282 >>h1(300,,,1000,-0.6,0.6)","summary->flags.pulser == 5","colz")

anita4->Draw("FFTtools::wrap(summary->peak[0][0].phi-summary->wais.phi,360,0):summary->wais.phi ","summary->flags.pulser == 1","")
anita4->Draw("FFTtools::wrap(summary->peak[1][0].phi-summary->wais.phi,360,0):summary->wais.phi ","summary->flags.pulser == 5","")



gStyle->SetOptFit();
TF1 *fitCos = new TF1("fitCos","[0]*cos((x - [1])/180.0*3.141) + [2]", 0, 360);
fitCos->SetParName(0,"tilt");
fitCos->SetParName(1,"azimuth");
fitCos->SetParName(2,"constant");
fitCos->SetParameter(0, -1);
fitCos->SetParameter(1, 0);
fitCos->SetParameter(2, 0);

anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0):summary->wais.phi >>h1(100,,,300,,)","summary->flags.pulser == 1","colz")
h1->GetXaxis()->SetTitle("wais.phi/[degree]")
h1->GetYaxis()->SetTitle("deltaTheta/[degree]")

h1->FitSlicesY();
h1_1->Draw("")
h1_1->Fit("fitCos");
h1_1->GetXaxis()->SetTitle("wais.phi/[degree]")
h1_1->GetYaxis()->SetTitle("deltaTheta/[degree]")

anita4->Draw("FFTtools::wrap(summary->peak[1][0].theta-summary->wais.theta,360,0):summary->wais.phi >>h2(100,,,300,,)","summary->flags.pulser == 5","colz")
h2->GetXaxis()->SetTitle("wais.phi[degree]")
h2->GetYaxis()->SetTitle("deltaTheta/[degree]")

h2->FitSlicesY();
h2_1->Draw("")
h2_1->Fit("fitCos");
h2_1->GetXaxis()->SetTitle("wais.phi/[degree]")
h2_1->GetYaxis()->SetTitle("deltaTheta/[degree]")


gStyle->SetOptFit();
TF1 *fitCos = new TF1("fitCos","[0]*cos((x - [1])/180.0*3.141) + [2]", 0, 360);
fitCos->SetParName(0,"tilt");
fitCos->SetParName(1,"azimuth");
fitCos->SetParName(2,"constant");
fitCos->SetParameter(0, -1);
fitCos->SetParameter(1, 0);
fitCos->SetParameter(2, 0);
thetaPhi->FitSlicesY();
thetaPhi_1->Draw("");
thetaPhi_1->Fit("fitCos");
thetaPhi_1->GetXaxis()->SetTitle("wais.phi/[degree]");
thetaPhi_1->GetYaxis()->SetTitle("deltaTheta/[degree]");

sun cut : 
anita4->Draw("FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0):FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0) >>h1(300,,,300,,)","summary->flags.pulser == 0 && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)>-6) && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)<4) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)>-2) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)<2)","colz")

h1->GetXaxis()->SetTitle("deltaTheta/[degree]")
h1->GetYaxis()->SetTitle("deltaPhi/[degree]")

anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->sun.theta,360,0):FFTtools::wrap(summary->sun.phi,360,180) >>h1(100,,,100,,)","summary->flags.pulser == 0 && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)>-6) && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)<4) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)>-2) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)<2)","colz")

anita4->Draw("FFTtools::wrap(summary->peak[1][0].theta-summary->sun.theta,360,0):FFTtools::wrap(summary->sun.phi,360,180) >>h2(100,,,100,,)","summary->flags.pulser == 0 && (FFTtools::wrap(summary->peak[1][0].phi - summary->sun.phi,360,0)>-6) && (FFTtools::wrap(summary->peak[1][0].phi - summary->sun.phi,360,0)<4) && (FFTtools::wrap(-1*summary->peak[1][0].theta + summary->sun.theta,360,0)>-2) && (FFTtools::wrap(-1*summary->peak[1][0].theta + summary->sun.theta,360,0)<2)","colz")



anita4->Draw("FFTtools::wrap(summary->peak[0][0].theta - summary->sun.theta,360,0): FFTtools::wrap(summary->sun.phi,360,180)>>h3(100,,,300,,)","summary->flags.pulser == 0 && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)>-6) && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)<4) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)>-2) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)<2)","")



anita4->Draw("summary->channels[0][11].rms:summary->realTime>>h(3000,,,300,,)","summary->flags.pulser == 0 && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)>-6) && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)<4) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)>-2) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)<2)","")

anita4->Draw("summary->channels[0][11].avgPower: FFTtools::wrap(summary->sun.phi - (summary->channels[0][11].getPhi() - 2)*22.5,360,0)>>hist4(100,,,300,,)","summary->flags.pulser == 0 && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)>-6) && (FFTtools::wrap(summary->peak[0][0].phi - summary->sun.phi,360,0)<4) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)>-2) && (FFTtools::wrap(-1*summary->peak[0][0].theta + summary->sun.theta,360,0)<2)","col")

anita4->Draw("summary->channels[0][11].avgPower: FFTtools::wrap(summary->wais.phi - (summary->channels[0][11].getPhi() - 2)*22.5,360,0)>>h5(100,,,300,,)","summary->flags.pulser == 1 ","")

h->ProjectionX()->Draw()

anita4->Draw("FFTtools::wrap(summary->sun.phi - (summary->channels[0][11].getPhi() - 2)*22.5,360,0)","","")
anita4->Draw("pat->heading","","")


wana projection?
h2->ProjectionX()->Draw()
h->ProjectionY()->Draw()


anita4->Draw("summary->channels[0][].avgPower : 10 - summary->wais.theta: FFTtools::wrap(summary->wais.phi - (summary->channels[0][].getPhi() - 2)*22.5,360,0) >>dTheta_dPhi_Power_H(100,,,100,,,100,,)","summary->flags.pulser == 1 ","SURF3")
dTheta_dPhi_Power_H->FitSlicesZ()
dTheta_dPhi_Power_H_1->Draw("colz")

pl->Project3DProfile("zy")

TF1 *fit20 = new TF1("fit20","0.93*pow(cos((x - [1])/360.0*3.141),18)", -100, 100);


gStyle->SetOptFit();
TF1 *fit2 = new TF1("fit2","exp(-1*pow((x - [0])/[1], 2)/2.0)", -10, 10);
fit2->SetParName(0,"amplitude");
fit2->SetParName(1,"azimuth");
fit2->SetParameter(0, 0);
fit2->SetParameter(1, 20);
Phi_Power_H->GetXaxis()->SetTitle("dphi/[degree]")
Phi_Power_H->GetYaxis()->SetTitle("power")
Phi_Power_H_1->Fit("fit2");



Phi_Power_H->FitSlicesY();
Phi_Power_H_1->Draw("same")
Phi_Power_H_1->Fit("fit3");
Phi_Power_H_1->GetXaxis()->SetTitle("dphi/[degree]")
Phi_Power_H_1->GetYaxis()->SetTitle("power")



corrTree->Draw("expectedDeltaT[] - maxCorTimes[]: FFTtools::wrap((phiWave - (fAntPhi[][0] + fAntPhi[][1])/2.0)/3.14*180, 360 , 0)", " abs(expectedDeltaT[] - maxCorTimes[]) < 0.2  && abs(FFTtools::wrap((phiWave - (fAntPhi[][0] + fAntPhi[][1])/2.0)/3.14*180, 360 , 0))<100" , "colz")  


wais->Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0):summary->wais.phi >>h1(300,,,300,,)","summary->flags.pulser == 1","colz")
wais->Draw("FFTtools::wrap(summary->peak[1][0].theta-summary->wais.theta,360,0):summary->wais.phi >>h1(300,,,300,,)","summary->flags.pulser == 5","colz")



wais->Draw("FFTtools::wrap(summary->peak[0][2].theta,360,0):summary->wais.phi >>h1(300,,,300,,)","summary->flags.pulser == 1","colz")


wais->Draw("summary->mostImpulsivePol(1)","","colz")

TestTree->Draw("expectedDeltaT - MLPBFGS")

TrainTree->Draw("expectedDeltaT - maxCorTimes")
TrainTree->Draw("MLPBFGS - maxCorTimes","firstAnt == 1 && secondAnt == 16")

corrTree->Draw("expectedDeltaT-maxCorTimes","pol == 0 && firstAnt == 1 && secondAnt == 16 && abs(phiWave - fAntPhi[][])<22.5 * TMath::DegToRad() && (maxCorTimes - expectedDeltaT) * (maxCorTimes - expectedDeltaT) < 1","")

corrTree->Draw("expectedDeltaT-maxCorTimes","firstAnt == 1"）

TrainTree->Draw("MLPBFGS - maxCorTimes_0_>>h(10000)","firstAnt_0_ ==0 && secondAnt_0_ == 16")
TestTree->Draw("MLPBFGS - maxCorTimes")
TestTree->Draw("expectedDeltaT - maxCorTimes")

TrainTree->Draw("MLPBFGS - maxCorTimes")
TrainTree->Draw("expectedDeltaT - maxCorTimes")
TrainTree->Draw("MLPBFGS - maxCorTimes : thetaWave*180/3.14 >>h(300,,,300,,)", "", "colz")
TrainTree->Draw("expectedDeltaT - maxCorTimes : thetaWave*180/3.14 >>h(300,,,300,,)", "", "colz")

TrainTree->Draw("MLPBFGS - maxCorTimes : phiWave*180/3.14 >>h(300,,,300,,)", "", "colz")
TrainTree->Draw("expectedDeltaT - maxCorTimes : phiWave*180/3.14 >>h(300,,,300,,)", "", "colz")

TrainTree->Draw("thetaWave*60:phiWave*60")

pairCorTree->Draw("expectedDeltaT - maxCorTimes:ant2:ant1>>h(48,,,48,,,1000,,)","pol==0 && isValid","")
h->FitSlicesZ()
h_1->Draw("colz")


pairCorTree->Draw("expectedDeltaT - maxCorTimes: FFTtools::wrap(phiWave*180/3.14 - ((ant1+ant2)/2.0 - 2)*22.5, 180,0)>>hist(1000,,,1000,,)","pol==1&&isValid","colz")

dTheta vs time:

wais->Draw("FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0):realTime>>h1(1000,,,1000,-2,2)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0))<2","colz")
wais->Draw("FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0):realTime>>h1(1000,,,1000,-2,2)","flags.pulser == 1 && abs(FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0))<2","colz")
wais->Draw("FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0):realTime>>h1(1000,,,1000,-2,2)","flags.pulser == 5 && abs(FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0))<2","colz")
wais->Draw("FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0):realTime>>h1(1000,,,1000,-2,2)","flags.pulser == 5 && abs(FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0))<2","colz")

dPhi vs time:
wais->Draw("FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0):realTime>>h1(1000,,,300,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")
wais->Draw("FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0):realTime>>h1(1000,,,300,-5,5)","flags.pulser == 1&& abs(FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0))<5","colz")
wais->Draw("FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0):realTime>>h1(1000,,,300,-5,5)","flags.pulser == 5&& abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")
wais->Draw("FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0):realTime>>h1(1000,,,300,-5,5)","flags.pulser == 5 && abs(FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0))<5","colz")


wais->Draw("FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0):realTime>>h1(1000,,,10000,,)","flags.pulser == 1 && mostImpulsivePol() == 1","colz")
wais->Draw("peak[1][0].theta:realTime>>h1(1000,,,1000,,)","flags.pulser == 1 ","colz")



corrTree->Draw("thetaWave*180/3.14:phiWave*180/3.14","","")
wais->Draw("FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0):FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0)>>h(300,,,300,,)","flags.pulser == 1 ","colz")
wais->Draw("FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0):FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0)>>h(300,,,300,,)","flags.pulser == 5 ","colz")
wais->Draw("FFTtools::wrap(wais.theta, 360, 0):realTime>>h1(1000,,,10000,,)","flags.pulser == 1","colz")

testAllPair->Draw("thetaWave*180/3.14:phiWave*180/3.14:expectedDeltaT - MLP>>t1()","ant1 == 0 && ant2 == 32 && pol == 0","colz")
testAllPair->Draw("thetaWave*180/3.14:phiWave*180/3.14:expectedDeltaT>>t2()","ant1 == 0 && ant2 == 32 && pol == 0","colz")
testAllPair->Draw("thetaWave*180/3.14:phiWave*180/3.14:MLP>>t2()","ant1 == 0 && ant2 == 32 && pol == 0","colz")





pairCorTree->Draw("maxCorTimes - expectedDeltaT: FFTtools::wrap(phiWave*180/3.14 - ((ant1+ant2)/2.0 - 2)*22.5, 180,0)>>g(300,,,1000,,)","pol==0&& isValid && abs(expectedDeltaT - maxCorTimes)<0.6","colz")
g->FitSlicesY()
g_1->Draw("same")


gStyle->SetOptFit();
TF1 *ploySym = new TF1("ploySym","[0] + [1]*pow(x , 2) + [2]*pow(x, 4) + [3]*pow(x,6)", -44.5, 44.5);
ploySym->SetParName(0,"0");
ploySym->SetParName(1,"2");
ploySym->SetParName(2,"3");
ploySym->SetParName(3,"4");
ploySym->SetParameter(0, 0);
ploySym->SetParameter(1, 0);
ploySym->SetParameter(2, 0);
ploySym->SetParameter(3, 1);




pairCorTree->Draw("maxCorTimes - expectedDeltaT: FFTtools::wrap(phiWave*180/3.14 - ((ant1+ant2)/2.0 - 2)*22.5, 180,0)>>g1(300,,,1000,,)","pol==0 && abs(expectedDeltaT - maxCorTimes)<0.6","colz")
g1->FitSlicesY()
g1_1->Draw("same")



pairCorTree->Draw("maxCorTimes - expectedDeltaT: FFTtools::wrap(phiWave*180/3.14 - ((ant1+ant2)/2.0 - 2)*22.5, 180,0)>>g2(300,,,1000,,)","pol==1 && isValid && abs(expectedDeltaT - maxCorTimes)<0.6","colz")
g2->FitSlicesY()
g2_1->Draw("same")


pairCorTree->Draw("maxCorTimes - expectedDeltaT: FFTtools::wrap(phiWave*180/3.14 - ((ant1+ant2)/2.0 - 2)*22.5, 180,0)>>g3(300,,,1000,,)","pol==1 && abs(expectedDeltaT - maxCorTimes)<0.6","colz")
g3->FitSlicesY()
g3_1->Draw("same")

pairCorTree->Draw(" pow(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0),2) + pow(thetaWave*180/3.14 -10, 2) : pow(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0),2) + pow(thetaWave*180/3.14 -10, 2):expectedDeltaT - maxCorTimes>>yy()","pol==0 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")


pairCorTree->Draw(" pow(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0),2) : pow(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0),2):(expectedDeltaT - maxCorTimes)>>yy0()","pol==0 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")

pairCorTree->Draw(" pow(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0),2) - pow(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0),2):(expectedDeltaT - maxCorTimes)>>yy1()","pol==0 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")

pairCorTree->Draw(" (expectedDeltaT - maxCorTimes):pow(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0),2) - pow(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0),2)>>yy2()","pol==0 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")


pairCorTree->Draw(" (expectedDeltaT - maxCorTimes):abs(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0)) - abs(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0))>>yy3()","pol==0 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")
yy3->FitSlicesY()
yy3_1->Draw("same")


pairCorTree->Draw(" (expectedDeltaT - maxCorTimes):abs(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0)) - abs(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0))>>yy4()","pol==1 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")
yy4->FitSlicesY()
yy4_1->Draw("same")


pairCorTree->Draw(" (expectedDeltaT - maxCorTimes):abs(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0)) - abs(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0))>>yy4()","pol==1 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")


pairCorTree->Draw(" theta:abs(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0)) - abs(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0))>>yy4()","pol==1 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")


pairCorTree->Draw(" (expectedDeltaT - maxCorTimes):abs(FFTtools::wrap(phiWave*180/3.14 - (ant1 - 2)*22.5, 360,0)) - abs(FFTtools::wrap(phiWave*180/3.14 - (ant2 - 2)*22.5, 360,0))>>yy4()","ant1==0&& ant2 == 2 &&pol==0 && isValid && abs(expectedDeltaT - maxCorTimes)<0.4","colz")



pairCorTree->Draw("maxCorTimes - expectedDeltaT: phiWave*180/3.14 >>g(300,,,1000,,)","ant1<16 && ant2 < 16 &&pol==0&& isValid && abs(expectedDeltaT - maxCorTimes)<0.6","colz")

corrTree->Draw("expectedDeltaT[] - maxCorTimes[]: FFTtools::wrap((phiWave - (fAntPhi[][0] + fAntPhi[][1])/2.0)/3.14*180, 360 , 0)", " abs(expectedDeltaT[] - maxCorTimes[]) < 0.2  && abs(FFTtools::wrap((phiWave - (fAntPhi[][0] + fAntPhi[][1])/2.0)/3.14*180, 360 , 0))<100" , "colz")


anita4->Draw("-1*summary->peak[0][0].theta:FFTtools::wrap(summary->peak[0][0].phi-pat.heading)>>h(1000,,,1000,,)","(summary->flags.isPayloadBlast != 1)","colz")  
anita4->Draw("summary->peak[0][0].longitude:summary->peak[0][0].latitude>>skymap(1000,,,1000,,)","abs(summary->peak[0][0].longitude)<360 && abs(summary->peak[0][0].latitude) < 360  && flags.pulser == 1","colz");

corrTree->Draw("expectedDeltaT[] - maxCorTimes[]: FFTtools::wrap((phiWave - (fAntPhi[][0] + fAntPhi[][1])/2.0)/3.14*180, 360 , 0)", " abs(expectedDeltaT[] - maxCorTimes[]) < 0.2  && abs(FFTtools::wrap((phiWave - (fAntPhi[][0] + fAntPhi[][1])/2.0)/3.14*180, 360 , 0))<100" , "colz")  

corrTree->Draw("expectedDeltaT[] - maxCorTimes[]", "abs(expectedDeltaT[] - maxCorTimes[])<1" , "colz")

corrTree->Draw("expectedDeltaT[] - maxCorTimes[]:phiWave >> h(1000,,,1000,,)", "abs(expectedDeltaT[] - maxCorTimes[])<1 && pol == 0" , "colz")

corrTree->Draw("expectedDeltaT[] - maxCorTimes[]:phiWave >> h(1000,,,1000,,)", "firstAnt == 0 && secondAnt ==32  &&firstAnt ! = 45 && secondAnt != 45 && abs(expectedDeltaT[] - maxCorTimes[])<1 && pol == 0" , "colz")

TF1 *offAxisFit2 = new TF1("offAxisFit2","[0]*(x+11.25)^2 + [1]*(x+11.25)^4 + [2]*(x+11.25)^6 + [3]*(x+11.25)^8 -[0]*(x-11.25)^2 - [1]*(x-11.25)^4 - [2]*(x-11.25)^6 - [3]*(x-11.25)^8",-45,45);
TF1 *offAxisFit1 = new TF1("offAxisFit1","[0]*(x+22.5)^2 + [1]*(x+22.5)^4 + [2]*(x+22.5)^6 + [3]*(x+22.5)^8 -[0]*(x-22.5)^2 - [1]*(x-22.5)^4 - [2]*(x-22.5)^6 - [3]*(x-22.5)^8",-33.75,33.75);



OffAxis fit 1 phi sectors:


GroupDelayvsPhiAll->Draw("colz")
GroupDelayvsPhiAll_1->Draw("same")
gStyle->SetOptFit();
TF1 *offAxisFit3 = new TF1("offAxisFit3","[0]*(x+11.25)^2 + [1]*(x+11.25)^4 -[0]*(x-11.25)^2 - [1]*(x-11.25)^4",-45,45);
offAxisFit3->SetParNames("c2","c4");
offAxisFit3->SetParameters(0,0);

GroupDelayvsPhiAll->Fit("offAxisFit3")
GroupDelayvsPhiAll->GetYaxis()->SetTitle("#Delta#Deltat/[ns]");

OffAxis fit 2 phi sectors:

GroupDelayvsPhiAll->Draw("colz")
GroupDelayvsPhiAll_1->Draw("same")
gStyle->SetOptFit();
TF1 *offAxisFit0 = new TF1("offAxisFit0","[0]*(x+22.5)^2 + [1]*(x+22.5)^4 -[0]*(x-22.5)^2 - [1]*(x-22.5)^4 ",-33.75,33.75);
offAxisFit0->SetParNames("c2","c4");
offAxisFit0->SetParameters(0,0);

GroupDelayvsPhiAll->Fit("offAxisFit0")
GroupDelayvsPhiAll->GetYaxis()->SetTitle("#Delta#Deltat/[ns]");


wais->Draw("FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0):run","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")


wais->Draw("realTime","flags.pulser == 1","colz")


wais->Draw("coherent_filtered[0][0].spectrumIntercept:run","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")

wais->Draw("FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0):notchFreq >>h4(100,0.2,1.2,100,,)","flags.pulser == 5","colz")
wais->Draw("FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0):notchFreq >>h2(100,0.2,1.2,100,,)","flags.pulser == 1","colz")


wais->Draw("FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0):realTime","flags.pulser == 1","")

wais->Draw("FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0):coherent[0][0].linearPolAngle()>>h(300,,,3000,,)","flags.pulser==1","colz")
wais->Draw("FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0):coherent[1][0].linearPolAngle()>>h(300,,,3000,,)","flags.pulser==5","colz")

wais->Draw("FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0):wais.distance>>h(500,,,3000,,)","flags.pulser==1","colz")
wais->Draw("FFTtools::wrap(peak[1][0].phi - wais.phi, 360, 0):wais.distance>>h(500,,,3000,,)","flags.pulser==5","colz")




noOffAxis 30013:

wais->Draw("FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0):peak[0][0].snr>>Hpol_noOffAxisDelay(30,,,100,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0))<5","colz")
Hpol_noOffAxisDelay->FitSlicesY(0,0,-1,10,"QNRSG5")
Hpol_noOffAxisDelay_2->Draw()
Hpol_noOffAxisDelay_2->GetYaxis()->SetTitle("Theta Angular Resolution(degrees)");
Hpol_noOffAxisDelay_2->GetXaxis()->SetTitle("SNR");


wais->Draw("FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0):peak[1][0].snr>>Vpol_noOffAxisDelay(30,,,100,-5,5)","flags.pulser == 5 && abs(FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0))<5","colz")
Vpol_noOffAxisDelay->FitSlicesY(0,0,-1,10,"QNRSG5")
Vpol_noOffAxisDelay_2->Draw()
Vpol_noOffAxisDelay_2->GetYaxis()->SetTitle("Theta Angular Resolution(degrees)");
Vpol_noOffAxisDelay_2->GetXaxis()->SetTitle("SNR");


withOffAxis 30015:

wais->Draw("FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0):peak[0][0].snr>>Hpol_withOffAxisDelay(30,,,100,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].theta - wais.theta, 360, 0))<5","colz")
Hpol_withOffAxisDelay->FitSlicesY(0,0,-1,10,"QNRSG5")
Hpol_withOffAxisDelay_2->Draw()
Hpol_withOffAxisDelay_2->GetYaxis()->SetTitle("Theta Angular Resolution(degrees)");
Hpol_withOffAxisDelay_2->GetXaxis()->SetTitle("SNR");


wais->Draw("FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0):peak[1][0].snr>>Vpol_withOffAxisDelay(30,,,100,-5,5)","flags.pulser == 5 && abs(FFTtools::wrap(peak[1][0].theta - wais.theta, 360, 0))<5","colz")
Vpol_withOffAxisDelay->FitSlicesY(0,0,-1,10,"QNRSG5")
Vpol_withOffAxisDelay_2->Draw()
Vpol_withOffAxisDelay_2->GetYaxis()->SetTitle("Theta Angular Resolution(degrees)");
Vpol_withOffAxisDelay_2->GetXaxis()->SetTitle("SNR");


#include "macro/cuts.C"

wais->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].impulsivityMeasure>>none_blastH(300,,,300,,)",blastCut,"colz")
wais->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].impulsivityMeasure>>blastH(300,,,300,,)",!notBlast,"")
none_blastH->GetXaxis()->SetTitle("impulsivityMeasure")
none_blastH->GetYaxis()->SetTitle("mapPeak")

wais->Draw("summary->peak[1][0].value:summary->coherent_filtered[1][0].impulsivityMeasure>>none_blastV(300,,,300,,)",blastCut,"colz")
wais->Draw("summary->peak[1][0].value:summary->coherent_filtered[1][0].impulsivityMeasure>>blastV(300,,,300,,)",!notBlast,"")
none_blastV->GetXaxis()->SetTitle("impulsivityMeasure")
none_blastV->GetYaxis()->SetTitle("mapPeak")

wais->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].impulsivityMeasure>>all_h(300,,,300,,)","flags.pulser == 1","colz")
wais->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].impulsivityMeasure>>run45deg_h(300,,,300,,)","flags.pulser == 1 && run>128 && run<133","colz")
wais->Draw("summary->peak[0][0].value:summary->coherent_filtered[0][0].impulsivityMeasure>>norun45deg_h(300,,,300,,)","flags.pulser == 1 && (run<129 || run>132)","colz")


wais->Draw("summary->peak[1][0].value:summary->coherent_filtered[1][0].impulsivityMeasure>>all_v(300,,,300,,)","flags.pulser == 5","colz")
wais->Draw("summary->peak[1][0].value:summary->coherent_filtered[1][0].impulsivityMeasure>>run45deg_v(300,,,300,,)","flags.pulser == 5 && run>128 && run<133","colz")
wais->Draw("summary->peak[1][0].value:summary->coherent_filtered[1][0].impulsivityMeasure>>norun45deg_v(300,,,300,,)","flags.pulser == 5 && (run<129 || run>132)","colz")




wais->Draw("coherent_filtered[0][0].linearPolFrac():realTime>>h1(1000,,,300,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")
wais->Draw("mostImpulsiveCoherentFiltered().linearPolFrac():realTime>>h1(1000,,,300,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")
wais->Draw("mostImpulsivePol(2):realTime>>h1(1000,,,300,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")
wais->Draw("mostImpulsiveInd():realTime>>h1(1000,,,300,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")
wais->Draw("altitude:run>>h1(1000,,,300,-5,5)","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")



wais->Draw("deconvolved_filtered[0][0].peakTime - deconvolved_filtered[1][0].peakTime","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5","colz")


TestTree->Draw("Fisher","classID==1","")
TestTree->Draw("Fisher","classID==0","same")


wais->Draw("coherent_filtered[0][0].linearPolFrac()","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5 && (run<129 || run>132)","colz")
wais->Draw("-1*mostImpulsivePeak().theta","flags.pulser == 1 && abs(FFTtools::wrap(peak[0][0].phi - wais.phi, 360, 0))<5 && (run<129 || run>132)","colz")


#include "macro/cuts.C"
wais->Draw("coherent_filtered[0][0].bandwidth[0]>>Thermal(100,,)",thermal_sample,"colz")
wais->Draw("coherent_filtered[0][0].bandwidth[0]>>Wais(100,,)",isWaisV,"colz")
wais->Draw("coherent_filtered[0][0].bandwidth[0]>>Blast(100,,)",!notBlast,"colz")

wais->Draw("peak[1][0].value>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("peak[1][0].value>>Wais(300,,)",isWaisV,"colz")
wais->Draw("peak[1][0].value>>Blast(300,,)",!notBlast,"colz")


wais->Draw("coherent_filtered[0][0].peakHilbert>>Thermal(100,,)",thermal_sample,"colz")
wais->Draw("coherent_filtered[0][0].peakHilbert>>Wais(100,,)",isWaisV,"colz")
wais->Draw("coherent_filtered[0][0].peakHilbert>>Blast(100,,)",!notBlast,"colz")



wais->Draw("peak[1][0].value/coherent_filtered[0][0].peakHilbert>>Thermal(100,,)",thermal_sample,"colz")
wais->Draw("peak[1][0].value/coherent_filtered[0][0].peakHilbert>>Wais(100,,)",isWaisV,"colz")
wais->Draw("peak[1][0].value/coherent_filtered[0][0].peakHilbert>>Blast(100,,)",!notBlast,"colz")


wais->Draw("channels[0][].peakHilbert>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("channels[0][].peakHilbert>>Wais(300,,)",isWaisV,"colz")
wais->Draw("channels[0][].peakHilbert>>Blast(300,,)",!notBlast,"colz")


wais->Draw("eventNumber:Iteration$","(flags.maxBottomToTopRatio[0]>2.7) && (flags.maxBottomToTopRatio[1]>2.7) ","")

wais->Scan("eventNumber:Iteration$",flags.meanPower[3]/flags.meanPower[1] && "flags.meanPower[3]/flags.meanPower[1]>2")
blastCut original:
wais->Scan("eventNumber:Iteration$","flags.meanPowerFiltered[0]>1000000 || flags.medianPowerFiltered[0]>1000000 || flags.maxBottomToTopRatio[0]>5|| flags.maxBottomToTopRatio[1]>5|| flags.maxBottomToTopRatio>2.7|| flags.maxBottomToTopRatio<1 ","")

peng's blast cut:
wais->Scan("eventNumber:Iteration$","flags.meanPowerFiltered[0]>1000000 || flags.medianPowerFiltered[0]>1000000 || flags.maxBottomToTopRatio>2.7|| flags.maxBottomToTopRatio<0.9 ","")

wais->Scan("eventNumber:run:flags.meanPower[3]/flags.meanPower[1]:flags.maxBottomToTopRatio[0]:flags.maxBottomToTopRatio[1]",isWais &&"flags.meanPower[3]/flags.meanPower[1]>1.5 ","")
wais->Draw("summary->peak[1][0].value:summary->coherent_filtered[0][0].peakHilbert>>h2(300,,,300,,)","flags.pulser==1 && flags.meanPower[3]/flags.meanPower[1]>1.5  && flags.hasGlitch == 0","colz")



wais->Draw("flags.meanPowerFiltered[0]>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("flags.meanPowerFiltered[0]>>Wais(300,,)",isWaisV,"colz")
wais->Draw("flags.meanPowerFiltered[0]>>Blast(3000,,)",!notBlast,"colz")

wais->Draw("flags.medianPowerFiltered[0]>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("flags.medianPowerFiltered[0]>>Wais(300,,)",isWaisV,"colz")
wais->Draw("flags.medianPowerFiltered[0]>>Blast(3000,,)",!notBlast,"colz")


wais->Draw("flags.meanPowerFiltered[2]>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("flags.meanPowerFiltered[2]>>Wais(300,,)",isWaisV,"colz")
wais->Draw("flags.meanPowerFiltered[2]>>Blast(3000,,)",!notBlast,"colz")

wais->Draw("flags.medianPowerFiltered[2]>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("flags.medianPowerFiltered[2]>>Wais(300,,)",isWaisV,"colz")
wais->Draw("flags.medianPowerFiltered[2]>>Blast(3000,,)",!notBlast,"colz")




wais->Draw("flags.maxBottomToTopRatio[0]>>Thermal(300,0,14)",thermal_sample,"colz")
wais->Draw("flags.maxBottomToTopRatio[0]>>Wais(300,0,14)",isWaisV,"colz")
wais->Draw("flags.maxBottomToTopRatio[0]>>Blast(300,0,14)",!notBlast,"colz")

wais->Draw("flags.maxBottomToTopRatio[1]>>Thermal(300,0,14)",thermal_sample,"colz")
wais->Draw("flags.maxBottomToTopRatio[1]>>Wais(300,0,14)",isWaisV,"colz")
wais->Draw("flags.maxBottomToTopRatio[1]>>Blast(300,0,14)",!notBlast,"colz")

wais->Draw("flags.maxBottomToTopRatio>>Thermal(300,0,14)",thermal_sample,"colz")
wais->Draw("flags.maxBottomToTopRatio>>Wais(300,0,14)",isWaisV,"colz")
wais->Draw("flags.maxBottomToTopRatio>>Blast(300,0,14)",!notBlast,"colz")





wais->Draw("flags.nSectorsWhereBottomExceedsTop>>Thermal(300,0,14)",thermal_sample,"colz")
wais->Draw("flags.nSectorsWhereBottomExceedsTop>>Wais(300,0,14)",isWaisV,"colz")
wais->Draw("flags.nSectorsWhereBottomExceedsTop>>Blast(300,0,14)",!notBlast,"colz")




wais->Draw("flags.meanPowerFiltered[3]/flags.meanPowerFiltered[1]>>Thermal(1000,0,40)",thermal_sample,"colz")
wais->Draw("flags.meanPowerFiltered[3]/flags.meanPowerFiltered[1]>>Wais(1000,0,40)",isWaisV,"colz")
wais->Draw("flags.meanPowerFiltered[3]/flags.meanPowerFiltered[1]>>Blast(1000,0,40)",!notBlast,"colz")


wais->Draw("flags.meanPower[0]>>Blast(300,,)",!notBlast,"colz")
wais->Draw("flags.meanPower[1]>>Blast1(300,,)",!notBlast,"colz")
wais->Draw("flags.meanPower[2]>>Blast2(300,,)",!notBlast,"colz")
wais->Draw("flags.meanPower[3]>>Blast3(300,,)",!notBlast,"colz")


wais->Draw("peak[1][0].value/peak[1][0].mapRMS>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("peak[1][0].value/peak[1][0].mapRMS>>Wais(300,,)",isWaisV,"colz")
wais->Draw("peak[1][0].value/peak[1][0].mapRMS>>Blast(300,,)",!notBlast,"colz")


wais->Draw("peak[1][0].snr>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("peak[1][0].snr>>Wais(300,,)",isWaisV,"colz")
wais->Draw("peak[1][0].snr>>Blast(300,,)",!notBlast,"colz")



wais->Draw("value[0][0]/value[1][0]*mapRMS[1][0]/mapRMS[0][1]+value[1][0]/value[0][0]*mapRMS[0][0]/mapRMS[1][0]>>Thermal(300,,)",isReal && brightestPeak && !isWais && !aboveHorizon && isntSidelobe && badReconstructionCut && blastCut && triggered && notMasked,"colz")
wais->Draw("value[0][0]/value[1][0]*mapRMS[1][0]/mapRMS[0][1]+value[1][0]/value[0][0]*mapRMS[0][0]/mapRMS[1][0]>>Wais(300,,)",isWaisV,"colz")
wais->Draw("value[0][0]/value[1][0]*mapRMS[1][0]/mapRMS[0][1]+value[1][0]/value[0][0]*mapRMS[0][0]/mapRMS[1][0]>>Blast(300,,)",!notBlast,"colz")


wais->Draw("peak[1][0].value/peak[1][0].mapRMS>>Thermal(300,,)",thermal_sample,"colz")
wais->Draw("peak[1][0].value/peak[1][0].mapRMS>>Wais(300,,)",isWaisV,"colz")
wais->Draw("peak[1][0].value/peak[1][0].mapRMS>>Blast(300,,)",!notBlast,"colz")


wais->Draw("channels[][].rms:eventNumber>>Thermal(1000,,,300,,)",thermal_sample,"colz")
wais->Draw("channels[][].rms:eventNumber>>Wais(300,,,300,,)",isWaisV,"colz")
wais->Draw("channels[][].rms:eventNumber>>Blast(1000,,,300,,)",!notBlast,"colz")



wais->Draw("countChannelAboveThreshold(20)>>Thermal(100,,)",thermal_sample,"colz")
wais->Draw("countChannelAboveThreshold(20)>>Wais(100,,)",isWaisV,"colz")
wais->Draw("countChannelAboveThreshold(20)>>Blast(100,,)",!notBlast,"colz")

wais->Scan("run:eventNumber:countChannelAboveThreshold(30)",thermal_sample &&"countChannelAboveThreshold(30)>15","colz")


wais->Draw("run","flags.meanPower[3]/flags.meanPower[1]>0.2 && flags.meanPower[3]/flags.meanPower[1]<1.7 ","colz")
wais->Draw("run","flags.maxBottomToTopRatio > 2.6 || flags.maxBottomToTopRatio < 0.9","colz")
wais->Draw("run","flags.maxBottomToTopRatio > 2.6 || flags.maxBottomToTopRatio < 0.9 || flags.meanPower[3]/flags.meanPower[1]<0.2 || flags.meanPower[3]/flags.meanPower[1]>1.7 ","colz")
wais->Draw("run","countChannelAboveThreshold()>30","colz")

wais->Draw("mostImpulsiveInd():eventNumber>>waisV(1000,,,10,,)",isWaisV,"colz")
wais->Draw("FFTtools::wrap(peak[1][0].phi-wais.phi,360,0):FFTtools::wrap(peak[1][0].theta-wais.theta,360,0)>>waisV(300,,,300,,)",isWaisV && "mostImpulsiveInd()>0","colz")
wais->Scan("run:eventNumber",isWaisV && "mostImpulsiveInd()>0 && run == 140","colz")


wais->Draw("peak[1][0].snr>>Thermal(100,,)",thermal_sample && aboveHorizon,"colz")


map_unweighted->segmentationScheme()->Draw("colz",map_unweighted->getProbSums(false));
map_unweighted->segmentationScheme()->Draw("colz",map_unweighted->getOccludedFractionSum());
map_unweighted->segmentationScheme()->Draw("colz",map_unweighted->getUniformPS());
map_unweighted->segmentationScheme()->Draw("colz",map_unweighted->getBaseWeightedUniformPS());
map_unweighted->segmentationScheme()->Draw("colz",map_unweighted->getUniformPSwithBase());
map_unweighted->segmentationScheme()->Draw("colz",map_unweighted->getUniformPSwithoutBase());

map_unweighted->getProbSumsIntegral(true)
tmp_copy->Integral()
drawGroupings(map_unweighted)
map_unweighted->showClusters(1,0,"mapcolz")

trend->Draw("avgNSinglets:numOfFiles>>h(20,,,20,,)","","*")
trend->Draw("N_singlets:percentOfData>>h(100,0,1,100,,)","","*")

TestTree->Draw("LD0>>h1(300,-2,3)","classID==1","")
TestTree->Draw("LD0>>h2(300,-2,3)","(classID==0)*(weight)","same")

#include "RampdemReader.h"
RampdemReader::getMap(RampdemReader::surface, 10)->Draw("colz"); 


TProfile2D * bg = RampdemReader::getMap(RampdemReader::surface, 10); 


#include "TGraphAntarctica.h"
#include "BaseList.h"
base =  BaseList::getBase(2)
TGraphAntarctica* plot = new TGraphAntarctica(base)
plot->Draw()

BaseList::getPath(1).getPosition(0)

TGraphAntarctica::makeGpsGraph(140,140)->Draw()


#include "BaseList.h"
#include "UsefulAdu5Pat.h"
.x drawPathsCloseToAnita.C


TestTree->Scan("run:evnum1:evnum2:Fisher2:mapPeak:deconvHilbertPeak:deconvImpulsivity:deconvLinearPolFraction:secondPeakRatio","Fisher2>1.5 && classID==1","")


anita4->Draw("coherent[1][0].totalPower - coherent_filtered[1][0].totalPower","","")
anita4->Draw("deconvolved[1][0].totalPowerXpol:deconvolved_filtered[1][0].totalPowerXpol>>h(300,,,300,,)","deconvolved[1][0].totalPowerXpol<100000 && deconvolved_filtered[1][0].totalPowerXpol<100000","colz")
anita4->Draw("coherent[0][0].peakPower:coherent_filtered[0][0].peakPower>>h(300,,,300,,)","coherent[0][0].peakPower<100000 && coherent_filtered[0][0].peakPower<100000","colz")
anita4->Draw("coherent_filtered[1][0].bandwidth:coherent_filtered[1][0].peakFrequency>>h(300,,,300,,)","","colz")

anita4->Draw("coherent[1][0].peakHilbert:coherent_filtered[1][0].peakHilbert>>h(300,,,300,,)","","colz")

anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure>>h(300)","flags.isPayloadBlast == 1","colz")




#include "macro/cuts.C"
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure>>ThermalAbove(100,,)",thermal_sample && aboveHorizon,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure>>ThermalBelow(100,,)",thermal_sample && !aboveHorizon,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure>>Blast(100,,)",!notBlast,"colz")

anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>blast(300,0,1,300,0,1)",!notBlast,"colz")

#include "macro/cuts.C"
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h0v0(300,0,1,300,0,1)",thermal_sample && h0 && v0,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h0v1(300,0,1,300,0,1)",thermal_sample && h0 && v1,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h0v2(300,0,1,300,0,1)",thermal_sample && h0 && v2,"colz")

anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h1v0(300,0,1,300,0,1)",thermal_sample && h1 && v0,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h1v1(300,0,1,300,0,1)",thermal_sample && h1 && v1,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h1v2(300,0,1,300,0,1)",thermal_sample && h1 && v2,"colz")

anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h2v0(300,0,1,300,0,1)",thermal_sample && h2 && v0,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h2v1(300,0,1,300,0,1)",thermal_sample && h2 && v1,"colz")
anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>h2v2(300,0,1,300,0,1)",thermal_sample && h2 && v2,"colz")

#include "macro/cuts.C"
TChain anita4("anita4");
anita4.Add("*30001*.root");
anita4.Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>belowHorizon(300,0,1,300,0,1)",thermal_sample && belowHorizon,"colz")
anita4.Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>aboveHorizontal(300,0,1,300,0,1)",thermal_sample && aboveHorizontal,"colz")
anita4.Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>betweenHorizonAndHorizontal(300,0,1,300,0,1)",thermal_sample && betweenHorizonAndHorizontal,"colz")





anita4->Scan("run:eventNumber",thermal_sample && aboveHorizon && "deconvolved_filtered[0][0].impulsivityMeasure>0.8","colz")
anita4->Scan("run:eventNumber",thermal_sample && aboveHorizon,"colz")

anita4->Draw("-1*peak.theta[1][0]:-1*peak.theta[0][0]>>theta(500,,,500,,)","","colz")


anita4->Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure>>thermal_abovehoriz(300,0,1,300,0,1)",thermal_sample && "peak[1][0].theta > 6 && peak[0][0].theta > 6","colz")


anita4->Scan("run:eventNumber",thermal_sample && h1 && v1 && "deconvolved_filtered[0][0].impulsivityMeasure>0.7","colz")
anita4->Scan("run:eventNumber",thermal_sample && h0 && v0 && "deconvolved_filtered[0][0].impulsivityMeasure<0.75 && deconvolved_filtered[0][0].impulsivityMeasure>0.6 && deconvolved_filtered[1][0].impulsivityMeasure<0.75 && deconvolved_filtered[1][0].impulsivityMeasure>0.6","colz")



#include "/Users/sylarcp/anitaNeutrino/anitaBuildTool/components/UCorrelator/macro/cuts.C"
TChain chain3("simulation");
chain3.Add("*1000*.root");

skymap
Vpol:
anita4->Draw("-1*peak[1][0].theta:FFTtools::wrap(peak[1][0].phi - pat.Adu5Pat.heading,360,180)>>h(300,,,300,,)","","colz")
anita4->Draw("-1*peak[0][0].theta:FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(300,,,300,,)","abs(FFTtools::wrap(peak[0][0].phi-sun.phi,360,0))<10","colz")
anita4->Draw("pat.Adu5Pat.longitude:FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(300,,,300,,)","peak[0][0].theta>6","colz")
anita4->Draw("pat.Adu5Pat.longitude:FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(300,,,300,,)","peak[0][0].theta<0","colz")
anita4->Draw("pat.Adu5Pat.longitude:FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(300,,,300,,)","abs(FFTtools::wrap(peak[0][0].phi-sun.phi,360,0))<10","colz")


anita4.Draw("sqrt(TMath::Power(mostImpulsiveDeconvolvedFiltered(2).Q,2) + TMath::Power(mostImpulsiveDeconvolvedFiltered(2).U,2))/mostImpulsiveDeconvolvedFiltered(2).I:run>>h(300,,,300,,)","flags.pulser==1","colz")

anita4.Draw("sqrt(TMath::Power(deconvolved_filtered[0][0].Q,2) + TMath::Power(deconvolved_filtered[0][0].U,2))/deconvolved_filtered[0][0].I:run>>h(300,,,300,,)","flags.pulser==1","colz")

anita4.Draw("sqrt(TMath::Power(deconvolved_filtered[1][0].Q,2) + TMath::Power(deconvolved_filtered[1][0].U,2))/deconvolved_filtered[1][0].I:run>>h(300,,,300,,)","flags.pulser==1","colz")

anita4.Draw("deconvolved_filtered[0][0].V/deconvolved_filtered[0][0].I:realTime>>h(300,,,300,,)","flags.pulser==1 &&run>120&&run<155","colz")

anita4.Draw("deconvolved_filtered[1][0].V/deconvolved_filtered[1][0].I:realTime>>h(300,,,300,,)","flags.pulser==1 &&run>120&&run<155","colz")


anita4.Draw("deconvolved_filtered[0][0].I:realTime>>h(300,,,300,,)","flags.pulser==1 &&run>120&&run<155","colz")
anita4.Draw("deconvolved_filtered[1][0].I:realTime>>h(300,,,300,,)","flags.pulser==1 &&run>120&&run<155","colz")
anita4.Draw("FFTtools::wrap(wais.phi - pat.Adu5Pat.heading,360,0):FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(200,,,200,,)","flags.pulser==1","colz")

anita4.Draw("FFTtools::wrap(wais.phi - pat.Adu5Pat.heading,360,0):FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(200,,,200,,)","abs(FFTtools::wrap(peak[0][0].phi-sun.phi,360,0))>10","colz")


anita4.Draw("pat.Adu5Pat.longitude:FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h0(300,,,300,,)","(((flags.pulser==1)*100+ 1)/100) * (abs(FFTtools::wrap(peak[0][0].phi-sun.phi,360,0))>10)","colz")

anita4.Draw("pat.Adu5Pat.longitude:FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(300,,,300,,)","flags.pulser==1","same")


anita4.Draw("pat.Adu5Pat.longitude:FFTtools::wrap(peak[0][0].phi - pat.Adu5Pat.heading,360,180)>>h(300,,,300,,)","","colz")
 anita4.Draw("FFTtools::wrap(pat.Adu5Pat.longitude-(peak[0][0].phi - pat.Adu5Pat.heading),360,180):realTime","abs(FFTtools::wrap(peak[0][0].phi-sun.phi,360,0))>10","colz")

#include "/Users/sylarcp/anitaNeutrino/anitaBuildTool/components/UCorrelator/macro/cuts.C"
wais.Draw("FFTtools::wrap(summary->peak[0][0].theta-summary->wais.theta,360,0): FFTtools::wrap(summary->peak[0][0].phi - summary->wais.phi,360,0)>>h(300,,,300,,)","wais_sample && isWaisH","colz")
w.Draw("FFTtools::wrap(mostImpulsivePeak().theta-summary->wais.theta,360,0): FFTtools::wrap(mostImpulsivePeak().phi - summary->wais.phi,360,0)>>h(300,,,300,,)",wais_sample,"colz")

"

wais.Draw("deconvolved_filtered[1][0].impulsivityMeasure:deconvolved_filtered[0][0].impulsivityMeasure","run>200 && (flags.pulser == 1)","colz")

anita4.Draw("FFTtools::wrap(pat.Adu5Pat.longitude-(peak[0][0].phi - pat.Adu5Pat.heading),360,180):run>>HpolImpulsiveEvents(366,,,300,,)","deconvolved_filtered[0][0].impulsivityMeasure>0.75","colz")
anita4.Draw("FFTtools::wrap(pat.Adu5Pat.longitude-(peak[0][0].phi - pat.Adu5Pat.heading),360,180):run>>HpolThermalEvents(366,,,300,,)","deconvolved_filtered[0][0].impulsivityMeasure<0.75","colz")
anita4.Draw("FFTtools::wrap(pat.Adu5Pat.longitude-(peak[0][0].phi - pat.Adu5Pat.heading),360,180):run>>HpolImpulsiveEvents(368,0,368,300,,)","","colz")


icemc data:
eventTree->Draw("fVolts[][]>>mcEvents(300,,)","","colz")


gStyle->SetOptFit();
TF1 *snrFit = new TF1("snrFit","[0]/x^[1]+[2]",1,35);
snrFit->SetParNames("c1","c2","c3");
snrFit->SetParameters(0.5,0.2,0.1);

dPhivsSNR_2->Fit("snrFit")

dThetavsSNR_2->Fit("snrFit")


anita4.Draw("deconvolved_filtered[0][0].impulsivityMeasure:coherent_filtered[0][0].impulsivityMeasure>>hpol(300,0,1,300,0,1)","(-1*summary->peak[0][0].theta + summary->sun.theta>-2)&&(-1*summary->peak[0][0].theta + summary->sun.theta<2)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)>-6)&&(FFTtools::wrap(summary->peak[0][0].phi-summary->sun.phi,360,0)<6)&&(summary->flags.isPayloadBlast != 1)","colz")

hpol->GetXaxis()->SetTitle("coherent Impulsivity")
hpol->GetYaxis()->SetTitle("deconvolved impulsivity")


//when swtich to hical2 branch , need to restart the terminal to get it work!!!!!!
anita4->Draw("realTime","Hical2::isHical(eventNumber, FFTtools::wrap(anitaLocation.heading - peak[0][0].phi, 360))!=1","colz")

TChain a4("anita4")
a4.Add("207*")

 a4.Draw("Hical2::isHical(summary):eventNumber>>h(1000,,,3,,)","","colz")

 a4.Draw("Hical2::isHical(eventNumber,header->triggerTime,FFTtools::wrap(anitaLocation.heading - peak[0][0].phi, 360)):eventNumber>>h(1000,,,3,,)","","colz")
 a4.Draw("Hical2::isHical(eventNumber,FFTtools::wrap(anitaLocation.heading - peak[0][0].phi, 360)):eventNumber>>h(1000,,,3,,)","","colz")
 a4.Draw("Hical2::isHical2(summary,header.triggerTime,header.triggerTimeNs,pat.longitude,pat.latitude):eventNumber>>h(1000,,,3,,)","","colz")

 a4.Draw("peak[0][0].theta:peak[0][0].phi","deconvolved_filtered[0][0].impulsivityMeasure > 0.75","colz")

anita4->Draw("deconvolved_filtered[0][0].linearPolAngle():deconvolved_filtered[1][0].linearPolAngle()","","colz")
anita4->Draw("deconvolved_filtered[0][0].I:deconvolved_filtered[1][0].I","","colz")
anita4->Draw("deconvolved_filtered[0][0].linearPolFrac():deconvolved_filtered[1][0].linearPolFrac()","","colz")
anita4->Draw("deconvolved_filtered[0][0].impulsivityMeasure:deconvolved_filtered[1][0].impulsivityMeasure","","colz")
anita4->Draw("deconvolved_filtered[0][0].V:deconvolved_filtered[1][0].V","","colz")


//gaussian fit
gStyle->SetOptFit();
TF1 *gaussion = new TF1("gaussion","[0]*exp(-x^2/(2*[1]^2))",-5,5);
gaussion->SetParNames("A","sigma");
gaussion->SetParameters(100,0.5);
dThetavsdPhi_px->Fit("gaussion")
dThetavsdPhi_py->Fit("gaussion")

//exponential fit:
gStyle->SetOptFit();
TF1 *exponential1 = new TF1("exponential1","[0]*exp(-[1]*abs(x))",-5,5);
exponential1->SetParNames("A","sigma");
exponential1->SetParameters(100,0.5);
dThetavsdPhi_px->Fit("exponential1")
dThetavsdPhi_py->Fit("exponential1")


//gaussian + exponential fit:
gStyle->SetOptFit();
TF1 *mixedFit1 = new TF1("mixedFit1","[0]*exp(-x^2/(2*[1]^2)) + [0]*exp(-[2]*abs(x))",-5,5);
mixedFit1->SetParNames("A","sigma","k");
mixedFit1->SetParameters(100,0.5,1);
dThetavsdPhi_px->Fit("mixedFit1")
dThetavsdPhi_py->Fit("mixedFit1")



//trend of singlets
TCanvas * dists = new TCanvas("c1","c1"); 
trend->Draw("N_singlets:percentOfData >> h0(100,0,1.1,100,0,10)","", "C*"); 
trend->Draw("N_singlets:percentOfData >> h1(100,0,1.1,100,0,10)","", "sameC*"); 
trend->Draw("N_singlets_nearbase:percentOfData >> h2(100,0,1.1,100,0,10)","", "sameC*"); 
trend->Draw("N_singlets_notnearbase:percentOfData >> h3(100,0,1.1,100,0,10)","", "sameC*"); 

auto legend = new TLegend(0.1,0.7,0.48,0.9);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
legend->AddEntry("h1","N_singlets","l");
legend->AddEntry("h2","N_singlets_nearbase","l");
legend->AddEntry("h3","N_singlets_notnearbase","l");
legend->Draw();

TChain a4("anita4");
a4.Add("*30001*");
#include "/Users/sylarcp/anitaNeutrino/anitaBuildTool/components/UCorrelator/macro/cuts.C";
a4.Draw("mostImpulsiveCoherent(2).snr/mostImpulsiveDeconvolved(2).snr>>h(300)",thermal_sample,"colz");
h->GetXaxis()->SetTitle("coherent_snr/deconvolved_snr")



///check hical is on or not:
double a 
double b
Hical2::angleToHical(68160718, &a, &b)
if a or b not = -999, means hical a or b is on.
And the peak phi - heading should be very close to this a or b.
anita4->Scan("FFTtools::wrap(anitaLocation.heading - peak[0][0].phi,360)","eventNumber == 68160718")




TChain a4("anita4")
a4.Add("*30002*")
a4.Draw("theta:deconvImpulsivity>>h(300,,,300,,)","pol==1 && deconvImpulsivity<0.7","colz")
h->FitSlicesX()
h_2->Draw("colz")
a4.Draw("powerV - powerH:impulsivityV- impulsivityH>>h(300,,,300,,)","theta < -6 && impulsivity>0.7","colz")


cluster->Scan("pol:n:base:noBase:H:V:avgLinearPolFrac:avgLinearPolAngle:longitude:latitude:powerH:powerV","","colsize=7 precision=7 col=::::::20.3:20.3:10.10:10.10:7.3:7.3")
cluster->Scan("pol:n:base:noBase:H:V:avgLinearPolFrac:avgLinearPolAngle:longitude:latitude:powerH:powerV","n==1","colsize=7 precision=7 col=::::::20.3:20.3:10.10:10.10:7.3:7.3")

events->Scan("run:event:pol:nsegs:NOverlapedBases:impulsivityV:impulsivityH:powerV:powerH:linearPolFrac:linearPolAngle:longitude:latitude","indexOfCluster==6")
 
convert -density 2000 Clusters.eps -flatten Clusters.png;


auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("h1","x=0.001","l");
   legend->AddEntry("h2","x=1","l");
   legend->AddEntry("h3","MinBias Energy222 MC","l");
   legend->AddEntry("h4","Wais data","l");
   legend->Draw();

near Base
x=2
cluster->Draw("n>>h1","n<10 && base==0")
x=3.5
cluster->Draw("n>>h2","n<10 && base==0","same")
x=4.5
cluster->Draw("n>>h3","n<10 && base==0","same")
h1->GetXaxis()->SetTitle("clusterSize")

auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("h1","x=2, noBase","l");
   legend->AddEntry("h2","x=3.5, noBase","l");
   legend->AddEntry("h3","x=4.5, noBase","l");
legend->Draw();


x=2
cluster->Draw("n>>h4","n<10 && base!=0")
x=3.5
cluster->Draw("n>>h5","n<10 && base!=0")
x=4.5
cluster->Draw("n>>h6","n<10 && base!=0")
h4->GetXaxis()->SetTitle("clusterSize")

h4->Draw()
h5->Draw("same")
h6->Draw("same")
auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("h4","x=2, nearBase","l");
   legend->AddEntry("h5","x=3.5, nearBase","l");
   legend->AddEntry("h6","x=4.5, nearBase","l");
   legend->Draw();




auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("h6","x=0.0001, minbias MC Kotera Max","l");
   legend->AddEntry("h7","x=0.0001, below horizon thermal events passing quality cuts","l");
   legend->Draw();

  h1->GetXaxis()->SetTitle("clusterSize")
  fit


gStyle->SetOptFit();
Double_t poissonf(Double_t*x,Double_t*par)                                         
{                                                                              
  return par[0];
}


cluster->Draw("n>>clusterSize","n<100 && base!=0")
gStyle->SetOptFit();
TF1 *fitPoisson = new TF1("fitPoisson","[0]*TMath::Power([1],x)/TMath::Gamma(x+1)/TMath::Exp([1])", 0, 100)
fitPoisson->SetParName(0,"p3");
fitPoisson->SetParName(1,"p4");
fitPoisson->SetParameter(0, 5);
fitPoisson->SetParameter(1, 3);
clusterSize->Fit("fitPoisson");
clusterSize->GetXaxis()->SetTitle("clusterSize")




cluster->Draw("n>>clusterSize(9,0,9)","n<6 && base!=0")
gStyle->SetOptFit();
TF1 *fitMaxwell3 = new TF1("fitMaxwell","[0]*(x-0.5)^2*TMath::Exp(-1*[1]*(x-0.5)^2)", 0, 6)
fitMaxwell3->SetParName(0,"c2");
fitMaxwell3->SetParName(1,"a2");
fitMaxwell3->SetParameter(0, 5);
fitMaxwell3->SetParameter(1, 3);
clusterSize->Fit("fitMaxwell");
clusterSize->GetXaxis()->SetTitle("clusterSize")





cluster->Draw("n>>clusterSize","n<100 ")
gStyle->SetOptFit();
TF1 *fitPow = new TF1("fitPow","[0]*TMath::Power(x,[1])", 0, 100)
fitPow->SetParName(0,"c1");
fitPow->SetParName(1,"a1");
fitPow->SetParameter(0, 5);
fitPow->SetParameter(1, -1);
clusterSize->Fit("fitPow");
clusterSize->GetXaxis()->SetTitle("clusterSize")





cluster->Draw("n>>noBaseClusterSize(9,0,9)","n<6 && base==0")
gStyle->SetOptFit();
TF1 *fitCombined = new TF1("fitCombined","[0]*TMath::Power(x,-2.3) + [1]*(x-0.5)^2*TMath::Exp(-0.21*(x-0.5)^2)", 0, 6)
fitCombined->SetParName(0,"c1");
fitCombined->SetParName(1,"c2");
fitCombined->SetParameter(0, 1);
fitCombined->SetParameter(1, 1);
noBaseClusterSize->Fit("fitCombined");
noBaseClusterSize->GetXaxis()->SetTitle("noBaseClusterSize")




cluster->Draw("n>>ClusterSize(9,0,9)","n>1 && n<6")
gStyle->SetOptFit();
TF1 *fitCombined = new TF1("fitCombined","[1]*x^2*TMath::Exp(-0.32*x^2)", 0, 10)
fitCombined->SetParName(0,"c1");
fitCombined->SetParName(1,"c2");
fitCombined->SetParameter(0, 1);
fitCombined->SetParameter(1, 1);
ClusterSize->Fit("fitCombined");
ClusterSize->GetXaxis()->SetTitle("ClusterSize")
